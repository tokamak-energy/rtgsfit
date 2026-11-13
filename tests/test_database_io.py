"""Local files and generated C: round trips and validation."""

from __future__ import annotations

import copy
import re
import subprocess
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit import library
from rtgsfit.database_readers.local import read_preshot
from rtgsfit.database_readers.local import read_results
from rtgsfit.database_readers.rtgsfit_mdsplus import read_mdsplus_nodes
from rtgsfit.database_writers.c_source import preshot_checksum
from rtgsfit.database_writers.c_source import write_preshot_file
from rtgsfit.database_writers.local import LocalWriter
from rtgsfit.database_writers.local import write_npz
from rtgsfit.library import Rtgsfit
from rtgsfit.library import build_directory
from rtgsfit.library import build_library
from rtgsfit.library import c_directory
from rtgsfit.library import cache_directory
from rtgsfit.library import run_make
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import PreshotError
from rtgsfit.replay import build_preshot
from rtgsfit.results import Results
from rtgsfit.settings import load_settings


def _assert_same(actual: object, expected: object, key: str) -> None:
    if isinstance(expected, (str, list)):
        assert actual == expected, key
    else:
        assert_array_equal(actual, expected, err_msg=key)
        assert np.asarray(actual).dtype == np.asarray(expected).dtype, key


def test_results_round_trip(analytic_results: Results, tmp_path: Path) -> None:
    path = write_npz(tmp_path / "results.npz", analytic_results.nested)
    again = read_results(path).flat
    expected = analytic_results.flat
    assert set(again) == set(expected)
    for key, value in expected.items():
        _assert_same(again[key], value, key)


def test_preshot_round_trip(analytic_preshot: Preshot, tmp_path: Path) -> None:
    path = write_npz(tmp_path / "preshot.npz", analytic_preshot.to_nested())
    again = read_preshot(path)
    assert again.content_hash() == analytic_preshot.content_hash()
    assert again["FIT:TARGET_NAMES"] == analytic_preshot["FIT:TARGET_NAMES"]


def test_cache_directory_is_repository_local(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RTGSFIT_CACHE_DIR", raising=False)
    monkeypatch.delenv("RTGSFIT_C_DIR", raising=False)
    monkeypatch.chdir(tmp_path)
    assert cache_directory() == Path(__file__).resolve().parents[1] / "c" / "cache"


def test_cache_directory_follows_explicit_c_sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RTGSFIT_CACHE_DIR", raising=False)
    monkeypatch.setenv("RTGSFIT_C_DIR", str(tmp_path / "c"))
    assert cache_directory() == tmp_path / "c" / "cache"


def test_cache_directory_can_be_explicitly_overridden(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RTGSFIT_CACHE_DIR", str(tmp_path / "export"))
    assert cache_directory() == tmp_path / "export"


@pytest.fixture
def installed_package(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """A normal (non-editable) install: no checkout, and read-only C sources inside the package."""
    packaged = tmp_path / "site-packages" / "rtgsfit" / "_c"
    (packaged / "src").mkdir(parents=True)
    (packaged / "Makefile").write_text("all:\n")
    (packaged / "src" / "rtgsfit.c").write_text("int x;\n")
    paths = [packaged / "Makefile", packaged / "src" / "rtgsfit.c", packaged / "src", packaged]
    for path in paths:
        path.chmod(0o555 if path.is_dir() else 0o444)
    monkeypatch.setattr(library, "_checkout_c_directory", lambda: None)
    monkeypatch.setattr(library, "PACKAGED_C_DIRECTORY", packaged)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("RTGSFIT_C_DIR", raising=False)
    monkeypatch.delenv("RTGSFIT_CACHE_DIR", raising=False)
    yield packaged
    for path in paths:
        path.chmod(0o755 if path.is_dir() else 0o644)


def test_installed_package_builds_from_a_writable_copy(installed_package: Path, tmp_path: Path) -> None:
    c = c_directory()
    assert c.parent == tmp_path / "xdg" / "rtgsfit"
    assert (c / "src" / "rtgsfit.c").read_text() == "int x;\n"
    (c / "src" / "preshot_data.c").write_text("")  # make writes the active preshot here
    (c / "build").mkdir()
    assert c_directory() == c
    assert sorted(path.name for path in installed_package.rglob("*")) == ["Makefile", "rtgsfit.c", "src"]
    assert cache_directory() == tmp_path / "xdg" / "rtgsfit" / "outputs"


def test_installed_package_with_other_sources_builds_from_another_copy(installed_package: Path) -> None:
    before = c_directory()
    source = installed_package / "src" / "rtgsfit.c"
    for path in (installed_package, source.parent, source):
        path.chmod(0o755 if path.is_dir() else 0o644)
    source.write_text("int y;\n")
    after = c_directory()
    assert after != before
    assert (after / "src" / "rtgsfit.c").read_text() == "int y;\n"


def test_build_writes_the_active_preshot_and_reuses_it(analytic_preshot: Preshot) -> None:
    library = build_library(analytic_preshot)
    assert library == build_directory() / "librtgsfit.so"
    files = [library, c_directory() / "src" / "preshot_sizes.h", c_directory() / "src" / "preshot_data.c"]
    mtimes = [path.stat().st_mtime_ns for path in files]
    assert build_library(analytic_preshot) == library
    assert [path.stat().st_mtime_ns for path in files] == mtimes


def test_rebuilding_for_another_preshot_reloads_the_library(analytic_preshot: Preshot) -> None:
    def call(rtgsfit: Rtgsfit, preshot: Preshot) -> int:
        return rtgsfit(np.zeros(preshot.sizes["N_MEAS"]), np.zeros(preshot.sizes["N_COIL"])).status

    coarse = build_preshot(load_settings("analytic").copy_with({"analytic/grid/n_r": 33, "analytic/grid/n_z": 65}))
    before = Rtgsfit.build(analytic_preshot)
    after = Rtgsfit.build(coarse)
    assert after.sizes == coarse.sizes
    call(after, coarse)
    with pytest.raises(RuntimeError, match="rebuilt"):
        call(before, analytic_preshot)
    assert Rtgsfit.build(analytic_preshot).sizes == analytic_preshot.sizes


def test_library_of_another_preshot_is_rejected(analytic_preshot: Preshot) -> None:
    library = build_library(analytic_preshot)
    tolerance = analytic_preshot["SETTINGS:TOLERANCE"]
    other = Preshot(_modified(analytic_preshot, "SETTINGS:TOLERANCE", tolerance * 2))
    assert other.sizes == analytic_preshot.sizes
    with pytest.raises(RuntimeError, match="checksum"):
        Rtgsfit(library, other)


@pytest.mark.parametrize("absolute_paths", [False, True])
def test_local_writers_resolve_paths(
    analytic_preshot: Preshot,
    analytic_results: Results,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    absolute_paths: bool,
) -> None:
    cache = tmp_path / "cache"
    working_directory = tmp_path / "working"
    working_directory.mkdir()
    monkeypatch.setenv("RTGSFIT_CACHE_DIR", str(cache))
    monkeypatch.chdir(working_directory)
    settings = load_settings("analytic")
    if absolute_paths:
        output = tmp_path / "export"
        settings = settings.copy_with({
            "writers/local/path": str(output / "rtgsfit_analytic.npz"),
            "writers/local/preshot_path": str(output / "rtgsfit_preshot.npz"),
        })  # fmt: skip
    else:
        output = cache
    writer = LocalWriter()
    writer.write_preshot(analytic_preshot, settings)
    writer.write_results(analytic_results, settings)
    assert read_preshot(output / "rtgsfit_preshot.npz").content_hash() == analytic_preshot.content_hash()
    assert_array_equal(read_results(output / "rtgsfit_analytic.npz")["GLOBAL.IP"], analytic_results["GLOBAL.IP"])
    assert not list(working_directory.iterdir())


def _modified(preshot: Preshot, key: str, value: object) -> dict:
    groups = copy.deepcopy(preshot.groups)
    group, _, name = key.partition(":")
    groups[group][name] = value
    return groups


@pytest.mark.parametrize(
    ("key", "change"),
    [
        ("GRID:R", lambda value: value[:-1]),
        ("GRID:G_COIL", lambda value: value.astype(np.float32)),
        ("FIT:WEIGHT", lambda value: np.where(np.arange(len(value)) == 0, np.nan, value)),
        ("GRID:MASK_LIMITER", lambda value: value * 2),
        ("POISSON:PERMUTATION", lambda value: np.zeros_like(value)),
        ("LIMITER:GRID_INDEX", lambda value: value - 1_000_000),
        ("CONSTRAINT:NAMES", lambda value: value[::-1]),
    ],
)
def test_preshot_validation_rejects(analytic_preshot: Preshot, key: str, change) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(PreshotError):
        Preshot(_modified(analytic_preshot, key, change(analytic_preshot[key])))


def test_preshot_validation_rejects_unknown_node(analytic_preshot: Preshot) -> None:
    with pytest.raises(PreshotError):
        Preshot(_modified(analytic_preshot, "GRID:EXTRA", np.zeros(1)))


# ----------------------------------------------------------------------------------------------
# The C generator, c/tools/write_preshot.c, run into a temporary directory (not the active c/src)
# ----------------------------------------------------------------------------------------------
def _generator(target: str = "tools", name: str = "write_preshot") -> Path:
    run_make([target])
    return build_directory() / name


def _generate(preshot: Preshot, directory: Path) -> str:
    preshot_file = write_preshot_file(preshot, directory / "preshot.bin")
    command = [str(_generator()), "--file", str(preshot_file), str(directory)]
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


@pytest.fixture(scope="module")
def generated(analytic_preshot: Preshot, tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("generated")
    _generate(analytic_preshot, directory)
    return directory


def test_generator_lists_the_data_dictionary() -> None:
    """The generator's own table of nodes (c/tools/write_preshot.c) is `data_dictionary.C_NODES`."""
    listing = subprocess.run([str(_generator()), "--list"], check=True, capture_output=True, text=True).stdout
    expected = [" ".join([f"{node.path[1]}:{node.path[2]}", node.dtype, *node.shape]) for node in dd.C_NODES]
    assert listing.splitlines() == expected


def test_c_sizes_are_enum_constants(analytic_preshot: Preshot, generated: Path) -> None:
    """Issue #97: array sizes must be integer constant expressions, so they are enum members."""
    header = (generated / "preshot_sizes.h").read_text()
    entries = re.findall(r"^ +(\w+) += +(\d+), /\* (.*) \*/$", header, flags=re.MULTILINE)
    assert {name: int(value) for name, value, _ in entries} == analytic_preshot.sizes
    assert {name: description for name, _, description in entries} == dd.SIZE_DESCRIPTIONS


def test_c_arrays_have_the_dictionary_lengths(analytic_preshot: Preshot, generated: Path) -> None:
    source = (generated / "preshot_data.c").read_text()
    lengths = {name: int(length) for name, length in re.findall(r"const \w+ (\w+)\[(\d+)\] = \{", source)}
    sizes = analytic_preshot.symbolic_sizes()
    for node in dd.C_NODES:
        if node.shape:
            assert lengths[node.c_name] == int(np.prod([sizes[dim] for dim in node.shape])), node.c_name


def test_c_checksum_is_the_python_checksum(analytic_preshot: Preshot, generated: Path) -> None:
    match = re.search(r"PRESHOT_CHECKSUM = 0x([0-9A-F]{8})u;", (generated / "preshot_data.c").read_text())
    assert match is not None
    assert int(match.group(1), 16) == preshot_checksum(analytic_preshot)


def test_c_source_is_only_rewritten_when_it_changes(analytic_preshot: Preshot, tmp_path: Path) -> None:
    assert _generate(analytic_preshot, tmp_path).count("Wrote") == 2
    files = [tmp_path / "preshot_sizes.h", tmp_path / "preshot_data.c"]
    mtimes = [path.stat().st_mtime_ns for path in files]
    assert _generate(analytic_preshot, tmp_path).count("is unchanged") == 2
    assert [path.stat().st_mtime_ns for path in files] == mtimes
    assert sorted(path.name for path in tmp_path.iterdir()) == ["preshot.bin", "preshot_data.c", "preshot_sizes.h"]


@pytest.mark.te_network
def test_mdsplus_generator_reads_what_python_reads(tmp_path: Path) -> None:
    """`make SHOT=... RUN_NAME=...` (the PCS, no Python) generates the same C as the preshot file route.

    Reads the private RTGSFIT 52000230 TEST03 fixture, also used by the ST40 example.
    This test never writes to MDSplus or reads the shared PCS preshot.
    """
    pulse, run_name = 52_000_230, "TEST03"
    nodes = tuple(node for node in dd.PRESHOT_NODES if node.c_name is not None or node.path[1:] == ("INFO", "SOURCE"))
    groups = read_mdsplus_nodes(pulse, run_name, nodes)["PRESHOT"]

    class CNodes:
        def __getitem__(self, key: str) -> object:
            group, _, name = key.partition(":")
            return groups[group][name]

    from_file = tmp_path / "file"
    from_file.mkdir()
    _generate(CNodes(), from_file)  # type: ignore[arg-type]
    from_mdsplus = tmp_path / "mdsplus"
    from_mdsplus.mkdir()
    generator = _generator("tools-mdsplus", "write_preshot_mdsplus")
    subprocess.run([str(generator), "--mdsplus", "smaug", str(pulse), run_name, str(from_mdsplus)], check=True)
    for name in ("preshot_sizes.h", "preshot_data.c"):
        assert (from_mdsplus / name).read_text() == (from_file / name).read_text(), name
