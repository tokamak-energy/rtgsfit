"""The C core's own tests (`make -C c test`): finite differences, nulls, boundary, Poisson, LAPACK."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest
from rtgsfit.library import c_directory
from rtgsfit.library import run_make


@pytest.mark.skipif(shutil.which("make") is None, reason="needs make")
def test_c_unit_tests() -> None:
    completed = subprocess.run(["make", "-C", str(c_directory()), "test"], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.parametrize("existing_flags", [False, True])
def test_missing_cmp_fails_without_rewriting_flags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, existing_flags: bool
) -> None:
    binaries = tmp_path / "bin"
    binaries.mkdir()
    for command in ("make", "mkdir"):
        executable = shutil.which(command)
        assert executable is not None, f"{command} is required for the C build tests"
        (binaries / command).symlink_to(executable)
    monkeypatch.setenv("PATH", str(binaries))

    flags = tmp_path / "obj" / "flags"
    if existing_flags:
        flags.parent.mkdir()
        flags.write_text("previous flags\n")

    with pytest.raises(RuntimeError, match="cmp not found; install diffutils"):
        run_make([f"OBJ_DIR={flags.parent}", "ARCH_FLAGS=", str(flags)])

    if existing_flags:
        assert flags.read_text() == "previous flags\n"
    else:
        assert not flags.exists()


def test_every_c_array_parameter_has_an_adjacent_checked_length() -> None:
    """Include private helpers, not just the exported API."""
    functions = re.compile(r"\b(\w+)\s*\(([^{};]*)\)\s*\{")
    arrays = re.compile(r"\b(\w+)\s*\[\s*\]")
    paths = [
        *c_directory().glob("src/*.[ch]"),
        *c_directory().glob("tests/*.[ch]"),
    ]
    for path in sorted(paths):
        text = path.read_text()
        for function in functions.finditer(text):
            depth = 1
            end = function.end()
            while depth and end < len(text):
                depth += (text[end] == "{") - (text[end] == "}")
                end += 1
            body = text[function.end() : end]
            parameters = function[2]
            for array in arrays.finditer(parameters):
                location = f"{path.name}:{function[1]}:{array[1]}"
                length = re.match(r"\s*,\s*(?:const\s+)?int32_t\s+(\w+)", parameters[array.end() :])
                assert length, f"{location}: array must be followed immediately by its length"
                assert re.search(rf"RTGSFIT_CHECK_(?:SIZE|MIN_SIZE|ARRAY)\([^;]*\b{length[1]}\b", body), (
                    f"{location}: missing length validation"
                )
