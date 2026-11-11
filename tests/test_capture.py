"""The same native NetCDF contract for replay and a future PCS logger."""

from __future__ import annotations

import copy
from dataclasses import fields
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit import capture as capture_module
from rtgsfit import data_dictionary as dd
from rtgsfit.__main__ import main
from rtgsfit.capture import OUTPUTS
from rtgsfit.capture import SCHEMA
from rtgsfit.capture import Capture
from rtgsfit.capture import read_capture
from rtgsfit.capture import write_capture
from rtgsfit.database_readers.local import read_results
from rtgsfit.library import BuildOptions
from rtgsfit.library import Rtgsfit
from rtgsfit.library import SliceResult
from rtgsfit.preshot import Preshot
from rtgsfit.replay import read_signals
from rtgsfit.replay import replay
from rtgsfit.settings import load_settings
from rtgsfit.settings import parse_settings


@pytest.fixture(scope="module")
def native_capture(analytic_preshot: Preshot) -> Capture:
    settings = load_settings("analytic").copy_with({"replay/time_end": 0.002, "replay/repeat": 1})
    signals = read_signals(settings, analytic_preshot)
    signals.meas[-1] = 0.0
    signals.coil_signal[-1] = 0.0
    solver = Rtgsfit.build(analytic_preshot, BuildOptions(timing=True))
    coils = signals.coil_current(analytic_preshot)
    slices = [solver(signals.meas[i], coils[i]) for i in range(len(signals.time))]
    assert slices[0].status == 0
    assert slices[-1].status != 0
    return Capture(analytic_preshot, signals, slices, settings.json_text, "PCS logger", "test-logger-version")


def assert_same_results(actual: dict, expected: dict) -> None:
    assert actual.keys() == expected.keys()
    for key in expected:
        assert_array_equal(actual[key], expected[key], err_msg=key)


def test_contract_covers_every_c_output() -> None:
    assert set(OUTPUTS) | {"timing"} == {field.name for field in fields(SliceResult)}


@pytest.mark.parametrize("timed", [False, True])
def test_native_round_trip(native_capture: Capture, tmp_path: Path, timed: bool) -> None:
    captured = copy.deepcopy(native_capture)
    assert captured.slices[0].timing is not None
    if not timed:
        captured.slices = [replace(item, timing=None) for item in captured.slices]
    path = write_capture(tmp_path / "capture.nc", captured)
    restored = read_capture(path)
    assert restored.code_version == "test-logger-version"
    assert restored.settings_json == captured.settings_json
    assert restored.preshot_source == captured.preshot_source
    assert restored.preshot.content_hash() == captured.preshot.content_hash()
    for node in dd.PRESHOT_NODES:
        key = ":".join(node.path[1:])
        assert_array_equal(restored.preshot[key], captured.preshot[key], err_msg=key)
    for name in ("time", "meas", "coil_signal"):
        assert_array_equal(getattr(restored.signals, name), getattr(captured.signals, name))
    for actual, expected in zip(restored.slices, captured.slices, strict=True):
        for field in fields(SliceResult):
            assert_array_equal(getattr(actual, field.name), getattr(expected, field.name), err_msg=field.name)
    assert_same_results(restored.to_results().flat, captured.to_results().flat)
    with Dataset(path) as file:
        assert file.schema == SCHEMA
        assert file.flux_convention == "native"
        assert file.data_model == "NETCDF4"
        assert int(file.solver_checksum) == capture_module.solver_preshot_checksum(captured.preshot)
        assert int(file.preshot_checksum) == capture_module._full_preshot_checksum(captured.preshot)
        assert file["OUTPUT/psi"].dimensions == ("N_TIME", "N_Z", "N_R")
        assert file["OUTPUT/mask"].dtype == np.dtype("int32")
        assert "TARGET" not in file["OUTPUT"].variables
    assert list(tmp_path.iterdir()) == [path]


def test_processing_native_units_once(native_capture: Capture, tmp_path: Path) -> None:
    captured = copy.deepcopy(native_capture)
    captured.preshot.groups["INFO"]["PSI_TO_WB"] = np.float64(2.0 * np.pi)
    native_psi = captured.slices[0].psi.copy()
    restored = read_capture(write_capture(tmp_path / "native.nc", captured))
    first = restored.to_results()
    second = restored.to_results()
    assert_allclose(first["PROFILES_2D.R_Z.PSI"][0].ravel(), native_psi * 2.0 * np.pi)
    assert_array_equal(restored.slices[0].psi, native_psi)
    assert_same_results(first.flat, second.flat)
    assert np.isnan(first["PROFILES_2D.R_Z.B_FIELD_R"][-1]).all()


def test_replay_capture_no_passives_and_early_stop(tmp_path: Path) -> None:
    path = tmp_path / "no_passives.nc"
    settings = load_settings("analytic").copy_with(
        {
            "analytic/passives/r": [],
            "analytic/passives/z": [],
            "build/timing": False,
            "replay/stop_on_status": True,
            "analytic/plasma/ip": 0.0,
            "capture/path": str(path),
        }
    )
    result = replay(settings, verbose=False)
    captured = read_capture(path)
    assert len(captured.slices) == len(captured.signals.time) == 1
    assert captured.slices[0].status != 0
    assert captured.slices[0].timing is None
    assert captured.preshot.sizes["N_PASSIVE"] == captured.preshot.sizes["N_REG"] == 0
    assert captured.preshot["GRID:G_PASSIVE"].shape == (captured.preshot.sizes["N_GRID"], 0)
    assert_same_results(captured.to_results().flat, result.flat)


def test_offline_cli_does_not_run_solver_or_read_databases(
    native_capture: Capture,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = write_capture(tmp_path / "pcs.nc", native_capture)

    def forbidden(*args: object, **kwargs: object) -> None:
        pytest.fail("post-processing must not run the solver, build, or read a database")

    monkeypatch.setattr(Rtgsfit, "__call__", forbidden)
    monkeypatch.setattr(Rtgsfit, "build", forbidden)
    monkeypatch.setattr("rtgsfit.replay.read_preshot", forbidden)
    monkeypatch.setattr("rtgsfit.replay.read_signals", forbidden)
    main(["post-process", str(path)])
    output = read_results(tmp_path / "pcs_processed.npz")
    assert_same_results(output.flat, native_capture.to_results().flat)
    assert output["INPUT.CODE_VERSION"] == "test-logger-version"
    custom = tmp_path / "custom.npz"
    main(["post-process", str(path), "--set", f"writers/local/path={custom}"])
    assert custom.is_file()
    with pytest.raises(ValueError, match="writers/"):
        main(["post-process", str(path), "--set", "replay/pulse=42"])


def test_offline_cli_does_not_enable_archived_database_writers(
    native_capture: Capture,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    settings = parse_settings(native_capture.settings_json).copy_with(
        {
            "writers/names": ["rtgsfit_mdsplus"],
            "writers/rtgsfit_mdsplus/enabled": True,
        }
    )
    path = write_capture(tmp_path / "pcs.nc", replace(native_capture, settings_json=settings.json_text))

    def forbidden(*args: object, **kwargs: object) -> None:
        pytest.fail("archived database writers must not be run implicitly")

    monkeypatch.setattr("rtgsfit.database_writers.rtgsfit_mdsplus.RtgsfitMdsplusWriter.write_results", forbidden)
    main(["post-process", str(path)])
    assert (tmp_path / "pcs_processed.npz").is_file()


def test_capture_written_before_postprocessing(
    native_capture: Capture,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "before_processing.nc"
    settings = parse_settings(native_capture.settings_json).copy_with({"capture/path": str(path)})

    def fail_processing(*args: object, **kwargs: object) -> None:
        raise RuntimeError("offline processing failed")

    monkeypatch.setattr(Capture, "to_results", fail_processing)
    with pytest.raises(RuntimeError, match="offline processing failed"):
        replay(settings, preshot=native_capture.preshot, signals=native_capture.signals, verbose=False)
    captured = read_capture(path)
    assert len(captured.slices) == len(native_capture.slices)


@pytest.mark.parametrize(
    "change",
    [
        "schema",
        "flux_convention",
        "solver_checksum",
        "preshot_checksum",
        "preshot_value",
        "coil_current",
        "extra_preshot",
    ],
)
def test_bad_captures_are_rejected(native_capture: Capture, tmp_path: Path, change: str) -> None:
    path = write_capture(tmp_path / "invalid.nc", native_capture)
    with Dataset(path, "r+") as file:
        if change in ("schema", "flux_convention"):
            file.setncattr(change, "wrong")
        elif change in ("solver_checksum", "preshot_checksum"):
            file.setncattr(change, 0)
        elif change == "preshot_value":
            file["PRESHOT/INFO/SOURCE"][...] = "modified"
        elif change == "coil_current":
            file["INPUT/coil_current"][0, 0] += 100.0
        else:
            file["PRESHOT/GRID"].createVariable("unexpected", "float64")[...] = 1.0
    with pytest.raises(ValueError):
        read_capture(path)


def test_failed_write_is_atomic(native_capture: Capture, tmp_path: Path) -> None:
    target = tmp_path / "existing-directory"
    target.mkdir()
    with pytest.raises(OSError, match="cannot replace"):
        write_capture(target, native_capture)
    assert list(tmp_path.iterdir()) == [target]


def test_failed_nonfinite_values_are_retained(native_capture: Capture, tmp_path: Path) -> None:
    captured = copy.deepcopy(native_capture)
    captured.signals.meas[-1, 0] = np.nan
    captured.slices[-1].j_phi[:] = np.nan
    captured.slices[-1].status = 1024
    restored = read_capture(write_capture(tmp_path / "nonfinite.nc", captured))
    assert np.isnan(restored.signals.meas[-1, 0])
    assert np.isnan(restored.slices[-1].j_phi).all()
    assert restored.slices[-1].status == 1024
    results = restored.to_results()
    assert np.isnan(results["PROFILES_2D.R_Z.B_FIELD_R"][-1]).all()


def test_legacy_pcs_file_is_not_silently_reinterpreted(tmp_path: Path) -> None:
    path = tmp_path / "legacy.nc"
    with Dataset(path, "w", format="NETCDF3_CLASSIC") as file:
        file.createDimension("n_time", 1)
        file.createVariable("time", "float64", ("n_time",))[:] = [0.0]
    with pytest.raises(ValueError, match="legacy PCS"):
        read_capture(path)


def test_invalid_capture_does_not_overwrite_existing_file(native_capture: Capture, tmp_path: Path) -> None:
    path = write_capture(tmp_path / "capture.nc", native_capture)
    original = path.read_bytes()
    invalid = replace(native_capture, slices=native_capture.slices[:-1])
    with pytest.raises(ValueError, match="one output slice"):
        write_capture(path, invalid)
    assert path.read_bytes() == original


def test_relative_capture_path_uses_cache(
    tmp_path: Path,
    analytic_preshot: Preshot,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("RTGSFIT_CACHE_DIR", str(tmp_path))
    settings = load_settings("analytic").copy_with({"capture/path": "native.nc", "replay/time_end": 0.0})
    replay(settings, preshot=analytic_preshot, verbose=False)
    assert (tmp_path / "native.nc").is_file()
