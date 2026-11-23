"""The PCS-like ST40 mode, `st40_tepcs`: the PCS's own inputs, from the raw TEPCS counts.

The `st40_tepcs` reader calibrates the TEPCS digitiser counts as the PCS does, with the pulse's PSS
calibration and the sensor-replacement matrix. The PSS and TEPCS tests at the top use small synthetic
records; the rest need the TE network (--te-network).
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mdsthin.exceptions import TreeNNF
from mdsthin.exceptions import TreeNODATA
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.database_readers import InputDataError
from rtgsfit.database_readers.st40_tepcs import pcs_signals
from rtgsfit.replay import build_preshot
from rtgsfit.replay import read_preshot
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import Settings
from rtgsfit.settings import load_settings

PRESHOT = {"MEAS:NAMES": ["MEAS"], "COIL:SIGNALS": ["COIL"]}

CASES = {
    (11089, 0.030): 0,
    (11089, 0.100): 0,
    (14211, 0.030): 0,
    (14211, 0.100): 0,
    (14816, 0.030): 0,
    (14816, 0.100): 1,
}
"""(pulse, time [s]) -> expected topology: 0 wall-limited, 1 diverted."""


class Connection:
    def __init__(self, *, missing_mag: bool = False, partial: bool = False, short: bool = False) -> None:
        self.missing_mag = missing_mag
        self.partial = partial
        self.samples = np.arange(23000, dtype=float)
        self.short = short

    def openTree(self, tree: str, pulse: int) -> None:
        pass

    def closeTree(self, tree: str, pulse: int) -> None:
        pass

    def get(self, expression: str) -> SimpleNamespace:
        values = {
            "\\PSS::TOP.AI.LIST": np.array(["MEAS", "COIL", "COMP"]),
            "\\PSS::TOP.AI.OFFSETCOUNTS": np.array([1, 2, 3]),
            "\\PSS::TOP.AI.COUNT2MKS": np.array([2, 3, 4]),
            "\\PSS::TOP.TIMES.T_OFFSET": -1.1,
            "\\PSS::TOP.TIMES.DT_CYCLE": 5e-5,
            "\\PSS::TOP.MAG.SUB_NAMES": np.array(["MEAS"]),
            "\\PSS::TOP.MAG.COMP_MATRIX": np.array([[0.5]]),
            "\\PSS::TOP.MAG.INDEX_COMP": np.array([3]),
        }
        if expression.startswith("\\PSS::TOP.MAG.") and (
            self.missing_mag or (self.partial and expression.endswith("INDEX_COMP"))
        ):
            raise TreeNNF()
        if expression.startswith("raw_of("):
            value = self.samples[:-1] if self.short else self.samples
        elif expression.startswith("dim_of("):
            value = -1.1 + 5e-5 * self.samples
        else:
            value = values[expression]
        return SimpleNamespace(data=lambda: value)


def _settings(pulse: int) -> Settings:
    return load_settings("st40_tepcs").copy_with(
        {"replay/pulse": pulse, "replay/time_start": 0.0055, "replay/time_end": 0.0055}
    )


def test_pcs_calibration_and_drift_compensation() -> None:
    signals = pcs_signals(Connection(), _settings(14928), PRESHOT)
    cycle = 22110
    avg1 = 2 * (np.arange(1, 601) + 1).mean()
    avg2 = 2 * (np.arange(9401, 10001) + 1).mean()
    expected = 2 * (cycle + 1) - (avg2 + (avg2 - avg1) / 0.47 * 5e-5 * (cycle - 10000))
    expected += 0.5 * 4 * (cycle + 3)
    assert_allclose(signals.meas[:, 0], [expected])
    assert "compensation" in signals.input_quality["correction"]


def test_missing_compensation_is_flagged_but_partial_metadata_is_rejected() -> None:
    settings = _settings(11089)
    signals = pcs_signals(Connection(missing_mag=True), settings, PRESHOT)
    assert_allclose(signals.meas[:, 0], [2 * (22110 + 1)])
    assert "unknown" in signals.input_quality["correction"]
    with pytest.raises(InputDataError, match="incomplete"):
        pcs_signals(Connection(partial=True), settings, PRESHOT)
    with pytest.raises(InputDataError, match="inconsistent counts and time"):
        pcs_signals(Connection(short=True), settings, PRESHOT)


def test_missing_tepcs_record_identifies_signal_and_channel() -> None:
    class MissingChannel(Connection):
        def get(self, expression: str) -> SimpleNamespace:
            if expression.startswith("raw_of("):
                raise TreeNODATA()
            return super().get(expression)

    with pytest.raises(InputDataError, match=r"14894.*MEAS.*ACQ2106_032:CH001.*TreeNODATA"):
        pcs_signals(MissingChannel(), _settings(14894), PRESHOT)


@pytest.mark.te_network
def test_mdsplus_preshot_is_the_gsfit_build() -> None:
    """The PRESHOT in MDSplus (`preshot/reader`) equals a fresh build from GSFit (`preshot/builder`).

    Floating-point nodes are compared to round-off rather than by `content_hash`: the Poisson LU
    bands come from threaded LAPACK, whose rounding depends on the BLAS thread count. The archived
    `INFO:SETTINGS` must agree in the sections that make the preshot; replay and writer settings may
    change without publishing a new preshot.
    """
    settings = load_settings("st40_tepcs")
    stored, built = read_preshot(settings), build_preshot(settings)
    assert built.groups.keys() == stored.groups.keys()
    for group, nodes in stored.groups.items():
        assert built.groups[group].keys() == nodes.keys(), group
        for name, value in nodes.items():
            if f"{group}:{name}" == "INFO:SETTINGS":
                archived, current = json.loads(str(value)), json.loads(str(built.groups[group][name]))
                for section in ("preshot", "gsfit"):
                    assert archived[section] == current[section], f"INFO:SETTINGS {section}"
                continue
            expected, actual = np.asarray(value), np.asarray(built.groups[group][name])
            if expected.dtype.kind == "f":
                scale = float(np.abs(expected).max(initial=0.0))
                assert_allclose(actual, expected, rtol=1e-12, atol=1e-12 * scale, err_msg=f"{group}:{name}")
            else:
                assert_array_equal(actual, expected, err_msg=f"{group}:{name}")


@pytest.mark.te_network
def test_default_replay() -> None:
    results = replay(load_settings("st40_tepcs").copy_with({"replay/time_end": 0.1}), verbose=False)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    flat_top = results["TIME"] > 0.03
    assert_array_equal(results["GLOBAL.STATUS"][flat_top], 0)
    assert np.all(np.abs(results["GLOBAL.IP"][flat_top]) > 1e5)
    assert_allclose(
        results["FIT.CHI_SQ"][flat_top].sum(axis=1),
        results["CONSTRAINTS.CHI_SQ_TOTAL"][flat_top],
        rtol=1e-8,
        atol=1e-8,
    )
    assert np.any(results["CONSTRAINTS.CHI_SQ_REG"][flat_top] > 0)


@pytest.fixture(scope="module", params=list(CASES), ids=[f"{pulse}@{time}" for pulse, time in CASES])
def case(request: pytest.FixtureRequest) -> tuple[Settings, Results]:
    pulse, time = request.param
    settings = load_settings("st40_tepcs").copy_with({
        "replay/pulse": pulse,
        "replay/time_start": time,
        "replay/time_end": time,
        "replay/repeat": 10,
        "preshot/rtgsfit_mdsplus/run_name": "PRESHOT05" if pulse < 12551 else "PRESHOT04",
    })  # fmt: skip
    return settings, replay(settings, verbose=False)


@pytest.mark.te_network
def test_fixed_time_replay(case: tuple[Settings, Results]) -> None:
    """A physical plasma of the expected topology, in both coil eras, without a GSFit input-equivalence claim."""
    settings, results = case
    assert_array_equal(results["GLOBAL.STATUS"], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    assert abs(results["GLOBAL.IP"][-1]) > 1e5
    assert results["BOUNDARY.OUTLINE.N"][-1] > 3
    expected = CASES[settings.get_path("replay/pulse"), settings.get_path("replay/time_start")]
    assert results["GLOBAL.XPT_DIVERTED"][-1] == expected


@pytest.mark.te_network
@pytest.mark.figures
def test_figures(case: tuple[Settings, Results], figure_dir: Path) -> None:
    """TEPCS reconstructions and residuals against their own fixed-time measurements, not GSFit's inputs."""
    import plotting

    settings, results = case
    figures = {
        "psi": plotting.plot_reconstruction(results, "PSI"),
        "j_phi": plotting.plot_reconstruction(results, "J_PHI"),
        "constraints": plotting.plot_constraints(results),
    }
    pulse, time = settings.get_path("replay/pulse"), settings.get_path("replay/time_start")
    paths = plotting.save(figures, figure_dir / "st40_tepcs" / f"{pulse}@{time}")
    assert all(path.stat().st_size > 0 for path in paths)
