"""The PCS-like ST40 mode, `st40_tepcs`: the PCS's own inputs, from the raw TEPCS counts.

The `st40_tepcs` reader calibrates the TEPCS digitiser counts as the PCS does, with the pulse's PSS
calibration and the sensor-replacement matrix. The PSS and TEPCS tests at the top use small synthetic
records; the rest need the TE network (--te-network).
"""

from __future__ import annotations

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
    def __init__(
        self,
        *,
        missing_mag: bool = False,
        partial: bool = False,
        short: bool = False,
        per_cycle: int = 1,
        comp_clock: float = 1.0,
    ) -> None:
        self.missing_mag = missing_mag
        self.partial = partial
        self.per_cycle = per_cycle
        self.comp_clock = comp_clock  # stored-time rate of the compensation channel, CH003
        self.samples = np.arange(23000 * per_cycle, dtype=float)
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
            clock = self.comp_clock if "CH003" in expression else 1.0
            value = -1.1 + 5e-5 / self.per_cycle * clock * self.samples
        else:
            value = values[expression]
        return SimpleNamespace(data=lambda: value)


def _settings(pulse: int) -> Settings:
    return load_settings("st40_tepcs").copy_with(
        {"replay/pulse": pulse, "replay/time_start": 0.0055, "replay/time_end": 0.0055}
    )


@pytest.mark.parametrize(("per_cycle", "comp_clock"), [(1, 1.0), (5, 1 - 1e-4)])
def test_pcs_calibration_and_drift_compensation(per_cycle: int, comp_clock: float) -> None:
    """Each cycle reads the last sample before it by count, even where the stored times run fast."""
    signals = pcs_signals(Connection(per_cycle=per_cycle, comp_clock=comp_clock), _settings(14928), PRESHOT)
    cycle = 22110
    sample = per_cycle * cycle - 1
    avg1 = 2 * per_cycle * np.arange(1, 601).mean()
    avg2 = 2 * per_cycle * np.arange(9401, 10001).mean()
    expected = 2 * (sample + 1) - (avg2 + (avg2 - avg1) / 0.47 * 5e-5 * (cycle - 10000))
    expected += 0.5 * 4 * (sample + 3)
    assert_allclose(signals.meas[:, 0], [expected])
    assert_allclose(signals.coil_signal[:, 0], [3 * (sample + 2)])
    assert "compensation" in signals.input_quality["correction"]


def test_missing_compensation_is_flagged_but_partial_metadata_is_rejected() -> None:
    settings = _settings(11089)
    signals = pcs_signals(Connection(missing_mag=True), settings, PRESHOT)
    assert_allclose(signals.meas[:, 0], [2 * 22110])
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
    from stored_preshot import assert_mdsplus_preshot_is_the_gsfit_build

    assert_mdsplus_preshot_is_the_gsfit_build(load_settings("st40_tepcs"))


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
        "preshot/rtgsfit_mdsplus/run_name": "PRESHOT08" if pulse < 12551 else "PRESHOT07",
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
