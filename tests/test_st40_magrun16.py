"""The default ST40 mode, `st40_magrun16`: RT-GSFit fed the inputs GSFit's post-shot reconstructions fit.

GSFIT RUN04 fits the processed MAG RUN16 signals and the PSU2COIL coil currents. The `st40_magrun16`
reader gives RT-GSFit exactly those inputs, and its preshot is built by GSFit from the packaged
settings `st40_magrun16/gsfit`. These are RUN04's settings with one change: p' and FF' have one degree
of freedom each, as in RT-GSFit. GSFit is run here with those settings on their 81x161 grid; RT-GSFit
uses its own 33x65 grid. So the two codes differ only by their grid and algorithm.

The reader tests at the top use a fake MDSplus connection; the rest need the TE network (--te-network).
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from mdsthin.exceptions import TreeNNF
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.database_readers import InputDataError
from rtgsfit.database_readers.st40_magrun16 import St40Magrun16Reader
from rtgsfit.database_readers.st40_magrun16 import mag_node
from rtgsfit.database_readers.st40_magrun16 import magrun16_signals
from rtgsfit.database_readers.st40_magrun16 import psu2coil_channels
from rtgsfit.replay import read_preshot
from rtgsfit.replay import read_signals
from rtgsfit.replay import replay
from rtgsfit.settings import Settings
from rtgsfit.settings import load_settings

# The reader, with a fake MDSplus connection.

GROUPS = {"SOL": ["SOL"], "MCT": ["MCT"], "MCB": ["MCB"], "DIV": ["DIVT", "DIVB"]}
FAKE_PRESHOT = {
    "CONSTRAINT:MAG_SIGNALS": ["PSI_FLOOP_007", "I_ROG_INIVC000"],
    "COIL:COIL_SIGNALS": ["I_SOL_PSU2COIL", "I_MCT_PSU2COIL", "I_MCB_PSU2COIL", "I_DIV_PSU2COIL"],
}
OLD_ERA_CHANNELS = np.array([["DIVB", "DIVT"], ["MCB", "MCT"], ["SOL", ""]])
"""PSU2COIL `PF.ALL:COILS` before pulse 12551: MCB and MCT are in series, in one channel."""


class Connection:
    """MAG RUN16 is recorded at 0, 10 and 20 ms with values 1, 2 and 3; PSU2COIL at `coil_time`, with
    the current of channel `k` equal to `(k + 1) * 100` times the sample number plus one."""

    def __init__(self, coil_time: tuple[float, ...] = (0.0, 0.01, 0.02), missing: str = "") -> None:
        self.coil_time = np.array(coil_time)
        self.missing = missing
        self.paths: list[str] = []

    def openTree(self, tree: str, pulse: int) -> None:
        pass

    def closeTree(self, tree: str, pulse: int) -> None:
        pass

    def get(self, expression: str) -> SimpleNamespace:
        self.paths.append(expression)
        if self.missing and self.missing in expression:
            raise TreeNNF()
        if "\\MAG::" in expression:
            value = np.array([0.0, 0.01, 0.02] if expression.startswith("dim_of(") else [1.0, 2.0, 3.0])
        elif expression.endswith("PF.ALL:COILS"):
            value = OLD_ERA_CHANNELS
        elif expression.endswith(":TIME"):
            value = self.coil_time
        elif expression.endswith("PF.ALL:I"):
            value = np.outer(np.arange(1, len(self.coil_time) + 1), 100.0 * np.arange(1, 4))
        else:
            assert expression == "\\PSU2COIL::TOP.BEST:BEST_RUN"
            value = np.array("RUN02")
        return SimpleNamespace(data=lambda: value)


def _fake_settings(time_end: float = 0.02) -> Settings:
    return load_settings("st40_magrun16").copy_with({
        "replay/time_start": 0.0,
        "replay/time_end": time_end,
        "replay/time_step": 0.005,
        "gsfit/coils/groups": GROUPS,
    })  # fmt: skip


def test_mag_node_is_the_processed_signal_gsfit_fits() -> None:
    assert mag_node("PSI_FLOOP_007") == "FLOOP.L007:PSI_PROCESS"
    assert mag_node("B_BPPROBE_101") == "BPPROBE.P101:B_PROCESS"
    assert mag_node("I_ROG_INIVC000") == "ROG.INIVC000:I_PROCESS"
    with pytest.raises(ValueError, match="I_SOL_PSU"):
        mag_node("I_SOL_PSU")


def test_psu2coil_channels() -> None:
    signals = FAKE_PRESHOT["COIL:COIL_SIGNALS"]
    assert psu2coil_channels(signals, GROUPS, OLD_ERA_CHANNELS) == [2, 1, 1, 0]
    with pytest.raises(InputDataError, match=r"does not have: \['PSH'\]"):
        psu2coil_channels(signals, GROUPS, np.vstack([OLD_ERA_CHANNELS, [["PSH", ""]]]))
    with pytest.raises(InputDataError, match=r"no current for the PF coils \['DIVB'\] of DIV"):
        psu2coil_channels(signals, GROUPS, np.array([["DIVT", ""], ["MCB", "MCT"], ["SOL", ""]]))
    with pytest.raises(InputDataError, match="DIV are in different PSU2COIL channels"):
        psu2coil_channels(signals, GROUPS, np.array([["DIVT", "MCB"], ["DIVB", "MCT"], ["SOL", ""]]))


def test_reader_interpolates_and_needs_only_mdsthin(monkeypatch: pytest.MonkeyPatch) -> None:
    import mdsthin

    connection = Connection()
    monkeypatch.setitem(sys.modules, "MDSplus", None)
    monkeypatch.setitem(sys.modules, "st40_database", None)
    monkeypatch.setattr(mdsthin, "Connection", _context(connection))
    signals = St40Magrun16Reader().read_signals(_fake_settings(), FAKE_PRESHOT)  # type: ignore[arg-type]
    assert_allclose(signals.time, [0.0, 0.005, 0.01, 0.015, 0.02])
    assert_allclose(signals.mag_signal, np.repeat([[1.0], [1.5], [2.0], [2.5], [3.0]], 2, axis=1))
    sample = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    assert_allclose(signals.coil_signal, np.column_stack([300 * sample, 200 * sample, 200 * sample, 100 * sample]))
    assert signals.input_quality == {"source": "MAG RUN16 processed signals and PSU2COIL BEST (RUN02)"}
    assert "data(\\MAG::TOP.RUN16.FLOOP.L007:PSI_PROCESS)" in connection.paths
    assert "\\PSU2COIL::TOP.BEST.PF.ALL:I" in connection.paths


def _context(connection: Connection) -> Any:
    class Context:
        def __init__(self, server: str) -> None:
            assert server == "smaug"

        def __enter__(self) -> Connection:
            return connection

        def __exit__(self, *args: object) -> None:
            pass

    return Context


def test_reader_drops_times_a_signal_did_not_record() -> None:
    signals = magrun16_signals(Connection(coil_time=(0.0, 0.01)), _fake_settings(), FAKE_PRESHOT)  # type: ignore[arg-type]
    assert_allclose(signals.time, [0.0, 0.005, 0.01])
    with pytest.raises(InputDataError, match="no replay times recorded"):
        magrun16_signals(Connection(coil_time=(0.05, 0.06)), _fake_settings(), FAKE_PRESHOT)  # type: ignore[arg-type]


def test_missing_node_is_an_input_error() -> None:
    connection = Connection(missing="ROG.INIVC000:I_PROCESS")
    with pytest.raises(InputDataError, match=r"pulse 14816 .*ROG\.INIVC000:I_PROCESS.*TreeNNF"):
        magrun16_signals(connection, _fake_settings(), FAKE_PRESHOT)  # type: ignore[arg-type]


# RT-GSFit against GSFit, on the TE network.

CASES = {
    (11089, 0.030): 0,
    (11089, 0.100): 0,
    (14211, 0.030): 0,
    (14211, 0.100): 0,
    (14816, 0.030): 0,
    (14816, 0.100): 1,
}
"""(pulse, time [s]) -> expected topology: 0 wall-limited, 1 diverted."""

# Tolerances of RT-GSFit's last iteration against GSFit, 2-3 times the largest difference of the six cases.
# Measured (largest of the six): Ip 5.0e-5 relative; the magnetic axis 7.6e-4 m in R and 4.7e-4 m in Z; psi at the
# axis 2.3e-3 Wb and at the boundary 1.2e-3 Wb; flux loops 9.2e-4 Wb, Bp probes 3.0e-3 T, Rogowski coils 1.3 A;
# IVC eigenmode currents 0.95 A.
RTOL = {"GLOBAL.IP": 1.5e-4}
"""Relative tolerances against GSFit."""
ATOL = {"GLOBAL.MAG_AXIS.R": 2e-3, "GLOBAL.MAG_AXIS.Z": 1.2e-3, "GLOBAL.PSI_MAG_AXIS": 6e-3, "BOUNDARY.PSI": 3e-3}
"""Absolute tolerances against GSFit [m, Wb]; psi at the boundary can be near 0, so not relative."""
# (group, GSFit sensor quantity, atol vs the measurements, atol vs GSFit's reconstruction) [Wb, T, A]. The atol vs
# the measurements bounds RT-GSFit's fit, whose largest residuals are 1.7e-2 Wb, 9.9e-2 T and 351 A.
SENSORS = {
    "flux_loop": ("FLUX_LOOP", "psi", 4.5e-2, 2.5e-3),
    "bp_probe": ("BP_PROBE", "b", 0.2, 6e-3),
    "rogowski": ("ROGOWSKI", "i", 750.0, 3.0),
}
IVC_ATOL = 2.5
"""[A] for each IVC eigenmode current."""


def _gsfit(settings: Settings, solve: bool = True) -> Any:
    """GSFit with the `gsfit/settings_path` settings, on their own grid, at the replay times."""
    from gsfit import Gsfit
    from rtgsfit.database_readers.gsfit import gsfit_settings_path
    from rtgsfit.settings import replay_times

    controller = Gsfit(
        pulseNo=int(settings["replay"]["pulse"]),
        run_name="TEST",
        settings_path=gsfit_settings_path(settings),
        write_to_mds=False,
    )
    code_settings = controller.settings["GSFIT_code_settings.json"]
    code_settings["timeslices"]["method"] = "user_defined"
    code_settings["timeslices"]["user_defined"] = [float(time) for time in replay_times(settings)]
    controller.set_environment_variables()
    controller.setup_timeslices()
    controller.setup_objects()
    if solve:
        controller.calculate_greens()
        controller.inverse_solver_rust()
    return controller


@pytest.mark.te_network
@pytest.mark.parametrize("pulse", [11089, 14816])
def test_reader_gives_rt_gsfit_the_inputs_gsfit_fits(pulse: int) -> None:
    """Exactly: 11089 has MCB and MCT in one PSU2COIL channel, 14816 in two."""
    from rtgsfit.database_readers.gsfit import coil_circuits
    from rtgsfit.database_readers.gsfit import gsfit_coil_current
    from rtgsfit.database_readers.gsfit import gsfit_signals

    settings = load_settings("st40_magrun16").copy_with({
        "replay/pulse": pulse, "replay/time_start": 0.03, "replay/time_end": 0.1, "replay/time_step": 0.01,
    })  # fmt: skip
    preshot = read_preshot(settings)
    controller = _gsfit(settings, solve=False)
    signals = St40Magrun16Reader().read_signals(settings, preshot)
    assert_array_equal(signals.mag_signal, gsfit_signals(controller, settings, preshot).mag_signal)
    circuits = coil_circuits(controller, settings)
    assert_array_equal(signals.coil_signal, gsfit_coil_current(controller, circuits, signals.time))


@pytest.mark.te_network
def test_mdsplus_preshot_is_the_gsfit_build() -> None:
    """`PRESHOT09` is GSFit's build, which is the same for every pulse in both coil eras."""
    from stored_preshot import assert_mdsplus_preshot_is_the_gsfit_build

    assert_mdsplus_preshot_is_the_gsfit_build(load_settings("st40_magrun16"))


@pytest.fixture(scope="module", params=list(CASES), ids=[f"{pulse}@{time}" for pulse, time in CASES])
def case(request: pytest.FixtureRequest):  # type: ignore[no-untyped-def]
    pulse, time = request.param
    settings = load_settings("st40_magrun16").copy_with({
        "replay/pulse": pulse,
        "replay/time_start": time,
        "replay/time_end": time,
        "replay/repeat": 10,
    })  # fmt: skip
    preshot = read_preshot(settings)
    results = replay(settings, preshot=preshot, signals=read_signals(settings, preshot), verbose=False)
    return settings, results, _gsfit(settings)


def _gsfit_equilibrium(controller) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    """GSFit's reconstruction of the time slice, by RT-GSFit's results node."""
    from gsfit_rs.imas import equilibrium_paths

    time_slice = equilibrium_paths.time_slice[0]
    profiles_2d, quantities = time_slice.profiles_2d[0], time_slice.global_quantities

    def array(path) -> np.ndarray:  # type: ignore[no-untyped-def]
        return np.asarray(controller.plasma.get(path))

    def scalar(path) -> float:  # type: ignore[no-untyped-def]
        return float(np.ravel(array(path))[0])

    return {
        "PROFILES_2D.R_Z.R": array(profiles_2d.grid.dim1),
        "PROFILES_2D.R_Z.Z": array(profiles_2d.grid.dim2),
        "PROFILES_2D.R_Z.PSI": array(profiles_2d.psi),
        "PROFILES_2D.R_Z.J_PHI": array(profiles_2d.j_phi),
        "GLOBAL.IP": scalar(quantities.ip),
        "GLOBAL.MAG_AXIS.R": scalar(quantities.magnetic_axis.r),
        "GLOBAL.MAG_AXIS.Z": scalar(quantities.magnetic_axis.z),
        "GLOBAL.PSI_MAG_AXIS": scalar(quantities.psi_magnetic_axis),
        "GLOBAL.XPT_DIVERTED": scalar(time_slice.boundary.type),
        "BOUNDARY.PSI": scalar(time_slice.boundary.psi),
        "BOUNDARY.OUTLINE.R": array(time_slice.boundary.outline.r),
        "BOUNDARY.OUTLINE.Z": array(time_slice.boundary.outline.z),
    }


def _gsfit_sensors(settings, results, controller, kind: str) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[no-untyped-def]
    """GSFit's measured and reconstructed values of the `kind` sensors, in RT-GSFit's order of `CONSTRAINTS`."""
    from rtgsfit.database_readers.gsfit import constraint_sensors
    from rtgsfit.database_readers.gsfit import sensor_name

    group, quantity, _, _ = SENSORS[kind]
    sensors = {
        "flux_loop": controller.flux_loops,
        "bp_probe": controller.bp_probes,
        "rogowski": controller.rogowski_coils,
    }[kind]
    names = list(results[f"CONSTRAINTS.{group}.ALL.NAMES"])
    measured, gsfit = np.zeros(len(names)), np.zeros(len(names))
    for sensor in (sensor for sensor in constraint_sensors(controller, settings) if sensor.kind == kind):
        index = names.index(sensor_name(sensor, "st40_pcs"))
        measured[index] = sensors.get_array1([sensor.name, quantity, "measured", "value"])[0]
        gsfit[index] = sensors.get_array1([sensor.name, quantity, "calculated", "value"])[0]
    return measured, gsfit


def _ivc_eigenmodes(results, controller) -> tuple[list[str], np.ndarray]:  # type: ignore[no-untyped-def]
    """The names of RT-GSFit's IVC eigenmode DOFs, and GSFit's values of them."""
    names = [str(name) for name in results["FIT.DOF_NAMES"] if str(name).startswith("IVC_EIG")]
    gsfit = np.array([
        controller.passives.get_array1(["IVC", "dof", name.removeprefix("IVC_").lower(), "calculated"])[0]
        for name in names
    ])  # fmt: skip
    return names, gsfit


@pytest.mark.te_network
def test_succeeds(case) -> None:  # type: ignore[no-untyped-def]
    _, results, _ = case
    assert_array_equal(results["GLOBAL.STATUS"], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    assert abs(results["GLOBAL.IP"][-1]) > 1e5
    assert results["BOUNDARY.OUTLINE.N"][-1] > 3


@pytest.mark.te_network
def test_topology(case) -> None:  # type: ignore[no-untyped-def]
    settings, results, controller = case
    expected = CASES[settings.get_path("replay/pulse"), settings.get_path("replay/time_start")]
    assert results["GLOBAL.XPT_DIVERTED"][-1] == expected
    assert _gsfit_equilibrium(controller)["GLOBAL.XPT_DIVERTED"] == expected


@pytest.mark.te_network
def test_matches_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    gsfit = _gsfit_equilibrium(controller)
    for name, rtol in RTOL.items():
        assert results[name][-1] == pytest.approx(gsfit[name], rel=rtol, abs=0.0), name
    for name, atol in ATOL.items():
        assert results[name][-1] == pytest.approx(gsfit[name], rel=0.0, abs=atol), name


@pytest.mark.te_network
@pytest.mark.parametrize("kind", SENSORS)
def test_sensors_match_gsfit(case, kind: str) -> None:  # type: ignore[no-untyped-def]
    settings, results, controller = case
    group, _, vs_measured, vs_gsfit = SENSORS[kind]
    reconstructed = results[f"CONSTRAINTS.{group}.ALL.RECONSTRUCT"][-1]
    measured, gsfit = _gsfit_sensors(settings, results, controller, kind)
    # The same measurements; GSFit's interpolation rounds differently.
    scale = float(np.abs(measured).max())
    assert_allclose(results[f"CONSTRAINTS.{group}.ALL.MEASURED"][-1], measured, rtol=0.0, atol=1e-14 * scale)
    assert_allclose(reconstructed, measured, rtol=0.0, atol=vs_measured)
    assert_allclose(reconstructed, gsfit, rtol=0.0, atol=vs_gsfit)


@pytest.mark.te_network
def test_ivc_eigenmode_currents_match_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    names, gsfit = _ivc_eigenmodes(results, controller)
    assert names
    all_names = list(results["FIT.DOF_NAMES"])
    rtgsfit = results["FIT.DOF_VALUES"][-1][[all_names.index(name) for name in names]]
    assert_allclose(rtgsfit, gsfit, rtol=0.0, atol=IVC_ATOL)


@pytest.mark.te_network
def test_ppac_replay() -> None:
    """The PPAC post-shot replay of the default pulse, without writing to MDSplus."""
    from rtgsfit.ppac import RtgsfitReplay

    controller = RtgsfitReplay(14816, write_to_mds=False, overrides={"replay/time_end": 0.1}, verbose=False)
    results = controller.run()
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


@pytest.mark.te_network
@pytest.mark.figures
def test_figures(case, figure_dir: Path) -> None:  # type: ignore[no-untyped-def]
    """What the tests above compare with GSFit, as figures in `<figure_dir>/st40_magrun16/<pulse>@<time>/`."""
    import plotting

    settings, results, controller = case
    gsfit = _gsfit_equilibrium(controller)
    constraints = {}
    for kind, (group, _, vs_measured, vs_gsfit) in SENSORS.items():
        measured, reconstructed = _gsfit_sensors(settings, results, controller, kind)
        constraints[group] = [
            plotting.Expected("GSFit reconstruction", reconstructed, atol=vs_gsfit),
            plotting.Expected("GSFit measured", measured, atol=vs_measured),
        ]
    reference = plotting.Reference(
        name="GSFit",
        r=gsfit["PROFILES_2D.R_Z.R"],
        z=gsfit["PROFILES_2D.R_Z.Z"],
        psi=gsfit["PROFILES_2D.R_Z.PSI"],
        j_phi=gsfit["PROFILES_2D.R_Z.J_PHI"],
        ip=gsfit["GLOBAL.IP"],
        axis_r=gsfit["GLOBAL.MAG_AXIS.R"],
        axis_z=gsfit["GLOBAL.MAG_AXIS.Z"],
        axis_psi=gsfit["GLOBAL.PSI_MAG_AXIS"],
        boundary_psi=gsfit["BOUNDARY.PSI"],
        outline=(gsfit["BOUNDARY.OUTLINE.R"], gsfit["BOUNDARY.OUTLINE.Z"]),
        sensors=plotting.gsfit_sensor_positions(controller, settings),
        constraints=constraints,
    )
    names, ivc = _ivc_eigenmodes(results, controller)
    filaments, currents = plotting.gsfit_passive_currents(results, controller)
    figures = {
        "psi": plotting.plot_map(results, reference, "PSI"),
        "j_phi": plotting.plot_map(results, reference, "J_PHI"),
        "constraints": plotting.plot_constraints(results, reference),
        "ivc_eigenmodes": plotting.plot_dofs(
            results, names, plotting.Expected("GSFit", ivc, atol=IVC_ATOL), "IVC eigenmode currents"
        ),
        "passive_currents": plotting.plot_passive_currents(results, filaments, currents, "GSFit"),
        "convergence": plotting.plot_convergence(results, reference, RTOL, ATOL),
    }
    pulse, time = settings.get_path("replay/pulse"), settings.get_path("replay/time_start")
    paths = plotting.save(figures, figure_dir / "st40_magrun16" / f"{pulse}@{time}")
    assert all(path.stat().st_size > 0 for path in paths)
