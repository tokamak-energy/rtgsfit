"""ST40 replays against GSFit and the MDSplus preshot: need the TE network (run with --te-network).

The vs-GSFit cases are the pulses and times of the v1 `rtgsfit_vs_gsfit` suite. RT-GSFit replays
the same time slice ten times (so it converges at a fixed time) with the preshot built by GSFit's
ST40 readers, and is compared with GSFit's own reconstruction of that time slice. The measurements
come from the PCS-equivalent MDSplus nodes (`st40_mdsplus`), as on the PCS, but the coil currents
are GSFit's, because the PCS PSU signals differ from them by up to 1 % (v1 used PSU2COIL here).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.replay import build_preshot
from rtgsfit.replay import read_preshot
from rtgsfit.replay import read_signals
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

pytestmark = pytest.mark.te_network

CASES = [(pulse, time) for pulse in (13343, 13345, 13346) for time in (0.030, 0.100)]

# Tolerances of RT-GSFit's last iteration against GSFit, 2-3 times the largest difference over the six cases. The two
# codes fit the same measurements with the same preshot and coil currents, so they differ only by their algorithms.
# Measured (largest of the six): psi 4.7e-3 Wb (at 100 ms, when RT-GSFit's psi changes by 2e-6 Wb per iteration, so
# not from convergence); j_phi 1.2e4 A/m^2 (peak 3.1e6); Ip 6.8e-5 relative; the magnetic axis 5.1e-4 m in R and
# 2.4e-4 m in Z (the cell is 0.03 m); psi at the axis 7.5e-4 Wb and at the boundary 6.3e-4 Wb; flux loops
# 2.9e-4 Wb, Bp probes 9.0e-4 T, Rogowski coils 2.4 A; IVC eigenmode currents 1.1 A.
RTOL = {"GLOBAL.IP": 2e-4}
"""Relative tolerances against GSFit."""
ATOL = {"GLOBAL.MAG_AXIS.R": 1.5e-3, "GLOBAL.MAG_AXIS.Z": 1.5e-3, "GLOBAL.PSI_MAG_AXIS": 2e-3, "BOUNDARY.PSI": 2e-3}
"""Absolute tolerances against GSFit [m, Wb]; psi at the boundary can be near 0, so not relative."""
PSI_ATOL = 1e-2
"""[Wb] at every grid point."""
J_PHI_ATOL = 3e4
"""[A/m^2] at every grid point."""
IVC_ATOL = 3.0
"""[A] for each IVC eigenmode current."""
# (group, GSFit sensor quantity, atol vs the measurements, atol vs GSFit's reconstruction) [Wb, T, A]. RT-GSFit's
# residuals are GSFit's to within the atol vs GSFit (|RT-GSFit - measured| <= |GSFit - measured| + atol), so the
# atol vs the measurements only bounds the fit of both codes: 2-3 times the largest residual of the six cases,
# flux loops 2.4e-2 Wb, Bp probes 6.2e-2 T and Rogowski coils 64 A.
SENSORS = {
    "flux_loop": ("FLUX_LOOP", "psi", 5e-2, 1e-3),
    "bp_probe": ("BP_PROBE", "b", 0.15, 2e-3),
    "rogowski": ("ROGOWSKI", "i", 150.0, 5.0),
}


def test_mdsplus_preshot_is_the_gsfit_build() -> None:
    """The PRESHOT in MDSplus (`preshot/reader`) equals a fresh build from GSFit (`preshot/builder`).

    Floating-point nodes are compared to round-off rather than by `content_hash`: the Poisson LU
    bands come from threaded LAPACK, whose rounding depends on the BLAS thread count.
    """
    settings = load_settings("st40")
    stored, built = read_preshot(settings), build_preshot(settings)
    assert built.groups.keys() == stored.groups.keys()
    for group, nodes in stored.groups.items():
        assert built.groups[group].keys() == nodes.keys(), group
        for name, value in nodes.items():
            expected, actual = np.asarray(value), np.asarray(built.groups[group][name])
            if expected.dtype.kind == "f":
                scale = float(np.abs(expected).max(initial=0.0))
                assert_allclose(actual, expected, rtol=1e-12, atol=1e-12 * scale, err_msg=f"{group}:{name}")
            else:
                assert_array_equal(actual, expected, err_msg=f"{group}:{name}")


def test_default_replay() -> None:
    settings = load_settings("st40").copy_with({"replay/time_end": 0.1})
    results = replay(settings, verbose=False)
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


@pytest.fixture(scope="module", params=CASES, ids=[f"{pulse}@{time}" for pulse, time in CASES])
def case(request: pytest.FixtureRequest):  # type: ignore[no-untyped-def]
    from rtgsfit.database_readers.gsfit import coil_circuits
    from rtgsfit.database_readers.gsfit import gsfit_coil_current
    from rtgsfit.database_readers.gsfit import gsfit_controller

    pulse, time = request.param
    settings = load_settings("st40").copy_with({
        "preshot/reader": "gsfit",
        "replay/pulse": pulse,
        "replay/time_start": time,
        "replay/time_end": time,
        "replay/repeat": 10,
    })  # fmt: skip
    controller = gsfit_controller(settings)
    preshot = read_preshot(settings)
    signals = read_signals(settings, preshot)
    # The PCS PSU signals differ from GSFit's coil currents by up to 1 %, which alone moves the
    # Bp probes by more than the tolerance. Use GSFit's currents so only the algorithms differ.
    coil_current = gsfit_coil_current(controller, coil_circuits(controller, settings), signals.time)
    signals = replace(signals, coil_signal=coil_current @ np.linalg.pinv(preshot["COIL:MATRIX"]).T)
    results = replay(settings, preshot=preshot, signals=signals, verbose=False)
    controller.inverse_solver_rust()
    return settings, results, controller


def _gsfit_equilibrium(controller) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    """GSFit's reconstruction of the time slice, by RT-GSFit's results node (on the same grid)."""
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


def test_succeeds(case) -> None:  # type: ignore[no-untyped-def]
    _, results, _ = case
    assert_array_equal(results["GLOBAL.STATUS"], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)


def test_matches_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    gsfit = _gsfit_equilibrium(controller)
    for name, rtol in RTOL.items():
        assert results[name][-1] == pytest.approx(gsfit[name], rel=rtol, abs=0.0), name
    for name, atol in ATOL.items():
        assert results[name][-1] == pytest.approx(gsfit[name], rel=0.0, abs=atol), name


def test_psi_and_j_phi_match_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    gsfit = _gsfit_equilibrium(controller)
    for name in ("PROFILES_2D.R_Z.R", "PROFILES_2D.R_Z.Z"):
        assert_allclose(results[name], gsfit[name], rtol=0.0, atol=1e-12, err_msg=name)
    assert_allclose(results["PROFILES_2D.R_Z.PSI"][-1], gsfit["PROFILES_2D.R_Z.PSI"], rtol=0.0, atol=PSI_ATOL)
    assert_allclose(results["PROFILES_2D.R_Z.J_PHI"][-1], gsfit["PROFILES_2D.R_Z.J_PHI"], rtol=0.0, atol=J_PHI_ATOL)


@pytest.mark.parametrize("kind", SENSORS)
def test_sensors_match_gsfit(case, kind: str) -> None:  # type: ignore[no-untyped-def]
    settings, results, controller = case
    group, _, vs_measured, vs_gsfit = SENSORS[kind]
    reconstructed = results[f"CONSTRAINTS.{group}.ALL.RECONSTRUCT"][-1]
    measured, gsfit = _gsfit_sensors(settings, results, controller, kind)
    assert_allclose(reconstructed, measured, rtol=0.0, atol=vs_measured)
    assert_allclose(reconstructed, gsfit, rtol=0.0, atol=vs_gsfit)


def test_ivc_eigenmode_currents_match_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    names, gsfit = _ivc_eigenmodes(results, controller)
    assert names
    all_names = list(results["FIT.DOF_NAMES"])
    rtgsfit = results["FIT.DOF_VALUES"][-1][[all_names.index(name) for name in names]]
    assert_allclose(rtgsfit, gsfit, rtol=0.0, atol=IVC_ATOL)


def test_results_are_complete(case) -> None:  # type: ignore[no-untyped-def]
    """Everything is stored: the settings, the preshot and the raw PCS inputs, as in MDSplus."""
    settings, results, _ = case
    assert isinstance(results, Results)
    assert results["INPUT.SETTINGS"] == settings.json_text
    assert results["INPUT.MEAS"].shape[1] == results["PRESHOT.SIZES.N_MEAS"]


@pytest.mark.figures
def test_figures(case, figure_dir: Path) -> None:  # type: ignore[no-untyped-def]
    """What the tests above compare with GSFit, as figures in `<figure_dir>/st40/<pulse>@<time>/`."""
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
        "psi": plotting.plot_map(results, reference, "PSI", atol=PSI_ATOL),
        "j_phi": plotting.plot_map(results, reference, "J_PHI", atol=J_PHI_ATOL),
        "constraints": plotting.plot_constraints(results, reference),
        "ivc_eigenmodes": plotting.plot_dofs(
            results, names, plotting.Expected("GSFit", ivc, atol=IVC_ATOL), "IVC eigenmode currents"
        ),
        "passive_currents": plotting.plot_passive_currents(results, filaments, currents, "GSFit"),
        "convergence": plotting.plot_convergence(results, reference, RTOL, ATOL),
    }
    pulse, time = settings.get_path("replay/pulse"), settings.get_path("replay/time_start")
    paths = plotting.save(figures, figure_dir / "st40" / f"{pulse}@{time}")
    assert all(path.stat().st_size > 0 for path in paths)
