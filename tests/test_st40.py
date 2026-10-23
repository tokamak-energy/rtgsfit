"""ST40 replays against GSFit and the MDSplus preshot: need the TE network (run with --te-network).

The vs-GSFit cases are the pulses and times of the v1 `rtgsfit_vs_gsfit` suite. RT-GSFit replays
the same time slice ten times (so it converges at a fixed time) with the preshot built by GSFit's
ST40 readers, and is compared with GSFit's own reconstruction of that time slice. The measurements
come from the PCS-equivalent MDSplus nodes (`st40_mdsplus`), as on the PCS, but the coil currents
are GSFit's, because the PCS PSU signals differ from them by up to 1 % (v1 used PSU2COIL here).
"""

from __future__ import annotations

from dataclasses import replace

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

# (GSFit sensor quantity, vs measured rtol/atol, vs GSFit's reconstruction rtol/atol); from v1
SENSORS = {
    "flux_loop": ("FLUX_LOOP", "psi", (5e-2, 5e-2), (1e-3, 1e-3)),
    "bp_probe": ("BP_PROBE", "b", (0.1, 0.5), (1e-3, 1e-3)),
    "rogowski": ("ROGOWSKI", "i", (5e-2, 5e2), (5e-2, 10.0)),
}


def test_mdsplus_preshot_is_the_gsfit_build() -> None:
    """The PRESHOT in MDSplus (`preshot/reader`) equals a fresh build from GSFit (`preshot/builder`)."""
    settings = load_settings("st40")
    assert read_preshot(settings).content_hash() == build_preshot(settings).content_hash()


def test_default_replay() -> None:
    settings = load_settings("st40").copy_with({"replay/time_end": 0.1})
    results = replay(settings, verbose=False)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    flat_top = results["TIME"] > 0.03
    assert_array_equal(results["GLOBAL.STATUS"][flat_top], 0)
    assert np.all(np.abs(results["GLOBAL.IP"][flat_top]) > 1e5)


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


def test_succeeds(case) -> None:  # type: ignore[no-untyped-def]
    _, results, _ = case
    assert_array_equal(results["GLOBAL.STATUS"], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)


def test_psi_matches_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    from gsfit_rs.imas import equilibrium_paths

    _, results, controller = case
    psi_gsfit = np.asarray(controller.plasma.get(equilibrium_paths.time_slice[0].profiles_2d[0].psi))
    assert_allclose(results["PROFILES_2D.R_Z.PSI"][-1], psi_gsfit, rtol=0.33, atol=5e-3)


@pytest.mark.parametrize("kind", SENSORS)
def test_sensors_match_gsfit(case, kind: str) -> None:  # type: ignore[no-untyped-def]
    from rtgsfit.database_readers.gsfit import constraint_sensors
    from rtgsfit.database_readers.gsfit import sensor_name

    settings, results, controller = case
    group, quantity, vs_measured, vs_gsfit = SENSORS[kind]
    sensors = {
        "flux_loop": controller.flux_loops,
        "bp_probe": controller.bp_probes,
        "rogowski": controller.rogowski_coils,
    }[kind]
    names = list(results[f"CONSTRAINTS.{group}.ALL.NAMES"])
    reconstructed = results[f"CONSTRAINTS.{group}.ALL.RECONSTRUCT"][-1]
    measured, gsfit = np.zeros(len(names)), np.zeros(len(names))
    for sensor in (sensor for sensor in constraint_sensors(controller, settings) if sensor.kind == kind):
        index = names.index(sensor_name(sensor, "st40_pcs"))
        measured[index] = sensors.get_array1([sensor.name, quantity, "measured", "value"])[0]
        gsfit[index] = sensors.get_array1([sensor.name, quantity, "calculated", "value"])[0]
    assert_allclose(reconstructed, measured, rtol=vs_measured[0], atol=vs_measured[1])
    assert_allclose(reconstructed, gsfit, rtol=vs_gsfit[0], atol=vs_gsfit[1])


def test_ivc_eigenmode_currents_match_gsfit(case) -> None:  # type: ignore[no-untyped-def]
    _, results, controller = case
    names = list(results["FIT.X_VEC_NAMES"])
    x_vec = results["FIT.X_VEC"][-1]
    eigenmodes = [(index, name) for index, name in enumerate(names) if name.startswith("IVC_EIG")]
    assert eigenmodes
    rtgsfit = np.array([x_vec[index] for index, _ in eigenmodes])
    gsfit = np.array([
        controller.passives.get_array1(["IVC", "dof", name.removeprefix("IVC_").lower(), "calculated"])[0]
        for _, name in eigenmodes
    ])  # fmt: skip
    assert_allclose(rtgsfit, gsfit, rtol=1e-3, atol=3.0)


def test_results_are_complete(case) -> None:  # type: ignore[no-untyped-def]
    """Everything is stored: the settings, the preshot and the raw PCS inputs, as in MDSplus."""
    settings, results, _ = case
    assert isinstance(results, Results)
    assert results["INPUT.SETTINGS"] == settings.xml_text
    assert results["INPUT.MEAS"].shape[1] == results["PRESHOT.SIZES.N_MEAS"]
