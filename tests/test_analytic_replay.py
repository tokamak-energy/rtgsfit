"""The analytic replay: a circular plasma with a known exact solution (no external data needed).

The same measurements are replayed at every time slice, so the slices are iterations that must
converge on the exact solution in `rtgsfit.database_readers.analytic.solution`. The v1 reference
(`data/analytic_v1_reference.npz`) is the last iteration of the v1 code on the same problem, with
the v1 names translated to the v2 node names (see the migration table in the manual). It was made
with `tests/rtgsfit_verify_analytic` at commit 1139238, the last commit with the v1 code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.constants import MU_0
from rtgsfit.database_readers import analytic
from rtgsfit.database_readers.analytic import analytic_plasma
from rtgsfit.database_readers.analytic.solution import b_theta_from_psi
from rtgsfit.library import STATUS_FLAGS
from rtgsfit.library import Rtgsfit
from rtgsfit.preshot import Preshot
from rtgsfit.replay import read_signals
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

DATA_DIR = Path(__file__).parent / "data"
CONSTRAINT_GROUPS = ("FLUX_LOOP", "BP_PROBE", "ROGOWSKI")
PSI_RTOL = 1e-2
J_PHI_ATOL = 1e4  # the peak current density is about 3e6 A/m^2


def _grid(results: Results) -> tuple[np.ndarray, np.ndarray]:
    return np.meshgrid(results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"])


def _constraints(results: Results, signal: str) -> np.ndarray:
    """`fit_target` (MEASURED) or `fit_prediction` (RECONSTRUCT) of the last slice, in the order of the fit."""
    return np.concatenate([results[f"CONSTRAINTS.{group}.ALL.{signal}"][-1] for group in CONSTRAINT_GROUPS])


@pytest.fixture(scope="module")
def plasma():  # type: ignore[no-untyped-def]
    return analytic_plasma(load_settings("analytic"))


# ----------------------------------------------------------------------------------------------
# Against the exact solution
# ----------------------------------------------------------------------------------------------
def test_every_slice_succeeds(analytic_results: Results) -> None:
    assert_array_equal(analytic_results["GLOBAL.STATUS"], 0)
    assert_array_equal(analytic_results["FIT.LAPACK_INFO"], 0)


def test_psi_matches_exact_solution(analytic_results: Results, plasma) -> None:  # type: ignore[no-untyped-def]
    r_2d, z_2d = _grid(analytic_results)
    assert_allclose(analytic_results["PROFILES_2D.R_Z.PSI"][-1], plasma.psi(r_2d, z_2d), rtol=PSI_RTOL)


def test_b_theta_matches_exact_solution(analytic_results: Results, plasma) -> None:  # type: ignore[no-untyped-def]
    r_2d, z_2d = _grid(analytic_results)
    b_theta = b_theta_from_psi(
        analytic_results["PROFILES_2D.R_Z.PSI"][-1],
        analytic_results["PROFILES_2D.R_Z.R"],
        analytic_results["PROFILES_2D.R_Z.Z"],
        plasma.r_0,
        plasma.z_0,
    )
    assert_allclose(b_theta, plasma.b_theta(r_2d, z_2d), atol=1e-2)


def test_j_phi_matches_exact_solution(analytic_results: Results, plasma) -> None:  # type: ignore[no-untyped-def]
    r_2d, z_2d = _grid(analytic_results)
    assert_allclose(analytic_results["PROFILES_2D.R_Z.J_PHI"][-1], plasma.j_phi(r_2d, z_2d), atol=J_PHI_ATOL)


def test_global_quantities_match_exact_solution(analytic_results: Results, plasma) -> None:  # type: ignore[no-untyped-def]
    assert analytic_results["GLOBAL.IP"][-1] == pytest.approx(plasma.ip, rel=1e-6)
    assert analytic_results["GLOBAL.MAG_AXIS.R"][-1] == pytest.approx(plasma.r_0, abs=1e-3)
    assert analytic_results["GLOBAL.MAG_AXIS.Z"][-1] == pytest.approx(plasma.z_0, abs=1e-3)
    assert analytic_results["BOUNDARY.PSI"][-1] == pytest.approx(plasma.boundary_psi, rel=1e-3)
    assert analytic_results["GLOBAL.PSI_MAG_AXIS"][-1] == pytest.approx(plasma.axis_psi, rel=1e-3)
    assert analytic_results["GLOBAL.XPT_DIVERTED"][-1] == 0


def test_boundary_shape_is_a_circle(analytic_results: Results, plasma) -> None:  # type: ignore[no-untyped-def]
    assert analytic_results["BOUNDARY.GEO_AXIS.R"][-1] == pytest.approx(plasma.r_0, abs=1e-3)
    assert analytic_results["BOUNDARY.GEO_AXIS.Z"][-1] == pytest.approx(plasma.z_0, abs=1e-3)
    assert analytic_results["BOUNDARY.MINOR_RADIUS"][-1] == pytest.approx(plasma.rho_b, rel=1e-3)
    assert analytic_results["BOUNDARY.ELONGATION"][-1] == pytest.approx(1.0, abs=1e-3)
    assert analytic_results["BOUNDARY.TRIANG"][-1] == pytest.approx(0.0, abs=1e-3)
    n_outline = int(analytic_results["BOUNDARY.OUTLINE.N"][-1])
    r_outline = analytic_results["BOUNDARY.OUTLINE.R"][-1][:n_outline]
    z_outline = analytic_results["BOUNDARY.OUTLINE.Z"][-1][:n_outline]
    assert_allclose(np.hypot(r_outline - plasma.r_0, z_outline - plasma.z_0), plasma.rho_b, rtol=2e-3)


def test_iterations_converge(analytic_results: Results) -> None:
    chi_sq = analytic_results["CONSTRAINTS.CHI_SQ_MAG"]
    assert chi_sq[-1] <= chi_sq[0]
    psi = analytic_results["PROFILES_2D.R_Z.PSI"]
    steps = np.max(np.abs(np.diff(psi, axis=0)), axis=(1, 2))
    # linear convergence: every step after the first few is smaller than the one before
    assert np.all(np.diff(steps[2:]) < 0)
    assert steps[-1] < 1e-4 * np.max(np.abs(psi[-1]))


def test_replay_is_deterministic(analytic_results: Results, analytic_preshot: Preshot) -> None:
    """A second replay with the same inputs must give bit-identical results."""
    again = replay(load_settings("analytic"), preshot=analytic_preshot, verbose=False)
    assert_array_equal(again["PROFILES_2D.R_Z.PSI"], analytic_results["PROFILES_2D.R_Z.PSI"])
    assert_array_equal(again["FIT.DOF_VALUES"], analytic_results["FIT.DOF_VALUES"])


@pytest.mark.figures
def test_figures(analytic_results: Results, plasma, figure_dir: Path) -> None:  # type: ignore[no-untyped-def]
    """Diagnostic figures (`pytest --figures`): the iterations against the exact solution."""
    import plotting

    settings = load_settings("analytic")
    r_2d, z_2d = _grid(analytic_results)
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    loop_r, loop_z = analytic.flux_loop_geometry(settings)
    probe_r, probe_z, _ = analytic.probe_geometry(settings)
    exact = plotting.Reference(
        name="exact",
        r=analytic_results["PROFILES_2D.R_Z.R"],
        z=analytic_results["PROFILES_2D.R_Z.Z"],
        psi=plasma.psi(r_2d, z_2d),
        j_phi=plasma.j_phi(r_2d, z_2d),
        ip=plasma.ip,
        axis_r=plasma.r_0,
        axis_z=plasma.z_0,
        axis_psi=plasma.axis_psi,
        boundary_psi=plasma.boundary_psi,
        outline=(plasma.r_0 + plasma.rho_b * np.cos(theta), plasma.z_0 + plasma.rho_b * np.sin(theta)),
        sensors={"flux loops": (loop_r, loop_z), "Bp probes": (probe_r, probe_z)},
    )
    figures = {
        "psi": plotting.plot_map(analytic_results, exact, "PSI", rtol=PSI_RTOL),
        "j_phi": plotting.plot_map(analytic_results, exact, "J_PHI", atol=J_PHI_ATOL),
        "midplane": plotting.plot_midplane(analytic_results, plasma),
        "constraints": plotting.plot_constraints(analytic_results),
        "convergence": plotting.plot_convergence(analytic_results, exact),
    }
    paths = plotting.save(figures, figure_dir / "analytic")
    assert all(path.stat().st_size > 0 for path in paths)


# ----------------------------------------------------------------------------------------------
# Against v1
# ----------------------------------------------------------------------------------------------
V1_KEYS = (
    "PROFILES_2D.R_Z.PSI",
    "PROFILES_2D.R_Z.PSI_NORM",
    "PROFILES_2D.R_Z.MASK",
    "FIT.DOF_VALUES",
    "BOUNDARY.PSI",
    "GLOBAL.PSI_MAG_AXIS",
    "GLOBAL.MAG_AXIS.R",
    "GLOBAL.MAG_AXIS.Z",
    "GLOBAL.IP",
    "GLOBAL.CURRENT_CENT.R",
    "GLOBAL.CURRENT_CENT.Z",
    "CONSTRAINTS.CHI_SQ_MAG",
    "XPT.ALL.N",
    "GLOBAL.XPT_DIVERTED",
    "GLOBAL.STATUS",
    "FIT.LAPACK_INFO",
)
# Names in the v1 reference file that have since been renamed
V1_OLD_KEYS = {"FIT.DOF_VALUES": "FIT.X_VEC", "XPT.ALL.N": "XPOINTS.ALL.N"}


def _v1_results(preshot: Preshot, ip: float | None) -> Results:
    """Replay v1's exact inputs: v1 read the sensor positions back from text files written with
    six decimals and computed the measurements from those rounded positions (its Green's tables
    used the exact positions). The fit is ill-conditioned along the P_PRIME / DELTA_Z directions,
    so the ~1e-7 change in the measurements moves DOF_VALUES by up to 1 % while psi moves by only 1e-7."""
    settings = load_settings("analytic")
    if ip is not None:
        settings = settings.copy_with({"analytic/plasma/ip": ip})

    def rounded(function):  # type: ignore[no-untyped-def]
        return lambda settings: tuple(np.round(array, 6) for array in function(settings))

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(analytic, "flux_loop_geometry", rounded(analytic.flux_loop_geometry))
        patch.setattr(analytic, "probe_geometry", rounded(analytic.probe_geometry))
        signals = read_signals(settings, preshot)
    return replay(settings, preshot=preshot, signals=signals, verbose=False)


def _assert_matches_v1(results: Results, case: str) -> None:
    reference = np.load(DATA_DIR / "analytic_v1_reference.npz")
    actual = {key: np.ravel(results[key][-1]) for key in V1_KEYS}
    actual["B_PRED"] = _constraints(results, "RECONSTRUCT")
    actual["B_VEC_MEAS"] = _constraints(results, "MEASURED")
    for key, value in actual.items():
        old_key = V1_OLD_KEYS.get(key, key)
        expected = np.ravel(reference[f"{case}/{old_key}"]).copy()
        if key == "FIT.DOF_VALUES":
            expected[1] /= MU_0
        if key in ("PROFILES_2D.R_Z.PSI", "BOUNDARY.PSI", "GLOBAL.PSI_MAG_AXIS"):
            expected *= 2 * np.pi
        if key in ("B_PRED", "B_VEC_MEAS"):
            expected[: results["PRESHOT.SIZES.N_FLUX_LOOP"]] *= 2 * np.pi
        rtol = 1e-4
        atol = rtol * float(np.max(np.abs(expected))) if expected.size else 0.0
        assert_allclose(value, expected, rtol=rtol, atol=atol, err_msg=key)


def test_matches_v1(analytic_preshot: Preshot) -> None:
    _assert_matches_v1(_v1_results(analytic_preshot, None), "main")


def test_zero_current_fails_cleanly(analytic_zero_results: Results) -> None:
    """No plasma current: no magnetic axis, so every slice reports NO_AXIS and zero current."""
    no_axis = next(bit for bit, name in STATUS_FLAGS.items() if name == "NO_AXIS")
    assert_array_equal(analytic_zero_results["GLOBAL.STATUS"], no_axis)
    assert_array_equal(analytic_zero_results["GLOBAL.IP"], 0.0)
    assert np.all(np.isnan(analytic_zero_results["BOUNDARY.MINOR_RADIUS"]))


def test_zero_current_matches_v1(analytic_preshot: Preshot) -> None:
    _assert_matches_v1(_v1_results(analytic_preshot, 0.0), "zero")


def test_wrong_input_shape_is_rejected(analytic_preshot: Preshot) -> None:
    """A 2D array with the right length but no data must not reach C, which trusts the lengths."""
    rtgsfit = Rtgsfit.build(analytic_preshot)
    n_meas, n_coil = analytic_preshot.sizes["N_MEAS"], analytic_preshot.sizes["N_COIL"]
    with pytest.raises(ValueError, match="shape"):
        rtgsfit(np.empty((n_meas, 0)), np.zeros(n_coil))
    with pytest.raises(ValueError, match="shape"):
        rtgsfit(np.zeros(n_meas), np.zeros(n_coil + 1))
