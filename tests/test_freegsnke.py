"""Replay a FreeGSNKE equilibrium (MAST-U-like machine): no TE network, but needs GSFit and FreeGSNKE.

FreeGSNKE solves a diverted equilibrium, GSFit turns its machine description and synthetic
measurements into a preshot and signals, and RT-GSFit must reconstruct the equilibrium. FreeGSNKE's
p' and FF' are linear in `psi_norm`, as in RT-GSFit's current model, on the same grid, so the
reconstruction can be compared closely.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

pytestmark = pytest.mark.freegsnke

# Measured with the bundled settings: Ip 1e-4, psi at the axis 2.3e-3 and at the boundary 1.6e-3 relative;
# the axis 0.025 of a cell; psi 1.8e-3 Wb; j_phi 3.1e3 A/m^2 (peak 8.8e5); P_PRIME 0.4 %, FF_PRIME_MU0 0.2 %
RTOL = {"GLOBAL.IP": 1e-3, "GLOBAL.PSI_MAG_AXIS": 1e-2, "BOUNDARY.PSI": 5e-3}
"""Relative tolerances of RT-GSFit's last iteration against FreeGSNKE."""
AXIS_CELLS = 0.1
"""The magnetic axis must be within this fraction of a grid cell of FreeGSNKE's, in R and in Z."""
PSI_ATOL = 5e-3
"""[Wb] at every grid point."""
J_PHI_ATOL = 1e4
"""[A/m^2] at every grid point."""
DOF_RTOL = 1e-2
"""P_PRIME and FF_PRIME_MU0 against FreeGSNKE's coefficients of the same basis functions."""


@pytest.fixture(scope="module")
def settings():  # type: ignore[no-untyped-def]
    settings = load_settings("freegsnke")
    section = settings["freegsnke"]
    if (
        not Path(section["machine_dir"]).expanduser().is_dir()
        or not Path(section["currents_path"]).expanduser().is_file()
    ):
        pytest.skip("needs the FreeGSNKE repository (machine pickles and example currents)")
    return settings


@pytest.fixture(scope="module")
def results(settings) -> Results:  # type: ignore[no-untyped-def]
    return replay(settings, verbose=False)


@pytest.fixture(scope="module")
def equilibrium(settings, results):  # type: ignore[no-untyped-def]
    from rtgsfit.database_readers.freegsnke import freegsnke_equilibrium

    return freegsnke_equilibrium(settings)


@pytest.fixture(scope="module")
def profiles(settings, results):  # type: ignore[no-untyped-def]
    from rtgsfit.database_readers.freegsnke import freegsnke_profiles

    return freegsnke_profiles(settings)


def test_converges_without_errors(results: Results) -> None:
    assert_array_equal(results["GLOBAL.STATUS"][-5:], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    psi = results["PROFILES_2D.R_Z.PSI"]
    assert np.max(np.abs(psi[-1] - psi[-2])) < 1e-6 * np.max(np.abs(psi[-1]))


def _grid_cell(results: Results) -> tuple[float, float]:
    return float(np.diff(results["PROFILES_2D.R_Z.R"])[0]), float(np.diff(results["PROFILES_2D.R_Z.Z"])[0])


def _atol(results: Results) -> dict[str, float]:
    """Absolute tolerances of RT-GSFit's last iteration against FreeGSNKE: the magnetic axis [m]."""
    d_r, d_z = _grid_cell(results)
    return {"GLOBAL.MAG_AXIS.R": AXIS_CELLS * d_r, "GLOBAL.MAG_AXIS.Z": AXIS_CELLS * d_z}


def test_matches_freegsnke(results: Results, equilibrium) -> None:  # type: ignore[no-untyped-def]
    # FreeGSNKE's psi is per radian
    axis_r, axis_z = equilibrium.magneticAxis()[:2]
    expected = {
        "GLOBAL.IP": equilibrium.plasmaCurrent(),
        "GLOBAL.PSI_MAG_AXIS": 2.0 * np.pi * equilibrium.psi_axis,
        "BOUNDARY.PSI": 2.0 * np.pi * equilibrium.psi_bndry,
        "GLOBAL.MAG_AXIS.R": axis_r,
        "GLOBAL.MAG_AXIS.Z": axis_z,
    }
    for name, rtol in RTOL.items():
        assert results[name][-1] == pytest.approx(expected[name], rel=rtol, abs=0.0), name
    for name, atol in _atol(results).items():
        assert results[name][-1] == pytest.approx(expected[name], rel=0.0, abs=atol), name
    assert results["GLOBAL.XPT_DIVERTED"][-1] == 1


def test_psi_and_j_phi_match_freegsnke(results: Results, equilibrium, profiles) -> None:  # type: ignore[no-untyped-def]
    # FreeGSNKE's grids are [n_r, n_z], RT-GSFit's [n_z, n_r]; the grids are the same
    assert_allclose(results["PROFILES_2D.R_Z.R"], equilibrium.R[:, 0], rtol=0.0, atol=1e-12)
    assert_allclose(results["PROFILES_2D.R_Z.Z"], equilibrium.Z[0, :], rtol=0.0, atol=1e-12)
    assert_allclose(results["PROFILES_2D.R_Z.PSI"][-1], 2.0 * np.pi * equilibrium.psi().T, rtol=0.0, atol=PSI_ATOL)
    assert_allclose(results["PROFILES_2D.R_Z.J_PHI"][-1], profiles.jtor.T, rtol=0.0, atol=J_PHI_ATOL)


def test_profile_matches_freegsnke(settings, results: Results, profiles) -> None:  # type: ignore[no-untyped-def]
    """FreeGSNKE's j_phi = L (Beta0 R / Raxis + (1 - Beta0) Raxis / R) (1 - psi_norm) is in RT-GSFit's basis, so the
    fitted coefficients of (1 - psi_norm) R and (1 - psi_norm) / R are FreeGSNKE's."""
    profile = settings["freegsnke"]["profile"]
    assert profile["alpha_m"] == profile["alpha_n"] == 1.0, "p' and FF' must be linear in psi_norm"
    dofs = dict(zip(results["FIT.DOF_NAMES"], results["FIT.DOF_VALUES"][-1], strict=True))
    assert dofs["P_PRIME"] == pytest.approx(profiles.L * profiles.Beta0 / profiles.Raxis, rel=DOF_RTOL)
    assert dofs["FF_PRIME_MU0"] == pytest.approx(profiles.L * (1.0 - profiles.Beta0) * profiles.Raxis, rel=DOF_RTOL)


@pytest.mark.figures
def test_figures(settings, results: Results, equilibrium, profiles, figure_dir: Path) -> None:  # type: ignore[no-untyped-def]
    """What the tests above compare with FreeGSNKE, as figures in `<figure_dir>/freegsnke/`."""
    import plotting
    from rtgsfit.database_readers.freegsnke import _controller

    # FreeGSNKE's grids are [n_r, n_z] and its psi is per radian
    axis_r, axis_z = equilibrium.magneticAxis()[:2]
    outline = np.asarray(equilibrium.separatrix())
    outline = np.vstack([outline, outline[:1]])
    reference = plotting.Reference(
        name="FreeGSNKE",
        r=equilibrium.R[:, 0],
        z=equilibrium.Z[0, :],
        psi=2.0 * np.pi * equilibrium.psi().T,
        j_phi=np.asarray(profiles.jtor).T,
        ip=float(equilibrium.plasmaCurrent()),
        axis_r=float(axis_r),
        axis_z=float(axis_z),
        axis_psi=2.0 * np.pi * float(equilibrium.psi_axis),
        boundary_psi=2.0 * np.pi * float(equilibrium.psi_bndry),
        outline=(outline[:, 0], outline[:, 1]),
        xpt=(equilibrium.xpt[:, 0], equilibrium.xpt[:, 1]),
        sensors=plotting.gsfit_sensor_positions(_controller(settings), settings),
    )
    figures = {
        "psi": plotting.plot_map(results, reference, "PSI", atol=PSI_ATOL),
        "j_phi": plotting.plot_map(results, reference, "J_PHI", atol=J_PHI_ATOL),
        "constraints": plotting.plot_constraints(results),
        "convergence": plotting.plot_convergence(results, reference, RTOL, _atol(results)),
    }
    paths = plotting.save(figures, figure_dir / "freegsnke")
    assert all(path.stat().st_size > 0 for path in paths)
