import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose

from rtgsfit_verify_analytic import analytic_soln, cnst

@pytest.fixture(scope="session")
def output_dict():
    output_file = Path(cnst.DATA_DIR) / "output_dict.npy"
    if not output_file.exists():
        pytest.fail(f"{output_file} does not exist. Run test_run_pipeline first.")
    return np.load(output_file, allow_pickle=True).item()

def test_psi_matches_expected(output_dict):
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    expected_psi = analytic_soln.analytic_psi(
        r_grid, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB, 
        cnst.ANALYTIC_PLASMA_CURRENT
    )
    actual_psi = output_dict['flux_total'][cnst.N_ITERS - 1].reshape(cnst.N_Z, cnst.N_R)
    assert_allclose(actual_psi, expected_psi, rtol=1e-2)

def test_b_theta_matches_expected(output_dict):
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)

    expected_b_theta = analytic_soln.analytic_b_theta(
        r_vec, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB,
        cnst.ANALYTIC_PLASMA_CURRENT
    )
    
    psi_num = output_dict['flux_total'][cnst.N_ITERS - 1].reshape(cnst.N_Z, cnst.N_R)
    br = -np.gradient(psi_num, z_vec, axis=0) / r_grid
    bz = np.gradient(psi_num, r_vec, axis=1) / r_grid
    theta = np.arctan2((z_grid - cnst.ANALYTIC_ZO), (r_grid - cnst.ANALYTIC_RO))
    actual_b_theta = -br * np.sin(theta) + bz * np.cos(theta)

    # Peak |B_theta| is about 0.5 T
    assert_allclose(actual_b_theta, expected_b_theta, atol=1e-2)

def test_current_density_matches_expected(output_dict):
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)

    expected_current_density = analytic_soln.analytic_j_phi(
        r_grid, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB,
        cnst.ANALYTIC_PLASMA_CURRENT
    )

    coef = output_dict['coef'][cnst.N_ITERS - 1]
    flux_norm = output_dict['flux_norm'][cnst.N_ITERS - 1].reshape(cnst.N_Z, cnst.N_R)
    flux_norm_z = np.gradient(flux_norm, z_vec, axis=0)

    current0 = coef[0] * (1 - flux_norm) * r_grid
    current1 = coef[1] * (1 - flux_norm) / (r_grid * cnst.MU_0)
    current2 = coef[2] * flux_norm_z
    actual_current_density = current0 + current1 + current2

    # Peak current density is about 3e6 A/m^2
    assert_allclose(actual_current_density, expected_current_density, atol=1e4)


