"""
Tests for utility/poisson_stencil.py, which recovers the POISSON_A, POISSON_B and POISSON_C
arrays of constants.c from the banded LU factors that GSFit stores in the PRESHOT node
(the path used by utility/mdsplus_const_to_file.py in the ST40 PCS build).
"""

import importlib.util
import os

import numpy as np
import pytest
from scipy.linalg import lu

from rtgsfit_verify_analytic import cnst, poisson_matrix


def load_poisson_stencil_module():
    path = os.path.join(cnst.RTGSFIT_PATH, "utility", "poisson_stencil.py")
    spec = importlib.util.spec_from_file_location("poisson_stencil", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def lup_bands_as_gsfit_writes_them(poiss_matrix, n_r):
    """
    LOWER_BAND, UPPER_BAND and PERM_IDX exactly as GSFit's RTGSFIT writer produces them
    (gsfit/database_writers/rtgsfit_mdsplus/poisson_matrix.py, poisson_matrix_lup_bands).
    """
    perm, lower, upper = lu(poiss_matrix)
    perm = perm.T
    lower_pad = np.pad(lower, ((0, 0), (n_r, 0)), mode="constant", constant_values=0)
    upper_pad = np.pad(upper, ((0, 0), (0, n_r + 1)), mode="constant", constant_values=0)
    n_grid = np.shape(poiss_matrix)[0]
    lower_band = np.zeros((n_grid, n_r), dtype=np.float64)
    upper_band = np.zeros((n_grid, n_r + 2), dtype=np.float64)
    for ii in range(n_grid):
        lower_band[ii, :] = lower_pad[ii, ii : ii + n_r]
        upper_band[ii, :] = upper_pad[ii, ii : ii + n_r + 2]
    upper_band[:, 0] = 1 / upper_band[:, 0]
    _, perm_idx = np.nonzero(perm)
    return lower_band.flatten(), upper_band.flatten(), perm_idx


def dense_operator(a, b, c, n_r, n_z):
    """Dense operator with identity boundary rows and the column-dependent stencil a, b, c."""
    n_grid = n_r * n_z
    mat = np.zeros((n_grid, n_grid), dtype=np.float64)
    for i in range(n_z):
        for j in range(n_r):
            idx = i * n_r + j
            if i == 0 or i == n_z - 1 or j == 0 or j == n_r - 1:
                mat[idx, idx] = 1.0
            else:
                mat[idx, idx] = -c[j]
                mat[idx, idx + 1] = a[j]
                mat[idx, idx - 1] = b[j]
                mat[idx, idx + n_r] = 1.0
                mat[idx, idx - n_r] = 1.0
    return mat


def gsfit_stencil(r_vec, z_vec):
    """The LIUQE eq. (45) coefficients GSFit's writer uses for the ST40 constants."""
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]
    a = (dz / dr) ** 2 * r_vec / (r_vec + dr / 2)
    b = (dz / dr) ** 2 * r_vec / (r_vec - dr / 2)
    return a, b, 2 + a + b


def test_recovers_the_st40_operator_written_by_gsfit():
    """ST40 grid, radius-dependent eq. (45) coefficients; partial pivoting permutes rows."""
    ps = load_poisson_stencil_module()
    r_vec = np.linspace(110.000e-3, 1.070e0, 33)
    z_vec = np.linspace(-960.000e-3, 960.000e-3, 65)
    a, b, c = gsfit_stencil(r_vec, z_vec)
    lower_band, upper_band, perm_idx = lup_bands_as_gsfit_writes_them(
        dense_operator(a, b, c, r_vec.size, z_vec.size), r_vec.size
    )
    assert not np.array_equal(perm_idx, np.arange(perm_idx.size)), "expected a nontrivial permutation"

    poisson_a, poisson_b, poisson_c = ps.stencil_from_lu_bands(
        lower_band, upper_band, perm_idx, r_vec.size, z_vec.size
    )

    np.testing.assert_allclose(poisson_a[1:-1], a[1:-1], rtol=1e-13, atol=0)
    np.testing.assert_allclose(poisson_b[1:-1], b[1:-1], rtol=1e-13, atol=0)
    np.testing.assert_allclose(poisson_c[1:-1], c[1:-1], rtol=1e-13, atol=0)
    for arr in (poisson_a, poisson_b, poisson_c):
        assert arr[0] == 0.0 and arr[-1] == 0.0


def test_matches_the_analytic_generator():
    """The recovered stencil equals what the analytic test writes to constants.c directly."""
    ps = load_poisson_stencil_module()
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, 17)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, 33)
    lower_band, upper_band, perm_idx = lup_bands_as_gsfit_writes_them(
        poisson_matrix.poisson_matrix(r_vec, z_vec), r_vec.size
    )

    recovered = ps.stencil_from_lu_bands(lower_band, upper_band, perm_idx, r_vec.size, z_vec.size)
    written = poisson_matrix.poisson_stencil(r_vec, z_vec)

    for got, expected in zip(recovered, written):
        np.testing.assert_allclose(got[1:-1], expected[1:-1], rtol=1e-13, atol=0)


def test_rejects_an_operator_the_separable_solver_cannot_represent():
    ps = load_poisson_stencil_module()
    r_vec = np.linspace(0.2, 1.0, 9)
    z_vec = np.linspace(-0.5, 0.5, 11)
    a, b, c = gsfit_stencil(r_vec, z_vec)

    # Vertical coupling that is not one on one interior row
    mat = dense_operator(a, b, c, r_vec.size, z_vec.size)
    row = 5 * r_vec.size + 4
    mat[row, row + r_vec.size] = 1.3
    with pytest.raises(ValueError, match="vertical coupling"):
        ps.stencil_from_lu_bands(*lup_bands_as_gsfit_writes_them(mat, r_vec.size), r_vec.size, z_vec.size)

    # Radial coefficient that differs between two rows of the same column
    mat = dense_operator(a, b, c, r_vec.size, z_vec.size)
    mat[row, row + 1] *= 1.01
    with pytest.raises(ValueError, match="differ between rows"):
        ps.stencil_from_lu_bands(*lup_bands_as_gsfit_writes_them(mat, r_vec.size), r_vec.size, z_vec.size)

    # An entry outside the five-point stencil
    mat = dense_operator(a, b, c, r_vec.size, z_vec.size)
    mat[row, row + 2] = 0.5
    with pytest.raises(ValueError, match="outside the five-point stencil"):
        ps.stencil_from_lu_bands(*lup_bands_as_gsfit_writes_them(mat, r_vec.size), r_vec.size, z_vec.size)
