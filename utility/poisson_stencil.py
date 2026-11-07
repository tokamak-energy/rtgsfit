"""
Poisson stencil coefficients for constants.c, recovered from the banded LU
factors that GSFit stores in the RTGSFIT PRESHOT node.

GSFit factorises the LIUQE eq. (45) operator A as A = P L U and stores the
strictly lower band of L (LOWER_BAND, N_GRID x N_R), the upper band of U with
the reciprocal of the diagonal in its first column (UPPER_BAND,
N_GRID x (N_R + 2)) and the permutation (PERM_IDX, with row i of L U equal to
row PERM_IDX[i] of A). Multiplying the factors back together recovers A
exactly up to rounding, so constants.c can carry the three stencil arrays
POISSON_A, POISSON_B and POISSON_C without any change to the PRESHOT node or
to the GSFit writer.

Only numpy is needed, so this runs in the same environment as
mdsplus_const_to_file.py.
"""

import numpy as np


def stencil_from_lu_bands(lower_band, upper_band, perm_idx, n_r, n_z, rtol=1e-9):
    """
    Recover the stencil coefficients of the Poisson operator from its banded LU factors.

    The operator must have identity rows on the grid boundary and, for every
    interior point (i, j) (i: vertical index, j: radial index),

        psi[i+1][j] + psi[i-1][j] + a_j psi[i][j+1] + b_j psi[i][j-1] - c_j psi[i][j]

    with a_j, b_j and c_j depending on the column j only. Every row is checked
    for that structure (entries outside the stencil zero, vertical coupling
    one, the same coefficients on every row of a column) to a tolerance of
    rtol times the largest entry of the row, and a ValueError is raised if a
    check fails, so a constants.c is never written for an operator that the
    separable solver in src/poisson_fast.c cannot represent.

    Parameters:
        lower_band: LOWER_BAND as stored by GSFit, N_GRID * N_R values.
        upper_band: UPPER_BAND as stored by GSFit, N_GRID * (N_R + 2) values.
        perm_idx: PERM_IDX as stored by GSFit, N_GRID values.
        n_r: number of radial grid points (N_R).
        n_z: number of vertical grid points (N_Z).
        rtol: relative tolerance of the structure checks.

    Returns:
        (poisson_a, poisson_b, poisson_c), each a float64 array of length n_r
        indexed by grid column. The boundary columns (j = 0 and j = n_r - 1)
        carry identity rows, are not read by the solver and are set to zero.
    """
    n_grid = n_r * n_z
    lower = np.asarray(lower_band, dtype=np.float64).reshape(n_grid, n_r)
    upper = np.array(upper_band, dtype=np.float64).reshape(n_grid, n_r + 2)
    perm = np.asarray(perm_idx, dtype=np.int64).reshape(n_grid)
    if not np.array_equal(np.sort(perm), np.arange(n_grid)):
        raise ValueError("PERM_IDX is not a permutation of 0 .. N_GRID - 1")
    if np.any(upper[:, 0] == 0.0):
        raise ValueError("UPPER_BAND has a zero reciprocal diagonal entry")
    upper[:, 0] = 1.0 / upper[:, 0]

    # Row i of L U over the columns i - n_r .. i + n_r + 1, which hold all of its
    # possible nonzeros: L[i][i - n_r + k] for k = 0 .. n_r (unit diagonal at
    # k = n_r) times row i - n_r + k of U, whose band starts at its diagonal.
    width = 2 * n_r + 2
    lu_rows = np.zeros((n_grid, width), dtype=np.float64)
    l_rows = np.concatenate([lower, np.ones((n_grid, 1))], axis=1)
    rows = np.arange(n_grid)
    for k in range(n_r + 1):
        u_row = rows - n_r + k
        valid = u_row >= 0
        lu_rows[valid, k:k + n_r + 2] += l_rows[valid, k:k + 1] * upper[u_row[valid], :]

    # Column window of each row of L U, and the operator row it belongs to.
    col = rows[:, None] - n_r + np.arange(width)[None, :]
    a_row = perm
    i_z, j_r = np.divmod(a_row, n_r)
    is_boundary = (i_z == 0) | (i_z == n_z - 1) | (j_r == 0) | (j_r == n_r - 1)
    row_scale = np.max(np.abs(lu_rows), axis=1)
    if np.any(row_scale == 0.0):
        raise ValueError("an operator row recovered from the LU factors is zero")
    tol = rtol * row_scale

    offsets = {"down": -n_r, "left": -1, "centre": 0, "right": 1, "up": n_r}
    entries = {}
    in_stencil = np.zeros_like(lu_rows, dtype=bool)
    for name, offset in offsets.items():
        target = a_row + offset
        hit = col == target[:, None]
        found = hit.any(axis=1)
        entries[name] = np.where(found, lu_rows[rows, np.argmax(hit, axis=1)], 0.0)
        if name == "centre":
            in_stencil |= hit
        else:
            in_stencil |= hit & ~is_boundary[:, None]

    stray = np.where(in_stencil, 0.0, np.abs(lu_rows))
    if np.any(stray.max(axis=1) > tol):
        raise ValueError("the operator has nonzero entries outside the five-point stencil")

    if np.any(np.abs(entries["centre"][is_boundary] - 1.0) > tol[is_boundary]):
        raise ValueError("a boundary row of the operator is not an identity row")

    interior = ~is_boundary
    if np.any(np.abs(entries["down"][interior] - 1.0) > tol[interior]) or np.any(
        np.abs(entries["up"][interior] - 1.0) > tol[interior]
    ):
        raise ValueError("the vertical coupling of the operator is not one on every interior row")

    poisson_a = np.zeros(n_r, dtype=np.float64)
    poisson_b = np.zeros(n_r, dtype=np.float64)
    poisson_c = np.zeros(n_r, dtype=np.float64)
    for j in range(1, n_r - 1):
        sel = interior & (j_r == j)
        a_j = entries["right"][sel]
        b_j = entries["left"][sel]
        c_j = -entries["centre"][sel]
        col_tol = tol[sel]
        # Use the middle row of the column; every other row must agree with it.
        ref = int(np.argsort(i_z[sel])[a_j.size // 2])
        for values in (a_j, b_j, c_j):
            if np.any(np.abs(values - values[ref]) > col_tol):
                raise ValueError(f"the stencil coefficients of column {j} differ between rows")
        poisson_a[j] = a_j[ref]
        poisson_b[j] = b_j[ref]
        poisson_c[j] = c_j[ref]

    return poisson_a, poisson_b, poisson_c
