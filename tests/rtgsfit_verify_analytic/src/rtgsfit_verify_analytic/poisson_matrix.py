"""
This module discretises the left-hand side of the Grad-Shafranov equation (the
modified Laplacian applied to the poloidal flux) and provides the stencil
coefficients that RTGSFIT's Poisson solver reads from constants.c.
"""

import numpy as np


def poisson_stencil(
    r_vec: np.ndarray,
    z_vec: np.ndarray
) -> tuple:
    """
    Stencil coefficients of the discretised Grad-Shafranov operator used by RTGSFIT.

    This implementation follows the finite difference scheme described in Equation (45) of:
    Moret, J.-M., et al. "Tokamak equilibrium reconstruction code LIUQE and its real time implementation."
    Fusion Engineering and Design 91 (2015): 1-15.

    The discretized equation has the form:

        ψ_{i+1,j} + ψ_{i-1,j} + a_j * ψ_{i,j+1} + b_j * ψ_{i,j-1} - c_j * ψ_{i,j}

    where:
        - ψ is the poloidal flux function,
        - a_j = (dz / dr)^2 * r_j / (r_j + dr / 2),
        - b_j = (dz / dr)^2 * r_j / (r_j - dr / 2),
        - c_j = 2 + a_j + b_j
    but we take the large major radius limit of the Grad-Shafranov equation, hence we simplify
    a_j and b_j to:
        - a_j = a = (dz / dr)^2
        - b_j = b = (dz / dr)^2

    Grid points on the boundary of the (r, z) box carry identity (Dirichlet) rows; the solver
    (src/poisson_fast.c) treats them that way, so only the interior coefficients matter.

    Parameters:
        r_vec (np.ndarray): 1D array of radial grid coordinates.
        z_vec (np.ndarray): 1D array of vertical grid coordinates.

    Returns:
        tuple: (a, b, c), each a 1D array of length len(r_vec) indexed by grid column.
    """

    n_r = len(r_vec)
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]

    a = np.full(n_r, (dz / dr) ** 2, dtype=np.float64)
    b = np.full(n_r, (dz / dr) ** 2, dtype=np.float64)
    # a = (dz / dr) ** 2 * r_vec / (r_vec + dr / 2)
    # b = (dz / dr) ** 2 * r_vec / (r_vec - dr / 2)
    c = 2 + a + b

    return a, b, c


def poisson_matrix(
    r_vec: np.ndarray,
    z_vec: np.ndarray
) -> np.ndarray:
    """
    Dense (n_grid x n_grid) form of the operator described in poisson_stencil(), for reference
    and testing. Not used by the RTGSFIT build.
    """

    n_r = len(r_vec)
    n_z = len(z_vec)
    n_grid = n_r * n_z
    a, b, c = poisson_stencil(r_vec, z_vec)

    poiss_matrix = np.zeros((n_grid, n_grid), dtype=np.float64)

    idx = -1
    for i in range(n_z):
        for j in range(n_r):
            idx += 1
            if i == 0 or i == n_z - 1 or j == 0 or j == n_r - 1:
                poiss_matrix[idx, idx] = 1  # Dirichlet BC
            else:
                poiss_matrix[idx, idx] = -c[j]  # coefficient of ψ_{i,j}
                poiss_matrix[idx, idx + 1] = a[j]  # coefficient of ψ_{i,j+1}
                poiss_matrix[idx, idx - 1] = b[j]  # coefficient of ψ_{i,j-1}
                poiss_matrix[idx, idx + n_r] = 1  # coefficient of ψ_{i+1,j}
                poiss_matrix[idx, idx - n_r] = 1  # coefficient of ψ_{i-1,j}

    return poiss_matrix


if __name__ == "__main__":
    # Example usage
    r_vec = np.linspace(110.000e-3, 1.070e0, 33)
    z_vec = np.linspace(-960.000e-3, 960.000e-3, 65)
    a, b, c = poisson_stencil(r_vec, z_vec)
