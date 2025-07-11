"""
This module contains functions to generate Green's matrices for RTGSFIT, namely the
G_GRID_MEAS_WEIGHT, G_COEF_MEAS_WEIGHT, G_LTRB, G_MEAS_COIL, G_GRID_COIL and G_GRID_VESSEL
matrices.

Note that store the Green's matrices on column major order, as this is default in matlab and LAPACK
so it saves us from having to transpose the matrices when we use them in RTGSFIT.
"""
import numpy as np
from rtgsfit_verify_analytic import mutual_inductance


def g_grid_meas_weight(
    r_vec: np.ndarray,
    z_vec: np.ndarray,
    r_flux_loops: np.ndarray,
    z_flux_loops: np.ndarray,
    r_b_probes: np.ndarray,
    z_b_probes: np.ndarray,
    angle_b_probes: np.ndarray,
    n_rog: int,
    n_meas: int,
    weights: np.ndarray) -> np.ndarray:
    """
    Generate the G_GRID_MEAS_WEIGHT matrix for RTGSFIT.

    This is effectively  M_{fy}, B_{my} terms in eqn (61) of the  Moret et al. (2015)
    LIUQE paper  transposed and multiplied by ΔR * Δz and the weights.

    G_GRID_MEAS_WEIGHT has shape (N_GRID, N_MEAS)

    Note that:

        g_pls_grid * G_GRID_MEAS_WEIGHT = G_COEF_MEAS_WEIGHT[:, :N_PLS]

    also note that:

        g_pls_grid has shape (N_PLS, N_GRID) and is effecitely the 
        T_{yg} matrix in eqn. (61) of the Moret et al. (2015) LIUQE paper
        transposed without the ΔR * Δz factor.
    
    in the RTGSFIT code and is also stored in stored in column major order.
    """

    g_grid_meas_weight = np.zeros((len(r_vec) * len(z_vec), n_meas), dtype=np.float64)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    delta_r = r_vec[1] - r_vec[0]
    delta_z = z_vec[1] - z_vec[0]

    for i, (r_g, z_g) in enumerate(zip(r_grid.flatten(), z_grid.flatten())):
        
        j = -1

        # Calculate flux loop terms
        # This is effectively the M_{fy} term in eqn. (61) of the Moret et al. (2015) LIUQE paper
        # multiplied by the weights and ΔR * Δz (we do the ΔR * Δz
        # multiplication at the end, as we need to do it for all terms)
        for r_meas, z_meas in zip(r_flux_loops, z_flux_loops):
            j += 1
            psi = mutual_inductance.mutual_inductance_psi(r_meas, z_meas, r_g, z_g)
            g_grid_meas_weight[i, j] = psi * weights[j]

        # Calculate B probe terms
        # This is effectively the B_{my} term in eqn. (61) of the Moret et al. (2015) LIUQE paper
        # multiplied by the weights and ΔR * Δz
        for r_meas, z_meas, angle_meas in zip(r_b_probes, z_b_probes, angle_b_probes):
            j += 1
            br = mutual_inductance.mutual_inductance_br(r_meas, z_meas, r_g, z_g)
            bz = mutual_inductance.mutual_inductance_bz(r_meas, z_meas, r_g, z_g)
            b_meas = br * np.cos(angle_meas) + bz * np.sin(angle_meas)
            g_grid_meas_weight[i, j] = b_meas * weights[j]

        # Calculate Rogovski coil terms
        # We assume that there are two Rogovski coils
        # The first one goes around the plasma, we set it ΔR * Δz
        # as we need to create the T_{pg} term in the eqn. (61) of the Moret et al. (2015) LIUQE paper
        # when g_pls_grid is multiplied by G_GRID_MEAS_WEIGHT we will create the sum over
        # T_{yg} terms as required.
        # The second Rogovski coil goes just around the vessel and so we set it to 0.
        if n_rog != 2:
            raise ValueError("n_rog must be 2, as we assume there are two Rogovski coils.")
        j += 1
        g_grid_meas_weight[i, j] = weights[j]
        j += 1
        g_grid_meas_weight[i, j] = 0

    return g_grid_meas_weight.flatten() * (delta_r * delta_z)

def g_coef_meas_weight(
    r_vessel: np.ndarray,
    z_vessel: np.ndarray,
    r_flux_loops: np.ndarray,
    z_flux_loops: np.ndarray,
    r_b_probes: np.ndarray,
    z_b_probes: np.ndarray,
    angle_b_probes: np.ndarray,
    n_rog: int,
    n_pls: int,
    n_coef: int,
    n_meas: int,
    weights: np.ndarray) -> np.ndarray:
    """

    Generate the G_COEF_MEAS_WEIGHT matrix for RTGSFIT.

    This is effectively the matrix in eqn. (61) of the Moret et al. 2015 LIUQE paper.

    G_COEF_MEAS_WEIGHT has shape (N_MEAS, N_COEF) (but stored in column major order).
    
    Entries G_COEF_MEAS_WEIGHT[:, :N_PLS] are all 0, but could
    be any values as they are overwritten by RT-GSFIT.

    G_COEF_MEAS_WEIGHT satisfies:

        G_COEF_MEAS_WEIGHT * coef = meas_no_coil ∘ WEIGHT

    (assuming RTGSFIT is able to perfectly match the measurements)
    where coef is the degree of freedom row vector of length N_COEF.

    The above equation is solved by RTGSFIT using LAPACKE_dgelss, this is
    why G_COEF_MEAS_WEIGHT is stored in column major order as LAPACK is a
    Fortran library and uses column major order, so we save computational time
    by not having to transpose the matrix when we use it in RTGSFIT.

    coef[:N_PLS] = [p' coefficent, ff' coefficent, asymmetry coefficent]
    coef[N_PLS:] = [vessel_element 1 current, vessel_element 2 current, ...]

    Vessel here refers to passive elements in the whole tokamak, such as the vacuum vessel.

    """

    g_coef_meas_w = np.zeros(n_coef * n_meas, dtype=np.float64)

    for j, (r_v, z_v) in enumerate(zip(r_vessel, z_vessel)):

        i = -1

        # Calculate flux loop terms
        # This is effectively the M_{fs} term in eqn. (61) of the Moret et al. (2015) LIUQE paper
        # multiplied by the weights
        for r_meas, z_meas in zip(r_flux_loops, z_flux_loops):
            i += 1
            psi = mutual_inductance.mutual_inductance_psi(r_meas, z_meas, r_v, z_v)
            g_coef_meas_w[i + n_meas * (n_pls + j)] = psi * weights[i]

        # Calculate B probe terms
        # This is effectively the B_{ms} term in eqn. (61) of the Moret et al. (2015) LIUQE paper
        # multiplied by the weights
        for r_meas, z_meas, angle_meas in zip(r_b_probes, z_b_probes, angle_b_probes):
            i += 1
            br = mutual_inductance.mutual_inductance_br(r_meas, z_meas, r_v, z_v)
            bz = mutual_inductance.mutual_inductance_bz(r_meas, z_meas, r_v, z_v)
            b_meas = br * np.cos(angle_meas) + bz * np.sin(angle_meas)
            g_coef_meas_w[i + n_meas * (n_pls + j)] = b_meas * weights[i]

        # Calculate rogovski coil terms
        # We assume that there are two Rogovski coils
        # The first one goes around the plasma, hence the vessel currents do not affect it
        # The second one goes around just the vessel.
        # This is effectively the 1_s term in eqn. (61) of the Moret et al. (2015) LIUQE paper
        # multiplied by the weights and by the delta R and delta Z
        # check n_rog == 2
        if n_rog != 2:
            raise ValueError("n_rog must be 2, as we assume there are two Rogovski coils.")
        i += 2
        g_coef_meas_w[i + n_meas * (n_pls + j)] = 1 * weights[i]

    return g_coef_meas_w

def g_ltrb(
    r_ltrb: np.ndarray,
    z_ltrb: np.ndarray) -> np.ndarray:
    """
    Generate the G_LTRB matrix for RTGSFIT.

    This is effectively the discretised version of the M(r, z, r', z') function
    in the second line of eqn. (51) of the Moret et al. (2015) LIUQE paper.

    G_LTRB has shape (N_LTRB, N_LTRB)

    Note that G_LTRB satisfies:

        G_LTRB * dpsi_ltrb * INV_R_LTRB_MU0 = psi_ltrb

    where:
        - dpsi_ltrb is derivate of psi normal to the boundaries.
        - INV_R_LTRB_MU0 is the inverse of the radial coordinate multiplied by mu_0.

    Parameters:
        r_ltrb (np.ndarray): Radial coordinates of the left, top, right, and bottom boundaries.
        z_ltrb (np.ndarray): Axial coordinates of the left, top, right, and bottom boundaries.
        dr (float): Radial grid spacing.
        dz (float): Vertical grid spacing.
    Returns:
        np.ndarray: The G_LTRB matrix, which contains the mutual inductance values
                    between the left, top, right, and bottom boundaries.
    """

    n_ltrb = len(r_ltrb)

    delta = np.zeros(n_ltrb, dtype=np.float64)
    delta[:-1] = np.sqrt((r_ltrb[1:] - r_ltrb[:-1])**2 + (z_ltrb[1:] - z_ltrb[:-1])**2)
    delta[-1] = np.sqrt((r_ltrb[0] - r_ltrb[-1])**2 + (z_ltrb[0] - z_ltrb[-1])**2)

    g_ltrb = np.zeros((n_ltrb, n_ltrb), dtype=np.float64)
    for i, (r_i, z_i, delta) in enumerate(zip(r_ltrb, z_ltrb, delta)):
        for j, (r_j, z_j) in enumerate(zip(r_ltrb, z_ltrb)):
            if i == j:
                g_ltrb[i, i] = mutual_inductance.self_inductance_linear_avg(r_i, z_i, delta)
            else:
                g_ltrb[i, j] = mutual_inductance.mutual_inductance_psi(r_i, z_i, r_j, z_j)

    return g_ltrb.flatten()

def g_grid_vessel(
    r_vessel: np.ndarray,
    z_vessel: np.ndarray,
    r_vec: np.ndarray,
    z_vec: np.ndarray) -> np.ndarray:
    """
    Generate the G_GRID_VESSEL matrix for RTGSFIT.

    It has shape (N_GRID, N_VESSEL).

    The matrix satisfies:

        G_GRID_VESSEL * coef[N_PLS:N_PLS+N_NESSEL] = flux_vessel

    where:

        - coef[N_PLS:N_PLS+N_NESSEL] is the current for each of the vessel elements
        - flux_vessel is ψ on the grid generated by the vessel currents.

    Parameters:
        r_vessel (np.ndarray): Radial coordinates of the vessel elements.
        z_vessel (np.ndarray): Vertical coordinates of the vessel elements.
        r_vec (np.ndarray): Radial coordinates of the grid points.
        z_vec (np.ndarray): Vertical coordinates of the grid points.

    Returns:

        np.ndarray: The G_GRID_VESSEL matrix, which contains the mutual inductance values
                    between the grid points and the vessel elements.
    """

    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    n_grid = len(r_vec) * len(z_vec)
    n_vessel = len(r_vessel)
    g_grid_vessel = np.zeros((n_grid, n_vessel), dtype=np.float64)

    for i, (r_g, z_g) in enumerate(zip(r_grid.flatten(), z_grid.flatten())):
        for j, (r_v, z_v) in enumerate(zip(r_vessel, z_vessel)):
            g_grid_vessel[i, j] = mutual_inductance.mutual_inductance_psi(r_g, z_g, r_v, z_v)

    return g_grid_vessel.flatten()
