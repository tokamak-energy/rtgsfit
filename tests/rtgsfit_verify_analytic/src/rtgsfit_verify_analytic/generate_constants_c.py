"""
This module genreates the constants.c file used by RTGSFIT.
It also outputs the flux loop coordinates and the BP probe coordinates.
"""

import os
import re

import numpy as np
import matplotlib.pyplot as plt

from rtgsfit_verify_analytic import greens_matrices, poisson_matrix, cnst

def compute_limit_idx_and_weights(r, z, lim_r, lim_z, n_intrp, n_lim):
    """
    Computes the indices and weights for interpolation at specified limiter points.
    Parameters:
    - r: 1D array of radial grid points.
    - z: 1D array of vertical grid points.
    - lim_r: 1D array of radial limiter points.
    - lim_z: 1D array of vertical limiter points.
    - n_intrp: Number of interpolation points per limiter point.
    - n_lim: Number of limiter points.
    Returns:
    - limit_idx: 1D array of indices for the interpolation points.
    - limit_w: 1D array of weights for the interpolation points.
    """

    n_r = len(r)
    n_lim = len(lim_r)
    limit_idx = np.zeros(n_lim * n_intrp, dtype=int)
    limit_w = np.zeros(n_lim * n_intrp, dtype=float)

    for i, (lr, lz) in enumerate(zip(lim_r, lim_z)):
        if lr < r[0] or lr > r[-1] or lz < z[0] or lz > z[-1]:
            raise ValueError(f"Limiter point ({lr}, {lz}) is out of bounds of the grid.")
        r_idx = np.searchsorted(r, lr) - 1
        z_idx = np.searchsorted(z, lz) - 1

        limit_idx[n_intrp * i + 0] = n_r * z_idx + r_idx            # (r_idx, z_idx)
        limit_idx[n_intrp * i + 1] = n_r * z_idx + r_idx + 1        # (r_idx + 1, z_idx)
        limit_idx[n_intrp * i + 2] = n_r * (z_idx + 1) + r_idx      # (r_idx, z_idx + 1)
        limit_idx[n_intrp * i + 3] = n_r * (z_idx + 1) + r_idx + 1  # (r_idx + 1, z_idx + 1)

        r0, r1 = r[r_idx], r[r_idx + 1]
        z0, z1 = z[z_idx], z[z_idx + 1]
        dr = r1 - r0
        dz = z1 - z0
        limit_w[n_intrp * i + 0] = (r1 - lr) * (z1 - lz) / (dr * dz)  # (r_idx, z_idx)
        limit_w[n_intrp * i + 1] = (lr - r0) * (z1 - lz) / (dr * dz)  # (r_idx + 1, z_idx)
        limit_w[n_intrp * i + 2] = (r1 - lr) * (lz - z0) / (dr * dz)  # (r_idx, z_idx + 1)
        limit_w[n_intrp * i + 3] = (lr - r0) * (lz - z0) / (dr * dz)  # (r_idx + 1, z_idx + 1)

    return limit_idx, limit_w

def generate_data_dictionary():
    """
    Generates a dictionary called `data_dictionary` that contains
    the constant name and values that will be used by the 
    `const_to_file` function that will generate the contants.c file
    used by RTGSFIT.

    This function also generates the flux loop and BP probe
    coordinates and saves them to files in the same directory as this script.
    The flux loop coordinates are saved in `flux_loops.txt` and the BP probe
    coordinates are saved in `bp_probes.txt`. These are later used by the measurements
    module to generate the measurements for RTGSFIT at these coordinates.

    Finally, it generates a plot of the grid, limiter coordinates, flux loop coordinates,
    and BP probe coordinates and saves it in the `plots` directory.
    """

    # Define the grid related quantities
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    r_flat = r_grid.flatten()
    z_flat = z_grid.flatten()
    r_ltrb = np.concatenate(
        (
            [r_vec[0]],  # (bottom, left)
            r_vec[0] * np.ones(cnst.N_Z - 2),  # traverse (bottom, left) to (top, left) (excluding corners)
            [r_vec[0]],  # (top, left)
            r_vec[1:-1],  # traverse (top, left) to (top, right) (excluding corners)
            [r_vec[-1]],  # (top, right)
            r_vec[-1] * np.ones(cnst.N_Z - 2),  # traverse (top, right) to (bottom, right) (excluding corners)
            [r_vec[-1]],  # (bottom, right)
            np.flip(r_vec[1:-1]),  # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    z_ltrb = np.concatenate(
        (
            [z_vec[0]],               # (bottom, left)
            z_vec[1:-1],              # traverse (bottom, left) to (top, left) (excluding corners)
            [z_vec[-1]],              # (top, left)
            z_vec[-1] * np.ones(cnst.N_R - 2), # traverse (top, left) to (top, right) (excluding corners)
            [z_vec[-1]],              # (top, right)
            np.flip(z_vec[1:-1]),     # traverse (top, right) to (bottom, right) (excluding corners)
            [z_vec[0]],               # (bottom, right)
            z_vec[0] * np.ones(cnst.N_R - 2),  # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    n_ltrb = 2 * cnst.N_R + 2 * cnst.N_Z - 4
    if len(r_ltrb) != n_ltrb:
        raise ValueError("r_ltrb must have length 2 * n_r + 2 * n_z - 4.")
    if len(z_ltrb) != n_ltrb:
        raise ValueError("z_ltrb must have length 2 * n_r + 2 * n_z - 4.")
    inv_r_ltrb_mu0 = 1 / (r_ltrb * cnst.MU_0)
    inv_r_mu0 = 1 / (r_flat * cnst.MU_0)

    # Calculate limiter quantities
    # Define the limiter values as a set of points along a rectangle with
    # corners at (r_lim_min, z_lim_min), (r_lim_min, z_lim_max),
    # (r_lim_max, z_lim_max), and (r_lim_max, z_lim_min)
    r_lim_unique = np.linspace(cnst.R_LIM_MIN, cnst.R_LIM_MAX, cnst.N_R_LIM)
    z_lim_unique = np.linspace(cnst.Z_LIM_MIN, cnst.Z_LIM_MAX, cnst.N_Z_LIM)
    r_lim = np.concatenate(
        (
            [r_lim_unique[0]],                   # (bottom, left)
            r_lim_unique[0] * np.ones(cnst.N_Z_LIM - 2),  # traverse (bottom, left) to (top, left) (excluding corners)
            [r_lim_unique[0]],                   # (top, left)
            r_lim_unique[1:-1],                  # traverse (top, left) to (top, right) (excluding corners)
            [r_lim_unique[-1]],                  # (top, right)
            r_lim_unique[-1] * np.ones(cnst.N_Z_LIM - 2), # traverse (top, right) to (bottom, right) (excluding corners)
            [r_lim_unique[-1]],                             # (bottom, right)
            np.flip(r_lim_unique[1:-1]),         # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    z_lim = np.concatenate(
        (
            [z_lim_unique[0]],                   # (bottom, left)
            z_lim_unique[1:-1],                  # traverse (bottom, left) to (top, left) (excluding corners)
            [z_lim_unique[-1]],                  # (top, left)
            z_lim_unique[-1] * np.ones(cnst.N_R_LIM - 2), # traverse (top, left) to (top, right) (excluding corners)
            [z_lim_unique[-1]],                  # (top, right)
            np.flip(z_lim_unique[1:-1]),         # traverse (top, right) to (bottom, right) (excluding corners)
            [z_lim_unique[0]],                   # (bottom, right)
            z_lim_unique[0] * np.ones(cnst.N_R_LIM - 2),  # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    n_limit = 2 * cnst.N_R_LIM + 2 * cnst.N_Z_LIM - 4
    if len(r_lim) != n_limit:
        print(f"len(r_lim): {len(r_lim)}, n_limit: {n_limit}")
        raise ValueError("r_lim must have length n_limit.")
    if len(z_lim) != n_limit:
        raise ValueError("z_lim must have length n_limit.")
    limit_idx, limit_weight = compute_limit_idx_and_weights(
        r_vec, z_vec, r_lim, z_lim, cnst.N_INTRP, n_limit
    )
    mask_lim_2d = (r_grid >= cnst.R_LIM_MIN - cnst.D_R) & (r_grid <= cnst.R_LIM_MAX + cnst.D_R) & \
                  (z_grid >= cnst.Z_LIM_MIN - cnst.D_Z) & (z_grid <= cnst.Z_LIM_MAX + cnst.D_Z)
    mask_lim = mask_lim_2d.flatten()
    
    # Calculate poisson matrix quantities
    lower_band, upper_band, perm_idx = poisson_matrix.compute_lup_bands(r_vec, z_vec)

    # Generate the flux loop locations
    r_fl_vec = np.linspace(cnst.R_FL_MIN, cnst.R_FL_MAX, cnst.N_FL_R)
    z_fl_vec = np.linspace(cnst.Z_FL_MIN, cnst.Z_FL_MAX, cnst.N_FL_Z)
    r_fl = np.concatenate(
        (
            [cnst.R_FL_MIN],                     # (bottom, left)
            cnst.R_FL_MIN * np.ones(cnst.N_FL_Z - 2), # traverse (bottom, left) to (top, left) (excluding corners)
            [cnst.R_FL_MIN],                     # (top, left)
            r_fl_vec[1:-1],                 # traverse (top, left) to (top, right) (excluding corners)
            [cnst.R_FL_MAX],                     # (top, right)
            cnst.R_FL_MAX * np.ones(cnst.N_FL_Z - 2), # traverse (top, right) to (bottom, right) (excluding corners)
            [cnst.R_FL_MAX],                     # (bottom, right)
            np.flip(r_fl_vec[1:-1]),        # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    z_fl = np.concatenate(
        (
            [cnst.Z_FL_MIN],                     # (bottom, left)
            z_fl_vec[1:-1],                 # traverse (bottom, left) to (top, left) (excluding corners)
            [cnst.Z_FL_MAX],                     # (top, left)
            cnst.Z_FL_MAX * np.ones(cnst.N_FL_R - 2), # traverse (top, left) to (top, right) (excluding corners)
            [cnst.Z_FL_MAX],                     # (top, right)
            np.flip(z_fl_vec[1:-1]),        # traverse (top, right) to (bottom, right) (excluding corners)
            [cnst.Z_FL_MIN],                     # (bottom, right)
            cnst.Z_FL_MIN * np.ones(cnst.N_FL_R - 2), # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    n_flux_loops = 2 * cnst.N_FL_R + 2 * cnst.N_FL_Z - 4
    if len(r_fl) != n_flux_loops:
        raise ValueError("r_fl must have length 2 * N_FL_R + 2 * N_FL_Z - 4.")
    if len(z_fl) != n_flux_loops:
        raise ValueError("z_fl must have length 2 * N_FL_R + 2 * N_FL_Z - 4.")
    # Save flux loop r and z coordinates to a file
    flux_loop_file = os.path.join(cnst.DATA_DIR, "flux_loops.txt")
    with open(flux_loop_file, 'w') as f:
        f.write("r_fl\tz_fl\n")
        for r, z in zip(r_fl, z_fl):
            f.write(f"{r:.6f}\t{z:.6f}\n")

    # Generate BP probe locations
    r_bp_vec = np.linspace(cnst.R_BP_MIN, cnst.R_BP_MAX, cnst.N_BP_R)
    z_bp_vec = np.linspace(cnst.Z_BP_MIN, cnst.Z_BP_MAX, cnst.N_BP_Z)
    r_bp_vec = r_bp_vec[1:-1]  # Exclude corners for BP probe locations
    z_bp_vec = z_bp_vec[1:-1]  # Exclude corners for BP probe locations
    r_bp = np.concatenate(
        (
            cnst.R_BP_MIN * np.ones(cnst.N_BP_Z - 2), # (bottom, left) to (top, left)
            r_bp_vec,                       # (top, left) to (top, right)
            cnst.R_BP_MAX * np.ones(cnst.N_BP_Z - 2), # (top, right) to (bottom, right)
            np.flip(r_bp_vec),              # (bottom, right) to (bottom, left)       
        )
    )
    z_bp = np.concatenate(
        (
            z_bp_vec,                       # (bottom, left) to (top, left)
            cnst.Z_BP_MAX * np.ones(cnst.N_BP_R - 2), # (top, left) to (top, right)
            np.flip(z_bp_vec),              # (top, right) to (bottom, right)
            cnst.Z_BP_MIN * np.ones(cnst.N_BP_R - 2), # (bottom, right) to (bottom, left)
        )
    )
    angle_bp = np.concatenate(
        (
            0.5 * np.pi * np.ones(cnst.N_BP_Z -2),   # Point up on the left side
            np.zeros(cnst.N_BP_R - 2),               # Point right on the top side
            -0.5 * np.pi * np.ones(cnst.N_BP_Z - 2), # Point down on the right side
            np.pi * np.ones(cnst.N_BP_R - 2),        # Point left on the bottom side
        )
    )
    # Save BP probe r, z, and angle coordinates to a file
    n_bp_probes = 2 * (cnst.N_BP_R - 2) + 2 * (cnst.N_BP_Z - 2)  # We don't have corner BP probes
    if len(r_bp) != n_bp_probes:
        raise ValueError("r_bp must have length 2 * (N_BP_R - 2) + 2 * (N_BP_Z - 2).")
    if len(z_bp) != n_bp_probes:
        raise ValueError("z_bp must have length 2 * (N_BP_R - 2) + 2 * (N_BP_Z - 2).")
    bp_probe_file = os.path.join(cnst.DATA_DIR, "bp_probes.txt")
    with open(bp_probe_file, 'w') as f:
        f.write("r_bp\tz_bp\tangle_bp\n")
        for r, z, angle in zip(r_bp, z_bp, angle_bp):
            f.write(f"{r:.6f}\t{z:.6f}\t{angle:.6f}\n")

    # Generate the weights for the measurements
    weights = np.ones(cnst.N_MEAS, dtype=np.float64)

    g_grid_meas_weight = greens_matrices.g_grid_meas_weight(
        r_vec, z_vec,
        r_fl, z_fl,
        r_bp, z_bp, angle_bp,
        cnst.N_ROGOWSKI_COILS,
        cnst.N_MEAS,
        weights
    )

    g_coef_meas_weight = greens_matrices.g_coef_meas_weight(
        cnst.R_VESSEL, cnst.Z_VESSEL,
        r_fl, z_fl,
        r_bp, z_bp, angle_bp,
        cnst.N_ROGOWSKI_COILS,
        cnst.N_PLS,
        cnst.N_COEF,
        cnst.N_MEAS,
        weights)
    
    g_ltrb = greens_matrices.g_ltrb(r_ltrb, z_ltrb)

    # We have one PF coil but we set its current to zero
    n_pf_coils = 1
    g_meas_coil = np.zeros((cnst.N_MEAS, n_pf_coils), dtype=np.float64)
    g_grid_coil = np.zeros((cnst.N_GRID, n_pf_coils), dtype=np.float64)

    g_grid_vessel = greens_matrices.g_grid_vessel(
        cnst.R_VESSEL, cnst.Z_VESSEL,
        r_vec, z_vec)

    fig, ax = plt.subplots()
    ax.plot(r_ltrb, z_ltrb, label="Grid Boundary")
    ax.plot(r_lim, z_lim, label="Limiter Boundary")
    ax.scatter(r_flat, z_flat, s=0.5, label="Grid Points", color='black')
    ax.scatter(r_lim, z_lim, s=1, label="Limiter Points", color='red')
    ax.scatter(r_flat[mask_lim], z_flat[mask_lim], s=0.5, label="MASK_LIM True Grid Points", color='blue')
    ax.scatter(r_bp, z_bp, s=1, label="BP Probe Points", color='green')
    ax.scatter(r_fl, z_fl, s=1, label="Flux Loop Points", color='orange')
    ax.scatter(cnst.R_VESSEL, cnst.Z_VESSEL, s=1, label="Vessel Points", color='purple')
    # Add small lines for bp_probes with angles
    delta_bp_line = 0.05
    for r, z, angle in zip(r_bp, z_bp, angle_bp):
        dx = delta_bp_line * np.cos(angle)
        dy = delta_bp_line * np.sin(angle)
        ax.plot([r, r + dx], [z, z + dy], color='green', linewidth=0.5)
    ax.set_xlabel("Radial Position (R)")
    ax.set_ylabel("Vertical Position (z)")
    ax.set_title("Grid and Limiter Points")
    # Legend box outside the plot
    ax.legend(loc='upper right', bbox_to_anchor=(2.5, 1))
    ax.set_aspect('equal')
    fig.savefig(os.path.join(cnst.PLOTS_DIR, "grid_and_limiter_points.png"),
                bbox_inches='tight', dpi=1000)
    
    data_dictionary = {}
    data_dictionary["dr"] = cnst.D_R
    data_dictionary["dz"] = cnst.D_Z
    data_dictionary["frac"] = cnst.FRAC 
    data_dictionary["inv_r_ltrb_mu0"] = inv_r_ltrb_mu0
    data_dictionary["inv_r_mu0"] = inv_r_mu0
    data_dictionary["limit_idx"] = limit_idx
    data_dictionary["limit_weight"] = limit_weight
    data_dictionary["lower_band"] = lower_band
    # convert mask lim from bool to int
    data_dictionary["mask_lim"] = mask_lim.astype(int)
    data_dictionary["n_bp_probes"] = n_bp_probes
    data_dictionary["n_coef"] = cnst.N_COEF
    data_dictionary["n_coil"] = n_pf_coils
    data_dictionary["n_flux_loops"] = n_flux_loops
    data_dictionary["n_grid"] = cnst.N_GRID
    data_dictionary["n_intrp"] = cnst.N_INTRP
    data_dictionary["n_lcfs_max"] = cnst.N_LCFS_MAX
    data_dictionary["n_limit"] = n_limit
    data_dictionary["n_ltrb"] = n_ltrb
    data_dictionary["n_meas"] = cnst.N_MEAS
    data_dictionary["n_pls"] = cnst.N_PLS
    data_dictionary["n_r"] = cnst.N_R
    data_dictionary["n_rogowski_coils"] = cnst.N_ROGOWSKI_COILS
    data_dictionary["n_vess"] = cnst.N_VESSEL
    data_dictionary["n_xpt_max"] = cnst.N_XPT_MAX
    data_dictionary["n_z"] = cnst.N_Z
    data_dictionary["perm_idx"] = perm_idx
    data_dictionary["r_grid"] = r_grid
    data_dictionary["r_mu0_dz2"] = r_flat * cnst.MU_0 * cnst.D_Z**2
    data_dictionary["r_vec"] = r_vec
    data_dictionary["thresh"] = cnst.THRESH
    data_dictionary["upper_band"] = upper_band
    data_dictionary["weight"] = weights
    data_dictionary["z_grid"] = z_grid
    data_dictionary["z_vec"] = z_vec
    data_dictionary["g_coef_meas_weight"] = g_coef_meas_weight
    data_dictionary["g_grid_coil"] = g_grid_coil
    data_dictionary["g_grid_meas_weight"] = g_grid_meas_weight
    data_dictionary["g_grid_vessel"] = g_grid_vessel
    data_dictionary["g_ltrb"] = g_ltrb
    data_dictionary["g_meas_coil"] = g_meas_coil

    return data_dictionary

def const_to_file(const_dict,
                  file_read=os.path.join(cnst.DATA_DIR, 'constants.h'),
                  file_write=os.path.join(cnst.DATA_DIR, 'constants.c')):

    with open(file_read, 'r') as file:
        file_contents = file.read()

    # Define a function to perform variable replacements
    def replace_variables(match):

        assert((match.lastindex == 4) or (match.lastindex == 5))
        
        var_name = match.group(4).lower()
        var_name_in_file = " ".join([match.group(xx) for xx in range(2, 5)])
        if match.lastindex == 5:
            var_name_in_file = var_name_in_file + match.group(match.lastindex)
        
        if var_name in const_dict:
            replacement = const_dict[var_name]
            # Check the data type of the replacement value
            if isinstance(replacement, bool):
                return f'{var_name_in_file} = {int(replacement)};'
            elif isinstance(replacement, (int, float)):
                return f'{var_name_in_file} = {replacement};'
            elif isinstance(replacement, np.ndarray) and (match.lastindex == 5):
                data_list = map(str, replacement.flatten().tolist())
                data_list = [x if ((ii+1) % 20) else x+'\n' for ii, x in enumerate(data_list)]
                array_values = ', '.join(data_list)
                return f'{var_name_in_file} = {{{array_values}}};'
            else:
                raise Exception(f'{type(replacement)} is not catered for')
        else:
            raise Exception(f"{var_name} is not in dictionary")
            
        return match.group(0)

    pattern = r'\b(\w+)\s*\b(\w+)\s*\b(\w+)\s*\b(\w+)(\[\s*\d*\s*\])?\s*;'
    
    # Replace variables in the C code
    new_c_code = re.sub(pattern, replace_variables, file_contents)

    with open(file_write, 'w') as file:
        file.write(new_c_code)

if __name__ == "__main__":
    data_dictionary = generate_data_dictionary()
    const_to_file(data_dictionary)