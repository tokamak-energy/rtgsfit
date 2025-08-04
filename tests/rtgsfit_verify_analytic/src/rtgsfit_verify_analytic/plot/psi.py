"""
Python module for plotting comapring the numerical and analytic values of ψ in Wb.

This module assumes the output_dict.npy file created by the replay_rtgsfit.py script
has already been created and contains the necessary data for plotting.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from rtgsfit_verify_analytic import analytic_soln, cnst, read_constants_c

def contour_ana_vs_num(output_dict):
    """
    Plot the contour of the analytic ψ and the numerical ψ from the output_dict.
    
    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """

    constants_c_path = os.path.join(cnst.DATA_DIR, 'constants.c')
    constants_c_dict = read_constants_c.constants_c_dict(constants_c_path,
                                                         read_g_grid_meas_weight=True)

    r_vec = constants_c_dict['r_vec']
    z_vec = constants_c_dict['z_vec']
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    r_flat = r_grid.flatten()
    z_flat = z_grid.flatten()
    psi_ana = analytic_soln.analytic_psi(
        r_grid, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB, 
        cnst.ANALYTIC_PLASMA_CURRENT
    )

    # Get the limiter points
    n_lim = constants_c_dict['n_limit']
    lim_idx = constants_c_dict['limit_idx']
    lim_weight = constants_c_dict['limit_weight'].reshape(n_lim, cnst.N_INTRP)
    lim_r = np.sum(r_flat[lim_idx].reshape(n_lim, cnst.N_INTRP) * lim_weight, axis=1)
    lim_z = np.sum(z_flat[lim_idx].reshape(n_lim, cnst.N_INTRP) * lim_weight, axis=1)

    # Get Green's matrices
    g_grid_meas_weight = constants_c_dict['g_grid_meas_weight'].reshape(cnst.N_GRID, cnst.N_MEAS)
    # Extract just the first N_FLUX_LOOPS columns for g_coef_meas_weight
    g_grid_meas_weight = g_grid_meas_weight[:, :cnst.N_FLUX_LOOPS]

    # Get the flux loop coordinates
    flux_loop_coords = np.loadtxt(os.path.join(cnst.DATA_DIR, 'flux_loops.txt'),
                                  dtype=np.float64, skiprows=1)
    fl_r_coords = flux_loop_coords[:, 0]
    fl_z_coords = flux_loop_coords[:, 1]

    for i_iter in range(cnst.N_ITERS + 1):
        psi_num = output_dict['flux_total'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm = output_dict['flux_norm'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm_flat = flux_norm.flatten()
        lcfs_r = output_dict['lcfs_r'][i_iter]
        lcfs_z = output_dict['lcfs_z'][i_iter]
        lcfs_n = output_dict['lcfs_n'][i_iter][0]
        coef = output_dict['coef'][i_iter]
        # Extract just the first N_PLS coefficients
        coef = coef[:cnst.N_PLS]
        coef = coef.reshape(cnst.N_PLS, 1)
        meas = output_dict['meas'][i_iter]
        meas = meas.reshape(cnst.N_MEAS, 1)

        # Calculate g_pls_grid
        g_pls_grid = np.zeros((cnst.N_PLS, cnst.N_GRID), dtype=np.float64)
        g_pls_grid[0, :] = (1 - flux_norm_flat) * r_flat
        g_pls_grid[1, :] = (1 - flux_norm_flat)  / (r_flat * cnst.MU_0)
        # Take derivative of flux_norm in the z-direction (along columns)
        dflux_norm_dz = np.gradient(flux_norm, z_vec, axis=0)
        g_pls_grid[2, :] = dflux_norm_dz.flatten()

        g_coef_meas_weight = g_pls_grid @ g_grid_meas_weight
        g_coef_meas_weight = g_coef_meas_weight.T

        meas_fit = g_coef_meas_weight @ coef

        flux_loop_meas_fit = meas_fit[:cnst.N_FLUX_LOOPS]
        flux_loop_meas = meas[:cnst.N_FLUX_LOOPS]

        # Create a figure with two subplots (vertical layout)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 14), gridspec_kw={'height_ratios': [2, 1]})

        # --- Top subplot: Contour plot ---
        ax1.contour(r_grid, z_grid, psi_ana, levels=20, colors='tab:blue')
        ax1.contour(r_grid, z_grid, psi_num, levels=20, colors='tab:orange', linestyles='dashed')

        # Limiter and LCFS points
        ax1.scatter(lim_r, lim_z, color='tab:red', label='Limiter Points')
        ax1.scatter(lcfs_r[:lcfs_n], lcfs_z[:lcfs_n], color='tab:orange', label='LCFS Points')

        # Plot the flux loop coordinates
        ax1.scatter(fl_r_coords, fl_z_coords, color='tab:green', label='Flux Loop Points', marker='o')
        # Add integer label to each flux loop point
        for idx, (r, z) in enumerate(zip(fl_r_coords, fl_z_coords)):
            ax1.text(r, z, str(idx), fontsize=12, color='black',
                     ha='center', va='center', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            
        # Plot the vessel filament coordiantes
        ax1.scatter(cnst.R_VESSEL, cnst.Z_VESSEL, color='black', label='Vessel Filament Points', marker='o')

        # Legend using proxy artists
        contour_proxy = plt.Line2D([0], [0], color='tab:blue', label='Analytic ψ')
        num_proxy = plt.Line2D([0], [0], color='tab:orange', linestyle='dashed', label='Numerical ψ')
        lim_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', label='Limiter Points', linestyle='None')
        lcfs_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:orange', label='LCFS Points', linestyle='None')
        floop_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:green', label='Flux Loop Points', linestyle='None')
        vessel_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Vessel Filament Points', linestyle='None')
        ax1.legend(handles=[contour_proxy, num_proxy, lim_proxy, lcfs_proxy, floop_proxy, vessel_proxy],
                   loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)

        ax1.set_title(f'ψ in Wb (Iteration {i_iter:02d})')
        ax1.set_xlabel('R [m]')
        ax1.set_ylabel('z [m]')
        ax1.grid(True)
        ax1.set_aspect('equal')

        # --- Bottom subplot: Flux loop measurements ---
        ax2.plot(flux_loop_meas, label='Analytic Flux Loop Values', marker='o', linestyle='None', color='tab:blue')
        ax2.plot(flux_loop_meas_fit, label='RTGSFIT Flux Loop values', marker='x', linestyle='None', color='tab:orange')
        # Make x labels all integers
        ax2.set_xticks(np.arange(cnst.N_FLUX_LOOPS))
        # ax2.set_title('Flux Loop Analytic Values vs RTGSFIT Values')
        ax2.set_xlabel('Flux Loop Index')
        ax2.set_ylabel('Flux Loop Measurement (Wb)')
        ax2.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), borderaxespad=0)
        ax2.grid(True)

        # Save the combined plot
        plot_dir = "combined_psi_and_fluxloop"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'combined_psi_fluxloop_iter_{i_iter:02d}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def midplane_ana_vs_num(output_dict):
    """
    Plot the midplane ψ values for both analytic and numerical solutions.
    
    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """
    
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    z_mid = z_vec[cnst.N_Z // 2]  # Midplane z-coordinate
    psi_ana = analytic_soln.analytic_psi(r_vec, z_mid,
                                         cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
                                         cnst.ANALYTIC_RHOB, 
                                         cnst.ANALYTIC_PLASMA_CURRENT)

    for i_iter in range(cnst.N_ITERS + 1):
        psi_num = output_dict['flux_total'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        psi_num_mid = psi_num[cnst.N_Z // 2, :]

        plt.figure(figsize=(10, 6))
        plt.plot(r_vec, psi_ana, label='Analytic ψ', color='tab:blue')
        plt.plot(r_vec, psi_num_mid, label='Numerical ψ', color='tab:orange', linestyle='dashed')
        # Add vertical lines a cnts.R_LIM_MIN and cnts.R_LIM_MAX
        plt.axvline(cnst.R_LIM_MIN, color='red', linestyle='--', label='R_LIM_MIN')
        plt.axvline(cnst.R_LIM_MAX, color='red', linestyle='--', label='R_LIM_MAX')
        plt.xlabel('Radial Coordinate (R) [m]')
        plt.ylabel('ψ (Wb)')
        plt.title(f'Midplane ψ Comparison (Iteration {i_iter:02d})')
        plt.legend()
        plt.grid(True)
        plot_dir = "midplane_psi_plots"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'midplane_psi_iter_{i_iter:02d}.png'
        plt.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name), dpi=300)
        plt.close()

def boundary_ana_vs_num(output_dict):
    """
    Plot the boundary ψ values for both analytic and numerical solutions.
    
    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """

    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    n_r = len(r_vec)
    n_z = len(z_vec)
    n_ltrb = 2 * n_r + 2 * n_z - 4  # Total number of points on the boundary (left, top, right, bottom)
    
    # Check value of psi at the boundary.
    # Should be about equal to mutual_inductance_psi(r_ltrb, z_ltrb, analytic_ro, analytic_zo)
    r_ltrb = np.concatenate(
        (
            [r_vec[0]],  # (bottom, left)
            r_vec[0] * np.ones(n_z - 2),  # traverse (bottom, left) to (top, left) (excluding corners)
            [r_vec[0]],  # (top, left)
            r_vec[1:-1],  # traverse (top, left) to (top, right) (excluding corners)
            [r_vec[-1]],  # (top, right)
            r_vec[-1] * np.ones(n_z - 2),  # traverse (top, right) to (bottom, right) (excluding corners)
            [r_vec[-1]],  # (bottom, right)
            np.flip(r_vec[1:-1]),  # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    z_ltrb = np.concatenate(
        (
            [z_vec[0]],               # (bottom, left)
            z_vec[1:-1],              # traverse (bottom, left) to (top, left) (excluding corners)
            [z_vec[-1]],              # (top, left)
            z_vec[-1] * np.ones(n_r - 2), # traverse (top, left) to (top, right) (excluding corners)
            [z_vec[-1]],              # (top, right)
            np.flip(z_vec[1:-1]),     # traverse (top, right) to (bottom, right) (excluding corners)
            [z_vec[0]],               # (bottom, right)
            z_vec[0] * np.ones(n_r - 2),  # traverse (bottom, right) to (bottom, left) (excluding corners)
        )
    )
    psi_ltrb_ana = np.zeros(n_ltrb)
    for i, (r, z) in enumerate(zip(r_ltrb, z_ltrb)):
        psi_ltrb_ana[i] = analytic_soln.analytic_psi(
            r, z, cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
            cnst.ANALYTIC_RHOB, cnst.ANALYTIC_PLASMA_CURRENT
        )

    # # Calculate psi at the boundary
    # psi_numerical = output_dict['flux_total'][-1].reshape(n_z, n_r)
    # psi_num_ltrb = np.zeros_like(r_ltrb)
    # psi_num_ltrb[:n_z] = psi_numerical[:, 0] # left side
    # psi_num_ltrb[n_z:n_z+n_r-1] = psi_numerical[-1,1:] # top side
    # psi_num_ltrb[n_z+n_r-1:2*n_z+n_r-2] = np.flip(psi_numerical[:-1, -1]) # right side
    # psi_num_ltrb[2*n_z+n_r-2:] = np.flip(psi_numerical[0, 1:-1]) # bottom side

    # psi_expected = mutual_inductance_psi(r_ltrb, z_ltrb, ANALYTIC_RO, ANALYTIC_ZO) * ANALYTIC_PLASMA_CURRENT

    # fig, ax = plt.subplots(figsize=(8, 4))
    # ax.plot(psi_num_ltrb, label='Numerical Psi at Boundary')
    # ax.plot(psi_expected, label='Expected Psi at Boundary', linestyle='--')
    # ax.set_title('Psi at Boundary')
    # ax.set_xlabel('Index')
    # ax.set_ylabel('Psi')
    # ax.legend()

    for i_iter in range(cnst.N_ITERS + 1):
        psi_num = output_dict['flux_total'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        psi_num_ltrb = np.zeros(n_ltrb)
        psi_num_ltrb[:n_z] = psi_num[:, 0]
        psi_num_ltrb[n_z:n_z+n_r-1] = psi_num[-1, 1:]
        psi_num_ltrb[n_z+n_r-1:2*n_z+n_r-2] = np.flip(psi_num[:-1, -1])
        psi_num_ltrb[2*n_z+n_r-2:] = np.flip(psi_num[0, 1:-1])

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(psi_ltrb_ana, label='Analytic ψ at Boundary', color='tab:blue')
        ax.plot(psi_num_ltrb, label='Numerical ψ at Boundary', color='tab:orange', linestyle='dashed')
        ax.legend()
        ax.set_title(f'Boundary ψ Comparison (Iteration {i_iter:02d})')
        ax.set_xlabel('Boundary Point Index')
        ax.set_ylabel('ψ (Wb)')
        ax.grid(True)
        plot_dir = "boundary_psi_plots"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'boundary_psi_iter_{i_iter:02d}.png'
        plt.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name), dpi=300)
        plt.close()

if __name__ == "__main__":
    # Load the output dictionary
    output_dict = np.load(os.path.join(cnst.DATA_DIR, 'output_dict.npy'), allow_pickle=True).item()

    # Call the contour plotting function

    boundary_ana_vs_num(output_dict)
    midplane_ana_vs_num(output_dict)
    contour_ana_vs_num(output_dict)