"""
Python module for plotting comapring the numerical and analytic values of B_θ in T.

This module assumes the output_dict.npy file created by the replay_rtgsfit.py script
has already been created and contains the necessary data for plotting.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

from rtgsfit_verify_analytic import analytic_soln, cnst, read_constants_c

def contour_ana_vs_num(output_dict):
    """
    Plot the contour of the analytic B_θ and the numerical B_θ from the output_dict.
    
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

    # Get the limiter points
    n_lim = constants_c_dict['n_limit']
    lim_idx = constants_c_dict['limit_idx']
    lim_weight = constants_c_dict['limit_weight'].reshape(n_lim, cnst.N_INTRP)
    lim_r = np.sum(r_flat[lim_idx].reshape(n_lim, cnst.N_INTRP) * lim_weight, axis=1)
    lim_z = np.sum(z_flat[lim_idx].reshape(n_lim, cnst.N_INTRP) * lim_weight, axis=1)

    # Get Green's matrices
    g_grid_meas_weight = constants_c_dict['g_grid_meas_weight'].reshape(cnst.N_GRID, cnst.N_MEAS)
    # Extract just the first BP_PROBE columns for g_coef_meas_weight
    g_grid_meas_weight = g_grid_meas_weight[:, cnst.N_FLUX_LOOPS:cnst.N_FLUX_LOOPS + cnst.N_BP_PROBES]

    # Get the bp_probe coordinates
    bp_probe_coords = np.loadtxt(os.path.join(cnst.DATA_DIR, 'bp_probes.txt'),
                                  dtype=np.float64, skiprows=1)
    bp_r_coords = bp_probe_coords[:, 0]
    bp_z_coords = bp_probe_coords[:, 1]

    b_theta_ana = analytic_soln.analytic_b_theta(
        r_grid, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB,
        cnst.ANALYTIC_PLASMA_CURRENT)

    for i_iter in range(cnst.N_ITERS + 1):
        psi_num = output_dict['flux_total'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm = output_dict['flux_norm'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm_flat = flux_norm.flatten()
        lcfs_r = output_dict['lcfs_r'][i_iter]
        lcfs_z = output_dict['lcfs_z'][i_iter]
        lcfs_n = output_dict['lcfs_n'][i_iter]
        coef = output_dict['coef'][i_iter]
        # Extract just the first N_PLS coefficients
        coef = coef[:cnst.N_PLS]
        coef = coef.reshape(cnst.N_PLS, 1)
        meas = output_dict['meas'][i_iter]

        # Calcualte the B_θ values
        br = -np.gradient(psi_num, z_vec, axis=0) / r_grid
        bz = np.gradient(psi_num, r_vec, axis=1) / r_grid
        theta = np.arctan2((z_grid - cnst.ANALYTIC_ZO), (r_grid - cnst.ANALYTIC_RO))
        b_theta = -br * np.sin(theta) + bz * np.cos(theta)

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
        meas_fit = meas_fit.flatten()

        bp_probe_meas_fit = meas_fit
        bp_probe_meas = meas[cnst.N_FLUX_LOOPS:cnst.N_FLUX_LOOPS + cnst.N_BP_PROBES]

        # Create a figure with two subplots (vertical layout)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 14), gridspec_kw={'height_ratios': [2, 1]})

        # --- Top subplot: Contour plot ---
        ax1.contour(r_grid, z_grid, b_theta_ana, levels=20, colors='tab:blue')
        ax1.contour(r_grid, z_grid, b_theta, levels=20, colors='tab:orange')

        # Limiter and LCFS points
        ax1.scatter(lim_r, lim_z, color='tab:red', label='Limiter Points')
        ax1.scatter(lcfs_r[:lcfs_n], lcfs_z[:lcfs_n], color='tab:orange', label='LCFS Points')

        # Plot the flux loop coordinates
        ax1.scatter(bp_r_coords, bp_z_coords, color='tab:green', label='BP Probe points', marker='o')
        # Add integer label to each flux loop point
        for idx, (r, z) in enumerate(zip(bp_r_coords, bp_z_coords)):
            ax1.text(r, z, str(idx), fontsize=12, color='black',
                     ha='center', va='center', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            
        # Plot the vessel filament coordiantes
        ax1.scatter(cnst.R_VESSEL, cnst.Z_VESSEL, color='black', label='Vessel Filament Points', marker='o')

        # Legend using proxy artists
        contour_proxy = plt.Line2D([0], [0], color='tab:blue', label='Analytic B_θ')
        num_proxy = plt.Line2D([0], [0], color='tab:orange', label='Numerical B_θ')
        lim_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', label='Limiter Points', linestyle='None')
        lcfs_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:orange', label='LCFS Points', linestyle='None')
        floop_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:green', label='BP Probe Points', linestyle='None')
        vessel_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Vessel Filament Points', linestyle='None')
        ax1.legend(handles=[contour_proxy, num_proxy, lim_proxy, lcfs_proxy, floop_proxy, vessel_proxy],
                   loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)

        ax1.set_title(f'B_θ in T (Iteration {i_iter:02d})')
        ax1.set_xlabel('R [m]')
        ax1.set_ylabel('z [m]')
        ax1.grid(True)
        ax1.set_aspect('equal')

        # --- Bottom subplot: Flux loop measurements ---
        ax2.plot(bp_probe_meas, label='Analytic Flux Loop Values', marker='o', linestyle='None', color='tab:blue')
        ax2.plot(bp_probe_meas_fit, label='RTGSFIT Flux Loop values', marker='x', linestyle='None', color='tab:orange')
        # Make x labels all integers
        ax2.set_xticks(np.arange(cnst.N_BP_PROBES))
        # ax2.set_title('Flux Loop Analytic Values vs RTGSFIT Values')
        ax2.set_xlabel('BP Probe Index')
        ax2.set_ylabel('BP Probe Measurement (T)')
        ax2.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), borderaxespad=0)
        ax2.grid(True)

        # Save the combined plot
        plot_dir = "combined_b_theta_and_bp_probes"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'combined_b_theta_and_bp_probes_iter_{i_iter:02d}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def midplane_ana_vs_num(output_dict):
    """
    Plot the midplane B_θ comparison between analytic and numerical values.
    
    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """
    r_vec = np.linspace(cnst.R_MIN, cnst.R_MAX, cnst.N_R)
    z_vec = np.linspace(cnst.Z_MIN, cnst.Z_MAX, cnst.N_Z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    z_mid = z_vec[cnst.N_Z // 2]

    for i_iter in range(cnst.N_ITERS + 1):
        psi_num = output_dict['flux_total'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        theta = np.arctan2(z_grid - cnst.ANALYTIC_ZO, r_grid - cnst.ANALYTIC_RO)
        br = -np.gradient(psi_num, z_vec, axis=0) / r_grid
        bz = np.gradient(psi_num, r_vec, axis=1) / r_grid
        b_theta_num = -br * np.sin(theta) + bz * np.cos(theta)

        b_theta_ana = analytic_soln.analytic_b_theta(
            r_vec, z_mid,
            cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
            cnst.ANALYTIC_RHOB,
            cnst.ANALYTIC_PLASMA_CURRENT)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(r_vec, b_theta_ana, label='Analytic B_θ', color='tab:blue')
        ax.plot(r_vec, b_theta_num[cnst.N_Z // 2, :], label='Numerical B_θ', color='tab:orange', linestyle='dashed')
        # Add vertical lines at cnst.R_LIM_MIN and cnst.R_LIM_MAX
        ax.axvline(cnst.R_LIM_MIN, color='red', linestyle='--', label='R_LIM_MIN')
        ax.axvline(cnst.R_LIM_MAX, color='red', linestyle='--', label='R_LIM_MAX')
        ax.set_title(f'Midplane B_θ Comparison (Iteration {i_iter:02d})')
        ax.set_xlabel('Radial Coordinate (R) [m]')
        ax.set_ylabel('B_θ [T]')
        ax.legend()
        ax.grid(True)
        plot_dir = "midplane_b_theta_plots"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'midplane_b_theta_iter_{i_iter:02d}.png'
        plt.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name), dpi=300)
        plt.close() 

if __name__ == "__main__":
    # Load the output dictionary
    output_dict = np.load(os.path.join(cnst.DATA_DIR, 'output_dict.npy'), allow_pickle=True).item()

    # Call the contour plotting function
    midplane_ana_vs_num(output_dict)
    contour_ana_vs_num(output_dict)