"""
This module contains function to plot the plasma currents profiles
and the vessel currents.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from rtgsfit_verify_analytic import analytic_soln, cnst, read_constants_c

def vessel_j_j0_j1_j2_j_total_j_ana(output_dict):
    """
    Plot the plasma current profiles and the vessel currents.

    This function makes a figure with 6 subplots, 2 rows and 3 columns.
    
    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """
    
    constants_c_path = os.path.join(cnst.DATA_DIR, 'constants.c')
    constants_c_dict = read_constants_c.constants_c_dict(constants_c_path,
                                                          read_g_grid_meas_weight=True)
      
    r_vec = constants_c_dict['r_vec']
    z_vec = constants_c_dict['z_vec']
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    rho = np.sqrt((r_grid - cnst.ANALYTIC_RO)**2 + (z_grid - cnst.ANALYTIC_ZO)**2)

    # psi_ana = analytic_soln.analytic_psi(
    #     r_grid,
    #     z_grid, 
    #     cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
    #     cnst.ANALYTIC_RHOB, 
    #     cnst.ANALYTIC_PLASMA_CURRENT
    # )
    # psi_n_ana = (cnst.PSI_O - psi_ana) / (cnst.PSI_O - cnst.PSI_B) * (rho <= cnst.ANALYTIC_RHOB) \
    #           + (rho > cnst.ANALYTIC_RHOB)
    # curr_dens_ana = cnst.A_RAW * (1 - psi_n_ana) / (r_grid * cnst.MU_0)
    curr_dens_ana = analytic_soln.analytic_j_phi(
        r_grid, z_grid,
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB,
        cnst.ANALYTIC_PLASMA_CURRENT
    )
    # Check sum of current density over the grid * D_R * D_Z is 1 MA
    curr_dens_sum = np.sum(curr_dens_ana) * cnst.D_R * cnst.D_Z
    print(f"Sum of analytic current density over the grid: {curr_dens_sum / 1e6:.6f} MA")

    for i_iter in range(cnst.N_ITERS + 1):

        coef = output_dict['coef'][i_iter]
        flux_norm = output_dict['flux_norm'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm_z = np.gradient(flux_norm, z_vec, axis=0)

        current0 = coef[0] * (1 - flux_norm) * r_grid
        current1 = coef[1] * (1 - flux_norm) / (r_grid * cnst.MU_0)
        current2 = coef[2] * flux_norm_z
        current_total = current0 + current1 + current2


        fig, ax = plt.subplots(2, 3, figsize=(18, 10))
        ax[0, 0].plot(coef[cnst.N_PLS:] / 1e3,
                      'o',
                      label = 'RTGSFIT')
        ax[0, 0].plot(np.zeros(cnst.N_VESSEL),
                      'x',
                      label = 'Analytic')
        ax[0, 0].legend()
        ax[0, 0].set_ylabel('Vessel Current (kA)')
        # ax[0, 0].set_xlabel('Vessel Index')
        ax[0, 0].set_title(f'Vessel Currents vs. Vessel Filamanet Index')
        # Set x-axis ticks to be integers
        ax[0, 0].set_xticks(np.arange(cnst.N_VESSEL))

        ax[0, 1].contourf(r_grid, z_grid, current0 / 1e6, levels=50)
        cbar = fig.colorbar(ax[0, 1].collections[0], ax=ax[0, 1], orientation='vertical')
        ax[0, 1].set_aspect('equal', adjustable='box')
        # ax[0, 1].set_xlabel('R [m]')
        ax[0, 1].set_ylabel('Z [m]')
        ax[0, 1].set_title("p' Current Density")
        cbar.set_label('Current Density (MA/m$^2$)')

        ax[0, 2].contourf(r_grid, z_grid, current1 / 1e6, levels=50)
        cbar = fig.colorbar(ax[0, 2].collections[0], ax=ax[0, 2], orientation='vertical')
        ax[0, 2].set_aspect('equal', adjustable='box')
        # ax[0, 2].set_xlabel('R [m]')
        ax[0, 2].set_ylabel('Z [m]')
        ax[0, 2].set_title("FF' Current Density")
        cbar.set_label('Current Density (MA/m$^2$)')

        ax[1, 0].contourf(r_grid, z_grid, current2 / 1e6, levels=50)
        cbar = fig.colorbar(ax[1, 0].collections[0], ax=ax[1, 0], orientation='vertical')
        ax[1, 0].set_aspect('equal', adjustable='box')
        ax[1, 0].set_xlabel('R [m]')
        ax[1, 0].set_ylabel('Z [m]')
        ax[1, 0].set_title("dψ/dz Current Density")
        cbar.set_label('Current Density (MA/m$^2$)')

        ax[1, 1].contourf(r_grid, z_grid, current_total / 1e6, levels=50)
        cbar = fig.colorbar(ax[1, 1].collections[0], ax=ax[1, 1], orientation='vertical')
        ax[1, 1].set_aspect('equal', adjustable='box')
        ax[1, 1].set_xlabel('R [m]')
        ax[1, 1].set_ylabel('Z [m]')
        ax[1, 1].set_title("Total Current Density")
        cbar.set_label('Current Density (MA/m$^2$)')

        c5 = ax[1, 2].contourf(r_grid, z_grid, curr_dens_ana / 1e6, levels=50)
        cbar = fig.colorbar(c5, ax=ax[1, 2], orientation='vertical')
        ax[1, 2].set_aspect('equal', adjustable='box')
        ax[1, 2].set_xlabel('R [m]')
        ax[1, 2].set_ylabel('Z [m]')
        ax[1, 2].set_title("Analytic Current Density")
        cbar.set_label('Current Density (MA/m$^2$)')

        fig.suptitle(f'Plasma Current Profiles and Vessel Currents at Iteration {i_iter:02d}')
        
        plot_dir = "vessel_j_j0_j1_j2_j_total_j_ana"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'vessel_j_j0_j1_j2_j_total_j_ana_{i_iter:02d}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def midplane_ana_vs_num(output_dict):
    """
    Plot the midplane current density for both analytic and numerical solutions.

    Parameters:
    output_dict (dict): Dictionary containing the output data from RTGSFIT.
    """
    
    constants_c_path = os.path.join(cnst.DATA_DIR, 'constants.c')
    constants_c_dict = read_constants_c.constants_c_dict(constants_c_path,
                                                          read_g_grid_meas_weight=True)
      
    r_vec = constants_c_dict['r_vec']
    z_vec = constants_c_dict['z_vec']
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    rho = np.sqrt((r_grid - cnst.ANALYTIC_RO)**2 + (z_grid - cnst.ANALYTIC_ZO)**2)

    midplane_z = z_vec[cnst.N_Z // 2]  # Midplane z-coordinate

    psi_ana = analytic_soln.analytic_psi(
        r_grid,
        z_grid, 
        cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_RHOB, 
        cnst.ANALYTIC_PLASMA_CURRENT
    )
    psi_n_ana = (cnst.PSI_O - psi_ana) / (cnst.PSI_O - cnst.PSI_B) * (rho <= cnst.ANALYTIC_RHOB) \
              + (rho > cnst.ANALYTIC_RHOB)
    curr_dens_ana = cnst.A_RAW * (1 - psi_n_ana) / (r_grid * cnst.MU_0)
    # Check sum of current density over the grid * D_R * D_Z is 1 MA
    curr_dens_sum = np.sum(curr_dens_ana) * cnst.D_R * cnst.D_Z
    print(f"Sum of analytic current density over the grid: {curr_dens_sum / 1e6:.6f} MA")

    for i_iter in range(cnst.N_ITERS + 1):

        coef = output_dict['coef'][i_iter]
        flux_norm = output_dict['flux_norm'][i_iter].reshape(cnst.N_Z, cnst.N_R)
        flux_norm_z = np.gradient(flux_norm, z_vec, axis=0)

        current0 = coef[0] * (1 - flux_norm) * r_grid
        current1 = coef[1] * (1 - flux_norm) / (r_grid * cnst.MU_0)
        current2 = coef[2] * flux_norm_z
        current_total = current0 + current1 + current2


        fig, ax = plt.subplots(2, 3, figsize=(18, 10))
        ax[0, 0].plot(coef[cnst.N_PLS:] / 1e3,
                      'o',
                      label = 'RTGSFIT')
        ax[0, 0].plot(np.zeros(cnst.N_VESSEL),
                      'x',
                      label = 'Analytic')
        ax[0, 0].legend()
        ax[0, 0].set_ylabel('Vessel Current (kA)')
        ax[0, 0].set_xlabel('Vessel Index')
        # ax[0, 0].set_title(f'Vessel Currents vs. Vessel Filamanet Index')
        # Set x-axis ticks to be integers
        ax[0, 0].set_xticks(np.arange(cnst.N_VESSEL))

        ax[0, 1].plot(r_vec, current0[cnst.N_Z // 2, :] / 1e6, label='RTGSFIT')
        # Add a horizontal line at 0
        ax[0, 1].axhline(0, color = 'tab:orange', linestyle='dashed', label="Analytic")
        ax[0, 1].legend()
        ax[0, 1].set_xlabel('R [m]')
        ax[0, 1].set_ylabel("p' Current Density [MA/m$^2$]")

        # ax[0, 2].plot(r_grid, z_grid, current1 / 1e6, levels=50)
        ax[0, 2].plot(r_vec, current1[cnst.N_Z // 2, :] / 1e6, label='RTGSFIT')
        ax[0, 2].plot(r_vec, curr_dens_ana[cnst.N_Z // 2, :] / 1e6, label='Analytic', linestyle='dashed')
        ax[0, 2].legend()
        ax[0, 2].set_xlabel('R [m]')
        ax[0, 2].set_ylabel("FF' Current Density [MA/m$^2$]")

        ax[1, 0].plot(r_vec, current2[cnst.N_Z // 2, :] / 1e6, label='RTGSFIT')
        # Add a horizontal line at 0
        ax[1, 0].axhline(0, color = 'tab:orange', linestyle='dashed', label="Analytic")
        ax[1, 0].legend()
        ax[1, 0].set_xlabel('R [m]')
        ax[1, 0].set_ylabel("dψ/dz Current Density [MA/m$^2$]")

        ax[1, 1].plot(r_vec, current_total[cnst.N_Z // 2, :] / 1e6, label='RTGSFIT')
        ax[1, 1].plot(r_vec, curr_dens_ana[cnst.N_Z // 2, :] / 1e6, label='Analytic', linestyle='dashed')
        ax[1, 1].legend()
        ax[1, 1].set_xlabel('R [m]')
        ax[1, 1].set_ylabel('Total Current Density [MA/m$^2$]')

        # Delete the last subplot
        ax[1, 2].remove()

        fig.suptitle('Plasma Current Profiles and Vessel Currents \n'
                     f'Iteration {i_iter:02d} and z = {midplane_z:.2f} m')
        
        plot_dir = "midplane_vessel_j_j0_j1_j2_j_total_j_ana"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'midplane_vessel_j_j0_j1_j2_j_total_j_ana_{i_iter:02d}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")


if __name__ == "__main__":
    # Load the output dictionary
    output_dict = np.load(os.path.join(cnst.DATA_DIR, 'output_dict.npy'), allow_pickle=True).item()

    # Call the contour plotting function
    midplane_ana_vs_num(output_dict)
    vessel_j_j0_j1_j2_j_total_j_ana(output_dict)