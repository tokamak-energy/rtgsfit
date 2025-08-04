"""
Python module for plotting comparing the numerical and analytic values of ψ in Wb.

This module assumes the rtgsfit_output_dict.npy file and gsfit_output_dict.npy file
have already been created and contain the necessary data for plotting.
"""

import os

import matplotlib.pyplot as plt
import mdsthin
import numpy as np

from rtgsfit_vs_gsfit import cnst

def contour_gsfit_vs_rtgsfit(rtgsfit_output_dict: dict, gsfit_output_dict: dict):
    """
    Plot the contour of the RTGSFIT and GSFIT values of ψ from the output_dict.
    
    Parameters:
    rtgsfit_output_dict (dict): Dictionary containing the output data from RTGSFIT.
    gsfit_output_dict (dict): Dictionary containing the output data from GSFIT.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:Z_VEC").data()
        n_limit = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_LIMIT").data()
        n_intrp = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_INTRP").data()
        mask_lim = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:MASK_LIM").data()
        limit_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_R").data()
        limit_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_Z").data()

    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    r_flat = r_grid.flatten()
    z_flat = z_grid.flatten()

    n_r = len(r_vec)
    n_z = len(z_vec)
    # mask_lim = mask_lim.reshape((n_z, n_r))
    plasma_idx = np.where(mask_lim == 1)[0]
    vacuum_idx = np.where(mask_lim == 0)[0]

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter_size = 2
    ax.scatter(limit_r, limit_z, label='RTGSFIT Limiter', s=scatter_size)
    ax.scatter(r_flat[plasma_idx], z_flat[plasma_idx], label='RTGSFIT Plasma', s=scatter_size)
    ax.scatter(r_flat[vacuum_idx], z_flat[vacuum_idx], label='RTGSFIT Vacuum', s=scatter_size)
    ax.set_aspect('equal')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    # Append legend to the plot but outside the plot area
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
    plot_dir = "limiter"
    os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
    plot_name = f'limiter.png'
    fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                dpi=300,
                bbox_inches='tight')
    plt.close("all")

    gsfit_r_grid = gsfit_output_dict["grid"]["r"]
    gsfit_z_grid = gsfit_output_dict["grid"]["z"]
    gsfit_psi = gsfit_output_dict["two_d"]["psi"]
    gsfit_psi_n = gsfit_output_dict["two_d"]["psi_n"]
    gsfit_pboundary_rbnd = gsfit_output_dict["p_boundary"]["rbnd"]
    gsfit_pboundary_zbnd = gsfit_output_dict["p_boundary"]["zbnd"]
    gsfit_pboundary_nbnd = gsfit_output_dict["p_boundary"]["nbnd"]
    gsfit_pboundary_rbnd = gsfit_pboundary_rbnd[:gsfit_pboundary_nbnd]
    gsfit_pboundary_zbnd = gsfit_pboundary_zbnd[:gsfit_pboundary_nbnd]

    # Plot GSFIT ψ
    fig, ax = plt.subplots(figsize=(10, 8))
    cs = ax.contour(gsfit_r_grid, gsfit_z_grid, gsfit_psi,
                    levels=20)
    ax.scatter(gsfit_pboundary_rbnd, gsfit_pboundary_zbnd,
               label='GSFIT LCFS',
               s=scatter_size)
    ax.scatter(limit_r, limit_z, label='RTGSFIT Limiter', s=scatter_size)
    # Add the colorbar
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label('ψ (Wb)')
    ax.set_title('GSFIT ψ')
    ax.set_aspect('equal')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    plot_dir = "gsfit_psi"
    os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
    plot_name = f'gsfit_psi.png'
    fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                dpi=300,
                bbox_inches='tight')
    plt.close("all")

    # Plot RTGSFIT ψ vs GSFIT ψ for each iteration
    for i_iter in range(cnst.N_ITERS + 1):
        rtgsfit_psi = rtgsfit_output_dict["flux_total"][i_iter, :].reshape((n_z, n_r)) * 2 * np.pi
        rtgsfit_lcfs_r = rtgsfit_output_dict["lcfs_r"][i_iter, :]
        rtgsfit_lcfs_z = rtgsfit_output_dict["lcfs_z"][i_iter, :]
        rtgsfit_lcfs_n = rtgsfit_output_dict["lcfs_n"][i_iter, 0]
        rtgsfit_lcfs_r = rtgsfit_lcfs_r[:rtgsfit_lcfs_n]
        rtgsfit_lcfs_z = rtgsfit_lcfs_z[:rtgsfit_lcfs_n]
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.contour(r_grid, z_grid, rtgsfit_psi,
                   levels=20,
                   colors='tab:blue')
        ax.contour(gsfit_r_grid, gsfit_z_grid, gsfit_psi,
                   levels=20,
                   colors='tab:orange')
        ax.scatter(rtgsfit_lcfs_r, rtgsfit_lcfs_z,
                   label='RTGSFIT LCFS',
                   c='tab:green',
                   s=scatter_size)
        ax.scatter(gsfit_pboundary_rbnd, gsfit_pboundary_zbnd,
                   label='GSFIT LCFS',
                   c='tab:red',
                   s=scatter_size)
        ax.scatter(limit_r, limit_z,
                   label='RTGSFIT Limiter',
                   s=scatter_size,
                   c='tab:purple')
        # Add legend but we need to make fake handles for the legend
        blue_patch = plt.Line2D([], [], color='tab:blue', label='RTGSFIT ψ')
        orange_patch = plt.Line2D([], [], color='tab:orange', label='GSFIT ψ')
        green_patch = plt.Line2D([], [], color='tab:green', label='RTGSFIT LCFS')
        red_patch = plt.Line2D([], [], color='tab:red', label='GSFIT LCFS')
        purple_patch = plt.Line2D([], [], color='tab:purple', label='RTGSFIT Limiter')
        ax.legend(handles=[blue_patch, orange_patch, green_patch, red_patch, purple_patch],
                   loc='upper right')
        ax.set_aspect('equal')
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        ax.set_title(f'RTGSFIT vs GSFIT ψ (Iteration {i_iter})')
        plot_dir = "rtgsfit_vs_gsfit_psi"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'rtgsfit_vs_gsfit_psi_iter_{i_iter}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

    # Plot just the RTGSFIT ψ on its own
    for i_iter in range(cnst.N_ITERS + 1):
        rtgsfit_psi = rtgsfit_output_dict["flux_total"][i_iter, :].reshape((n_z, n_r))  * 2 * np.pi
        rtgsfit_lcfs_r = rtgsfit_output_dict["lcfs_r"][i_iter, :]
        rtgsfit_lcfs_z = rtgsfit_output_dict["lcfs_z"][i_iter, :]
        rtgsfit_lcfs_n = rtgsfit_output_dict["lcfs_n"][i_iter, 0]
        rtgsfit_lcfs_r = rtgsfit_lcfs_r[:rtgsfit_lcfs_n]
        rtgsfit_lcfs_z = rtgsfit_lcfs_z[:rtgsfit_lcfs_n]
        fig, ax = plt.subplots(figsize=(10, 8))
        cs = ax.contour(r_grid, z_grid, rtgsfit_psi,
                        levels=20)
        cbar = fig.colorbar(cs, ax=ax)
        cbar.set_label('ψ (Wb)')
        ax.scatter(rtgsfit_lcfs_r, rtgsfit_lcfs_z,
                   label='RTGSFIT LCFS',
                   s=scatter_size)
        ax.scatter(limit_r, limit_z,
                   label='RTGSFIT Limiter',
                   s=scatter_size)
        # Add legend for scatter plots
        ax.legend(loc='upper right')
        ax.set_aspect('equal')
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        ax.set_title(f'RTGSFIT ψ (Iteration {i_iter})')
        plot_dir = "rtgsfit_psi"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'rtgsfit_psi_iter_{i_iter}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def plot_psi_ltrb(rtgsfit_output_dict: dict, gsfit_output_dict: dict):
    """
    Plot the ψ values for the left, top, right, and bottom boundaries.

    Parameters:
    rtgsfit_output_dict (dict): Dictionary containing the output data from RTGSFIT.
    gsfit_output_dict (dict): Dictionary containing the output data from GSFIT.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_ltrb = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_LTRB").data()
        # inv_r_ltrb_mu0 = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:INV_R_L_MU0").data()
        nr = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_R").data()
        nz = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_Z").data()
    # r_ltrb = (cnst.MU_0 / inv_r_ltrb_mu0)

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cnst.PULSE_NUM_WRITE)
        gsfit_flux_total = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.TWO_D:PSI").data()[0, :, :]

    gsfit_psi_ltrb = np.zeros(n_ltrb, dtype=np.float64)
    gsfit_psi_ltrb[:nz] = gsfit_flux_total[:nz, 0]
    gsfit_psi_ltrb[nz:nz+nr - 1] = gsfit_flux_total[-1, 1:]
    gsfit_psi_ltrb[nz+nr-1:2*nz+nr-2] = np.flip(gsfit_flux_total[1:, -1])
    gsfit_psi_ltrb[2*nz+nr-2:] = np.flip(gsfit_flux_total[0, 1:-1])

    for i_iter in range(cnst.N_ITERS + 1):

        flux_total = rtgsfit_output_dict["flux_total"][i_iter, :].reshape((nz, nr))

        rtgsfit_psi_ltrb = np.zeros(n_ltrb, dtype=np.float64)
        rtgsfit_psi_ltrb[:nz] = flux_total[:nz, 0]
        rtgsfit_psi_ltrb[nz:nz+nr - 1] = flux_total[-1, 1:]
        rtgsfit_psi_ltrb[nz+nr-1:2*nz+nr-2] = np.flip(flux_total[1:, -1])
        rtgsfit_psi_ltrb[2*nz+nr-2:] = np.flip(flux_total[0, 1:-1])

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(rtgsfit_psi_ltrb, label=f'RTGSFIT ψ LTRB Iter {i_iter}')
        ax.plot(gsfit_psi_ltrb, label='GSFIT ψ LTRB')
        ylim = ax.get_ylim()
        ax.vlines([0, nz, nz+nr-1, 2*nz+nr-2],
                ymin=-1e6,
                ymax=1e6,
                linestyles='dashed',
                colors="k",
                label='Boundary Indices')
        ax.set_ylim(ylim)
        ax.set_xlabel('Boundary Index (LTRB)')
        ax.set_ylabel('ψ (Wb)')
        ax.set_title('ψ LTRB Comparison')
        ax.set_xticks(range(0, n_ltrb, 10))
        ax.set_xticklabels([f'{i}' for i in range(0, n_ltrb, 10)])
        ax.set_xlim(0, n_ltrb - 1)
        ax.set_ylim(np.min(gsfit_psi_ltrb) - 0.1, np.max(gsfit_psi_ltrb) + 0.1)
        ax.grid(True)
        ax.legend()
        plot_dir = "psi_ltrb"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'psi_ltrb_iter_{i_iter}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")


if __name__ == "__main__":

    # Load the output dictionary
    rtgsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'rtgsfit_output_dict.npy'),
                                  allow_pickle=True).item()
    gsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy'),
                                allow_pickle=True).item()

    # contour_gsfit_vs_rtgsfit(rtgsfit_output_dict, gsfit_output_dict)
    plot_psi_ltrb(rtgsfit_output_dict, gsfit_output_dict)