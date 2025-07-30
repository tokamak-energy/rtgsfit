"""
Python module for plotting measured values from the sensors and fitted values from RTGSFIT and GSFIT.
"""

import os

import matplotlib.pyplot as plt
import mdsthin
import numpy as np
from scipy.interpolate import RectBivariateSpline

from rtgsfit_vs_gsfit import cnst, calc_meas

def plot_flux_loop_at_sensors(rtgsfit_output_dict: dict,
                              gsfit_output_dict: dict) -> None:
    """
    Plot the flux loop values at the sensors.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_flux_loops = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_F_LOOPS").data()
        n_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_R").data()
        n_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_Z").data()
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:Z_VEC").data()
        limit_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_R").data()
        limit_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_Z").data()

    rtgsfit_fl_meas = rtgsfit_output_dict['meas'][0, :n_flux_loops]

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cnst.PULSE_NUM_WRITE)
        fl_include = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.FLOOP:INCLUDE").data() == 1
        gsfit_fl_meas = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.FLOOP:MVALUE").data()[0, fl_include]
        gsfit_fl_pred = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.FLOOP:CVALUE").data()[0, fl_include]
        gsfit_flux_total = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.TWO_D:PSI").data()[0, :, :]
        gsfit_lcfs_n = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:NBND").data()[0]
        gsfit_lcfs_r = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:RBND").data()[0, :gsfit_lcfs_n]
        gsfit_lcfs_z = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:ZBND").data()[0, :gsfit_lcfs_n]

    flux_loop_r = gsfit_output_dict["flux_loops"]["r"]
    flux_loop_z = gsfit_output_dict["flux_loops"]["z"]
    flux_loop_names = gsfit_output_dict["flux_loops"]["names"]

    # Calculate and plot the RTGSFIT predictions at the sensors
    for i_iter in range(cnst.N_ITERS + 1):

        flux_norm = rtgsfit_output_dict["flux_norm"][i_iter, :]
        coef = rtgsfit_output_dict["coef"][i_iter, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][i_iter, :]
        pred_meas = calc_meas.calculate_predicted_measurements(flux_norm, coef, coil_curr)
        pred_meas_fl = pred_meas[:n_flux_loops]

        lcfs_n = rtgsfit_output_dict["lcfs_n"][i_iter, 0]
        lcfs_r = rtgsfit_output_dict["lcfs_r"][i_iter, :lcfs_n]
        lcfs_z = rtgsfit_output_dict["lcfs_z"][i_iter, :lcfs_n]

        flux_total = rtgsfit_output_dict["flux_total"][i_iter, :]
        flux_total = flux_total.reshape((n_z, n_r))
        spline = RectBivariateSpline(z_vec, r_vec, flux_total)
        flux_total_interp = spline.ev(flux_loop_z, flux_loop_r)

        fig, (ax_top, ax_bottom) = plt.subplots(2, 1,
                                                figsize=(10, 12),
                                                gridspec_kw={'height_ratios': [2, 1]})
        dot_size = 5
        ax_top.contour(r_vec, z_vec, flux_total * 2 * np.pi,
                       levels=20,
                       colors='tab:blue')
        ax_top.contour(r_vec, z_vec, gsfit_flux_total,
                       levels=20,
                       colors='tab:orange')
        ax_top.scatter(lcfs_r, lcfs_z,
                       label='RTGSFIT LCFS',
                       color='tab:blue',
                       s=dot_size)
        ax_top.scatter(gsfit_lcfs_r, gsfit_lcfs_z,
                       label='GSFIT LCFS',
                       color='tab:orange',
                       s=dot_size)
        ax_top.scatter(limit_r, limit_z,
                       label='RTGSFIT Limiter',
                       color='tab:green',
                       s=dot_size)
        blue_patch = plt.Line2D([], [], color='tab:blue', label='RTGSFIT ψ x 2π')
        orange_patch = plt.Line2D([], [], color='tab:orange', label='GSFIT ψ')
        blue_dots = plt.Line2D([], [], marker='o', color='w', label='RTGSFIT LCFS',
                               markerfacecolor='tab:blue', markersize=dot_size)
        orange_dots = plt.Line2D([], [], marker='o', color='w', label='GSFIT LCFS',
                                 markerfacecolor='tab:orange', markersize=dot_size)
        green_dots = plt.Line2D([], [], marker='o', color='w', label='RTGSFIT Limiter',
                                markerfacecolor='tab:green', markersize=dot_size)
        ax_top.legend(handles=[blue_patch, orange_patch, blue_dots, orange_dots, green_dots],
                      loc='upper right')
        ax_top.set_aspect('equal')
        ax_top.set_xlabel('R (m)')
        ax_top.set_ylabel('Z (m)')
        ax_top.set_title(f'RTGSFIT vs GSFIT ψ (Iteration {i_iter})' '\n'
                         f'(Pulse #{cnst.PULSE_NUM}, Time {cnst.TIME*1e3:.1f} ms)')
        for j, name in enumerate(flux_loop_names):
            ax_top.annotate(name, xy=(flux_loop_r[j], flux_loop_z[j]),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=8, color='black', ha='center')
            
        ax_bottom.plot(flux_loop_names, pred_meas_fl, label='RTGSFIT Predicted', linestyle='--')
        ax_bottom.plot(flux_loop_names, gsfit_fl_pred, label='GSFIT Predicted', linestyle='-.')
        ax_bottom.plot(flux_loop_names, rtgsfit_fl_meas, label='Measured')
        ax_bottom.plot(flux_loop_names, flux_total_interp * 2 * np.pi,
                       label='RTGSFIT Interpolated x 2π', linestyle=':')
        ax_bottom.set_ylabel('Flux Loop Value')
        ax_bottom.tick_params(axis='x', rotation=90)
        ax_bottom.legend()

        plot_dir = "flux_loop_contours"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'flux_loop_contours_iter_{i_iter}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def plot_b_at_sensors(rtgsfit_output_dict, gsfit_output_dict):
    """
    Plot the B field values at the sensors.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_flux_loops = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_F_LOOPS").data()
        n_bp_probes = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_BP_PROBES").data()
        n_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_R").data()
        n_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_Z").data()
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:Z_VEC").data()
        limit_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_R").data()
        limit_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:LIMIT_Z").data()
        weights = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:WEIGHT").data()

    r_grid, _ = np.meshgrid(r_vec, z_vec) 

    bp_probe_range = np.arange(n_flux_loops, n_flux_loops + n_bp_probes)
    rtgsfit_bp_meas = rtgsfit_output_dict['meas'][0, bp_probe_range]
    weights = weights[bp_probe_range]

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cnst.PULSE_NUM_WRITE)
        bp_include = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.BPPROBE:INCLUDE").data() == 1
        gsfit_bp_meas = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.BPPROBE:MVALUE").data()[0, bp_include]
        gsfit_bp_pred = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.BPPROBE:CVALUE").data()[0, bp_include]
        gsfit_br = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.TWO_D:BR").data()[0, :, :]
        gsfit_bz = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.TWO_D:BZ").data()[0, :, :]
        gsfit_lcfs_n = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:NBND").data()[0]
        gsfit_lcfs_r = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:RBND").data()[0, :gsfit_lcfs_n]
        gsfit_lcfs_z = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:ZBND").data()[0, :gsfit_lcfs_n]

    bp_probe_r = gsfit_output_dict["bp_probes"]["r"]
    bp_probe_z = gsfit_output_dict["bp_probes"]["z"]
    bp_probe_angle = gsfit_output_dict["bp_probes"]["angle"]
    bp_probe_names = gsfit_output_dict["bp_probes"]["names"]

    # Calculate and plot the RTGSFIT predictions at the sensors
    for i_iter in range(cnst.N_ITERS + 1):

        flux_norm = rtgsfit_output_dict["flux_norm"][i_iter, :]
        coef = rtgsfit_output_dict["coef"][i_iter, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][i_iter, :]
        pred_meas = calc_meas.calculate_predicted_measurements(flux_norm, coef, coil_curr)
        pred_meas_bp = pred_meas[bp_probe_range]

        lcfs_n = rtgsfit_output_dict["lcfs_n"][i_iter, 0]
        lcfs_r = rtgsfit_output_dict["lcfs_r"][i_iter, :lcfs_n]
        lcfs_z = rtgsfit_output_dict["lcfs_z"][i_iter, :lcfs_n]

        flux_total = rtgsfit_output_dict["flux_total"][i_iter, :]
        flux_total = flux_total.reshape((n_z, n_r))
        # br = -(dψ/dz) / r
        br = -np.gradient(flux_total, z_vec, axis=0) / r_grid
        # bz = dψ/dr / r
        bz = np.gradient(flux_total, r_vec, axis=1) / r_grid
        br_spline = RectBivariateSpline(z_vec, r_vec, br)
        bz_spline = RectBivariateSpline(z_vec, r_vec, bz)
        br_interp = br_spline.ev(bp_probe_z, bp_probe_r)
        bz_interp = bz_spline.ev(bp_probe_z, bp_probe_r)
        bp_interp = np.cos(bp_probe_angle) * br_interp + \
                    np.sin(bp_probe_angle) * bz_interp

        # Make a figure with 3 sublots two on top and one on the bottom
        # The top two subplots are contours of the B field
        # The top left is Br and the top right is Bz
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])
        ax_top_left = fig.add_subplot(gs[0, 0])
        ax_top_right = fig.add_subplot(gs[0, 1])
        ax_bottom = fig.add_subplot(gs[1, :])

        dot_size = 5

        # Top left subplot for Br
        ax_top_left.contour(r_vec, z_vec, br,
                            levels=20,
                            colors='tab:blue')
        ax_top_left.contour(r_vec, z_vec, gsfit_br,
                            levels=100,
                            colors='tab:orange')
        ax_top_left.scatter(lcfs_r, lcfs_z,
                            label='RTGSFIT LCFS',
                            color='tab:blue',
                            s=dot_size)
        ax_top_left.scatter(gsfit_lcfs_r, gsfit_lcfs_z,
                            label='GSFIT LCFS',
                            color='tab:orange',
                            s=dot_size)
        ax_top_left.scatter(limit_r, limit_z,
                            label='RTGSFIT Limiter',
                            color='tab:green',
                            s=dot_size)
        blue_patch = plt.Line2D([], [], color='tab:blue', label='RTGSFIT')
        orange_patch = plt.Line2D([], [], color='tab:orange', label='GSFIT')
        blue_dots = plt.Line2D([], [], marker='o', color='w', label='RTGSFIT LCFS',
                               markerfacecolor='tab:blue', markersize=dot_size)
        orange_dots = plt.Line2D([], [], marker='o', color='w', label='GSFIT LCFS',
                                 markerfacecolor='tab:orange', markersize=dot_size)
        green_dots = plt.Line2D([], [], marker='o', color='w', label='RTGSFIT Limiter',
                                markerfacecolor='tab:green', markersize=dot_size)
        ax_top_left.legend(handles=[blue_patch, orange_patch, blue_dots, orange_dots, green_dots],
                           loc='upper right',
                           bbox_to_anchor=(1.9, 1))
        ax_top_left.set_aspect('equal')
        ax_top_left.set_xlabel('R (m)')
        ax_top_left.set_ylabel('Z (m)')
        for j, name in enumerate(bp_probe_names):
            ax_top_left.annotate(name, xy=(bp_probe_r[j], bp_probe_z[j]),
                                 xytext=(5, 5), textcoords='offset points',
                                 fontsize=8, color='black', ha='center')
        ax_top_left.set_title(r"$B_r$")
            
        # Top right subplot for Bz
        ax_top_right.contour(r_vec, z_vec, bz,
                             levels=20,
                             colors='tab:blue')
        ax_top_right.contour(r_vec, z_vec, gsfit_bz,
                             levels=20,
                             colors='tab:orange')
        ax_top_right.scatter(lcfs_r, lcfs_z,
                             label='RTGSFIT LCFS',
                             color='tab:blue',
                             s=dot_size)
        ax_top_right.scatter(gsfit_lcfs_r, gsfit_lcfs_z,
                             label='GSFIT LCFS',
                             color='tab:orange',
                             s=dot_size)
        ax_top_right.scatter(limit_r, limit_z,
                             label='RTGSFIT Limiter',
                             color='tab:green',
                             s=dot_size)
        ax_top_right.set_aspect('equal')
        ax_top_right.set_xlabel('R (m)')
        ax_top_right.set_title(r"$B_z$")
        ax_top_right.set_ylabel('Z (m)')
        for j, name in enumerate(bp_probe_names):
            ax_top_right.annotate(name, xy=(bp_probe_r[j], bp_probe_z[j]),
                                  xytext=(5, 5), textcoords='offset points',
                                  fontsize=8, color='black', ha='center')
            
        # Bottom subplot for the B field values at the sensors
        ax_bottom.plot(bp_probe_names, pred_meas_bp, label='RTGSFIT Predicted', linestyle='--')
        ax_bottom.plot(bp_probe_names, gsfit_bp_pred, label='GSFIT Predicted', linestyle='-.')
        ax_bottom.plot(bp_probe_names, rtgsfit_bp_meas, label='Measured')
        ax_bottom.plot(bp_probe_names, bp_interp,
                       label='RTGSFIT Interpolated', linestyle=':')
        ax_bottom.set_ylabel('B Field Value [T]')
        ax_bottom.tick_params(axis='x', rotation=90)
        ax_bottom.legend()

        fig.suptitle(f'(Pulse #{cnst.PULSE_NUM}, Time {cnst.TIME*1e3:.1f} ms)',
                     y=0.95)

        plot_dir = "b_field_contours"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        plot_name = f'b_field_contours_iter_{i_iter}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")

def plot_j_at_sensors(rtgsfit_output_dict, gsfit_output_dict):
    """
    Plot vessel currents as well as the current
    and in a table below show the current values at the rogovski coils
    and the vessel current contribution from each of the eigenvectors.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_flux_loops = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_F_LOOPS").data()
        n_bp_probes = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_BP_PROBES").data()
        n_rogowski_coils = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_ROG_COILS").data()
        n_r = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_R").data()
        n_z = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_Z").data()
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:Z_VEC").data()
        weights = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:WEIGHT").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COEF").data()
        n_pls = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_PLS").data()
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:SENS_NAMES").data()
        rtgsfit_coil_names = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:COIL_NAMES").data()
        n_coil_rt = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COIL").data()

    rogowski_range = np.arange(n_flux_loops + n_bp_probes,
                               n_flux_loops + n_bp_probes + n_rogowski_coils)
    rtgsfit_rog_meas = rtgsfit_output_dict['meas'][0, rogowski_range]
    rtgsfit_rog_meas = [f'{current:.2e}' for current in rtgsfit_rog_meas * 1e-6]
    weights = weights[rogowski_range]
    rogovski_names_rt = sens_names[rogowski_range]
    # for table replace I_ROG_ part of string with ''
    rogovski_names_table = [name.replace("I_ROG_", "") for name in rogovski_names_rt]

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cnst.PULSE_NUM_WRITE)
        # rog_include = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:INCLUDE").data() == 1
        # gsfit_rog_meas = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:MVALUE").data()[0, rog_include]
        # gsfit_rog_pred = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:CVALUE").data()[0, rog_include]
        rogovski_names_gs = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:NAME").data()
        gsfit_rog_meas_raw = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:MVALUE").data()[0, :]
        gsfit_rog_pred_raw = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:CVALUE").data()[0, :]
        gsfit_rog_meas = np.zeros(n_rogowski_coils)
        gsfit_rog_pred = np.zeros(n_rogowski_coils)
        for i, rogovski_name_rt in enumerate(rogovski_names_rt):
            for j, rogovski_name_gs in enumerate(rogovski_names_gs):
                if rogovski_name_rt == rogovski_name_gs:
                    gsfit_rog_meas[i] = gsfit_rog_meas_raw[j]
                    gsfit_rog_pred[i] = gsfit_rog_pred_raw[j]
        gsfit_lcfs_n = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:NBND").data()[0]
        gsfit_lcfs_r = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:RBND").data()[0, :gsfit_lcfs_n]
        gsfit_lcfs_z = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.P_BOUNDARY:ZBND").data()[0, :gsfit_lcfs_n]
        gsfit_psi = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.TWO_D:PSI").data()[0, :, :]

    gsfit_rog_pred = [f'{current:.2e}' for current in gsfit_rog_pred * 1e-6]
    gsfit_rog_meas = [f'{current:.2e}' for current in gsfit_rog_meas * 1e-6]

    rogovski_coils_dict = gsfit_output_dict["rogowski_coils"]
    coils_dict = gsfit_output_dict["coils"]
    ivc_dict = gsfit_output_dict["IVC"]
    n_ves_dof = len(ivc_dict["dof"].keys())
    n_ves_fil = len(gsfit_output_dict["IVC"]["r"])
    j_ves_gsfit = np.zeros(n_ves_fil)
    ivc_coef_range = np.arange(n_coef - n_ves_dof - 1, n_coef - 1)
    for i_ivc_dof in range(n_ves_dof):
        current_distribution = ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["current_distribution"]
        calculated = ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["calculated"]
        j_ves_gsfit += current_distribution * calculated / ivc_dict["area"]

    gsfit_coil_names = coils_dict.keys()
    rtgsfit_coil_names_extended = []
    for i in range (len(gsfit_coil_names)):
        if i < len(rtgsfit_coil_names):
            rtgsfit_coil_names_extended.append(rtgsfit_coil_names[i])
        else:
            rtgsfit_coil_names_extended.append(" ")
    rtgsfit_coil_currents_extended = np.zeros(len(gsfit_coil_names))
    rtgsfit_coil_currents = rtgsfit_output_dict["coil_curr"][0, :]
    rtgsfit_coil_currents_extended[:n_coil_rt] = rtgsfit_coil_currents[:n_coil_rt]
    gsfit_coil_currents = np.zeros(len(gsfit_coil_names))
    for i, name_gs in enumerate(gsfit_coil_names):
        gsfit_coil_currents[i] = coils_dict[name_gs]["i"]
    gsfit_coil_currents_shorter = np.zeros(n_coil_rt)
    for i, name_rt in enumerate(rtgsfit_coil_names):
        for name_gs in gsfit_coil_names:
            if name_rt in name_gs:
                gsfit_coil_currents_shorter[i] = coils_dict[name_gs]["i"]
                break

    for i_iter in range(cnst.N_ITERS + 1):

        flux_norm = rtgsfit_output_dict["flux_norm"][i_iter, :]
        coef = rtgsfit_output_dict["coef"][i_iter, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][i_iter, :]
        pred_meas = calc_meas.calculate_predicted_measurements(flux_norm, coef, coil_curr)

        lcfs_n = rtgsfit_output_dict["lcfs_n"][i_iter, 0]
        lcfs_r = rtgsfit_output_dict["lcfs_r"][i_iter, :lcfs_n]
        lcfs_z = rtgsfit_output_dict["lcfs_z"][i_iter, :lcfs_n]

        flux_total = rtgsfit_output_dict["flux_total"][i_iter, :]
        flux_total = flux_total.reshape((n_z, n_r))

        j_ves_rtgsfit = np.zeros(n_ves_fil)
        for i_ivc_dof in range(n_ves_dof):
            current_distribution = ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["current_distribution"]
            calculated = coef[ivc_coef_range[i_ivc_dof]]
            j_ves_rtgsfit += current_distribution * calculated / ivc_dict["area"]

        def plot_vessel_rogowski_coils():
            """"
            Make a figure with 2 sublots two on top
            and a table below
            The top left is rt-gsfit_psi and the top right is gsfit_psi
            """

            import matplotlib.gridspec as gridspec
            fig = plt.figure(figsize=(12, 10))
            fig.subplots_adjust(hspace=0.4)
            gs = gridspec.GridSpec(2, 2, height_ratios=[2.2, 0.8])
            ax_top_left = fig.add_subplot(gs[0, 0])
            ax_top_right = fig.add_subplot(gs[0, 1])
            ax_bottom = fig.add_subplot(gs[1, :],
                                        frame_on=False,
                                        xticks=[],
                                        yticks=[])

            dot_size = 5
            ax_top_left.contour(r_vec, z_vec, flux_total * 2 * np.pi,
                                levels=20,
                                alpha=0.5,
                                colors='tab:blue')
            ax_top_left.scatter(lcfs_r, lcfs_z,
                                alpha=0.5,
                                label='LCFS',
                                edgecolors=None,
                                color='tab:blue',
                                s=dot_size)
            ax_top_left.set_aspect('equal')
            ax_top_left.set_xlabel('R (m)')
            ax_top_left.set_ylabel('Z (m)')
            ax_top_left.set_title('RTGSFIT')
            for name in rogovski_coils_dict.keys():
                if name == "INIVC000": continue  # Skip the INIVC000 coil
                r = rogovski_coils_dict[name]["r"]
                z = rogovski_coils_dict[name]["z"]
                r_centroid = np.sum(r) / len(r)
                z_centroid = np.sum(z) / len(z)
                ax_top_left.plot(r, z,
                                c = "tab:red",
                                label=name)
                ax_top_left.annotate(name, xy=(r_centroid, z_centroid),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='black', ha='center')
            for name in coils_dict.keys():
                r = coils_dict[name]["r"]
                z = coils_dict[name]["z"]
                r_centroid = np.sum(r) / len(r)
                z_centroid = np.sum(z) / len(z)
                ax_top_left.scatter(r, z,
                                    c = "tab:purple",
                                    s = 1,
                                    label=name)
                if name == "BVLT" or name == "BVLB": continue  # Skip the BVLT and BVTB coils
                ax_top_left.annotate(name, xy=(r_centroid, z_centroid),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='black', ha='center')
            blue_patch = plt.Line2D([], [], color='tab:blue', label='ψ')
            blue_dots = plt.Line2D([], [], marker='o', color='w', label='LCFS',
                                markerfacecolor='tab:blue', markersize=dot_size)
            red_patch = plt.Line2D([], [], color='tab:red', label='Rogowski Coils')
            purple_dots = plt.Line2D([], [], marker='o', color='w', label='PF Coils',
                                    markerfacecolor='tab:purple', markersize=dot_size)
            ax_top_left.legend(handles=[blue_patch, blue_dots, red_patch, purple_dots],
                            loc='upper right', fontsize=8)

            # Plot vessel currents as filled polygons
            for idx in range(len(ivc_dict["r"])):
                if i_iter == 0: continue
                r = ivc_dict["r"][idx]
                z = ivc_dict["z"][idx]
                dr = ivc_dict["d_r"][idx]
                dz = ivc_dict["d_z"][idx]
                current = j_ves_rtgsfit[idx]
                filament_r = [r - 0.5 * dr, r + 0.5 * dr, r + 0.5 * dr, r - 0.5 * dr]
                filament_z = [z - 0.5 * dz, z - 0.5 * dz, z + 0.5 * dz, z + 0.5 * dz]
                color = plt.cm.viridis(
                    (current - j_ves_rtgsfit.min()) /
                    (j_ves_rtgsfit.max() - j_ves_rtgsfit.min())
                    )
                ax_top_left.fill(filament_r, filament_z, color=color, edgecolor='none', linewidth=0.5)
            # Put colorbar labels in scientific notation
            sm = plt.cm.ScalarMappable(
                cmap=plt.cm.viridis,
                norm=plt.Normalize(vmin=j_ves_rtgsfit.min() / 1e6, vmax=j_ves_rtgsfit.max() / 1e6)
            )
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax_top_left)
            cbar.set_label('Vessel Current Density (MA/m²)')
            cbar.ax.tick_params(labelsize=8)
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()

            ax_top_right.contour(r_vec, z_vec, gsfit_psi,
                                levels=20,
                                alpha=0.5,
                                colors='tab:orange')
            ax_top_right.scatter(gsfit_lcfs_r, gsfit_lcfs_z,
                                alpha=0.5,
                                label='LCFS',
                                edgecolors=None,
                                color='tab:orange',
                                s=dot_size)
            ax_top_right.set_aspect('equal')
            ax_top_right.set_xlabel('R (m)')
            ax_top_right.set_ylabel('Z (m)')
            ax_top_right.set_title('GSFIT')
            for name in rogovski_coils_dict.keys():
                # if name == "INIVC000": continue  # Skip the INIVC000 coil
                # if GAS is in the name, skip it
                if "GAS" in name: continue
                r = rogovski_coils_dict[name]["r"]
                z = rogovski_coils_dict[name]["z"]
                r_centroid = np.sum(r) / len(r)
                z_centroid = np.sum(z) / len(z)
                ax_top_right.plot(r, z,
                                c = "tab:red",
                                label=name)
                ax_top_right.annotate(name, xy=(r_centroid, z_centroid),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='black', ha='center')
            for name in coils_dict.keys():
                r = coils_dict[name]["r"]
                z = coils_dict[name]["z"]
                r_centroid = np.sum(r) / len(r)
                z_centroid = np.sum(z) / len(z)
                ax_top_right.scatter(r, z,
                                    c = "tab:purple",
                                    s = 1,
                                    label=name)
                if name == "BVLT" or name == "BVLB": continue  # Skip the BVLT and BVTB coils
                ax_top_right.annotate(name, xy=(r_centroid, z_centroid),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='black', ha='center')
            orange_patch = plt.Line2D([], [], color='tab:orange', label='ψ')
            orange_dots = plt.Line2D([], [], marker='o', color='w', label='LCFS',
                                    markerfacecolor='tab:orange', markersize=dot_size)
            ax_top_right.legend(handles=[orange_patch, orange_dots, red_patch, purple_dots],
                                loc='upper right', fontsize=8)
            
            # Plot vessel currents as filled polygons
            for idx in range(len(ivc_dict["r"])):
                if i_iter == 0: continue
                r = ivc_dict["r"][idx]
                z = ivc_dict["z"][idx]
                dr = ivc_dict["d_r"][idx]
                dz = ivc_dict["d_z"][idx]
                current = j_ves_gsfit[idx]
                filament_r = [r - 0.5 * dr, r + 0.5 * dr, r + 0.5 * dr, r - 0.5 * dr]
                filament_z = [z - 0.5 * dz, z - 0.5 * dz, z + 0.5 * dz, z + 0.5 * dz]
                color = plt.cm.viridis(
                    (current - j_ves_gsfit.min()) /
                    (j_ves_gsfit.max() - j_ves_gsfit.min())
                    )
                ax_top_right.fill(filament_r, filament_z, color=color, edgecolor='none', linewidth=0.5)
            sm = plt.cm.ScalarMappable(
                cmap=plt.cm.viridis,
                norm=plt.Normalize(vmin=j_ves_gsfit.min() / 1e6, vmax=j_ves_gsfit.max() / 1e6)
            )
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax_top_right)
            cbar.set_label('Vessel Current Density (MA/m²)')

            # Make a bottom table that looks like this:
            # | Vessel Current | Total | Eig 1 | Eig 2 | ... | Eig n_ves_dof |
            # | RTGSFIT        | np.sum | calculated values for each eigenvector |
            # | GSFIT         | np.sum | calculated values for each eigenvector |
            # --------------------------------------------------------------------
            # | Rog      | Rog Name 1 | Rog Name 2 | ... | Rog Name n_rogowski_coils |
            # | Measured | Current 1 (MA) | Current 2 (MA) | ... | Current n_rogowski_coils (MA) |
            # | RTGSFIT  | Current 1 (MA) | Current 2 (MA) | ... | Current n_rogowski_coils (MA) |
            # | GSFIT    | Current 1 (MA) | Current 2 (MA) | ... | Current n_rogowski_coils (MA) |
            # --------------------------------------------------------------------
            # | PF Coil Currents | Coil Name 1 | Coil Name 2 | ... | Coil Name n_coils |
            # | RTGSFIT           | Current 1 (MA) | Current 2 (MA) | ... | Current n_coils (MA) |
            # | GSFIT             | Current 1 (MA) | Current 2 (MA) | ... | Current n_coils (MA) |
            rtgsfit_vessel_currents = np.zeros(n_ves_dof + 1)
            rtgsfit_vessel_currents[1:] = \
                [np.sum(np.abs(coef[ivc_coef_range[i_ivc_dof]] *
                        ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["current_distribution"]))
                for i_ivc_dof in range(n_ves_dof)]
            rtgsfit_vessel_currents[0] = np.sum(np.abs(rtgsfit_vessel_currents[1:]))
            rtgsfit_vessel_currents *= 1e-6
            gsfit_vessel_currents = np.zeros(n_ves_dof + 1)
            gsfit_vessel_currents[1:] = \
                [np.sum(np.abs(ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["calculated"] *
                        ivc_dict["dof"][f"eig_{i_ivc_dof+1:02d}"]["current_distribution"]))
                for i_ivc_dof in range(n_ves_dof)]
            gsfit_vessel_currents[0] = np.sum(np.abs(gsfit_vessel_currents[1:]))
            gsfit_vessel_currents *= 1e-6
            
            # Make a table with the vessel currents
            # Store data in {.1e} format
            rtgsfit_vessel_currents = [f'{current:.1e}' for current in rtgsfit_vessel_currents]
            gsfit_vessel_currents = [f'{current:.1e}' for current in gsfit_vessel_currents]
            rtgsfit_rog_pred = [f'{current:.2e}' for current in pred_meas[rogowski_range] * 1e-6]
            table_data1 = []
            table_data1.append([r"$\Sigma$|$I_{vessel}$|", "Total"] + [f'Eig {i+1}' for i in range(n_ves_dof)])
            table_data1.append(["RTGS [MA]"] + list(rtgsfit_vessel_currents))
            table_data1.append(["GSFIT [MA]"] + list(gsfit_vessel_currents))
            table_data2 = []
            table_data2.append(["Rog I [MA]"] + rogovski_names_table)
            table_data2.append(["Meas RT"] + list(rtgsfit_rog_meas))
            table_data2.append(["Meas GS"] + list(gsfit_rog_meas))
            table_data2.append(["RTGSFIT"] + list(rtgsfit_rog_pred))
            table_data2.append(["GSFIT"] + list(gsfit_rog_pred))
            table_data3 = []
            # table_data3.append(["PF I [MA]"] + list(rtgsfit_coil_names_extended))
            # table_data3.append(["RTGSFIT"] + [f'{current:.2e}' for current in rtgsfit_coil_currents_extended * 1e-6])
            # table_data3.append(["PF I [MA]"] + list(gsfit_coil_names))
            # table_data3.append(["GSFIT"] + [f'{current:.2e}' for current in gsfit_coil_currents * 1e-6])
            table_data3.append(["PF I [MA]"] + list(rtgsfit_coil_names))
            table_data3.append(["RTGSFIT"] + [f'{current:.2e}' for current in rtgsfit_coil_currents * 1e-6])
            table_data3.append(["GSFIT"] + [f'{current:.2e}' for current in gsfit_coil_currents_shorter * 1e-6])
            table1 = ax_bottom.table(cellText=table_data1,
                                    colLabels=None,
                                    cellLoc='center',
                                    loc='top')
            table1.auto_set_font_size(False)
            table1.set_fontsize(8)
            table1.scale(1.2, 1.2)
            table2 = ax_bottom.table(cellText=table_data2,
                                    colLabels=None,
                                    cellLoc='center',
                                    loc='center')
            table2.auto_set_font_size(False)
            table2.set_fontsize(8)
            table2.scale(1.2, 1.2)
            table3 = ax_bottom.table(cellText=table_data3,
                                    colLabels=None,
                                    cellLoc='center',
                                    loc='bottom')
            table3.auto_set_font_size(False)
            table3.set_fontsize(8)
            table3.scale(1.2, 1.2)

            plot_dir = "vessel_rogowski_coils"
            os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
            plot_name = f'vessel_rogowski_coils_iter_{i_iter}.png'
            fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, plot_name),
                        dpi=300,
                        bbox_inches='tight')
            plt.close("all")
        plot_vessel_rogowski_coils()

if __name__ == "__main__":

    # Load the output dictionary
    rtgsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'rtgsfit_output_dict.npy'),
                                  allow_pickle=True).item()
    
    gsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy'),
                                allow_pickle=True).item()

    # Call the contour plotting function
    plot_flux_loop_at_sensors(rtgsfit_output_dict, gsfit_output_dict)
    plot_b_at_sensors(rtgsfit_output_dict, gsfit_output_dict)
    plot_j_at_sensors(rtgsfit_output_dict, gsfit_output_dict)
