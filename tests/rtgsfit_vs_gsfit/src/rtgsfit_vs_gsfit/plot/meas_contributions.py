"""
Python module for plotting measured values from the sensors and fitted values from RTGSFIT and GSFIT.
"""

import os

import matplotlib.pyplot as plt
import mdsthin
import numpy as np

from rtgsfit_vs_gsfit import calc_meas, cnst

def plot_meas_contributions(rtgsfit_output_dict, gsfit_output_dict):
    """
    Function to plot the contributions to the measurements from each term in the
    coef vector.
    """

    def plot_contributions_iter_i_coef_i():
        """
        Plot the contributions to the measurements for a given iteration and coefficient.

        Make the figure consist of 2 rows and 2 columns

        The first axis will contain the flux loop contributions,
        the second axis will contain the bp probe contributions,
        the third axis will contain the rogowski coil contributions.
        the fourth axis will contain the regularization contributions.
        """

        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        # Increase the space between subplots
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        fig.suptitle(f"Iteration {i_iter}, Coefficient {coef_names[coef_i]}")

        # Plot flux loop contributions
        ax[0, 0].set_title("Flux Loop Contributions")
        ax[0, 0].set_xlabel("Flux Loop")
        ax[0, 0].set_ylabel("Contribution [Wb]")
        ax[0, 0].plot(flux_loop_names,
                   (rtgsfit_pred_i[flux_loop_idx_range] - pf_contrib[flux_loop_idx_range]) * 2 * np.pi,
                   label="i'th contribution")
        # ax[0].plot(flux_loop_names,
        #            rtgsfit_pred[flux_loop_idx_range] - pf_contrib[flux_loop_idx_range],
        #            label="RTGSFIT prediction")
        # ax[0].plot(flux_loop_names,
        #            gsfit_pred_fl - pf_contrib[flux_loop_idx_range],
        #            label="GSFIT prediction")
        # ax[0].plot(flux_loop_names,
        #            meas[flux_loop_idx_range] - pf_contrib[flux_loop_idx_range],
        #            label="Measured")
        ax[0, 0].legend()
        ax[0, 0].tick_params(axis='x', rotation=90)
        ax[0, 0].grid(True)

        ax[0, 1].plot(bp_probe_names,
                      rtgsfit_pred_i[bp_probe_idx_range] - pf_contrib[bp_probe_idx_range],
                      label="i'th contribution")
        ax[0, 1].set_title("BP Probe Contributions")
        ax[0, 1].set_xlabel("BP Probe")
        ax[0, 1].set_ylabel("Contribution [T]")
        ax[0, 1].legend()
        ax[0, 1].tick_params(axis='x', rotation=90)
        ax[0, 1].grid(True)

        ax[1, 0].plot(rogowski_names,
                     (rtgsfit_pred_i[rogowski_range] - pf_contrib[rogowski_range]),
                   label="i'th contribution")
        # ax[1, 0].plot(rogowski_names,
        #               pf_contrib[rogowski_range] * 1e-6,
        #               label="PF coil contribution")
        ax[1, 0].set_title("Rogowski Coil Contributions")
        ax[1, 0].set_xlabel("Rogowski Coil")
        ax[1, 0].set_ylabel("Contribution [MA]")
        ax[1, 0].legend()
        ax[1, 0].tick_params(axis='x', rotation=90)
        ax[1, 0].grid(True)

        ax[1, 1].set_title("Regularization Contributions")
        ax[1, 1].set_xlabel("Regularization Term")
        ax[1, 1].set_ylabel("Contribution")
        ax[1, 1].plot(reg_names,
                      rtgsfit_pred_i[reg_range] - pf_contrib[reg_range],
                        label="i'th contribution")
        ax[1, 1].plot(reg_names,
                      pf_contrib[reg_range],
                      label="PF coil contribution")
        ax[1, 1].plot(reg_names,
                      coef_i_array[-len(reg_range):],
                      '--',
                      label="Coefficient contribution")
        ax[1, 1].plot(reg_names,
                      meas[reg_range],
                      label="Measured")
        ax[1, 1].legend()
        ax[1, 1].tick_params(axis='x', rotation=90)
        ax[1, 1].grid(True)

        plot_dir = "meas_contributions"
        sub_dir = f"iter_{i_iter}"
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir), exist_ok=True)
        os.makedirs(os.path.join(cnst.PLOTS_DIR, plot_dir, sub_dir), exist_ok=True)
        plot_name = f'meas_contributions_iter_{i_iter}_coef_{coef_names[coef_i]}.png'
        fig.savefig(os.path.join(cnst.PLOTS_DIR, plot_dir, sub_dir, plot_name),
                    dpi=300,
                    bbox_inches='tight')
        plt.close("all")


    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_pls = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_PLS").data()
        n_coil = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COIL").data()
        n_meas = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_MEAS").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COEF").data()
        n_flux_loops = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_F_LOOPS").data()
        n_bp_probes = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_BP_PROBES").data()
        n_rog_coils = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_ROG_COILS").data()
        g_meas_coil = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT.GREENS:MEAS_COIL").data()
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:SENS_NAMES").data()
        g_meas_coil = g_meas_coil.reshape((n_meas, n_coil)).T

    ivc_dict = gsfit_output_dict["IVC"]
    n_eigs = len(ivc_dict["dof"])
    flux_loop_names = gsfit_output_dict["flux_loops"]["names"]
    coef_names = ["p_prime", "ff_prime", "delta_z"] + \
                 ['BVLBCASE', 'BVLTCASE', 'DIVPSRB', 'DIVPSRT', 'HFSPSRB', 'HFSPSRT'] + \
                 list(ivc_dict["dof"].keys()) \
               + ["OVC"]
    flux_loop_idx_range = np.arange(n_flux_loops)
    bp_probe_idx_range = np.arange(n_flux_loops, n_flux_loops + n_bp_probes)
    bp_probe_names = gsfit_output_dict["bp_probes"]["names"]
    rogowski_range = np.arange(n_flux_loops + n_bp_probes, n_flux_loops + n_bp_probes + n_rog_coils)

    rogowski_names = sens_names[rogowski_range]
    rogowski_names = [name.replace("I_ROG_", "") for name in rogowski_names]

    reg_names = list(ivc_dict["dof"].keys()) + ["OVC"]
    reg_range = np.arange(n_flux_loops + n_bp_probes + n_rog_coils, n_meas)

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cnst.PULSE_NUM_WRITE)
        fl_include = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.FLOOP:INCLUDE").data() == 1
        gsfit_pred_fl = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.FLOOP:CVALUE").data()[0, fl_include]
        bp_include = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.BPPROBE:INCLUDE").data() == 1
        gsfit_pred_bp = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.BPPROBE:CVALUE").data()[0, bp_include]
        gsfit_pred_rog = conn.get(f"\\GSFIT::TOP.{cnst.RUN_NAME}.CONSTRAINTS.ROG:CVALUE").data()[0, :]

    # for i_iter in range(cnst.N_ITERS + 1):
    for i_iter in [cnst.N_ITERS - 1]:
        meas = rtgsfit_output_dict["meas"][i_iter, :]
        coef  = rtgsfit_output_dict["coef"][i_iter, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][0, :]
        rtgsfit_pred = calc_meas.calculate_predicted_measurements(
            rtgsfit_output_dict["flux_norm"][i_iter, :],
            coef,
            coil_curr
        )
        pf_contrib = coil_curr @ g_meas_coil
        for coef_i in range(len(coef)):
            coef_i_array = np.zeros_like(coef)
            coef_i_array[coef_i] = coef[coef_i]
            rtgsfit_pred_i = calc_meas.calculate_predicted_measurements(
                rtgsfit_output_dict["flux_norm"][i_iter, :],
                coef_i_array,
                rtgsfit_output_dict["coil_curr"][i_iter, :]
            )
            plot_contributions_iter_i_coef_i()

if __name__ == "__main__":

    rtgsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'rtgsfit_output_dict.npy'), allow_pickle=True).item()
    gsfit_output_dict = np.load(os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy'), allow_pickle=True).item()

    plot_meas_contributions(rtgsfit_output_dict, gsfit_output_dict)