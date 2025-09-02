"""
Full integration test checking that the RTGSFIT results agree 
with the GSFIT results to within a specified tolerance.
"""

import logging
import os

import mdsthin
import numpy as np
import pytest

from rtgsfit_vs_gsfit import config_loader, \
                             replay_gsfit, replay_rtgsfit, rtgsfit_compile_setup
from rtgsfit_vs_gsfit.plot import examples
from rtgsfit_vs_gsfit.table import save_to_csv

test_cases = [
    (13_343, 0.030),
    (13_345, 0.030),
    (13_346, 0.030),
]

@pytest.mark.parametrize("pulse_num,time", test_cases)
def test_rtgsfit_vs_gsfit_consistency(pulse_num, time):
    """
    Runs RTGSFIT and GSFIT for a given pulse and time, compares output.
    """

    def check_psi(cfg: dict,
                  rtol: float,
                  atol: float) -> None:
        """
        Check that the psi values from RTGSFIT and GSFIT agree to within a tolerance.
        at the last iteration
        """
        iteration = cfg["n_iters"] - 1

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_write"])
            n_r = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:N_R").data()
            n_z = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:N_Z").data()
        rtgsfit_output_dict =np.load(cfg["rtgsfit_output_dict_path"],
                                     allow_pickle=True).item()
        psi_rtgsfit = rtgsfit_output_dict["flux_total"][iteration].reshape(n_z, n_r)

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            psi_gsfit = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.TWO_D:PSI").data()[0, :, :]

        np.testing.assert_allclose(psi_rtgsfit, psi_gsfit,
                                   rtol=rtol, atol=atol)
        
    def check_psi_meas(cfg: dict,
                       rtol_meas: float,
                       atol_meas: float,
                       rtol_pred: float,
                       atol_pred: float) -> None:
        """
        Check that measured values of psi at the flux loops from RTGSFIT and GSFIT agree to within a tolerance.
        at the last iteration
        """
        from rtgsfit_vs_gsfit import rtgsfit_pred_meas

        iteration = cfg["n_iters"] - 1

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_write"])
            sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:SENS_NAMES").data()
        flux_loop_indices = np.where(np.char.startswith(sens_names, "PSI_FLOOP_"))[0]
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                      allow_pickle=True).item()
        flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
        coef = rtgsfit_output_dict["coef"][iteration, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
        pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)
        pred_meas_rtgsfit = pred_meas_rtgsfit[flux_loop_indices]

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            fl_include = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:INCLUDE").data() == 1
            gsfit_fl_meas = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:MVALUE").data()[0, fl_include]
            gsfit_fl_pred = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:CVALUE").data()[0, fl_include]

        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_fl_meas,
                                   rtol=rtol_meas, atol=atol_meas)
        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_fl_pred,
                                   rtol=rtol_pred, atol=atol_pred)
        
    def check_bp_meas(cfg: dict,
                      rtol_meas: float,
                      atol_meas: float,
                      rtol_pred: float,
                      atol_pred: float) -> None:
        """
        Check that the measured values of the magnetic field at the
        BP probes from RTGSFIT and GSFIT agree to within a
        tolerance.
        """
        from rtgsfit_vs_gsfit import rtgsfit_pred_meas

        iteration = cfg["n_iters"] - 1

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_write"])
            sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:SENS_NAMES").data()
        bp_probe_indices = np.where(np.char.startswith(sens_names, "B_BPPROBE_"))[0]
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
        coef = rtgsfit_output_dict["coef"][iteration, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
        pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)
        pred_meas_rtgsfit = pred_meas_rtgsfit[bp_probe_indices]

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            bp_include = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:INCLUDE").data() == 1
            gsfit_bp_meas = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:MVALUE").data()[0, bp_include]
            gsfit_bp_pred = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:CVALUE").data()[0, bp_include]

        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_bp_meas,
                                   rtol=rtol_meas, atol=atol_meas)
        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_bp_pred,
                                   rtol=rtol_pred, atol=atol_pred)

    def check_ivc_eigs(cfg: dict,
                       rtol: float,
                       atol: float) -> None:
        """
        Check the IVC eigenvalues agree to within the specified tolerances.
        """

        iteration = cfg["n_iters"] - 1

        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        coef = rtgsfit_output_dict["coef"][iteration, :]
        with open(cfg["coef_names_path"], "r") as f:
            coef_names = [line.strip() for line in f]
        ivc_indices = [i for i, name in enumerate(coef_names) if name.startswith("eig_")]
        eigs_rtgsfit =  coef[ivc_indices]

        ivc_dict = np.load(cfg["ivc_dict_path"], allow_pickle=True).item()
        n_eigs = ivc_dict["current_distributions"].shape[0]
        eigs_gsfit = np.zeros(n_eigs)
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            for eig_num in range(1, n_eigs + 1):
                eigs_gsfit[eig_num - 1] = \
                    conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.PASSIVES.IVC.DOF:EIG_{eig_num:02d}")[0]
        np.testing.assert_allclose(eigs_rtgsfit, eigs_gsfit,
                                   rtol=rtol, atol=atol)

    def check_ovc_current(cfg: dict,
                          rtol: float,
                          atol: float) -> None:
        """
        Check the OVC current measurements agree to within the specified tolerances.
        """

        from rtgsfit_vs_gsfit.plot.components import ovc_current_rtgsfit, ovc_current_gsfit

        iteration = cfg["n_iters"] - 1

        np.testing.assert_allclose(ovc_current_rtgsfit(iteration, cfg),
                                   ovc_current_gsfit(cfg),
                                   rtol=rtol, atol=atol)

    def check_rog_meas(cfg: dict,
                       rtol_meas: float,
                       atol_meas: float,
                       rtol_pred: float,
                       atol_pred: float):
        """
        Check the Rogowski measurements agree to within the specified tolerances.
        """
        
        from rtgsfit_vs_gsfit import rtgsfit_pred_meas

        iteration = cfg["n_iters"] - 1

        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
        coef = rtgsfit_output_dict["coef"][iteration, :]
        coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
        pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)
        with open(cfg["meas_names_path"], "r") as f:
            meas_names = [line.strip() for line in f]
        rog_indices = np.zeros(len(cfg["rogowski_names"]), dtype=int)
        for i, rog_name in enumerate(cfg["rogowski_names"]):
            for j, meas_name in enumerate(meas_names):
                if rog_name == meas_name:
                    rog_indices[i] = j
        pred_meas_rtgsfit = pred_meas_rtgsfit[rog_indices]

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            rog_names_gsfit = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:NAME")
            gsfit_rog_meas = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:MVALUE")[0]
            gsfit_rog_pred = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:CVALUE")[0]
        rog_indices = np.zeros(len(cfg["rogowski_names"]), dtype=int)
        for i, rog_name in enumerate(cfg["rogowski_names"]):
            for j, rog_name_gsfit in enumerate(rog_names_gsfit):
                if rog_name in rog_name_gsfit:
                    rog_indices[i] = j
        gsfit_rog_pred = gsfit_rog_pred[rog_indices]
        gsfit_rog_meas = gsfit_rog_meas[rog_indices]

        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_rog_meas,
                                   rtol=rtol_meas, atol=atol_meas)
        np.testing.assert_allclose(pred_meas_rtgsfit, gsfit_rog_pred,
                                   rtol=rtol_pred, atol=atol_pred)

    def check_regression(cfg: dict,
                         rtol: float,
                         atol: float) -> None:
        """
        Verify that the final iteration results from RTGSFIT match the 
        reference values stored in the corresponding 
        `rtgsfit_output_dict_*_regression.npy` file.

        This function ensures that recent code changes have not caused 
        significant deviations in the RTGSFIT output.
        """

        iteration = cfg["n_iters"] - 1

        rtgsfit_dict = np.load(cfg["rtgsfit_output_dict_path"],
                               allow_pickle=True).item()
        regress_path = cfg["rtgsfit_output_dict_path"].replace(".npy", "_regression.npy")
        regress_dict = np.load(regress_path, allow_pickle=True).item()

        for key in rtgsfit_dict.keys():
            np.testing.assert_allclose(rtgsfit_dict[key], regress_dict[key],
                                       rtol=rtol, atol=atol,
                                       err_msg=f"Key: {key}")

    run_name = f"t{int(time*1e3):03d}ms"
    cfg = config_loader.load_and_prepare_config(
        run_name=run_name,
        pulse_num=pulse_num)
    cfg["time"] = time

    logging.info(f"Running RTGSFIT vs GSFIT consistency "
                 f"test for pulse {pulse_num} at time {time}s with "
                 f"run name {run_name}")
    
    # Replay RTGSFIT
    logging.info(f"Clearing RTGSFIT node...")
    rtgsfit_compile_setup.rtgsfit_mds_nodeclear(cfg)
    logging.info(f"RTGSFIT node cleared.")
    logging.info(f"Initialising RTGSFIT node...")
    rtgsfit_compile_setup.initialise_rtgsfit_node(cfg)
    logging.info(f"RTGSFIT node initialised.")
    logging.info(f"Compiling RTGSFIT...")
    rtgsfit_compile_setup.compile_rtgsfit(cfg)
    logging.info(f"RTGSFIT compiled.")
    logging.info(f"Replaying RTGSFIT...")
    replay_rtgsfit.replay_rtgsfit(cfg)
    logging.info(f"RTGSFIT replay completed.")

    # Replay GSFIT
    logging.info(f"Replaying GSFIT...")
    replay_gsfit.replay_gsfit(cfg)
    logging.info(f"GSFIT replay completed.")

    # Plot results
    iterations = [0, 1, 8]
    logging.info(f"Plotting results...")
    examples.psi_fluxloop_bp(cfg,
                             iterations=iterations)
    examples.psi_ivc_passives(cfg,
                              iterations=iterations)
    examples.eigenvalue_bar_chart(cfg,
                                  iterations=iterations)
    save_to_csv.save_dfs_to_csv(iterations, cfg)
    logging.info(f"Results plotted.")

    tolerances = {
        "psi": {"rtol": 1e-3, "atol": 1e-3},
        "psi_meas": {"rtol_meas": 5e-2, "atol_meas": 5e-2, "rtol_pred": 1e-3, "atol_pred": 1e-3},
        "bp_meas": {"rtol_meas": 5e-2, "atol_meas": 5e-2, "rtol_pred": 1e-3, "atol_pred": 1e-3},
        "ivc_eigs": {"rtol": 1e-3, "atol": 5e-1},
        "ovc_current": {"rtol": 2e-2, "atol": 1e-2},
        "rog_meas": {"rtol_meas": 5e-2, "atol_meas": 5e2, "rtol_pred": 5e-2, "atol_pred": 1e1},
        "regression": {"rtol": 1e-8, "atol": 1e-8}
    }
    logging.info(f"Checking psi...")
    check_psi(cfg, **tolerances["psi"])
    logging.info(f"Checking psi_meas...")
    check_psi_meas(cfg, **tolerances["psi_meas"])
    logging.info(f"Checking bp_meas...")
    check_bp_meas(cfg, **tolerances["bp_meas"])
    logging.info(f"Checking ivc_eigs...")
    check_ivc_eigs(cfg, **tolerances["ivc_eigs"])
    logging.info(f"Checking ovc_current...")
    check_ovc_current(cfg, **tolerances["ovc_current"])
    logging.info(f"Checking rog_meas...")
    check_rog_meas(cfg, **tolerances["rog_meas"])
    logging.info(f"Checking regression...")
    check_regression(cfg, **tolerances["regression"])
    logging.info(f"Results checked successfully.")
