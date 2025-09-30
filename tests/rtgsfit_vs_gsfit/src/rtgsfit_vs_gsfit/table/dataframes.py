"""
This module contains routines for making DataFrames from simulation data.
"""

import mdsthin
import numpy as np
import pandas as pd

def ivc_df(iteration: int, cfg: dict) -> pd.DataFrame:
    """
    Create a DataFrame for the IVC eigenvalues and total current data.
    """
    def ivc_rtgsfit_row(iteration: int,
                        cfg: dict) -> np.ndarray:
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        coef = rtgsfit_output_dict["coef"][iteration, :]
        with open(cfg["coef_names_path"], "r") as f:
            coef_names = [line.strip() for line in f]
        ivc_indices = [i for i, name in enumerate(coef_names) if name.startswith("eig_")]
        return coef[ivc_indices]
    
    def ivc_gsfit_row(cfg: dict) -> np.ndarray:
        ivc_dict = np.load(cfg["ivc_dict_path"], allow_pickle=True).item()
        n_eigs = ivc_dict["current_distributions"].shape[0]
        eigenvalues = np.zeros(n_eigs)
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            for eig_num in range(1, n_eigs + 1):
                eigenvalues[eig_num - 1] = \
                    conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.PASSIVES.IVC.DOF:EIG_{eig_num:02d}")[0]
        return eigenvalues
        
    rtgsfit_vals = ivc_rtgsfit_row(iteration, cfg)
    gsfit_vals = ivc_gsfit_row(cfg)

    n_eigs = len(rtgsfit_vals)
    cols = [f"eig_{i:02d}" for i in range(1, n_eigs + 1)]

    df = pd.DataFrame([rtgsfit_vals, gsfit_vals], 
                      index=["RTGSFIT", "GSFIT  "],
                      columns=cols)
    df.index.name = " " * len("RTGSFIT")
    pd.set_option("display.float_format", "{:.2e}".format)

    return df

def rog_df(iteration: int, cfg: dict) -> pd.DataFrame:
    """
    Create a DataFrame for the Rogowski coil measurements
    and predicted values from RTGSFIT and GSFIT.
    """

    from rtgsfit_vs_gsfit import rtgsfit_pred_meas

    def rog_meas_row(cfg: dict) -> np.ndarray:
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            rog_names_gsfit = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:NAME")
            meas_values = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:MVALUE")[0]
        rog_indices = np.zeros(len(cfg["rogowski_names"]), dtype=int)
        for i, rog_name in enumerate(cfg["rogowski_names"]):
            for j, rog_name_gsfit in enumerate(rog_names_gsfit):
                if rog_name in rog_name_gsfit:
                    rog_indices[i] = j
        return meas_values[rog_indices]

    def rog_rtgsfit_row(iteration: int, cfg: dict) -> np.ndarray:
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
        return pred_meas_rtgsfit[rog_indices]
    
    def rog_gsfit_row(cfg: dict) -> np.ndarray:
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            rog_names_gsfit = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:NAME")
            calc_values = \
                conn.get("\\GSFIT::TOP." + cfg["run_name"] + ".CONSTRAINTS.ROG:CVALUE")[0]
        rog_indices = np.zeros(len(cfg["rogowski_names"]), dtype=int)
        for i, rog_name in enumerate(cfg["rogowski_names"]):
            for j, rog_name_gsfit in enumerate(rog_names_gsfit):
                if rog_name in rog_name_gsfit:
                    rog_indices[i] = j
        return calc_values[rog_indices]

    meas_vals = rog_meas_row(cfg)
    rtgsfit_vals = rog_rtgsfit_row(iteration, cfg)
    gsfit_vals = rog_gsfit_row(cfg)

    n_eigs = len(rtgsfit_vals)
    cols = [f"eig_{i:02d}" for i in range(n_eigs)]

    df = pd.DataFrame([meas_vals, rtgsfit_vals, gsfit_vals],
                      index=["MEASURED",
                             "RTGSFIT ",
                             "GSFIT   "],
                      columns=cfg["rogowski_names"])
    df.index.name = " " * len("MEASURED")
    pd.set_option("display.float_format", "{:.2e}".format)

    return df

def pf_coil_rtgsfit_df(cfg: dict) -> pd.DataFrame:
    """
    Create a DataFrame for the PF coil currents in both 
    GSFIT and RTGSFIT.

    Doesn't change with iterations.
    """

    def pf_coil_rtgsfit_row(cfg: dict) -> np.ndarray:
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        return rtgsfit_output_dict["coil_curr"][0, :]

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        coil_names = conn.get("\\RTGSFIT::TOP." + cfg["run_name"] + ".PRESHOT:COIL_NAMES")

    df = pd.DataFrame([pf_coil_rtgsfit_row(cfg)], columns=coil_names)

    return df

def pf_coil_gsfit_df(cfg: dict) -> pd.DataFrame:
    """
    Create a DataFrame for the PF coil currents in GSFIT.
    
    Doesn't change with iterations.
    """
    import json

    # Read the JSON file saved at cfg["gsfit_pf_coils_path"]
    # Put into dataframe where each columns is the key and its just one row
    with open(cfg["gsfit_pf_coils_path"], "r") as f:
        coil_currents = json.load(f)

    df = pd.DataFrame([coil_currents], columns=coil_currents.keys())

    return df

def plasma_current_df(iteration: int, cfg: dict) -> pd.DataFrame:
    """
    Create a DataFrame for the plasma current from RTGSFIT and GSFIT.
    Just one row with two columns. Literally just the plasma current.
    so just a 1x2 dataframe.
    """

    def plasma_current_rtgsfit(iteration: int, cfg: dict) -> np.ndarray:
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        return rtgsfit_output_dict["plasma_current"][iteration, 0]
    
    def plasma_current_gsfit(cfg: dict) -> np.ndarray:
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            plasma_current = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.GLOBAL:IP")[0]
        return plasma_current

    rtgsfit_val = plasma_current_rtgsfit(iteration, cfg)
    gsfit_val = plasma_current_gsfit(cfg)

    df = pd.DataFrame([[rtgsfit_val, gsfit_val]],
                    columns=["RTGSFIT", "GSFIT"])
    pd.set_option("display.float_format", "{:.2e}".format)

    return df

if __name__ == "__main__":

    import time

    from rtgsfit_vs_gsfit import config_loader

    cfg = config_loader.load_and_prepare_config()

    start_time = time.time()
    df = ivc_df(8, cfg)
    end_time = time.time()
    print(f"Time taken to create IVC DataFrame: {end_time - start_time:.2e} seconds")