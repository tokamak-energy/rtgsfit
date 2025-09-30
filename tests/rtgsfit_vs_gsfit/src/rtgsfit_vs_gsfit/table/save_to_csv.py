"""
Combines Pandas Dataframes together into a CSV and formats them to make them more readable.
"""

import os

import numpy as np
import pandas as pd

def save_dfs_to_csv(iterations: np.ndarray, cfg: dict) -> None:
    """
    Save DataFrames to a CSV with titles between them and spaces.
    """
    from rtgsfit_vs_gsfit.table.dataframes import ivc_df, pf_coil_gsfit_df, \
        pf_coil_rtgsfit_df, plasma_current_df, rog_df

    this_plot_dir = os.path.join(cfg["plots_this_run_dir"], "dataframes")
    os.makedirs(this_plot_dir, exist_ok=True)

    df_pf_coil_rtgsfit = pf_coil_rtgsfit_df(cfg)
    df_pf_coil_rtgsfit.columns = [f"{col:>8}" for col in df_pf_coil_rtgsfit.columns]
    pf_coil_gsfit_df = pf_coil_gsfit_df(cfg)
    pf_coil_gsfit_df.columns = [f"{col:>8}" for col in pf_coil_gsfit_df.columns]

    for iteration in iterations:
        filename = f"dataframes_{cfg['pulse_num']}_{cfg['run_name']}_{iteration:02d}.csv"
        # Open CSV write title then save dataframes
        with open(os.path.join(this_plot_dir, filename), "w") as f:
            f.write(f"Pulse number: {cfg['pulse_num']}\n")
            f.write(f"Run name: {cfg['run_name']}\n")
            f.write(f"Time: {cfg['time']:.1e} s\n")
            f.write(f"Iteration: {iteration}\n")
            f.write("\n")
            f.write(f"IVC Eigenvalues:\n")
            df_ivc = ivc_df(iteration, cfg)
            df_ivc.columns = [f"{col:>8}" for col in df_ivc.columns]
            df_ivc.to_csv(f, index=True, float_format="%8.1e")
            f.write("\n")
            f.write(f"Rogowski Coil Measurements and Predictions [A]:\n")
            df_rog = rog_df(iteration, cfg)
            df_rog.columns = [f"{col:>8}" for col in df_rog.columns]
            df_rog.to_csv(f, index=True, float_format="%8.1e")
            f.write("\n")
            f.write(f"RTGSFIT PF Coil Currents [A]:\n")
            df_pf_coil_rtgsfit.to_csv(f, index=False, float_format="%8.1e")
            f.write("\n")
            f.write(f"GSFIT PF Coil Currents [A]:\n")
            pf_coil_gsfit_df.to_csv(f, index=False, float_format="%8.1e")
            f.write("\n")
            f.write(f"Plasma Current [A]:\n")
            df_plasma_current = plasma_current_df(iteration, cfg)
            df_plasma_current.columns = [f"{col:>8}" for col in df_plasma_current.columns]
            df_plasma_current.to_csv(f, index=False, float_format="%8.1e")

if __name__ == "__main__":

    import time
    from rtgsfit_vs_gsfit import config_loader

    cfg = config_loader.load_and_prepare_config()

    iterations = np.array([0, 1, 8])

    start_time = time.time()
    save_dfs_to_csv(iterations=iterations, cfg=cfg)
    end_time = time.time()
    print(f"Save_df_to_csv took {end_time - start_time:.2e} seconds")
