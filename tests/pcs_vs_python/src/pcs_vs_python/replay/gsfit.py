"""
This module contains the function to run GSFIT and save the results to MDS+.
"""

from gsfit import Gsfit
import json
import mdsthin
import numpy as np

def replay_gsfit(cfg: dict):
    """
    Run GSFIT, then save the results to a dictionary.

    Input:
        cfg (dict): Configuration dictionary containing key settings.

    Output:
        Saves the results to MDS+.
    """
   
    gsfit_controller = Gsfit(
        pulseNo=cfg["pulse_num"],
        run_name=cfg["run_name_replay"],
        run_description=cfg["run_description_replay"],
        write_to_mds=True,
        pulseNo_write=cfg["pulse_num_replay"],
        settings_path=cfg["settings_path"],
    )

    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["method"] = "linspace"
    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["linspace"]["time_start"] = cfg["t_min"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["linspace"]["time_end"] = cfg["t_max"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["linspace"]["n_time"] = cfg["n_t"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["database_writer"]["method"] = "tokamak_energy_mdsplus"
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_r"] = cfg["n_r"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_z"] = cfg["n_z"]

    print("Run name:", gsfit_controller.run_name)
    print("Pulse number:", gsfit_controller.pulseNo)
    print("Pulse number to write:", gsfit_controller.pulseNo_write)
    print("Time slice:", gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["user_defined"])
    gsfit_controller.run()

if __name__ == "__main__":

    from pcs_vs_python import config_loader

    cfg = config_loader.load_and_prepare_config()
    replay_gsfit(cfg)
