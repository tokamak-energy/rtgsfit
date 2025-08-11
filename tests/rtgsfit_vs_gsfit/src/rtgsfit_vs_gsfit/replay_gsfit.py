"""
This module contains the function to run GSFIT and save the results to MDS+.
"""

from gsfit import Gsfit

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
        run_name=cfg["run_name"],
        run_description=cfg["run_description"],
        write_to_mds=True,
        pulseNo_write=cfg["pulse_num_write"],
        settings_path=cfg["settings_path"],
    )

    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["method"] = "user_defined"
    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["user_defined"] = [cfg["time"]]
    gsfit_controller.settings["GSFIT_code_settings.json"]["database_writer"]["method"] = "tokamak_energy_mdsplus"
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_r"] = cfg["n_r"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_z"] = cfg["n_z"]

    print("Run name:", gsfit_controller.run_name)
    print("Pulse number:", gsfit_controller.pulseNo)
    print("Pulse number to write:", gsfit_controller.pulseNo_write)
    print("Time slice:", gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["user_defined"])
    gsfit_controller.run()

    cfg["gsfit_replayed"] = True

if __name__ == "__main__":

    from rtgsfit_vs_gsfit import config_loader

    cfg = config_loader.load_and_prepare_config()
    replay_gsfit(cfg)
