"""
Module for interfacing with MDS+ and building RTGSFIT.

This module provides utility functions to:
- `rtgsfit_mds_nodeclear`: Delete existing MDS+ node for RTGSFIT.
- `initialise_rtgsfit_node`: Generate and populate the MDS+ node with the required data.
- `compile_rtgsfit`: Compile the RTGSFIT code using its `Makefile`.
"""
import os
import subprocess

import MDSplus
import numpy as np

from gsfit import Gsfit

def write_mds_node(cfg: dict):
    """
    Initialise the RTGSFIT MDSplus PRESHOT node by running GSFIT with the specified settings.
    This function sets up the GSFIT controller, modifies the settings, and runs the analysis
    to generate the necessary data for RTGSFIT. 
    """

    # Construct the GSFit object; using the "st40_setup_for_rtgsfit" settings
    gsfit_controller = Gsfit(
        pulseNo=cfg['pulse_num'],
        run_name=cfg['run_name_preshot'],
        run_description=cfg['run_description_preshot'],
        write_to_mds=True,
        pulseNo_write=cfg['pulse_num_preshot'],
        settings_path=cfg['settings_path'],
    )

    # Change the analysis_name, so that GSFit writes into RT-GSFit MDSplus tree
    gsfit_controller.analysis_name = "RTGSFIT"

    # Add a list of signals to be read using PCS formatting
    gsfit_controller.results["PRESHOT"]["COIL_SIGNALS"] = np.array(cfg["coil_signals"])
    gsfit_controller.results["PRESHOT"]["COIL_MATRIX"] = np.array(cfg["coil_matrix"])

    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_r"] = cfg["n_r"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_z"] = cfg["n_z"]

    # Run
    gsfit_controller.run()

def compile_rtgsfit(cfg: dict):
    """
    Clean and then compile RTGSFIT.
    """

    # Clean the RTGSFIT source directory
    os.chdir(cfg['rtgsfit_src_path'])
    subprocess.run(["make", "clean"], cwd=cfg['rtgsfit_src_path'], check=True)
    # Check if the object files are removed
    object_files = [f for f in os.listdir(cfg['rtgsfit_src_path']) if f.endswith('.o')]
    assert not object_files, "Object files were not removed during cleaning."
    # Check constants.c was removed
    constants_c_path = os.path.join(cfg['rtgsfit_src_path'], 'constants.c')
    assert not os.path.exists(constants_c_path), f"constants.c was not removed: {constants_c_path}"

    # Compile RTGSFIT
    os.chdir(cfg['rtgsfit_src_path'])
    subprocess.run(
        f"make SHOT={cfg['pulse_num_preshot']} RUN_NAME={cfg['run_name_preshot']} DEBUG=1",
        cwd=cfg['rtgsfit_src_path'],
        check=True,
        shell=True
    )

if __name__ == "__main__":

    from replay_rtgsfit_st40 import config_loader

    cfg = config_loader.load_and_prepare_config()
    write_mds_node(cfg)
    compile_rtgsfit(cfg)