"""
This module contains the function that initialises the MDSplus node for RTGSFIT which the
Makefile in the RTGSFIT repository uses to generate the constants.c file
before compiling the RTGSFIT code.
"""

import numpy as np
from gsfit import Gsfit
from rtgsfit_vs_gsfit import cnst

def initialise_rtgsfit_node():
    """
    Initialise the RTGSFIT MDSplus node by running GSFIT with the specified settings.
    This function sets up the GSFIT controller, modifies the settings, and runs the analysis
    to generate the necessary data for RTGSFIT.
    """

    # Construct the GSFit object; using the "st40_setup_for_rtgsfit" settings
    gsfit_controller = Gsfit(
        pulseNo=cnst.PULSE_NUM,
        run_name=cnst.RUN_NAME,
        run_description= cnst.RUN_DESCRIPTION,
        write_to_mds=True,
        pulseNo_write=cnst.PULSE_NUM_WRITE,
        settings_path=cnst.settings_path,
    )

    # Change the analysis_name, so that GSFit writes into RT-GSFit MDSplus tree
    gsfit_controller.analysis_name = "RTGSFIT"

    # Add a list of signals to be read using PCS formatting
    gsfit_controller.results["PRESHOT"]["COIL_SIGNALS"] = np.array(
        [
            "I_BVL_PSU",
            "I_BVUB_PSU",
            "I_BVUT_PSU",
            "I_DIV_PSU",
            "I_MCVC_PSU",
            "I_PSH_PSU",
            "I_ROG_MCWIRE",
            "I_SOL_PSU",
        ]
    )
    # fmt: off
    coil_matrix = np.array(
        [
            # BVL_PSU, BVUB_PSU, BVUT_PSU, DIV_PSU, MCVC_PSU, PSH_PSU, ROG_MCWIRE, SOL_PSU
            [1.0,      0.0,      0.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVL coil
            [0.0,      1.0,      0.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVUB coil
            [0.0,      0.0,      1.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVUT coil
            [0.0,      0.0,      0.0,      1.0,     0.0,      0.0,     0.0,        0.0],  # DIV coil
            [0.0,      0.0,      0.0,      0.0,     1.0,      0.0,     1.0,        0.0],  # MCB coil
            [0.0,      0.0,      0.0,      0.0,     1.0,      0.0,     0.0,        0.0],  # MCT coil
            [0.0,      0.0,      0.0,      0.0,     0.0,      1.0,     0.0,        0.0],  # PSH coil
            [0.0,      0.0,      0.0,      0.0,     0.0,      0.0,     0.0,        1.0],  # SOL coil
        ]
    )
    # fmt: on
    gsfit_controller.results["PRESHOT"]["COIL_MATRIX"] = coil_matrix

    # Run
    gsfit_controller.run()

if __name__ == "__main__":
    initialise_rtgsfit_node()
