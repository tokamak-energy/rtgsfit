import os
import shutil
import subprocess

import pytest
import mdsthin
import numpy as np
import standard_utility

from rtgsfit_vs_gsfit import cnst, initialise_rtgsfit_node, replay_gsfit

@pytest.fixture(scope="session")
def run_gsfit_and_save_dict():
    """
    Run GSFIT, then save the results to a dictionary.
    """

    replay_gsfit.replay_gsfit()

@pytest.fixture(scope="session")
def rtgsfit_mds_nodegen():
    """
    Delete the existing MDSplus node and generate a new one.
    """

    conn = mdsthin.Connection('smaug')
    conn.tcl(f'edit RTGSFIT/shot={cnst.PULSE_NUM_WRITE}')
    conn.tcl(f'DELETE NODE \\\\RTGSFIT::TOP:{cnst.RUN_NAME} /CONFIRM')
    conn.tcl('write')

    standard_utility.create_script_nodes(
        script_name='RTGSFIT',
        pulseNo_write=cnst.PULSE_NUM_WRITE,
        run_name=cnst.RUN_NAME,
    )

    initialise_rtgsfit_node.initialise_rtgsfit_node()
