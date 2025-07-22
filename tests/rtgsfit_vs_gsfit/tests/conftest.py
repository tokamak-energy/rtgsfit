import os
import shutil
import subprocess

import pytest
import mdsthin
import numpy as np
import standard_utility

from rtgsfit_vs_gsfit import cnst, initialise_rtgsfit_node, replay_gsfit

# @pytest.fixture(scope="session")
# def clone_gsfit_repo():
#     """
#     1. Delete any existing GSFIT repository and clone a fresh copy.
#     This ensures that the tests always run against the latest version of GSFIT.

#     2. Clone the GSFIT repository.

#     3. Switch to the branch specified in src/cnst.py file.
#     """ 
#     if not os.path.exists(cnst.GSFIT_PATH):
#         subprocess.run(
#             ["git", "clone", cnst.GSFIT_URL, str(cnst.GSFIT_PATH)],
#             check=True
#         )
#     else:
#         print(f"GSFIT repository already exists at {cnst.GSFIT_PATH}. Deleting and re-cloning.")
#         shutil.rmtree(cnst.GSFIT_PATH)
#         subprocess.run(
#             ["git", "clone", cnst.GSFIT_URL, str(cnst.GSFIT_PATH)],
#             check=True
#         )

#     # Change to the GSFIT directory
#     os.chdir(cnst.GSFIT_PATH)

#     # Checkout the correct branch
#     subprocess.run(
#         ["git", "checkout", cnst.GSFIT_BRANCH],
#         check=True
#     )
#     print(f"GSFIT repository cloned and switched to branch {cnst.GSFIT_BRANCH} at {cnst.GSFIT_PATH}")


@pytest.fixture(scope="session")
def run_gsfit_and_save_dict():
# def run_gsfit_and_save_dict(clone_gsfit_repo):
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
