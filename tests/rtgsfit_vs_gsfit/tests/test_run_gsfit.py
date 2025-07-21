import os
import shutil
import subprocess

import pytest

from rtgsfit_vs_gsfit import cnst

def test_clone_gsfit_successfully(clone_gsfit_repo):
    """
    Test that the GSFIT repository is cloned successfully.
    """
    assert os.path.exists(cnst.GSFIT_PATH), f"GSFIT repository not found at {cnst.GSFIT_PATH}"

    # Check if the correct branch is checked out
    current_branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        cwd=cnst.GSFIT_PATH
    ).strip().decode('utf-8')
    
    assert current_branch == cnst.GSFIT_BRANCH, f"Expected branch {cnst.GSFIT_BRANCH}, but found {current_branch}"

@pytest.mark.usefixtures("clone_gsfit_repo")
def test_run_gsfit_and_save_to_mdsplus():
    """
    Run GSFIT, which saves to MDSplus, and verify the result.
    """

    from gsfit import Gsfit
    
    gsfit_controller = Gsfit(
        pulseNo=cnst.PULSE_NUM,
        run_name=cnst.RUN_NAME,
        run_description=cnst.RUN_DESCRIPTION,
        write_to_mds=True,
        pulseNo_write=cnst.PULSE_NUM_WRITE,
        settings_path="st40_setup_for_rtgsfit",
    )

    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["method"] = "arange"
    gsfit_controller.settings["GSFIT_code_settings.json"]["database_writer"]["method"] = "tokamak_energy_mdsplus"

    gsfit_controller.run()
