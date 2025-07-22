import os
import shutil
import subprocess

import numpy as np
import pytest

from rtgsfit_vs_gsfit import cnst

# def test_clone_gsfit_successfully(clone_gsfit_repo):
#     """
#     Test that the GSFIT repository is cloned successfully.
#     """
#     assert os.path.exists(cnst.GSFIT_PATH), f"GSFIT repository not found at {cnst.GSFIT_PATH}"

#     # Check if the correct branch is checked out
#     current_branch = subprocess.check_output(
#         ["git", "rev-parse", "--abbrev-ref", "HEAD"],
#         cwd=cnst.GSFIT_PATH
#     ).strip().decode('utf-8')
    
#     assert current_branch == cnst.GSFIT_BRANCH, f"Expected branch {cnst.GSFIT_BRANCH}, but found {current_branch}"

def test_run_gsfit_and_save_dict(run_gsfit_and_save_dict):
    """
    Test that GSFIT runs and saves the output dictionary correctly.
    """
    output_file = os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy')
    assert os.path.exists(output_file), f"Output file {output_file} not found."

    # Load the output dictionary
    output_dict = np.load(output_file, allow_pickle=True).item()
    
    # Check if the expected keys are present in the output dictionary
    expected_keys = ["grid", "two_d", "p_boundary"]
    for key in expected_keys:
        assert key in output_dict, f"Key '{key}' not found in the output dictionary."

# @pytest.mark.usefixtures("clone_gsfit_repo")
# def test_run_gsfit_and_save_to_mdsplus():
#     """
#     Run GSFIT, then save the results to a dictionary.
#     """

#     from gsfit import Gsfit
    
#     gsfit_controller = Gsfit(
#         pulseNo=cnst.PULSE_NUM,
#         run_name=cnst.RUN_NAME,
#         run_description=cnst.RUN_DESCRIPTION,
#         write_to_mds=False,
#         pulseNo_write=cnst.PULSE_NUM_WRITE,
#         settings_path="st40_setup_for_rtgsfit",
#     )

#     gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["method"] = "user_defined"
#     gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["user_defined"] = [100.0e-3]
#     gsfit_controller.settings["GSFIT_code_settings.json"]["database_writer"]["method"] = "tokamak_energy_mdsplus"

#     gsfit_controller.run()

#     plasma = gsfit_controller.plasma
    
#     output_dict = {}
#     output_dict["grid"] = {}
#     output_dict["grid"]["r"] = plasma.get_array1(["grid", "r"])
#     output_dict["grid"]["z"] = plasma.get_array1(["grid", "z"])
#     output_dict["two_d"] = {}
#     output_dict["two_d"]["psi"] = plasma.get_array3(["two_d", "psi"])[0, :, :]
#     output_dict["two_d"]["psi_n"] = plasma.get_array3(["two_d", "psi_n"])[0, :, :]
#     output_dict["p_boundary"] = {}
#     output_dict["p_boundary"]["rbnd"] = plasma.get_array2(["p_boundary", "rbnd"])[0, :]
#     output_dict["p_boundary"]["zbnd"] = plasma.get_array2(["p_boundary", "zbnd"])[0, :]
#     output_dict["p_boundary"]["nbnd"] = plasma.get_vec_usize(["p_boundary", "nbnd"])[0]
#     # Save the output dictionary to a file
#     output_file = os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy')
#     os.makedirs(cnst.DATA_DIR, exist_ok=True)
#     with open(output_file, 'wb') as f:
#         np.save(f, output_dict, allow_pickle=True)
    