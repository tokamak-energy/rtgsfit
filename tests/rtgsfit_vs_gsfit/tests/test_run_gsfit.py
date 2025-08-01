import os

import numpy as np
import pytest

from rtgsfit_vs_gsfit import cnst

@pytest.mark.usefixtures("run_gsfit_and_save_dict")
def test_run_gsfit_and_save_dict():
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