"""
This module contains a pytest fixture to run the full RTGSFIT pipeline.
In other words, it will generate the constants.c file, then compile RTGSFIT,
and finally run the replay_rtgsfit.py script to generate the output_dict.npy file.
"""

import numpy as np
from pathlib import Path
from rtgsfit_verify_analytic import cnst

def test_pipeline_runs_successfully(run_pipeline):
    constants_c = Path(cnst.DATA_DIR) / "constants.c"
    assert constants_c.exists(), f"{constants_c} not found"

    output_file = Path(cnst.DATA_DIR) / "output_dict.npy"
    assert output_file.exists(), f"{output_file} not found"

    output_dict = np.load(output_file, allow_pickle=True).item()
    assert isinstance(output_dict, dict), "Output is not a dictionary"
