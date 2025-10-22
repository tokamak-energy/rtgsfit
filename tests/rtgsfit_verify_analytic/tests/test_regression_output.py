import os

import pytest
import numpy as np
from numpy.testing import assert_allclose

from rtgsfit_verify_analytic import cnst

def test_output_dict_regression(output_dict):
    """
    Regression test for the final step only to save storage space.
    """    

    output_dict_final_iter_path = os.path.join(cnst.DATA_DIR, 'output_dict_final_iter.npy')
    output_dict_final_iter = np.load(output_dict_final_iter_path, allow_pickle=True).item()
    reference_file = os.path.join(cnst.DATA_DIR, 'reference_output_dict.npy')
    reference_dict = np.load(reference_file, allow_pickle=True).item()

    assert output_dict.keys() == reference_dict.keys(), "Mismatch in output_dict keys"

    for key in output_dict:
        rtol = 1e-4
        atol = np.max(np.abs(reference_dict[key])) * rtol
        assert_allclose(
            output_dict_final_iter[key],
            reference_dict[key],
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in key '{key}'"
        )
