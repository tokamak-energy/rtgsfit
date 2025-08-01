import os

import mdsthin
import numpy as np
import pytest

from rtgsfit_vs_gsfit import cnst, replay_rtgsfit

@pytest.mark.usefixtures("rtgsfit_mds_nodegen")
def test_rtgsfit_mds_nodegen():
    """
    Test that the RTGSFIT MDSplus node was created correctly.
    """

    conn = mdsthin.Connection('smaug')
    conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)

    # Check PRESHOT.N_MEAS exists and is float equal to an integer
    n_meas = conn.get(f"\\RTGSFIT::TOP:{cnst.RUN_NAME}.PRESHOT:N_MEAS").data()
    assert np.isclose(n_meas, int(n_meas)), "PRESHOT.N_MEAS is not equal to an integer."

    # Check PRESHOT.N_REG exists and is float equal to an integer
    # n_reg = conn.tcl(f'print \\\\RTGSFIT::TOP:{cnst.RUN_NAME}:PRESHOT.N_REG')
    # assert np.isclose(n_reg, int(n_reg)), "PRESHOT.N_REG is not equal to an integer."

@pytest.mark.usefixtures("compile_rtgsfit")
def test_compile_rtgsfit():
    """
    Test that RTGSFIT compiles successfully.
    """

    # Check if object files are created
    object_files = [f for f in os.listdir(cnst.RTGSFIT_SRC_PATH) if f.endswith('.o')]
    assert object_files, "Object files were not created during compilation."

    # Check if constants.c was created
    constants_c_path = os.path.join(cnst.RTGSFIT_SRC_PATH, 'constants.c')
    assert os.path.exists(constants_c_path), f"constants.c was not created: {constants_c_path}"

@pytest.mark.usefixtures("compile_rtgsfit")
def test_replay_rtgsfit():
    """
    Test that replay_rtgsfit runs without errors.
    """

    # Run the replay function
    replay_rtgsfit.replay_rtgsfit(cnst.TIME)

    # Check if the output file was created
    output_file = os.path.join(cnst.DATA_DIR, 'rtgsfit_output_dict.npy')
    assert os.path.exists(output_file), f"Output file was not created: {output_file}"
    