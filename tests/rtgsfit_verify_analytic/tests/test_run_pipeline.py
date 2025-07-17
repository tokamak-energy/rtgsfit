"""
This module contains a pytest fixture to run the full RTGSFIT pipeline.
In other words, it will generate the constants.c file, then compile RTGSFIT,
and finally run the replay_rtgsfit.py script to generate the output_dict.npy file.
"""

import subprocess
from pathlib import Path
import pytest
import shutil
import numpy as np

from rtgsfit_verify_analytic import cnst

@pytest.fixture(scope="session")
def run_pipeline(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("pipeline_run")

    # Step 1: Generate constants.c
    subprocess.run(
        ["python", "-m", "rtgsfit_verify_analytic.generate_constants_c"],
        check=True
    )

    # Step 2: Clean RTGSFIT/src
    rtgsfit_src = Path(cnst.RTGSFIT_PATH) / "src"
    subprocess.run(["make", "clean"], cwd=rtgsfit_src, check=True)

    # Step 3: Copy constants.c to src
    constants_c = Path(cnst.DATA_DIR) / "constants.c"
    shutil.copy(constants_c, rtgsfit_src / "constants.c")

    # Step 4: Compile RTGSFIT
    # subprocess.run(
    #     "make SHOT=0 RUN_NAME=no_mds LDFLAGS='-shared -Wl,-Ofast -pthread -llapacke -llapack -lblas -lm'",
    #     cwd=rtgsfit_src,
    #     check=True,
    #     shell=True
    # )
    # subprocess.run(
    #     "make SHOT=0 RUN_NAME=no_mds LDFLAGS='-shared -Wl,-Ofast -pthread -L/usr/lib/x86_64-linux-gnu -lopenblas -llapacke -llapack -lm'",
    #     cwd=rtgsfit_src,
    #     check=True,
    #     shell=True
    # )
    subprocess.run(
        "make SHOT=0 RUN_NAME=no_mds",
        cwd=rtgsfit_src,
        check=True,
        shell=True
    )

    # Step 5: Run the RTGSFIT pipeline
    subprocess.run(
        ["python", "-m", "rtgsfit_verify_analytic.replay_rtgsfit"],
        check=True
    )

    return tmp_dir


def test_pipeline_runs_successfully(run_pipeline):
    output_file = Path(cnst.DATA_DIR) / "output_dict.npy"
    assert output_file.exists(), f"{output_file} not found"
    output_dict = np.load(output_file, allow_pickle=True).item()
    assert isinstance(output_dict, dict), "Output is not a dictionary"
