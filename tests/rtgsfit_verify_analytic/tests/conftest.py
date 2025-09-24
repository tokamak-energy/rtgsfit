import subprocess
import numpy as np
from pathlib import Path
import pytest
import shutil

from rtgsfit_verify_analytic import cnst, generate_constants_c, replay_rtgsfit

@pytest.fixture(scope="session")
def run_pipeline():
    """
    Run the full RTGSFIT pipeline once per test session.
    """

    generate_constants_c.generate_constants_c()

    rtgsfit_src = Path(cnst.RTGSFIT_PATH) / "src"
    subprocess.run(["make", "clean"], cwd=rtgsfit_src, check=True)

    constants_c = Path(cnst.DATA_DIR) / "constants.c"
    shutil.copy(constants_c, rtgsfit_src / "constants.c")

    subprocess.run(
        "make SHOT=0 RUN_NAME=no_mds DEBUG=1",
        cwd=rtgsfit_src,
        check=True,
        shell=True
    )

    replay_rtgsfit.replay_rtgsfit()
    

@pytest.fixture(scope="session")
def output_dict(run_pipeline):
    output_file = Path(cnst.DATA_DIR) / "output_dict.npy"
    if not output_file.exists():
        pytest.fail(f"{output_file} does not exist. Run test_run_pipeline first.")
    return np.load(output_file, allow_pickle=True).item()
