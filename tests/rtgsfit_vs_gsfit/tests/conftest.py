import os
import shutil
import subprocess

import pytest

from rtgsfit_vs_gsfit import cnst

@pytest.fixture(scope="session")
def clone_gsfit_repo():
    """
    1. Delete any existing GSFIT repository and clone a fresh copy.
    This ensures that the tests always run against the latest version of GSFIT.

    2. Clone the GSFIT repository.

    3. Switch to the branch specified in src/cnst.py file.
    """ 
    if not os.path.exists(cnst.GSFIT_PATH):
        subprocess.run(
            ["git", "clone", cnst.GSFIT_URL, str(cnst.GSFIT_PATH)],
            check=True
        )
    else:
        print(f"GSFIT repository already exists at {cnst.GSFIT_PATH}. Deleting and re-cloning.")
        shutil.rmtree(cnst.GSFIT_PATH)
        subprocess.run(
            ["git", "clone", cnst.GSFIT_URL, str(cnst.GSFIT_PATH)],
            check=True
        )

    # Change to the GSFIT directory
    os.chdir(cnst.GSFIT_PATH)

    # Checkout the correct branch
    subprocess.run(
        ["git", "checkout", cnst.GSFIT_BRANCH],
        check=True
    )
    print(f"GSFIT repository cloned and switched to branch {cnst.GSFIT_BRANCH} at {cnst.GSFIT_PATH}")