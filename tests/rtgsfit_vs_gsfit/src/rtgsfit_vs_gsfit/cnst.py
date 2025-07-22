"""
Python module containing a list of constants used throughout this repository.
"""

import os

# Github constants
GSFIT_URL = "https://github.com/tokamak-energy/gsfit.git"
R_VS_G_SUB_REPO_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GSFIT_PATH = os.path.join(R_VS_G_SUB_REPO_PATH, "gsfit")
GSFIT_BRANCH = "rtgsfit_mdsplus"

# MDSplus constants
PULSE_NUM = 13343
PULSE_NUM_WRITE = PULSE_NUM + 52_000_000
RUN_NAME = "RT_V_G_1"
RUN_DESCRIPTION = "Using a single degree of freedom for p_prime and ff_prime."

# Directory parameters
REPO_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(REPO_PATH, 'data')
PLOTS_DIR = os.path.join(REPO_PATH, 'plots')
# RTGSFIT_PATH = os.path.join(REPO_PATH, '..', '..')
RTGSFIT_PATH = os.path.dirname(os.path.dirname(REPO_PATH))
