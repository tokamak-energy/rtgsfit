"""
Python module containing a list of constants used throughout this repository.
"""

import os

import numpy as np

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability

# Replay parameters
TIME = 30e-3 # Time in seconds for which GSFIT and RTGSFIT are run

# GSFit constants
settings_path = "st40_setup_for_rtgsfit"

# MDSplus constants
PULSE_NUM = 13_343
PULSE_NUM_WRITE = PULSE_NUM + 52_000_000
RUN_NAME = "RT_V_G_1"
RUN_DESCRIPTION = "Using a single degree of freedom for p_prime and ff_prime."

# Directory parameters
REPO_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(REPO_PATH, 'data')
PLOTS_DIR = os.path.join(REPO_PATH, 'plots')
RTGSFIT_PATH = os.path.dirname(os.path.dirname(REPO_PATH))
RTGSFIT_SRC_PATH = os.path.join(RTGSFIT_PATH, 'src')

# PSU2COIL parameters
if PULSE_NUM > 13_000:
    PSU2COIL_RUN_NAME = "RUN05"
else:
    PSU2COIL_RUN_NAME = "RUN01"

# Initial flux norm parameters
R_AXIS0 = 0.5
Z_AXIS0 = 0.0
RHO_BOUNDARY0 = 0.2

# Replay RTGSFIT parameters
N_ITERS = 10
