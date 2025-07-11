"""
Python module containing a list of constants used throughout this repository.
"""
import os

import numpy as np
import scipy

# Directory parameters
REPO_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(REPO_PATH, 'data')
PLOTS_DIR = os.path.join(REPO_PATH, 'plots')

# Replay parameters
N_ITERS = 10

# Physical constants
MU_0 = 4 * np.pi * 1e-7 # Vacuum permeability (H/m)

# RTGSFIT constants
FRAC = 0.99
N_LCFS_MAX = 500
N_PLS = 3
N_XPT_MAX = 50
N_INTRP = 4
THRESH = 1e-10 # numerical precision threshold

# Analytic solution parameters
ANALYTIC_PLASMA_CURRENT = 1e6
ANALYTIC_RO = 1e1
ANALYTIC_ZO = 0.0
ANALYTIC_RHOB = 0.5
A_STAR = scipy.special.jn_zeros(0, 1)[0] # = 2.4048255576957727
A_J1_A = A_STAR * scipy.special.j1(A_STAR) # = 1.2484591696955067
J1_A = scipy.special.j1(A_STAR) # = 0.5191474972894669
MU0_BAR = 2e-7  # mu_0 / 2pi
DIST = np.sqrt(ANALYTIC_RO**2 + ANALYTIC_ZO**2)  # Distance from the origin to the o-point
PSI_O = MU0_BAR * ANALYTIC_PLASMA_CURRENT * ANALYTIC_RO * (np.log(DIST / ANALYTIC_RHOB) + 1 / A_J1_A)
PSI_B = MU0_BAR * ANALYTIC_PLASMA_CURRENT * ANALYTIC_RO * np.log(DIST / ANALYTIC_RHOB)
A_RAW = MU0_BAR * A_STAR * ANALYTIC_PLASMA_CURRENT * ANALYTIC_RO / (J1_A * ANALYTIC_RHOB**2)

# Grid parameters
R_MIN = ANALYTIC_RO - 1.5 * ANALYTIC_RHOB
R_MAX = ANALYTIC_RO + 1.5 * ANALYTIC_RHOB
Z_MIN = ANALYTIC_ZO - 3 * ANALYTIC_RHOB
Z_MAX = ANALYTIC_ZO + 3 * ANALYTIC_RHOB
N_R = 65
N_Z = 129
N_GRID = N_R * N_Z
D_R = (R_MAX - R_MIN) / (N_R - 1)
D_Z = (Z_MAX - Z_MIN) / (N_Z - 1)

# Limiter parameters
R_LIM_MIN = ANALYTIC_RO - ANALYTIC_RHOB
R_LIM_MAX = ANALYTIC_RO + ANALYTIC_RHOB
Z_LIM_MIN = ANALYTIC_ZO - 2 * ANALYTIC_RHOB
Z_LIM_MAX = ANALYTIC_ZO + 2 * ANALYTIC_RHOB
N_R_LIM = 50
N_Z_LIM = 100

# Flux loop parameters
N_FL_R = 4
N_FL_Z = 8
N_FLUX_LOOPS = 2 * N_FL_R + 2 * N_FL_Z - 4
R_FL_MIN = ANALYTIC_RO - 1.2 * ANALYTIC_RHOB
R_FL_MAX = ANALYTIC_RO + 1.2 * ANALYTIC_RHOB
Z_FL_MIN = ANALYTIC_ZO - 2.4 * ANALYTIC_RHOB
Z_FL_MAX = ANALYTIC_ZO + 2.4 * ANALYTIC_RHOB

# BP probe parameters
N_BP_R = 5
N_BP_Z = 10
N_BP_PROBES = 2 * (N_BP_R - 2) + 2 * (N_BP_Z - 2)  # We don't have corner BP probes
R_BP_MIN = ANALYTIC_RO - 1.2 * ANALYTIC_RHOB
R_BP_MAX = ANALYTIC_RO + 1.2 * ANALYTIC_RHOB
Z_BP_MIN = ANALYTIC_ZO - 2.4 * ANALYTIC_RHOB
Z_BP_MAX = ANALYTIC_ZO + 2.4 * ANALYTIC_RHOB

# Number of Rogovski coils
N_ROGOWSKI_COILS = 2  # One around the plasma, one around the vessel only

# Number of measurements
N_MEAS = N_FLUX_LOOPS + N_BP_PROBES + N_ROGOWSKI_COILS

# Vessel filamament parameters
N_VESSEL = 4
R_VESSEL = np.array([ANALYTIC_RO - 1.1 * ANALYTIC_RHOB,
                     ANALYTIC_RO,
                     ANALYTIC_RO + 1.1 * ANALYTIC_RHOB,
                     ANALYTIC_RO])
Z_VESSEL = np.array([ANALYTIC_ZO,
                     ANALYTIC_ZO + 2.2 * ANALYTIC_RHOB,
                     ANALYTIC_ZO,
                     ANALYTIC_ZO - 2.2 * ANALYTIC_RHOB])

# Number of coefficients
N_COEF = N_PLS + N_VESSEL

