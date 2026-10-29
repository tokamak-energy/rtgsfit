"""Physical and numerical constants shared by the whole package."""

import numpy as np

MU_0 = 4.0e-7 * np.pi
"""Vacuum permeability [H/m]. The pre-2019 exact value, as used by GSFit and the analytic test."""

N_PLASMA = 3
"""Plasma degrees of freedom; fixed by the current basis in `c/src/current_basis.c`."""

PLASMA_DOF_NAMES = ("P_PRIME", "FF_PRIME_MU0", "DELTA_Z")
"""Names of the plasma degrees of freedom, the first `N_PLASMA` entries of `fit_dof`."""

N_INTERP = 4
"""Grid points in the bilinear interpolation stencil of each limiter point; `PLASMA_BOUNDARY_N_INTERP` in C."""

MDSPLUS_NODE_NAME_MAX = 12
"""MDSplus limits node names to 12 characters."""
