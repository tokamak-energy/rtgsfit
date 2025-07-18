"""
This file generates the measurements which RTGSFIT uses to fit the plasma equilibrium.

Since our analytic solution is poloidally symmetric, we can replicate the measurements by
assuming that there is a filamanet of current at the o-point (r₀, z₀) with current equal to
the total plasma current.
"""

import numpy as np
import os

from rtgsfit_verify_analytic import mutual_inductance

def generate_measurements(fl_coords, bp_probe_coords, r_o, z_o, total_current):

    meas = np.zeros(len(fl_coords) + len(bp_probe_coords) + 2, dtype=np.float64) # We have two rogovski coils

    i = -1

    # Flux loop measurements are first
    for r, z in fl_coords:
        i += 1
        meas[i] = mutual_inductance.mutual_inductance_psi(r, z, r_o, z_o) * total_current

    # BP probe measurements are second
    for r, z, alpha in bp_probe_coords:
        i += 1
        br = mutual_inductance.mutual_inductance_br(r, z, r_o, z_o) * total_current
        bz = mutual_inductance.mutual_inductance_bz(r, z, r_o, z_o) * total_current
        meas[i] = br * np.cos(alpha) + bz * np.sin(alpha)

    # Finally rogovski coil measurements
    # The first rogovski goes around the plasma
    i += 1
    meas[i] = total_current
    # The second rogovski goes around the just the vessel
    i += 1
    meas[i] = 0

    return meas
