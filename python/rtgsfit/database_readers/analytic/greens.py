"""Greens functions of toroidal filaments in the large aspect ratio limit.

These are only used by the analytic test case, whose analytic solution is derived in the same limit.
`(r, z)` is the observation point and `(r_source, z_source)` the filament carrying unit current. The
flux is chosen to vanish at (R, Z) = (0, 0), so it is only meaningful near the filament's major radius.

    psi(r, z)  = -(mu0 / 2 pi) ln(rho / d) r_source,   rho = |(r, z) - (r_source, z_source)|,
                                                      d = |(r_source, z_source)|
    B_R(r, z)  = -(d psi / d z) / r =  (mu0 / 2 pi) (z - z_source) / rho^2 r_source / r
    B_Z(r, z)  =  (d psi / d r) / r = -(mu0 / 2 pi) (r - r_source) / rho^2 r_source / r
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from rtgsfit.constants import MU_0

FloatArray = npt.NDArray[np.float64]
ArrayLike = FloatArray | float


def psi(r: ArrayLike, z: ArrayLike, r_source: ArrayLike, z_source: ArrayLike) -> FloatArray:
    rho = np.hypot(np.subtract(r, r_source), np.subtract(z, z_source))
    distance = np.hypot(r_source, z_source)
    return np.asarray(-(MU_0 / (2.0 * np.pi)) * np.log(rho / distance) * r_source)


def psi_self(r: ArrayLike, z: ArrayLike, length: ArrayLike) -> FloatArray:
    """psi of a filament on itself, averaged over a straight segment of `length` through (r, z)."""
    distance = np.hypot(r, z)
    return np.asarray(-(MU_0 / (2.0 * np.pi)) * (np.log(0.5 * np.asarray(length) / distance) - 1.0) * r)


def b_r(r: ArrayLike, z: ArrayLike, r_source: ArrayLike, z_source: ArrayLike) -> FloatArray:
    d_r, d_z = np.subtract(r, r_source), np.subtract(z, z_source)
    return np.asarray(MU_0 / (2.0 * np.pi) * d_z / (d_r**2 + d_z**2) * np.divide(r_source, r))


def b_z(r: ArrayLike, z: ArrayLike, r_source: ArrayLike, z_source: ArrayLike) -> FloatArray:
    d_r, d_z = np.subtract(r, r_source), np.subtract(z, z_source)
    return np.asarray(-MU_0 / (2.0 * np.pi) * d_r / (d_r**2 + d_z**2) * np.divide(r_source, r))


def b_probe(r: ArrayLike, z: ArrayLike, angle: ArrayLike, r_source: ArrayLike, z_source: ArrayLike) -> FloatArray:
    """Field component along `angle`, measured anticlockwise from the R direction."""
    return np.asarray(b_r(r, z, r_source, z_source) * np.cos(angle) + b_z(r, z, r_source, z_source) * np.sin(angle))
