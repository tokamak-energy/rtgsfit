"""Plasma boundary outline and shape, computed in Python after `rtgsfit()` (not in real time).

`rtgsfit()` finds `boundary_psi` and the plasma `mask` on the grid. The outline is the contour
`psi = boundary_psi` that encloses the magnetic axis; the shape parameters follow the IMAS
`equilibrium.time_slice[:].boundary` definitions, as stored by GSFit:

* `GEO_AXIS`: ((R_max + R_min) / 2, (Z_max + Z_min) / 2) of the outline;
* `MINOR_RADIUS`: a = (R_max - R_min) / 2;
* `ELONGATION`: (Z_max - Z_min) / (2 a);
* `TRIANG_U`, `TRIANG_L`: (R_geo - R(Z_max)) / a and (R_geo - R(Z_min)) / a; `TRIANG` is their mean.
"""

from __future__ import annotations

from dataclasses import dataclass

import contourpy
import numpy as np
import numpy.typing as npt
from scipy.ndimage import binary_dilation

FloatArray = npt.NDArray[np.float64]

MASK_DILATION = 2
"""Grid cells the mask is grown by, so that the contour between the last inside and the first outside
point is found, while flux surfaces far outside the plasma are not."""


def outline_size(n_r: int, n_z: int) -> int:
    """`N_OUTLINE`: maximum number of outline points stored; longer outlines are resampled."""
    return 4 * (n_r + n_z)


@dataclass
class BoundaryShape:
    outline_r: FloatArray  # [N_OUTLINE], NaN padded
    outline_z: FloatArray  # [N_OUTLINE], NaN padded
    outline_n: int
    r_geo: float = np.nan
    z_geo: float = np.nan
    minor_radius: float = np.nan
    elongation: float = np.nan
    triangularity_upper: float = np.nan
    triangularity_lower: float = np.nan

    @property
    def triangularity(self) -> float:
        return 0.5 * (self.triangularity_upper + self.triangularity_lower)

    @classmethod
    def empty(cls, n_outline: int) -> BoundaryShape:
        return cls(np.full(n_outline, np.nan), np.full(n_outline, np.nan), 0)


def _contains(r_line: FloatArray, z_line: FloatArray, r: float, z: float) -> bool:
    """Even-odd rule: is (r, z) inside the closed polygon?"""
    r_next, z_next = np.roll(r_line, -1), np.roll(z_line, -1)
    crosses = (z_line > z) != (z_next > z)
    with np.errstate(divide="ignore", invalid="ignore"):
        r_cross = r_line + (z - z_line) * (r_next - r_line) / (z_next - z_line)
    return bool(np.count_nonzero(crosses & (r < r_cross)) % 2)


def boundary_shape(
    r: FloatArray,
    z: FloatArray,
    psi: FloatArray,
    mask: npt.NDArray[np.int32],
    boundary_psi: float,
    axis_psi: float,
    axis_r: float,
    axis_z: float,
) -> BoundaryShape:
    """The boundary of one time slice; `psi` and `mask` are [N_Z, N_R]."""
    n_outline = outline_size(len(r), len(z))
    values = (psi, boundary_psi, axis_psi, axis_r, axis_z)
    if not np.any(mask) or not all(np.all(np.isfinite(value)) for value in values) or axis_psi == boundary_psi:
        return BoundaryShape.empty(n_outline)

    # Outside the grown mask, put psi well on the outside of boundary_psi, so that only the plasma
    # boundary is traced and not open flux surfaces or private-flux regions.
    inside = binary_dilation(mask.astype(bool), iterations=MASK_DILATION)
    outward = np.sign(boundary_psi - axis_psi)
    far_outside = boundary_psi + outward * (np.ptp(psi) + 1.0)
    field = np.where(inside, psi, far_outside)

    lines = contourpy.contour_generator(x=r, y=z, z=field).lines(boundary_psi)
    closed = [line for line in lines if len(line) > 3 and np.allclose(line[0], line[-1])]
    around_axis = [line for line in closed if _contains(line[:, 0], line[:, 1], axis_r, axis_z)]
    candidates = around_axis or closed
    if not candidates:
        return BoundaryShape.empty(n_outline)
    line = max(candidates, key=len)
    if len(line) > n_outline:
        line = line[np.round(np.linspace(0, len(line) - 1, n_outline)).astype(int)]

    shape = BoundaryShape.empty(n_outline)
    shape.outline_n = len(line)
    shape.outline_r[: len(line)] = line[:, 0]
    shape.outline_z[: len(line)] = line[:, 1]
    r_line, z_line = line[:, 0], line[:, 1]
    r_min, r_max, z_min, z_max = r_line.min(), r_line.max(), z_line.min(), z_line.max()
    shape.r_geo = 0.5 * (r_max + r_min)
    shape.z_geo = 0.5 * (z_max + z_min)
    shape.minor_radius = 0.5 * (r_max - r_min)
    if shape.minor_radius > 0.0:
        shape.elongation = (z_max - z_min) / (2.0 * shape.minor_radius)
        shape.triangularity_upper = (shape.r_geo - r_line[np.argmax(z_line)]) / shape.minor_radius
        shape.triangularity_lower = (shape.r_geo - r_line[np.argmin(z_line)]) / shape.minor_radius
    return shape
