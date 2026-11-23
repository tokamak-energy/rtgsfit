"""Divertor legs, secondary separatrix and bounding point, computed in Python after `rtgsfit()`.

These are GSFit's `BOUNDARY.BOUNDING`, `SOL` and `SEC_XPT` outputs, so the same plots work for both:

* `BOUNDARY.BOUNDING` is the point that bounds the plasma: the active x-point `XPT.ALL[0]` when
  diverted, otherwise the limiter point that set `BOUNDARY:PSI` (the same choice as the C code);
* `SOL.HFS` and `SOL.LFS` are the divertor legs: the two branches of psi = psi(active x-point) that
  leave the active x-point away from the magnetic axis, from the x-point to the vessel;
* `SEC_XPT.{HFS,LFS}_{LOWER,UPPER}` are the four branches of the secondary separatrix, through the
  first other x-point of `XPT.ALL` that is on the other side of the magnetic axis (above or below).

psi is interpolated with a bicubic spline. Each x-point is refined to grad(psi) = 0 on the spline.
Its branches leave along the asymptotes d^T H d = 0 of the Hessian H and are integrated along the
flux surface with RK4, projected back onto psi = psi(x-point) after every step. A branch ends where
it first crosses the vessel (the strike point), where it leaves the grid or where it returns to its
x-point; only the first has a strike point. The vessel is the leading closed polygon of
`PRESHOT.LIMITER` (its first point repeated); without one, no branches are traced. Limited slices have
a bounding point but no branches.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field

import numpy as np
import numpy.typing as npt
from scipy.interpolate import RectBivariateSpline

from rtgsfit.boundary import polygon_contains
from rtgsfit.data_dictionary import SEC_XPT_LEGS
from rtgsfit.data_dictionary import SOL_LEGS

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int32]

STEP_FRACTION = 0.5
"""Integration step along a branch, as a fraction of the smaller grid spacing."""

NEWTON_ITERATIONS = 20
NEWTON_TOLERANCE = 1e-10
"""Convergence of the x-point refinement, as a fraction of the grid spacing."""

VESSEL_CLOSURE = 1e-9
"""Distance [m] within which a later limiter point closes the vessel polygon."""


def leg_size(n_r: int, n_z: int) -> int:
    """`N_LEG`: maximum number of points stored per branch; longer branches are resampled."""
    return 2 * (n_r + n_z)


@dataclass
class Leg:
    """One separatrix branch, from its x-point to the vessel."""

    r: FloatArray  # [N_LEG], NaN padded
    z: FloatArray  # [N_LEG], NaN padded
    n: int = 0
    strike_r: float = np.nan
    strike_z: float = np.nan

    @classmethod
    def empty(cls, n_leg: int) -> Leg:
        return cls(np.full(n_leg, np.nan), np.full(n_leg, np.nan))

    @classmethod
    def from_points(cls, points: FloatArray, strike: FloatArray | None, n_leg: int) -> Leg:
        if len(points) > n_leg:
            points = points[np.round(np.linspace(0, len(points) - 1, n_leg)).astype(int)]
        leg = cls.empty(n_leg)
        leg.n = len(points)
        leg.r[: len(points)] = points[:, 0]
        leg.z[: len(points)] = points[:, 1]
        if strike is not None:
            leg.strike_r, leg.strike_z = float(strike[0]), float(strike[1])
        return leg


@dataclass
class Separatrix:
    bounding_r: float
    bounding_z: float
    sol: dict[str, Leg]  # keyed by SOL_LEGS
    sec_xpt: dict[str, Leg]  # keyed by SEC_XPT_LEGS

    @classmethod
    def empty(cls, n_leg: int) -> Separatrix:
        return cls(
            np.nan,
            np.nan,
            {name: Leg.empty(n_leg) for name in SOL_LEGS},
            {name: Leg.empty(n_leg) for name in SEC_XPT_LEGS},
        )


@dataclass
class Limiter:
    """`PRESHOT.LIMITER`; `vessel` is its leading closed polygon [N, 2], or None."""

    r: FloatArray
    z: FloatArray
    grid_index: IntArray  # [N_LIMITER, N_INTERP]
    grid_weight: FloatArray  # [N_LIMITER, N_INTERP]
    vessel: FloatArray | None = field(init=False)

    def __post_init__(self) -> None:
        self.vessel = vessel_polygon(self.r, self.z)


def vessel_polygon(limiter_r: FloatArray, limiter_z: FloatArray) -> FloatArray | None:
    """The limiter points up to the first one that repeats the first point, [N, 2]; None if there is
    no such closed polygon. Points after it, such as isolated tiles, are not part of the vessel."""
    r, z = np.asarray(limiter_r, dtype=np.float64), np.asarray(limiter_z, dtype=np.float64)
    closes = np.nonzero(np.hypot(r[1:] - r[0], z[1:] - z[0]) <= VESSEL_CLOSURE)[0] + 1
    if len(closes) == 0 or closes[0] < 3:
        return None
    return np.column_stack((r[: closes[0] + 1], z[: closes[0] + 1]))


def limiter_contact(
    psi: FloatArray,
    limiter: Limiter,
    axis_r: float,
    axis_z: float,
    xpt_r: FloatArray,
    xpt_z: FloatArray,
) -> tuple[float, float]:
    """The limiter point with the largest psi on the core side of every x-point, as in
    `plasma_boundary_limiter_psi()`; native psi is largest on the magnetic axis."""
    psi_limiter = np.sum(np.ravel(psi)[limiter.grid_index] * limiter.grid_weight, axis=1)
    xpt_r, xpt_z = np.asarray(xpt_r), np.asarray(xpt_z)
    dot = (limiter.r[:, None] - xpt_r) * (axis_r - xpt_r) + (limiter.z[:, None] - xpt_z) * (axis_z - xpt_z)
    core_side = np.all(dot >= 0.0, axis=1)
    if not np.any(core_side):
        return np.nan, np.nan
    i = int(np.argmax(np.where(core_side, psi_limiter, -np.inf)))
    return float(limiter.r[i]), float(limiter.z[i])


class _Flux:
    """Bicubic spline of psi [N_Z, N_R]."""

    def __init__(self, r: FloatArray, z: FloatArray, psi: FloatArray) -> None:
        self._spline = RectBivariateSpline(z, r, psi)
        self.r_range = (float(r[0]), float(r[-1]))
        self.z_range = (float(z[0]), float(z[-1]))
        self.spacing = min(float(r[1] - r[0]), float(z[1] - z[0]))

    def value(self, point: FloatArray) -> float:
        return float(self._spline.ev(point[1], point[0]))

    def gradient(self, point: FloatArray) -> FloatArray:
        """(d psi / dR, d psi / dZ)."""
        return np.array([self._spline.ev(point[1], point[0], dy=1), self._spline.ev(point[1], point[0], dx=1)])

    def hessian(self, point: FloatArray) -> FloatArray:
        psi_rr = self._spline.ev(point[1], point[0], dy=2)
        psi_zz = self._spline.ev(point[1], point[0], dx=2)
        psi_rz = self._spline.ev(point[1], point[0], dx=1, dy=1)
        return np.array([[psi_rr, psi_rz], [psi_rz, psi_zz]])

    def contains(self, point: FloatArray) -> bool:
        return self.r_range[0] <= point[0] <= self.r_range[1] and self.z_range[0] <= point[1] <= self.z_range[1]

    def tangent(self, point: FloatArray, heading: FloatArray) -> FloatArray:
        """Unit vector along the flux surface, on the side of `heading`."""
        psi_r, psi_z = self.gradient(point)
        tangent = np.array([-psi_z, psi_r])
        norm = np.hypot(*tangent)
        if norm == 0.0:
            return heading
        tangent /= norm
        return tangent if np.dot(tangent, heading) >= 0.0 else -tangent

    def project(self, point: FloatArray, psi_target: float) -> FloatArray:
        """Newton steps along grad(psi) back onto psi = psi_target."""
        for _ in range(2):
            gradient = self.gradient(point)
            norm_sq = float(np.dot(gradient, gradient))
            if norm_sq == 0.0:
                break
            point = point - (self.value(point) - psi_target) * gradient / norm_sq
        return point


def _refine_xpt(flux: _Flux, start: FloatArray) -> FloatArray | None:
    """The saddle point of the spline near `start`; None if Newton's method leaves the grid cell."""
    point = start.copy()
    for _ in range(NEWTON_ITERATIONS):
        hessian = flux.hessian(point)
        if np.linalg.det(hessian) >= 0.0:
            return None
        step = np.linalg.solve(hessian, -flux.gradient(point))
        point = point + step
        if np.hypot(*(point - start)) > flux.spacing or not flux.contains(point):
            return None
        if np.hypot(*step) < NEWTON_TOLERANCE * flux.spacing:
            return point
    return None


def _branch_directions(
    flux: _Flux, xpt: FloatArray, axis: FloatArray, axis_psi: float
) -> tuple[dict[str, FloatArray], str]:
    """Unit vectors of the four branches at the x-point, keyed by SEC_XPT_LEGS, and "UPPER" or "LOWER",
    the side of the two branches that leave away from the magnetic axis.

    Along one eigenvector of the Hessian psi moves towards axis_psi on both sides: towards the core and
    into the private-flux region. The two branches on the private-flux side are the divertor legs.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(flux.hessian(xpt))
    i_rise = 1 if axis_psi > flux.value(xpt) else 0
    rise, across = eigenvectors[:, i_rise], eigenvectors[:, 1 - i_rise]
    if np.dot(rise, xpt - axis) < 0.0:
        rise = -rise
    angle = np.arctan(np.sqrt(abs(eigenvalues[i_rise] / eigenvalues[1 - i_rise])))
    away = [rise * np.cos(angle) + sign * across * np.sin(angle) for sign in (1.0, -1.0)]
    toward = [-rise * np.cos(angle) + sign * across * np.sin(angle) for sign in (1.0, -1.0)]
    away_side, toward_side = ("UPPER", "LOWER") if xpt[1] > axis[1] else ("LOWER", "UPPER")
    directions: dict[str, FloatArray] = {}
    for side, pair in ((away_side, away), (toward_side, toward)):
        hfs, lfs = sorted(pair, key=lambda direction: direction[0])
        directions[f"HFS_{side}"], directions[f"LFS_{side}"] = hfs, lfs
    return directions, away_side


def _first_crossing(start: FloatArray, end: FloatArray, vessel: FloatArray) -> FloatArray | None:
    """The first point of the segment start -> end on the vessel polygon, or None."""
    a = vessel[:-1]
    edge = vessel[1:] - a
    step = end - start
    offset = a - start
    denominator = step[0] * edge[:, 1] - step[1] * edge[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (offset[:, 0] * edge[:, 1] - offset[:, 1] * edge[:, 0]) / denominator
        u = (offset[:, 0] * step[1] - offset[:, 1] * step[0]) / denominator
    hit = (denominator != 0.0) & (t >= 0.0) & (t <= 1.0) & (u >= 0.0) & (u <= 1.0)
    if not np.any(hit):
        return None
    return start + np.min(t[hit]) * step


def _trace(
    flux: _Flux, xpt: FloatArray, direction: FloatArray, vessel: FloatArray, max_steps: int
) -> tuple[FloatArray, FloatArray | None] | None:
    """Points of the branch leaving `xpt` along `direction`, and its strike point on the vessel.
    None if the branch neither reaches the vessel, nor leaves the grid, nor closes within `max_steps`."""
    psi_x = flux.value(xpt)
    step = STEP_FRACTION * flux.spacing
    points = [xpt]
    position, heading = xpt, direction
    for i_step in range(max_steps):
        if i_step == 0:
            candidate = xpt + step * direction
        else:
            k1 = flux.tangent(position, heading)
            k2 = flux.tangent(position + 0.5 * step * k1, heading)
            k3 = flux.tangent(position + 0.5 * step * k2, heading)
            k4 = flux.tangent(position + step * k3, heading)
            candidate = position + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        candidate = flux.project(candidate, psi_x)
        strike = _first_crossing(position, candidate, vessel)
        if strike is not None:
            points.append(strike)
            return np.array(points), strike
        if not flux.contains(candidate):
            return np.array(points), None
        if i_step > 2 and np.hypot(*(candidate - xpt)) < step:
            points.append(xpt)
            return np.array(points), None
        heading = (candidate - position) / np.hypot(*(candidate - position))
        position = candidate
        points.append(position)
    return None


def _trace_legs(
    flux: _Flux, xpt: FloatArray, directions: dict[str, FloatArray], vessel: FloatArray, n_leg: int
) -> dict[str, Leg]:
    legs = {}
    for name, direction in directions.items():
        traced = _trace(flux, xpt, direction, vessel, 2 * n_leg)
        legs[name] = Leg.empty(n_leg) if traced is None else Leg.from_points(*traced, n_leg)
    return legs


def separatrix(
    r: FloatArray,
    z: FloatArray,
    psi: FloatArray,
    limiter: Limiter,
    *,
    axis_r: float,
    axis_z: float,
    axis_psi: float,
    xpt_r: FloatArray,
    xpt_z: FloatArray,
    xpt_n: int,
    xpt_diverted: int,
) -> Separatrix:
    """The bounding point and separatrix branches of one successful time slice; `psi` is native [N_Z, N_R].

    `xpt_r`, `xpt_z` are the x-points of `rtgsfit()`, sorted with the innermost in flux first.
    """
    n_leg = leg_size(len(r), len(z))
    result = Separatrix.empty(n_leg)
    xpts = np.column_stack((xpt_r[:xpt_n], xpt_z[:xpt_n])).astype(np.float64)
    axis = np.array([axis_r, axis_z], dtype=np.float64)
    if not (np.all(np.isfinite(psi)) and np.all(np.isfinite(xpts)) and np.all(np.isfinite([*axis, axis_psi]))):
        return result
    if not xpt_diverted or xpt_n == 0:
        result.bounding_r, result.bounding_z = limiter_contact(psi, limiter, axis_r, axis_z, xpts[:, 0], xpts[:, 1])
        return result

    flux = _Flux(np.asarray(r, dtype=np.float64), np.asarray(z, dtype=np.float64), psi)
    vessel = limiter.vessel
    active = _refine_xpt(flux, xpts[0])
    result.bounding_r, result.bounding_z = (xpts[0] if active is None else active).tolist()
    if vessel is None:
        return result

    def inside_vessel(point: FloatArray) -> bool:
        return polygon_contains(vessel[:, 0], vessel[:, 1], float(point[0]), float(point[1]))

    if active is not None and inside_vessel(active):
        directions, away_side = _branch_directions(flux, active, axis, axis_psi)
        legs = {name: directions[f"{name}_{away_side}"] for name in SOL_LEGS}
        result.sol = _trace_legs(flux, active, legs, vessel, n_leg)

    for candidate in xpts[1:]:
        if (candidate[1] - axis_z) * (xpts[0, 1] - axis_z) >= 0.0:
            continue
        secondary = _refine_xpt(flux, candidate)
        if secondary is not None and inside_vessel(secondary):
            directions, _ = _branch_directions(flux, secondary, axis, axis_psi)
            result.sec_xpt = _trace_legs(flux, secondary, directions, vessel, n_leg)
            break
    return result
