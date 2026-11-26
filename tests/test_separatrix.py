"""Bounding point, divertor legs and secondary separatrix on an analytic double-null flux."""

from __future__ import annotations

import itertools

import numpy as np
import pytest
from rtgsfit.data_dictionary import SEC_XPT_LEGS
from rtgsfit.data_dictionary import SOL_LEGS
from rtgsfit.separatrix import Leg
from rtgsfit.separatrix import Limiter
from rtgsfit.separatrix import Separatrix
from rtgsfit.separatrix import separatrix

R0, C, EPS = 0.6, 0.5, 0.02
R = np.linspace(0.1, 1.1, 51)
Z = np.linspace(-1.0, 1.0, 101)
VESSEL = (0.15, 1.05, -0.9, 0.9)
TILE = (0.3, 0.0)

# Measured worst errors: x-point 6.7e-8 m; psi on the contour 1.0e-8, and 6.3e-6 at the strike
# point, which lies on the straight step across the wall (psi_axis - psi_x = 0.053).
ATOL_XPT = 2e-7
ATOL_PSI_CONTOUR = 3e-8
ATOL_PSI_STRIKE = 2e-5


def flux(r: np.ndarray | float, z: np.ndarray | float) -> np.ndarray:
    """Native psi, largest on the axis, with saddle points on R = R0 near Z = -C (active) and Z = +C."""
    return -((r - R0) ** 2) + (z**2 - C**2) ** 2 / (4 * C**2) - EPS * z


def critical_z() -> np.ndarray:
    """Z of the lower x-point, the axis and the upper x-point: roots of C^2 dpsi/dZ on R = R0."""
    return np.sort(np.roots([1.0, 0.0, -(C**2), -EPS * C**2]).real)


def rectangle(r_min: float, r_max: float, z_min: float, z_max: float, closed: bool = True) -> np.ndarray:
    corners = [(r_min, z_min), (r_max, z_min), (r_max, z_max), (r_min, z_max), (r_min, z_min)]
    sides = [np.linspace(a, b, 20, endpoint=False) for a, b in itertools.pairwise(corners)]
    points = np.concatenate(sides)
    return np.vstack((points, points[:1])) if closed else points


def limiter(points: np.ndarray) -> Limiter:
    """Limiter with bilinear stencils into psi[N_Z, N_R] flattened, as in the preshot."""
    ir = np.clip(np.searchsorted(R, points[:, 0]) - 1, 0, len(R) - 2)
    iz = np.clip(np.searchsorted(Z, points[:, 1]) - 1, 0, len(Z) - 2)
    fr = (points[:, 0] - R[ir]) / (R[ir + 1] - R[ir])
    fz = (points[:, 1] - Z[iz]) / (Z[iz + 1] - Z[iz])
    corner = iz * len(R) + ir
    index = np.column_stack((corner, corner + 1, corner + len(R), corner + len(R) + 1)).astype(np.int32)
    weight = np.column_stack(((1 - fr) * (1 - fz), fr * (1 - fz), (1 - fr) * fz, fr * fz))
    return Limiter(points[:, 0], points[:, 1], index, weight)


def diverted(vessel: np.ndarray) -> Separatrix:
    """The x-points are given off their true positions, as from the grid search of `rtgsfit()`."""
    z_lower, z_axis, z_upper = critical_z()
    rr, zz = np.meshgrid(R, Z)
    return separatrix(
        R,
        Z,
        flux(rr, zz),
        limiter(np.vstack((vessel, TILE))),
        axis_r=R0,
        axis_z=z_axis,
        axis_psi=float(flux(R0, z_axis)),
        xpt_r=np.array([R0 + 0.004, R0 - 0.004, np.nan]),
        xpt_z=np.array([z_lower - 0.004, z_upper + 0.004, np.nan]),
        xpt_n=2,
        xpt_diverted=1,
    )


def legs(result: Separatrix) -> dict[str, Leg]:
    assert result.sol.keys() == set(SOL_LEGS) and result.sec_xpt.keys() == set(SEC_XPT_LEGS)
    return {f"SOL.{name}": leg for name, leg in result.sol.items()} | {
        f"SEC_XPT.{name}": leg for name, leg in result.sec_xpt.items()
    }


def test_diverted_branches_follow_the_flux_surfaces_to_the_vessel() -> None:
    z_lower, _, z_upper = critical_z()
    result = diverted(rectangle(*VESSEL))
    assert (result.bounding_r, result.bounding_z) == pytest.approx((R0, z_lower), abs=ATOL_XPT)

    r_min, r_max, z_min, z_max = VESSEL
    for name, leg in legs(result).items():
        xpt_z = z_lower if name.startswith("SOL") else z_upper
        r, z = leg.r[: leg.n], leg.z[: leg.n]
        assert leg.n > 10, name
        assert np.all(np.isnan(leg.r[leg.n :])), name
        assert (r[0], z[0]) == pytest.approx((R0, xpt_z), abs=ATOL_XPT), name
        assert (r[-1], z[-1]) == (leg.strike_r, leg.strike_z), name
        assert flux(r[:-1], z[:-1]) == pytest.approx(flux(R0, xpt_z), abs=ATOL_PSI_CONTOUR), name
        assert flux(r[-1], z[-1]) == pytest.approx(flux(R0, xpt_z), abs=ATOL_PSI_STRIKE), name
        distance_to_wall = min(
            abs(leg.strike_r - r_min), abs(leg.strike_r - r_max), abs(leg.strike_z - z_min), abs(leg.strike_z - z_max)
        )
        assert distance_to_wall < 1e-12, name
        assert (leg.strike_r < R0) == ("HFS" in name), name
        assert (leg.strike_z > xpt_z) == ("UPPER" in name), name

    for name in SOL_LEGS:
        assert result.sol[name].strike_z < z_lower, name


@pytest.mark.parametrize(
    ("vessel", "traced"),
    [
        pytest.param(rectangle(*VESSEL, closed=False), False, id="limiter-not-closed"),
        pytest.param(rectangle(0.0, 1.2, -1.2, 1.2), True, id="vessel-beyond-grid"),
    ],
)
def test_branches_without_a_strike_point(vessel: np.ndarray, traced: bool) -> None:
    result = diverted(vessel)
    assert np.isfinite([result.bounding_r, result.bounding_z]).all()
    for name, leg in legs(result).items():
        assert (leg.n > 0) == traced, name
        assert np.isnan([leg.strike_r, leg.strike_z]).all(), name


def test_limited_slice_is_bounded_by_the_limiter_point_with_the_largest_flux() -> None:
    """The isolated tile is not part of the vessel, but it can still be the limiter point."""
    vessel = rectangle(*VESSEL)
    tiles = limiter(np.vstack((vessel, TILE)))
    rr, zz = np.meshgrid(R, Z)
    result = separatrix(
        R,
        Z,
        -((rr - R0) ** 2) - zz**2,
        tiles,
        axis_r=R0,
        axis_z=0.0,
        axis_psi=0.0,
        xpt_r=np.full(2, np.nan),
        xpt_z=np.full(2, np.nan),
        xpt_n=0,
        xpt_diverted=0,
    )
    assert tiles.vessel is not None and len(tiles.vessel) == len(vessel)
    assert (result.bounding_r, result.bounding_z) == TILE
    assert all(leg.n == 0 for leg in legs(result).values())
