"""Analytic test case: a circular large aspect ratio plasma with a known exact solution.

The machine has a rectangular limiter, flux loops and Bp probes on a rectangle around the plasma,
two Rogowski coils (one around the plasma, one around the passive structures only), four passive
filaments and one PF coil carrying no current. The measurements are those of a filament carrying
`ip` at the plasma centre, which is what the exact solution looks like from outside the plasma, so
the reconstruction should converge to `solution.AnalyticPlasma`. Needs nothing but numpy and scipy.

All Greens functions use the large aspect ratio limit (see `greens`), with total poloidal flux in Wb.
The derivations are in `documentation/analytic_solution` and `documentation/mutual_inductance`.
"""

from __future__ import annotations

import numpy as np

from rtgsfit.constants import MU_0
from rtgsfit.database_readers import PreshotReader
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.database_readers.analytic import greens
from rtgsfit.database_readers.analytic.solution import AnalyticPlasma
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import PreshotInputs
from rtgsfit.preshot import assemble_preshot
from rtgsfit.preshot import cosine_bump_psi_norm
from rtgsfit.preshot import rectangle_points
from rtgsfit.settings import Settings
from rtgsfit.settings import replay_times

__all__ = ["AnalyticPlasma", "AnalyticReader", "analytic_plasma", "probe_geometry"]


def analytic_plasma(settings: Settings) -> AnalyticPlasma:
    plasma = settings["analytic"]["plasma"]
    return AnalyticPlasma(r_0=plasma["r_0"], z_0=plasma["z_0"], rho_b=plasma["rho_b"], ip=plasma["ip"])


def _linspace(section: dict, name: str) -> np.ndarray:
    return np.linspace(section[f"{name}_min"], section[f"{name}_max"], section[f"n_{name}"])


def flux_loop_geometry(settings: Settings) -> tuple[np.ndarray, np.ndarray]:
    section = settings["analytic"]["flux_loops"]
    return rectangle_points(_linspace(section, "r"), _linspace(section, "z"))


def probe_geometry(settings: Settings) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bp probes on the sides of a rectangle (not the corners), pointing anticlockwise."""
    section = settings["analytic"]["bp_probes"]
    r, z = _linspace(section, "r"), _linspace(section, "z")
    r_inner, z_inner = r[1:-1], z[1:-1]
    n_r, n_z = len(r_inner), len(z_inner)
    probe_r = np.concatenate([np.full(n_z, r[0]), r_inner, np.full(n_z, r[-1]), r_inner[::-1]])
    probe_z = np.concatenate([z_inner, np.full(n_r, z[-1]), z_inner[::-1], np.full(n_r, z[0])])
    angle = np.concatenate([np.full(n_z, 0.5 * np.pi), np.zeros(n_r), np.full(n_z, -0.5 * np.pi), np.full(n_r, np.pi)])
    return probe_r, probe_z, angle


def _constraint_greens(settings: Settings, r_source: np.ndarray, z_source: np.ndarray) -> np.ndarray:
    """[N_FLUX_LOOP + N_BP_PROBE, N_SOURCE]: flux loop psi and Bp per ampere in each source filament."""
    loop_r, loop_z = flux_loop_geometry(settings)
    probe_r, probe_z, angle = probe_geometry(settings)
    g_loop = greens.psi(loop_r[:, None], loop_z[:, None], r_source[None, :], z_source[None, :])
    g_probe = greens.b_probe(probe_r[:, None], probe_z[:, None], angle[:, None], r_source[None, :], z_source[None, :])
    return np.vstack([g_loop, g_probe])


def constraint_names(settings: Settings) -> tuple[list[str], int, int, int]:
    n_loop = len(flux_loop_geometry(settings)[0])
    n_probe = len(probe_geometry(settings)[0])
    names = [f"FLUX_LOOP_{i:03d}" for i in range(n_loop)] + [f"BP_PROBE_{i:03d}" for i in range(n_probe)]
    names += ["ROGOWSKI_PLASMA", "ROGOWSKI_PASSIVE"]
    return names, n_loop, n_probe, 2


class AnalyticReader(PreshotReader, SignalsReader):
    def read_preshot(self, settings: Settings) -> Preshot:
        analytic = settings["analytic"]
        grid_r, grid_z = _linspace(analytic["grid"], "r"), _linspace(analytic["grid"], "z")
        n_r, n_z = len(grid_r), len(grid_z)
        d_r, d_z = grid_r[1] - grid_r[0], grid_z[1] - grid_z[0]
        r_2d, z_2d = np.meshgrid(grid_r, grid_z)
        r_flat, z_flat = r_2d.ravel(), z_2d.ravel()

        limiter = analytic["limiter"]
        limiter_r, limiter_z = rectangle_points(_linspace(limiter, "r"), _linspace(limiter, "z"))
        mask_limiter = (
            (r_flat >= limiter["r_min"] - d_r) & (r_flat <= limiter["r_max"] + d_r)
            & (z_flat >= limiter["z_min"] - d_z) & (z_flat <= limiter["z_max"] + d_z)
        )  # fmt: skip

        edge_r, edge_z = rectangle_points(grid_r, grid_z)
        segment = np.hypot(np.roll(edge_r, -1) - edge_r, np.roll(edge_z, -1) - edge_z)
        with np.errstate(divide="ignore"):  # the diagonal is replaced by the self term
            g_edge_edge = greens.psi(edge_r[:, None], edge_z[:, None], edge_r[None, :], edge_z[None, :])
        np.fill_diagonal(g_edge_edge, greens.psi_self(edge_r, edge_z, segment))
        # Hagenow's surface-current conversion uses 1/(mu0 R); put 1/(2 pi) here, as in GSFit.
        g_edge_edge /= 2.0 * np.pi

        passive = analytic["passives"]
        passive_r, passive_z = np.asarray(passive["r"], dtype=float), np.asarray(passive["z"], dtype=float)
        n_passive = len(passive_r)

        names, n_loop, n_probe, n_rogowski = constraint_names(settings)
        n_constraint = len(names)
        g_constraint_grid = np.zeros((n_constraint, n_r * n_z))
        g_constraint_grid[: n_loop + n_probe] = _constraint_greens(settings, r_flat, z_flat)
        g_constraint_grid[n_loop + n_probe] = 1.0  # ROGOWSKI_PLASMA encloses the whole grid
        g_constraint_passive = np.zeros((n_constraint, n_passive))
        g_constraint_passive[: n_loop + n_probe] = _constraint_greens(settings, passive_r, passive_z)
        g_constraint_passive[n_loop + n_probe + 1] = 1.0  # ROGOWSKI_PASSIVE encloses only the passives

        initial = analytic["initial"]
        boundary = settings["preshot"]["boundary"]
        inputs = PreshotInputs(
            mu_0=MU_0,
            grid_r=grid_r,
            grid_z=grid_z,
            poisson_operator="laplacian",
            poisson_rhs_scale=2.0 * np.pi * r_flat * MU_0 * d_z**2,
            mask_limiter=mask_limiter,
            limiter_r=limiter_r,
            limiter_z=limiter_z,
            g_edge_edge=g_edge_edge,
            coil_names=["COIL"],
            g_grid_coil=np.zeros((n_r * n_z, 1)),
            passive_names=[f"PASSIVE_{i}" for i in range(n_passive)],
            passive_dof_units=["ampere"] * n_passive,
            g_grid_passive=greens.psi(r_flat[:, None], z_flat[:, None], passive_r[None, :], passive_z[None, :]),
            constraint_names=names,
            n_flux_loop=n_loop,
            n_bp_probe=n_probe,
            n_rogowski=n_rogowski,
            g_constraint_grid=g_constraint_grid,
            g_constraint_coil=np.zeros((n_constraint, 1)),
            g_constraint_passive=g_constraint_passive,
            # Preserve the physical uncertainties of the original per-radian reference.
            constraint_weight=np.concatenate([np.full(n_loop, 1.0 / (2.0 * np.pi)), np.ones(n_constraint - n_loop)]),
            regularisation_names=[],
            regularisation=np.zeros((0, n_passive)),
            regularisation_weight=np.zeros(0),
            sensor_replacements=[],
            coil_signals=["COIL"],
            coil_matrix=np.ones((1, 1)),
            initial_psi_norm=cosine_bump_psi_norm(grid_r, grid_z, initial["r"], initial["z"], initial["radius"]),
            initial_mask=np.ones(n_r * n_z, dtype=np.int32),
            xpt_fraction=boundary["xpt_fraction"],
            tolerance=boundary["tolerance"],
            ip_min=boundary["ip_min"],
            n_xpt_max=boundary["n_xpt_max"],
            settings_json=settings.json_text,
            source="analytic",
            psi_to_weber=1.0,
        )
        return assemble_preshot(inputs)

    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        """The same measurements at every time: a filament carrying ip at the plasma centre."""
        plasma = analytic_plasma(settings)
        time = replay_times(settings)
        filament = _constraint_greens(settings, np.array([plasma.r_0]), np.array([plasma.z_0]))[:, 0] * plasma.ip
        meas = np.concatenate([filament, [plasma.ip, 0.0]])
        return Signals(
            time=time,
            meas=np.tile(meas, (len(time), 1)),
            coil_signal=np.zeros((len(time), 1)),
        )
