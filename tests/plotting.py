"""Diagnostic figures of the tests: opt-in with `pytest --figures`, never made in CI.

Each pytest run writes to a new directory, `tests/figures/<UTC time>/` (the `figure_dir` fixture in
conftest.py), so earlier figures are never overwritten. The figures show what the tests compare:
RT-GSFit at its first, second and last iterations (as in v1) against a `Reference` (the exact
solution, GSFit or FreeGSNKE), with the tests' tolerances, and the convergence over every iteration.

The figures are drawn on `matplotlib.figure.Figure` rather than with pyplot, so no GUI backend or
global state is involved. Each is saved as PNG and PDF.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.colors import SymLogNorm
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from rtgsfit.library import describe_status
from rtgsfit.results import Results
from scipy.interpolate import RegularGridInterpolator

DPI = 150
FORMATS = ("png", "pdf")
RTGSFIT = "tab:red"
REFERENCE = "tab:blue"
EXPECTED_STYLES = (("+", REFERENCE), ("s", "tab:green"), ("D", "tab:purple"))
SENSOR_STYLES = (("s", "tab:orange"), ("^", "tab:purple"), ("D", "tab:brown"))
CONSTRAINT_GROUPS = {
    "FLUX_LOOP": ("Flux loops", "psi", "Wb", 1.0),
    "BP_PROBE": ("Bp probes", "B", "mT", 1e3),
    "ROGOWSKI": ("Rogowski coils", "I", "kA", 1e-3),
}
MAPS = {"PSI": ("psi", "Wb", 1.0), "J_PHI": ("j_phi", "MA/m²", 1e-6)}


@dataclass
class Expected:
    """Values that a test compares RT-GSFit's last iteration with, and the test's tolerance:
    |RT-GSFit - values| <= atol + rtol |values|, as `numpy.testing.assert_allclose`."""

    label: str
    values: np.ndarray
    rtol: float = 0.0
    atol: float = 0.0

    def tolerance(self) -> np.ndarray:
        return self.atol + self.rtol * np.abs(self.values)


@dataclass
class Reference:
    """What RT-GSFit is compared with: the exact solution, GSFit or FreeGSNKE.

    The maps are [n_z, n_r] on the reference's own grid (`r`, `z`), and every value is in RT-GSFit's
    units: psi is the total poloidal flux [Wb], not per radian. Unknown values are NaN or None.
    """

    name: str
    r: np.ndarray
    z: np.ndarray
    psi: np.ndarray
    j_phi: np.ndarray | None = None
    ip: float = math.nan
    axis_r: float = math.nan
    axis_z: float = math.nan
    axis_psi: float = math.nan
    boundary_psi: float = math.nan
    outline: tuple[np.ndarray, np.ndarray] | None = None
    """(R, Z) of the boundary; without it the contour psi = boundary_psi is drawn instead."""
    xpt: tuple[np.ndarray, np.ndarray] | None = None
    sensors: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    """Label -> (R, Z) of the sensors, drawn on the maps."""
    constraints: dict[str, list[Expected]] = field(default_factory=dict)
    """`CONSTRAINTS` group (e.g. FLUX_LOOP) -> what the tests compare its reconstruction with, in RT-GSFit's order."""


@dataclass
class Filaments:
    """Filaments of a passive structure: centres and sizes [m]."""

    r: np.ndarray
    z: np.ndarray
    d_r: np.ndarray
    d_z: np.ndarray


def save(figures: Mapping[str, Figure], directory: Path) -> list[Path]:
    """Save every figure as `<directory>/<name>.png` and `.pdf`."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, figure in figures.items():
        for suffix in FORMATS:
            path = directory / f"{name}.{suffix}"
            figure.savefig(path, dpi=DPI)
            paths.append(path)
    return paths


# ----------------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------------
def iterations(results: Results) -> list[int]:
    """The first, second and last calls of `rtgsfit()`, as in v1."""
    n = len(results["TIME"])
    return sorted({0, min(1, n - 1), n - 1})


def _iteration(results: Results, i: int) -> str:
    title = f"iteration {i + 1} of {len(results['TIME'])}"
    status = int(results["GLOBAL.STATUS"][i])
    return title if status == 0 else f"{title}\n{describe_status(status)}"


def _shades(n: int) -> list[Any]:
    """RT-GSFit's iterations, from light (first) to dark (last)."""
    return [colormaps["Reds"](x) for x in np.linspace(0.4, 0.95, n)]


def _short_names(names: Sequence[str]) -> list[str]:
    """Drop the leading words, separated by underscores, that every name shares: `PSI_FLOOP_007` -> `007`."""
    words = [str(name).split("_") for name in names]
    n = 0
    while all(len(word) > n + 1 and word[n] == words[0][n] for word in words):
        n += 1
    return ["_".join(word[n:]) for word in words]


def _on_grid(results: Results, reference: Reference, values: np.ndarray) -> np.ndarray:
    """A map of the reference on RT-GSFit's grid: as it is if the grids are the same, otherwise
    interpolated linearly (NaN outside the reference's grid)."""
    r, z = results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"]
    same_shape = (len(reference.r), len(reference.z)) == (len(r), len(z))
    if same_shape and np.allclose(reference.r, r) and np.allclose(reference.z, z):
        return values
    interpolator = RegularGridInterpolator((reference.z, reference.r), values, bounds_error=False, fill_value=np.nan)
    z_2d, r_2d = np.meshgrid(z, r, indexing="ij")
    return np.asarray(interpolator((z_2d, r_2d)))


def _plot_limiter(ax: Axes, results: Results) -> None:
    """The limiter contours as lines and separate limiter points as dots.

    The ST40 limiter is a closed contour followed by two separate points, so it is split wherever the
    gap to the next point is much larger than the usual spacing (as in the ST40 example).
    """
    r = np.asarray(results["PRESHOT.LIMITER.R"], dtype=float)
    z = np.asarray(results["PRESHOT.LIMITER.Z"], dtype=float)
    step = np.hypot(np.diff(r), np.diff(z))
    max_step = 4.0 * np.median(step)
    gaps = np.flatnonzero(step > max_step) + 1
    for k, (r_p, z_p) in enumerate(zip(np.split(r, gaps), np.split(z, gaps), strict=True)):
        if len(r_p) == 1:
            ax.plot(r_p, z_p, "ko", ms=3)
            continue
        if len(r_p) > 2 and np.hypot(r_p[-1] - r_p[0], z_p[-1] - z_p[0]) <= max_step:
            r_p, z_p = np.r_[r_p, r_p[:1]], np.r_[z_p, z_p[:1]]
        ax.plot(r_p, z_p, "k-", lw=1.2, label="limiter" if k == 0 else None)


def _draw_machine(ax: Axes, results: Results, i: int, reference: Reference, sensors: bool = True) -> None:
    """Limiter; boundary, magnetic axis and x-points of RT-GSFit (solid) and the reference (dashed); sensors."""
    _plot_limiter(ax, results)
    n = int(results["BOUNDARY.OUTLINE.N"][i])
    if n:
        ax.plot(results["BOUNDARY.OUTLINE.R"][i, :n], results["BOUNDARY.OUTLINE.Z"][i, :n], "-", color=RTGSFIT, lw=2,
                label="RT-GSFit boundary")  # fmt: skip
    if reference.outline is not None:
        ax.plot(*reference.outline, "--", color=REFERENCE, lw=2, label=f"{reference.name} boundary")
    elif math.isfinite(reference.boundary_psi):
        ax.contour(reference.r, reference.z, reference.psi, levels=[reference.boundary_psi], colors=REFERENCE,
                   linewidths=2, linestyles="--")  # fmt: skip
        ax.plot([], [], "--", color=REFERENCE, lw=2, label=f"{reference.name} boundary")
    ax.plot(results["GLOBAL.MAG_AXIS.R"][i], results["GLOBAL.MAG_AXIS.Z"][i], "+", color=RTGSFIT, ms=12, mew=2,
            label="RT-GSFit magnetic axis")  # fmt: skip
    ax.plot(reference.axis_r, reference.axis_z, "+", color=REFERENCE, ms=12, mew=2,
            label=f"{reference.name} magnetic axis")  # fmt: skip
    n_xpt = int(results["XPT.ALL.N"][i])
    ax.plot(results["XPT.ALL.R"][i, :n_xpt], results["XPT.ALL.Z"][i, :n_xpt], "x", color=RTGSFIT, ms=8, mew=1.5,
            label="RT-GSFit x-points")  # fmt: skip
    if reference.xpt is not None:
        ax.plot(*reference.xpt, "x", color=REFERENCE, ms=8, mew=1.5, label=f"{reference.name} x-points")
    if sensors:
        for (label, (r, z)), (marker, colour) in zip(reference.sensors.items(), SENSOR_STYLES, strict=False):
            ax.plot(r, z, marker, color=colour, mfc="none", ms=4, label=label)


def _set_extent(axes: np.ndarray, results: Results, reference: Reference) -> None:
    """Cover RT-GSFit's grid and the sensors, with equal aspect."""
    r = np.concatenate([results["PROFILES_2D.R_Z.R"], *(np.ravel(rz[0]) for rz in reference.sensors.values())])
    z = np.concatenate([results["PROFILES_2D.R_Z.Z"], *(np.ravel(rz[1]) for rz in reference.sensors.values())])
    margin = 0.03 * max(np.ptp(r), np.ptp(z))
    for ax in axes.flat:
        ax.set_aspect("equal")
        ax.set_xlim(r.min() - margin, r.max() + margin)
        ax.set_ylim(z.min() - margin, z.max() + margin)
    for ax in axes[-1]:
        ax.set_xlabel("R [m]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Z [m]")


def _legend_below(fig: Figure, ax: Axes, ncols: int = 4) -> None:
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=True))
    fig.legend(unique.values(), unique.keys(), loc="outside lower center", ncols=ncols, fontsize="small")


def _reference_line(ax: Axes, value: float, reference: Reference) -> None:
    if math.isfinite(value):
        ax.axhline(value, color=REFERENCE, ls="--", label=reference.name)


def _shade_failures(ax: Axes, x: np.ndarray, results: Results) -> None:
    failed = results["GLOBAL.STATUS"] != 0
    if np.any(failed):
        ax.fill_between(x, 0, 1, where=failed, step="mid", color=RTGSFIT, alpha=0.12, lw=0,
                        transform=ax.get_xaxis_transform(), label="status != OK")  # fmt: skip


# ----------------------------------------------------------------------------------------------
# Figures for every test suite
# ----------------------------------------------------------------------------------------------
def plot_map(results: Results, reference: Reference, quantity: str, rtol: float = 0.0, atol: float = 0.0) -> Figure:
    """RT-GSFit's `PSI` or `J_PHI` (`quantity`) at the first, second and last iterations against the reference.

    Top: contours of RT-GSFit (solid) and the reference (dashed) at the same levels. Bottom: the difference
    RT-GSFit - reference; the points outside the test's tolerance, atol + rtol |reference|, are outlined in black.
    """
    label, unit, scale = MAPS[quantity]
    expected = getattr(reference, quantity.lower())
    if expected is None:
        raise ValueError(f"{reference.name} has no {label} map")
    r, z = results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"]
    shown = iterations(results)
    values = results[f"PROFILES_2D.R_Z.{quantity}"][shown]
    difference = values - _on_grid(results, reference, expected)
    outside = np.abs(difference) > atol + rtol * np.abs(_on_grid(results, reference, expected))
    both = np.concatenate([values.ravel(), np.ravel(expected)]) * scale
    levels = np.linspace(np.nanmin(both), np.nanmax(both), 27)[1:-1]
    tested = bool(rtol or atol)

    extent_r, extent_z = np.ptp(r), np.ptp(z)
    height = float(np.clip(3.6 * extent_z / extent_r, 2.5, 7.0))
    fig = Figure(figsize=(3.6 * len(shown) + 1.0, 2 * height + 2.5), layout="constrained")
    axes = fig.subplots(2, len(shown), sharex=True, sharey=True, squeeze=False)
    for column, i in enumerate(shown):
        ax = axes[0, column]
        ax.contour(r, z, values[column] * scale, levels=levels, colors=RTGSFIT, linewidths=0.7)
        ax.contour(reference.r, reference.z, expected * scale, levels=levels, colors=REFERENCE, linewidths=0.7,
                   linestyles="--")  # fmt: skip
        _draw_machine(ax, results, i, reference)
        ax.set_title(_iteration(results, i), fontsize="medium")

        # Each iteration has its own colour scale, so that the small differences of the last are visible
        ax = axes[1, column]
        limit = float(np.nanmax(np.abs(difference[column]))) * scale or 1.0
        image = ax.pcolormesh(r, z, difference[column] * scale, shading="nearest", cmap="RdBu_r", vmin=-limit,
                              vmax=limit)  # fmt: skip
        if tested and np.any(outside[column]):
            ax.contour(r, z, outside[column].astype(float), levels=[0.5], colors="k", linewidths=1.0)
        _draw_machine(ax, results, i, reference, sensors=False)
        colorbar = fig.colorbar(image, ax=ax, location="bottom", shrink=0.9,
                                label=f"RT-GSFit - {reference.name} [{unit}]")  # fmt: skip
        # Few ticks, with the power of ten as an offset, so that the labels of small differences do not overlap
        colorbar.locator = MaxNLocator(5, symmetric=True)
        colorbar.formatter.set_powerlimits((-2, 3))
        colorbar.update_ticks()
        title = f"max |difference| = {limit:.3g} {unit}"
        if tested:
            title += f"\n{np.count_nonzero(outside[column])} points outside the tolerance"
        ax.set_title(title, fontsize="medium")
    _set_extent(axes, results, reference)
    suptitle = f"{label} [{unit}]: RT-GSFit (solid) and {reference.name} (dashed)"
    if tested:
        suptitle += f"\ntest tolerance rtol = {rtol:g}, atol = {atol:g}; points outside it are outlined in black"
    fig.suptitle(suptitle)
    _legend_below(fig, axes[0, 0])
    return fig


def plot_constraints(results: Results, reference: Reference | None = None) -> Figure:
    """Every magnetic constraint.

    Top: the measurement that RT-GSFit fits, the values the tests compare with, and RT-GSFit's reconstruction at the
    first, second and last iterations. Below, one row for each comparison of the tests: RT-GSFit's last iteration
    minus the values, with the test's tolerance shaded and the points outside it circled; without any, minus the
    measurement.
    """
    shown = iterations(results)
    groups = [group for group in CONSTRAINT_GROUPS if len(results[f"CONSTRAINTS.{group}.ALL.NAMES"])]
    names = {group: _short_names(results[f"CONSTRAINTS.{group}.ALL.NAMES"]) for group in groups}
    expected = {group: reference.constraints.get(group, []) if reference is not None else [] for group in groups}
    n_rows = 1 + max(max(len(items), 1) for items in expected.values())
    # Columns as wide as their number of sensors, but wide enough for the titles of a group with a few sensors
    n_sensors = sum(len(names[group]) for group in groups)
    widths = [max(len(names[group]), 0.2 * n_sensors) for group in groups]
    fig = Figure(figsize=(17, 3.4 * n_rows + 1.5), layout="constrained")
    axes = fig.subplots(n_rows, len(groups), sharex="col", squeeze=False, gridspec_kw={"width_ratios": widths})
    for column, group in enumerate(groups):
        title, quantity, unit, scale = CONSTRAINT_GROUPS[group]
        node = f"CONSTRAINTS.{group}.ALL"
        x = np.arange(len(names[group]))
        measured = results[f"{node}.MEASURED"][-1]
        reconstructed = results[f"{node}.RECONSTRUCT"]

        ax = axes[0, column]
        ax.plot(x, measured * scale, "o", color="k", mfc="none", ms=7, label="measured (RT-GSFit's fit target)")
        for item, (marker, colour) in zip(expected[group], EXPECTED_STYLES, strict=False):
            ax.plot(x, item.values * scale, marker, color=colour, mfc="none", ms=7, label=item.label)
        for i, shade in zip(shown, _shades(len(shown)), strict=True):
            ax.plot(x, reconstructed[i] * scale, "x", color=shade, ms=6, label=f"RT-GSFit, iteration {i + 1}")
        chi_sq = results[f"CONSTRAINTS.{group}.CHI_SQ"][-1]
        ax.set(title=f"{title}\nchi_sq = {chi_sq:.3g}", ylabel=f"{quantity} [{unit}]")

        if expected[group]:
            comparisons = list(zip(expected[group], EXPECTED_STYLES, strict=False))
        else:
            comparisons = [(Expected("measured", measured), ("o", "k"))]
        for row, (item, (marker, colour)) in enumerate(comparisons, start=1):
            ax = axes[row, column]
            error = reconstructed[-1] - item.values
            tolerance = item.tolerance()
            ax.plot(x, error * scale, marker, color=colour, mfc="none", ms=7)
            title = f"RT-GSFit - {item.label}"
            if np.any(tolerance > 0):
                outside = np.abs(error) > tolerance
                band = f"test tolerance: rtol = {item.rtol:g}, atol = {item.atol * scale:g} {unit}"
                ax.fill_between(x, -tolerance * scale, tolerance * scale, step="mid", color=colour, alpha=0.15,
                                label=band)  # fmt: skip
                if np.any(outside):
                    ax.plot(x[outside], error[outside] * scale, "o", color="k", ms=11, mfc="none", mew=1.5,
                            label="outside the tolerance")  # fmt: skip
                ax.legend(fontsize="x-small")
                title += f"\n{np.count_nonzero(outside)} of {len(x)} outside the tolerance"
            ax.axhline(0.0, color="k", lw=0.5)
            ax.set(title=title, ylabel=f"difference [{unit}]")
        for ax in axes[len(comparisons) + 1 :, column]:
            ax.set_visible(False)
        axes[len(comparisons), column].set_xticks(x, names[group], rotation=90, fontsize="x-small")
        axes[len(comparisons), column].tick_params(labelbottom=True)
    fig.suptitle(
        f"Magnetic constraints; chi_sq of all magnetics = {results['CONSTRAINTS.CHI_SQ_MAG'][-1]:.3g}\n"
        "chi_sq and the differences are of RT-GSFit's last iteration"
    )
    _legend_below(fig, axes[0, 0])
    return fig


def plot_convergence(
    results: Results,
    reference: Reference,
    rtol: Mapping[str, float] | None = None,
    atol: Mapping[str, float] | None = None,
) -> Figure:
    """Every iteration: the change of psi and its difference from the reference, chi_sq, the relative difference of
    Ip and psi at the axis and boundary from the reference, the magnetic axis, and psi at the axis and boundary;
    the reference's values are dashed.

    The tests' tolerances are drawn when given, as maps from a results node to the tolerance of its last iteration:
    `rtol` relative, for `GLOBAL.IP`, `GLOBAL.PSI_MAG_AXIS` and `BOUNDARY.PSI`, as lines; `atol` absolute, for
    `GLOBAL.MAG_AXIS.R`, `GLOBAL.MAG_AXIS.Z`, `GLOBAL.PSI_MAG_AXIS` and `BOUNDARY.PSI`, as bands around the reference.
    """
    x = np.arange(1, len(results["TIME"]) + 1)
    psi = results["PROFILES_2D.R_Z.PSI"]
    fig = Figure(figsize=(15, 8), layout="constrained")
    axes = fig.subplots(2, 3, sharex=True)

    ax = axes[0, 0]
    ax.semilogy(x[1:], np.nanmax(np.abs(np.diff(psi, axis=0)), axis=(1, 2)), "o-", color=RTGSFIT,
                label="max |psi(i) - psi(i - 1)|")  # fmt: skip
    error = np.nanmax(np.abs(psi - _on_grid(results, reference, reference.psi)), axis=(1, 2))
    ax.semilogy(x, error, "s--", color=REFERENCE, mfc="none", label=f"max |psi - {reference.name}|")
    ax.set(title="psi [Wb]")

    ax = axes[0, 1]
    ax.semilogy(x, results["CONSTRAINTS.CHI_SQ_MAG"], "o-", color=RTGSFIT, label="magnetics")
    ax.semilogy(x, results["CONSTRAINTS.CHI_SQ_TOTAL"], "s--", color="k", mfc="none", label="with regularisation")
    ax.set(title="chi_sq")

    ax = axes[0, 2]
    relative = {"GLOBAL.IP": ("Ip", reference.ip), "GLOBAL.PSI_MAG_AXIS": ("psi at the axis", reference.axis_psi),
                "BOUNDARY.PSI": ("psi at the boundary", reference.boundary_psi)}  # fmt: skip
    for (key, (label, value)), style in zip(relative.items(), ("o-", "s-", "^-"), strict=True):
        if value:
            ax.semilogy(x, np.abs(results[key] / value - 1.0), style, color=RTGSFIT, mfc="none", label=label)
    tested: dict[float, list[str]] = {}
    for key, value in (rtol or {}).items():
        tested.setdefault(value, []).append(relative[key][0])
    for (value, labels), style in zip(tested.items(), (":", "-.", (0, (6, 2, 1, 2, 1, 2))), strict=False):
        ax.axhline(value, color="k", ls=style, label=f"test tolerance of {', '.join(labels)}")
    ax.set(title=f"|RT-GSFit / {reference.name} - 1|")

    panels = (
        (axes[1, 0], "Magnetic axis R [m]", {"GLOBAL.MAG_AXIS.R": ("RT-GSFit", reference.axis_r)}),
        (axes[1, 1], "Magnetic axis Z [m]", {"GLOBAL.MAG_AXIS.Z": ("RT-GSFit", reference.axis_z)}),
        (axes[1, 2], "psi [Wb]", {"GLOBAL.PSI_MAG_AXIS": ("RT-GSFit, axis", reference.axis_psi),
                                  "BOUNDARY.PSI": ("RT-GSFit, boundary", reference.boundary_psi)}),
    )  # fmt: skip
    for ax, title, curves in panels:
        for (key, (label, value)), style in zip(curves.items(), ("o-", "s-"), strict=False):
            ax.plot(x, results[key], style, color=RTGSFIT, mfc="none", label=label)
            _reference_line(ax, value, reference)
            if atol and key in atol:
                ax.axhspan(value - atol[key], value + atol[key], color=REFERENCE, alpha=0.15, label="test tolerance")
        ax.set(title=title)
    for ax in axes.flat:
        _shade_failures(ax, x, results)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles, strict=True))
        ax.legend(unique.values(), unique.keys(), fontsize="small")
    for ax in axes[1]:
        ax.set_xlabel("iteration (call of rtgsfit() with the same measurements)")
    fig.suptitle(f"Convergence of RT-GSFit, {reference.name} dashed")
    return fig


# ----------------------------------------------------------------------------------------------
# Figures for one test suite
# ----------------------------------------------------------------------------------------------
def plot_midplane(results: Results, plasma: Any) -> Figure:
    """psi, B_theta and j_phi along the grid row nearest the centre of the analytic plasma (an
    `AnalyticPlasma`) against the exact solution (top), and the difference (bottom), as in v1."""
    from rtgsfit.database_readers.analytic.solution import b_theta_from_psi

    r, z = results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"]
    row = int(np.argmin(np.abs(z - plasma.z_0)))
    r_fine = np.linspace(r[0], r[-1], 1000)
    shown = iterations(results)

    def rtgsfit(i: int) -> dict[str, np.ndarray]:
        psi = results["PROFILES_2D.R_Z.PSI"][i]
        b_theta = b_theta_from_psi(psi, r, z, plasma.r_0, plasma.z_0)
        return {"psi": psi[row], "b_theta": b_theta[row], "j_phi": results["PROFILES_2D.R_Z.J_PHI"][i][row]}

    panels = {"psi": ("psi", "Wb", 1.0), "b_theta": ("B_theta", "mT", 1e3), "j_phi": ("j_phi", "MA/m²", 1e-6)}
    fig = Figure(figsize=(15, 7.5), layout="constrained")
    axes = fig.subplots(2, 3, sharex=True)
    for column, (key, (label, unit, scale)) in enumerate(panels.items()):
        exact = getattr(plasma, key)
        ax = axes[0, column]
        ax.plot(r_fine, exact(r_fine, np.full_like(r_fine, z[row])) * scale, "k-", label="exact")
        for i, shade in zip(shown, _shades(len(shown)), strict=True):
            ax.plot(r, rtgsfit(i)[key] * scale, ".-", color=shade, label=f"RT-GSFit, iteration {i + 1}")
            axes[1, column].plot(r, (rtgsfit(i)[key] - exact(r, np.full_like(r, z[row]))) * scale, ".-",
                                 color=shade, label=f"iteration {i + 1}")  # fmt: skip
        for ax in axes[:, column]:
            ax.axvspan(plasma.r_0 - plasma.rho_b, plasma.r_0 + plasma.rho_b, color="0.9", zorder=0)
        axes[0, column].set(title=f"{label} [{unit}]")
        axes[1, column].set(title=f"RT-GSFit - exact [{unit}]", xlabel="R [m]")
        axes[1, column].axhline(0.0, color="k", lw=0.5)
    axes[0, 0].legend(fontsize="small")
    fig.suptitle(f"Midplane Z = {z[row]:.3g} m (grey: inside the plasma)")
    return fig


def plot_dofs(results: Results, names: Sequence[str], expected: Expected, title: str) -> Figure:
    """The fitted values of the DOFs `names` at the first, second and last iterations against `expected` (top), and
    the difference of the last iteration with the test's tolerance (bottom)."""
    all_names = list(results["FIT.DOF_NAMES"])
    columns = [all_names.index(name) for name in names]
    values = results["FIT.DOF_VALUES"][:, columns]
    unit = str(results["FIT.DOF_UNITS"][columns[0]])
    shown = iterations(results)
    x = np.arange(len(names))
    width = 0.8 / (len(shown) + 1)

    fig = Figure(figsize=(max(8.0, 0.6 * len(names) + 3), 7.5), layout="constrained")
    top, bottom = fig.subplots(2, 1, sharex=True)
    for k, (i, shade) in enumerate(zip(shown, _shades(len(shown)), strict=True)):
        top.bar(x + (k - len(shown) / 2) * width, values[i], width, color=shade, label=f"RT-GSFit, iteration {i + 1}")
    top.bar(x + len(shown) / 2 * width, expected.values, width, color=REFERENCE, label=expected.label)
    top.axhline(0.0, color="k", lw=0.5)
    top.set(title=title, ylabel=f"[{unit}]")
    top.legend(fontsize="small")

    error = values[-1] - expected.values
    tolerance = expected.tolerance()
    bottom.bar(x, error, 0.6, color=np.where(np.abs(error) > tolerance, "k", RTGSFIT),
               label=f"RT-GSFit (last iteration) - {expected.label}")  # fmt: skip
    tolerance_label = f"tolerance: rtol = {expected.rtol:g}, atol = {expected.atol:g}; black: outside"
    bottom.fill_between(x, -tolerance, tolerance, step="mid", color=REFERENCE, alpha=0.15, label=tolerance_label)
    bottom.axhline(0.0, color="k", lw=0.5)
    bottom.set(ylabel=f"difference [{unit}]")
    bottom.set_xticks(x, _short_names(names), rotation=90)
    bottom.legend(fontsize="small")
    return fig


def gsfit_passive_currents(results: Results, controller: Any) -> tuple[dict[str, Filaments], dict[str, Any]]:
    """The filaments of every GSFit passive structure, and their currents [A] from RT-GSFit's first, second and
    last iterations and from GSFit: title -> structure -> filament currents.

    Each passive DOF is a pattern of filament currents (GSFit's `current_distribution`) times its value.
    """
    passives = controller.passives
    names = list(results["FIT.DOF_NAMES"])
    shown = iterations(results)
    filaments: dict[str, Filaments] = {}
    currents: dict[str, dict[str, np.ndarray]] = {f"RT-GSFit, iteration {i + 1}": {} for i in shown}
    currents["GSFit"] = {}
    for passive in passives.keys():
        filaments[passive] = Filaments(
            *(np.ravel(passives.get_array1([passive, "geometry", key])) for key in ("r", "z", "d_r", "d_z"))
        )
        for current in currents.values():
            current[passive] = np.zeros(len(filaments[passive].r))
        for dof in passives.keys([passive, "dof"]):
            # The same names as the preshot reader (rtgsfit.database_readers.gsfit.preshot_inputs)
            name = passive if dof == "constant_current_density" else f"{passive}_{dof.upper()}"
            if name not in names:
                continue
            pattern = np.ravel(passives.get_array1([passive, "dof", dof, "current_distribution"]))
            value = results["FIT.DOF_VALUES"][:, names.index(name)]
            for i in shown:
                currents[f"RT-GSFit, iteration {i + 1}"][passive] += value[i] * pattern
            currents["GSFit"][passive] += passives.get_array1([passive, "dof", dof, "calculated"])[0] * pattern
    return filaments, currents


def plot_passive_currents(
    results: Results, filaments: Mapping[str, Filaments], currents: Mapping[str, Any], reference: str
) -> Figure:
    """The current density in every passive filament for each entry of `currents` (title -> structure -> filament
    currents [A]), and RT-GSFit's last iteration minus `reference` (one of the titles), with the net current of each
    structure.

    The current densities of the structures differ by orders of magnitude, so their colour scale is symmetric
    logarithmic, linear within 0.1 % of the largest.
    """
    last = f"RT-GSFit, iteration {iterations(results)[-1] + 1}"
    panels = dict(currents)
    difference_title = f"{last} - {reference}"
    panels[difference_title] = {name: currents[last][name] - currents[reference][name] for name in currents[last]}
    density = {
        title: {name: current / (filaments[name].d_r * filaments[name].d_z) * 1e-6 for name, current in panel.items()}
        for title, panel in panels.items()
    }

    def largest(titles: Sequence[str]) -> float:
        return max(float(np.max(np.abs(j))) for title in titles for j in density[title].values()) or 1.0

    j_max = largest(list(currents))
    norms = dict.fromkeys(currents, SymLogNorm(1e-3 * j_max, vmin=-j_max, vmax=j_max))
    norms[difference_title] = Normalize(-largest([difference_title]), largest([difference_title]))
    fig = Figure(figsize=(3.2 * len(panels) + 1.0, 8.0), layout="constrained")
    axes = fig.subplots(1, len(panels), sharex=True, sharey=True, squeeze=False)[0]
    for ax, (title, panel) in zip(axes, density.items(), strict=True):
        for name, j in panel.items():
            f = filaments[name]
            corners = zip(f.r - f.d_r / 2, f.z - f.d_z / 2, strict=True)
            rectangles = [Rectangle(corner, d_r, d_z) for corner, d_r, d_z in zip(corners, f.d_r, f.d_z, strict=True)]
            collection = PatchCollection(rectangles, cmap="RdBu_r", norm=norms[title])
            collection.set_array(j)
            ax.add_collection(collection)
        _plot_limiter(ax, results)
        net = "\n".join(f"{name:<9}{current.sum() * 1e-3:8.2f} kA" for name, current in panels[title].items())
        ax.text(0.02, 0.02, net, transform=ax.transAxes, fontsize="x-small", family="monospace", va="bottom",
                bbox={"facecolor": "w", "alpha": 0.8, "edgecolor": "0.7"})  # fmt: skip
        ax.autoscale_view()
        ax.set_aspect("equal")
        ax.set_title(title.replace(" - ", "\n- "), fontsize="medium")
        ax.set_xlabel("R [m]")
    axes[0].set_ylabel("Z [m]")
    fig.colorbar(ScalarMappable(norm=norms[last], cmap="RdBu_r"), ax=axes[:-1].tolist(), label="j [MA/m²]",
                 location="bottom", shrink=0.6)  # fmt: skip
    fig.colorbar(ScalarMappable(norm=norms[difference_title], cmap="RdBu_r"), ax=axes[-1], label="difference [MA/m²]",
                 location="bottom")  # fmt: skip
    fig.suptitle("Passive structures: current density and net current")
    return fig


def gsfit_sensor_positions(controller: Any, settings: Any) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """(R, Z) of the flux loops and Bp probes that RT-GSFit fits, from the GSFit controller."""
    from rtgsfit.database_readers.gsfit import constraint_sensors

    positions = {}
    for kind, label, objects in (("flux_loop", "flux loops", controller.flux_loops),
                                 ("bp_probe", "Bp probes", controller.bp_probes)):  # fmt: skip
        sensors = [sensor.name for sensor in constraint_sensors(controller, settings) if sensor.kind == kind]
        positions[label] = tuple(
            np.array([objects.get_array1([name, "geometry", key])[0] for name in sensors]) for key in ("r", "z")
        )
    return positions  # type: ignore[return-value]
