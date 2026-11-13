"""Helpers for `example_01_run_full_pipeline_st40.ipynb`: checks, summaries and plots.

The notebook calls the RT-GSFit pipeline (`rtgsfit.replay`) itself. This module holds only the code
that would distract from it: checking the environment and the MDSplus runs, collecting the machine
geometry from GSFit, and plotting.
"""

from __future__ import annotations

import importlib
import os
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import rtgsfit
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers.rtgsfit_mdsplus import DEFAULT_SERVER
from rtgsfit.database_writers.rtgsfit_mdsplus import check_standard_utility_csv
from rtgsfit.library import describe_status
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results
from rtgsfit.settings import Settings


# ----------------------------------------------------------------------------------------------
# Environment and MDSplus runs
# ----------------------------------------------------------------------------------------------
def check_environment(server: str = DEFAULT_SERVER) -> None:
    """Check the packages, the MDSplus server and the compiler that this example needs.

    Every problem is reported at once, in a single `RuntimeError`.
    """
    problems: list[str] = []
    print(f"rtgsfit {rtgsfit.__version__} from {Path(rtgsfit.__file__).parent}")
    # st40_database is needed by GSFit's ST40 reader, which builds the preshot
    for name in ("gsfit", "st40_database", "mdsthin", "standard_utility"):
        try:
            module = importlib.import_module(name)
        except ImportError as error:
            problems.append(f"cannot import {name} ({error}); install the TE extras: uv pip install -e '.[te]'")
            continue
        print(f"{name} {getattr(module, '__version__', '')} from {Path(str(module.__file__)).parent}")
        if name == "standard_utility":
            try:
                check_standard_utility_csv(module)
                print("standard_utility RTGSFIT.csv matches rtgsfit.data_dictionary")
            except RuntimeError as error:
                problems.append(str(error))
    try:
        import mdsthin

        with mdsthin.Connection(server) as connection:
            connection.get("1")
        print(f"MDSplus server {server} is reachable")
    except Exception as error:
        problems.append(f"cannot reach the MDSplus server {server}: {error}")
    for tool in ("make", os.environ.get("CC", "gcc")):
        if shutil.which(tool) is None:
            problems.append(f"{tool} is not on PATH; it is needed to compile librtgsfit.so")
    if problems:
        raise RuntimeError("This example cannot run:\n  " + "\n  ".join(problems))


def _text(value: Any) -> str:
    return (value.decode() if isinstance(value, bytes) else str(value)).strip()


def existing_runs(pulse: int, server: str = DEFAULT_SERVER) -> list[str]:
    """The runs under `\\RTGSFIT::TOP` at `pulse`; empty if there is no RTGSFIT tree for `pulse` yet."""
    import mdsthin

    with mdsthin.Connection(server) as connection:
        try:
            connection.openTree(dd.TREE_NAME, pulse)
            names = connection.get(f'getnci("\\\\{dd.TREE_NAME}::TOP.*", "NODE_NAME")').data()
        except (mdsthin.TreeFOPENR, mdsthin.TreeFILE_NOT_FOUND, mdsthin.TreeNNF):
            return []
    return sorted(_text(name) for name in np.atleast_1d(names))


def check_run_is_free(pulse: int, run_name: str, overwrite: bool = False, server: str = DEFAULT_SERVER) -> None:
    """Stop if `run_name` already exists at `pulse`, unless `overwrite` is true."""
    runs = existing_runs(pulse, server)
    if run_name not in runs:
        print(f"{dd.TREE_NAME} {pulse}: {run_name} is free (existing runs: {', '.join(runs) or 'none'})")
    elif overwrite:
        print(f"{dd.TREE_NAME} {pulse}: {run_name} exists and will be overwritten (OVERWRITE = True)")
    else:
        raise FileExistsError(
            f"{dd.TREE_NAME} {pulse} already has a run {run_name}. Choose another RUN_NAME, "
            "or set OVERWRITE = True to overwrite it."
        )


def run_description(pulse: int, run_name: str, server: str = DEFAULT_SERVER) -> str:
    """The run description, which standard_utility stores at `\\RTGSFIT::TOP.<run_name>:HELP`."""
    import mdsthin

    with mdsthin.Connection(server) as connection:
        connection.openTree(dd.TREE_NAME, pulse)
        return _text(connection.get(f"\\{dd.TREE_NAME}::TOP.{run_name}:HELP").data())


def _same(expected: Any, stored: Any) -> bool:
    expected, stored = np.asarray(expected), np.asarray(stored)
    if expected.shape != stored.shape:
        return False
    if expected.dtype.kind in "USO":
        return all(_text(a) == _text(b) for a, b in zip(expected.ravel(), stored.ravel(), strict=True))
    return bool(np.array_equal(expected, stored, equal_nan=True))


def compare_results(results: Results, stored: Results) -> None:
    """Check that every node of `stored` (read back from MDSplus) equals `results`."""
    expected = results.flat
    differences = [key for key, value in stored.flat.items() if key not in expected or not _same(expected[key], value)]
    missing = sorted(set(expected) - set(stored.flat))
    if differences or missing:
        raise AssertionError(f"nodes that differ: {differences}; nodes not read back: {missing}")
    print(f"All {len(expected)} nodes read back from MDSplus are identical to the results in memory")


# ----------------------------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------------------------
def _time_ranges(time_ms: np.ndarray, flags: np.ndarray) -> str:
    """`'7.6-12.6 ms, 280.6-299.6 ms'`: the contiguous time ranges where `flags` is true."""
    edges = np.flatnonzero(np.diff(np.r_[0, flags.astype(int), 0]))
    return ", ".join(f"{time_ms[a]:.1f}-{time_ms[b - 1]:.1f} ms" for a, b in zip(edges[::2], edges[1::2], strict=True))


def print_replay_summary(results: Results) -> None:
    """Status of every call, plasma current, fit quality and the time taken by each section of `rtgsfit()`."""
    time_ms = results["TIME"] * 1e3
    status = results["GLOBAL.STATUS"]
    ok = status == 0
    print(f"{results.n_time} calls of rtgsfit(), t = {time_ms[0]:.1f} to {time_ms[-1]:.1f} ms")
    for value in np.unique(status):
        flags = status == value
        print(f"  {describe_status(int(value)):<24}{np.count_nonzero(flags):5d} calls: {_time_ranges(time_ms, flags)}")

    ip = results["GLOBAL.IP"]
    peak = int(np.argmax(np.abs(np.where(ok, ip, 0.0))))
    print(f"Peak |Ip| = {abs(ip[peak]) * 1e-3:.1f} kA at t = {time_ms[peak]:.1f} ms")
    if np.any(ok):
        print(f"Median magnetic chi_sq of the successful calls: {np.median(results['CONSTRAINTS.CHI_SQ_MAG'][ok]):.3g}")
    print(f"LAPACK failures: {np.count_nonzero(results['FIT.LAPACK_INFO'])}")

    duration_us = results["TIMING.DURATION"][ok] * 1e6
    if not np.any(duration_us):
        print("No timing: the library was built without build/timing")
        return
    print(f"\nTime per call over the {np.count_nonzero(ok)} successful calls [µs]")
    print(f"{'section':<14}{'mean':>9}{'median':>9}{'p99':>9}{'max':>9}")
    for i, section in enumerate(results["TIMING.SECTIONS"]):
        column = duration_us[:, i]
        print(f"{section:<14}{column.mean():9.1f}{np.median(column):9.1f}"
              f"{np.percentile(column, 99):9.1f}{column.max():9.1f}")  # fmt: skip


# ----------------------------------------------------------------------------------------------
# Machine geometry from GSFit
# ----------------------------------------------------------------------------------------------
@dataclass
class Filaments:
    """Rectangular filaments of a PF coil or passive structure: centres and sizes [m]."""

    r: np.ndarray
    z: np.ndarray
    d_r: np.ndarray
    d_z: np.ndarray


@dataclass
class Geometry:
    """The geometry behind the preshot's Greens functions.

    The preshot holds only what `rtgsfit()` needs: the flux and field of each PF coil and passive degree of
    freedom (DOF), not the conductors that make them. This comes from the GSFit controller that built it.
    """

    coil_members: dict[str, list[str]]  # COIL:NAMES -> the GSFit PF coils in series in that circuit
    coils: dict[str, Filaments]  # GSFit PF coil -> filaments
    passives: dict[str, Filaments]  # GSFit passive structure -> filaments
    dof_passive: dict[str, str]  # passive DOF (FIT:DOF_NAMES) -> GSFit passive structure
    dof_current: dict[str, np.ndarray]  # passive DOF -> current in each filament of its structure when the DOF is 1
    dof_unit: dict[str, str]  # passive DOF -> its unit: A/m² for a constant current density, A for an eigenmode
    time_constant: dict[str, float]  # eigenmode DOF -> its L/R time constant [s]


def gsfit_geometry(settings: Settings) -> Geometry:
    """The PF coils, passive structures and passive DOFs of the GSFit controller that built the preshot.

    Each passive DOF is a fixed pattern of filament currents times an amplitude, the DOF. A structure with
    `constant_current_density` has one DOF, the current density [A/m²], so the pattern is the filament areas.
    The in-vessel components (IVC) instead have the slowest `eig` modes of their L/R circuit equation,
    `M x = tau R x`, with the filaments' mutual inductance matrix M and resistance matrix R. GSFit scales each
    eigenvector `x` to a 2-norm of 100, so the DOF is the amplitude of the mode in amperes.
    """
    from rtgsfit.database_readers.gsfit import coil_circuits
    from rtgsfit.database_readers.gsfit import gsfit_controller

    controller = gsfit_controller(settings)  # cached by build_preshot, so this does not set up GSFit again

    def filaments(objects: Any, path: list[str]) -> Filaments:
        return Filaments(*(np.ravel(objects.get_array1([*path, "geometry", key])) for key in ("r", "z", "d_r", "d_z")))

    passives = controller.passives
    geometry = Geometry(
        coil_members=coil_circuits(controller, settings).groups,
        coils={name: filaments(controller.coils, ["pf", name]) for name in controller.coils.keys(["pf"])},
        passives={name: filaments(passives, [name]) for name in passives.keys()},
        dof_passive={}, dof_current={}, dof_unit={}, time_constant={},
    )  # fmt: skip
    for passive in passives.keys():
        for dof in passives.keys([passive, "dof"]):
            # The same names as the preshot reader (rtgsfit.database_readers.gsfit.preshot_inputs)
            name = passive if dof == "constant_current_density" else f"{passive}_{dof.upper()}"
            geometry.dof_passive[name] = passive
            geometry.dof_current[name] = np.ravel(passives.get_array1([passive, "dof", dof, "current_distribution"]))
            geometry.dof_unit[name] = "A/m²" if dof == "constant_current_density" else "A"
            if "time_constant" in passives.keys([passive, "dof", dof]):
                geometry.time_constant[name] = float(passives.get_f64([passive, "dof", dof, "time_constant"]))
    return geometry


def up_down_parity(filaments: Filaments, current: np.ndarray) -> float:
    """+1 if the filament currents are up-down symmetric, -1 if antisymmetric (for an up-down symmetric structure)."""
    distance = np.hypot(filaments.r[:, None] - filaments.r, filaments.z[:, None] + filaments.z)
    mirror = np.argmin(distance, axis=1)
    return float(current @ current[mirror] / (current @ current))


# ----------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------
def _plot_limiter(ax: Axes, r: np.ndarray, z: np.ndarray, labels: bool = False) -> None:
    """The limiter contours as lines and the separate limiter points as dots.

    `rtgsfit()` uses only the limiter points. The ST40 limiter is a closed contour followed by two separate points
    on the MC coils, so a single line through every point would cross the vessel. The points are therefore split
    wherever the gap to the next point is much larger than the usual spacing.
    """
    r, z = np.asarray(r, dtype=float), np.asarray(z, dtype=float)
    step = np.hypot(np.diff(r), np.diff(z))
    max_step = 4.0 * np.median(step)
    gaps = np.flatnonzero(step > max_step) + 1
    pieces = list(zip(np.split(r, gaps), np.split(z, gaps), strict=True))
    contours = [(r_p, z_p) for r_p, z_p in pieces if len(r_p) > 1]
    points = [(r_p, z_p) for r_p, z_p in pieces if len(r_p) == 1]
    for k, (r_p, z_p) in enumerate(contours):
        if len(r_p) > 2 and np.hypot(r_p[-1] - r_p[0], z_p[-1] - z_p[0]) <= max_step:
            r_p, z_p = np.r_[r_p, r_p[:1]], np.r_[z_p, z_p[:1]]
        ax.plot(r_p, z_p, "k-", lw=1.5, label="limiter" if labels and k == 0 else None)
    if points:
        ax.plot(np.concatenate([p[0] for p in points]), np.concatenate([p[1] for p in points]), "ko", ms=4,
                label="separate limiter points" if labels else None)  # fmt: skip


def _short_name(name: str) -> str:
    """`PSI_FLOOP_007` -> `007`, `I_ROG_INIVC000` -> `INIVC000`."""
    return str(name).split("_", 2)[-1]


def _shade_failures(ax: Axes, time_ms: np.ndarray, ok: np.ndarray) -> None:
    ax.fill_between(time_ms, 0, 1, where=~ok, step="mid", color="tab:red", alpha=0.12, lw=0,
                    transform=ax.get_xaxis_transform(), label="status != OK")  # fmt: skip


def _grid_shape(preshot: Preshot) -> tuple[int, int]:
    return len(preshot["GRID:Z"]), len(preshot["GRID:R"])


def _si(value: float, unit: str) -> str:
    """`2.04e-7, "Wb"` -> `204 nWb`."""
    prefixes = (("", 1.0), ("m", 1e-3), ("µ", 1e-6), ("n", 1e-9), ("p", 1e-12))
    prefix, scale = next(((p, s) for p, s in prefixes if abs(value) >= s), prefixes[-1])
    return f"{value / scale:.3g} {prefix}{unit}"


def plot_preshot(preshot: Preshot) -> Figure:
    """The grid's allowed-plasma region and limiter points, and the initial plasma."""
    r, z = preshot["GRID:R"], preshot["GRID:Z"]
    shape = _grid_shape(preshot)
    r_2d, z_2d = np.meshgrid(r, z)
    allowed = preshot["GRID:MASK_LIMITER"].reshape(shape).astype(bool)

    fig, axes = plt.subplots(1, 2, figsize=(9, 6.5), sharey=True, constrained_layout=True)
    ax = axes[0]
    ax.plot(r_2d[~allowed], z_2d[~allowed], ".", color="0.7", ms=3, label="plasma current excluded")
    ax.plot(r_2d[allowed], z_2d[allowed], ".", color="tab:blue", ms=3, label="plasma current allowed")
    ax.plot(preshot["EDGE:R"], preshot["EDGE:Z"], "s", mfc="none", color="tab:orange", ms=4, label="edge point")
    ax.set_title(f"{len(r)} x {len(z)} grid")

    ax = axes[1]
    image = ax.contourf(r, z, preshot["INITIAL:PSI_NORM"].reshape(shape), levels=np.linspace(0.0, 1.0, 11))
    ax.contour(r, z, preshot["INITIAL:MASK"].reshape(shape), levels=[0.5], colors="w", linewidths=2)
    fig.colorbar(image, ax=ax, location="bottom", label="psi_norm (white: plasma mask)")
    ax.set_title("Initial psi_norm and mask")

    for ax in axes:
        _plot_limiter(ax, preshot["LIMITER:R"], preshot["LIMITER:Z"], labels=True)
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
    axes[0].set_ylabel("Z [m]")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize="small")
    return fig


def _plot_filaments(ax: Axes, filaments: Filaments, values: np.ndarray | None = None, **kwargs: Any) -> None:
    """Filaments as rectangles, in one colour or coloured by `values`."""
    corners = zip(filaments.r - filaments.d_r / 2, filaments.z - filaments.d_z / 2, strict=True)
    rectangles = [
        Rectangle(corner, d_r, d_z) for corner, d_r, d_z in zip(corners, filaments.d_r, filaments.d_z, strict=True)
    ]
    collection = PatchCollection(rectangles, **kwargs)
    if values is not None:
        collection.set_array(values)
    ax.add_collection(collection)


def _panels(n: int, n_columns: int, height: float) -> tuple[Figure, np.ndarray]:
    """Small multiples of R-Z plots; the panels after the first `n` are hidden."""
    n_rows = -(-n // n_columns)
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(3.2 * n_columns, height * n_rows), sharex=True, sharey=True,
                             squeeze=False, layout="compressed")  # fmt: skip
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("Z [m]")
    for ax in axes[-1]:
        ax.set_xlabel("R [m]")
    if n_rows > 1:
        # With sharex, only the bottom row has x tick labels; also label the panels above the hidden ones
        for ax in axes[-2, n % n_columns or n_columns :]:
            ax.xaxis.set_tick_params(labelbottom=True)
            ax.set_xlabel("R [m]")
    return fig, axes


def _set_limits(axes: np.ndarray, preshot: Preshot, filaments: Sequence[Filaments]) -> None:
    """Cover the grid and every filament drawn."""
    r = np.concatenate(
        [preshot["GRID:R"], *(f.r + f.d_r / 2 for f in filaments), *(f.r - f.d_r / 2 for f in filaments)]
    )
    z = np.concatenate(
        [preshot["GRID:Z"], *(f.z + f.d_z / 2 for f in filaments), *(f.z - f.d_z / 2 for f in filaments)]
    )
    margin = 0.03
    for ax in axes.flat:
        ax.set_aspect("equal")
        ax.set_xlim(min(r.min(), 0.0) - margin, r.max() + margin)
        ax.set_ylim(z.min() - margin, z.max() + margin)


def _circuit(preshot: Preshot, i_coil: int) -> str:
    """The current of coil `i_coil` from the PCS signals: row `i_coil` of COIL:MATRIX."""
    terms = []
    for coefficient, signal in zip(preshot["COIL:MATRIX"][i_coil], preshot["COIL:SIGNALS"], strict=True):
        if coefficient != 0.0:
            sign = "-" if coefficient < 0 else "+"
            size = "" if abs(coefficient) == 1.0 else f"{abs(coefficient):g} "
            terms.append(f"{sign} {size}{signal}")
    if not terms:
        return "I = 0: no signal"
    return "I = " + " ".join(terms).removeprefix("+ ")


def plot_coil_greens(preshot: Preshot, geometry: Geometry | None = None) -> Figure:
    """GRID:G_COIL, the poloidal flux per ampere of every PF coil circuit, with its coils in red."""
    r, z = preshot["GRID:R"], preshot["GRID:Z"]
    shape = _grid_shape(preshot)
    names = list(preshot["COIL:NAMES"])
    fig, axes = _panels(len(names), 4, 5.2)
    for i_coil, (ax, name) in enumerate(zip(axes.flat, names, strict=False)):
        image = ax.contourf(r, z, preshot["GRID:G_COIL"][:, i_coil].reshape(shape) * 1e6, levels=20, cmap="viridis")
        fig.colorbar(image, ax=ax, location="bottom", shrink=0.9, label="psi [µWb / A]", ticks=MaxNLocator(5))
        _plot_limiter(ax, preshot["LIMITER:R"], preshot["LIMITER:Z"])
        title = name
        if geometry is not None:
            members = geometry.coil_members[name]
            for coil, filaments in geometry.coils.items():
                in_circuit = coil in members
                _plot_filaments(ax, filaments, color="tab:red" if in_circuit else "0.75", zorder=3 if in_circuit else 2)
            if members != [name]:
                title += " = " + " + ".join(members)
        ax.set_title(f"{title}\n{_circuit(preshot, i_coil)}", fontsize="medium")
    _set_limits(axes, preshot, list(geometry.coils.values()) if geometry is not None else [])
    return fig


def plot_passive_greens(preshot: Preshot, geometry: Geometry, names: Sequence[str], n_columns: int = 5) -> Figure:
    """GRID:G_PASSIVE, the poloidal flux of each passive DOF in `names`, and the filament currents that make it.

    Both are divided by their largest magnitude in each panel, so they share one colour scale. The title gives
    that magnitude, for a DOF of 1.
    """
    r, z = preshot["GRID:R"], preshot["GRID:Z"]
    shape = _grid_shape(preshot)
    dof_names = list(preshot["FIT:DOF_NAMES"][int(preshot["SIZES:N_PLASMA"]) :])
    norm, cmap = Normalize(-1.0, 1.0), "RdBu_r"
    levels = np.linspace(-1.0, 1.0, 21)
    fig, axes = _panels(len(names), n_columns, 4.6)
    for ax, name in zip(axes.flat, names, strict=False):
        psi = preshot["GRID:G_PASSIVE"][:, dof_names.index(name)].reshape(shape)
        psi_max = np.abs(psi).max()
        ax.contourf(r, z, psi / psi_max, levels=levels, cmap=cmap, norm=norm)
        _plot_limiter(ax, preshot["LIMITER:R"], preshot["LIMITER:Z"])
        current = geometry.dof_current[name]
        _plot_filaments(ax, geometry.passives[geometry.dof_passive[name]], current / np.abs(current).max(),
                        cmap=cmap, norm=norm, edgecolor="k", linewidth=0.2, zorder=3)  # fmt: skip
        title = name
        if name in geometry.time_constant:
            title += f", tau = {geometry.time_constant[name] * 1e3:.1f} ms"
        ax.set_title(f"{title}\n|psi| ≤ {_si(psi_max, 'Wb')} per {geometry.dof_unit[name]}", fontsize="medium")
    _set_limits(axes, preshot, [geometry.passives[geometry.dof_passive[name]] for name in names])
    fig.colorbar(ScalarMappable(norm, cmap), ax=axes, location="bottom", shrink=0.5,
                 label="value / largest |value| in the panel: psi on the grid, current in the filaments")  # fmt: skip
    return fig


def _regularisation_weight(preshot: Preshot) -> np.ndarray:
    """The weight of the regularisation of each passive DOF: the size of its regularisation rows of FIT:G_PASSIVE."""
    n_constraint = int(preshot["SIZES:N_CONSTRAINT"])
    return np.linalg.norm(preshot["FIT:G_PASSIVE"][n_constraint:], axis=0)


def plot_eigenmodes(preshot: Preshot, geometry: Geometry, passive: str = "IVC") -> Figure:
    """For each eigenmode of `passive`: its L/R time constant, net current and regularisation weight."""
    dof_names = list(preshot["FIT:DOF_NAMES"][int(preshot["SIZES:N_PLASMA"]) :])
    names = [name for name in dof_names if geometry.dof_passive[name] == passive and name in geometry.time_constant]
    filaments = geometry.passives[passive]
    symmetric = np.array([up_down_parity(filaments, geometry.dof_current[name]) > 0 for name in names])
    colour = np.where(symmetric, "tab:blue", "tab:orange")
    weight = _regularisation_weight(preshot)
    panels = {
        "L/R time constant [ms]": [geometry.time_constant[name] * 1e3 for name in names],
        "Net current per ampere of amplitude": [geometry.dof_current[name].sum() for name in names],
        "Regularisation weight [1 / A]": [weight[dof_names.index(name)] for name in names],
    }
    labels = [name.removeprefix(f"{passive}_EIG_") for name in names]
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), sharex=True, constrained_layout=True)
    for ax, (title, values) in zip(axes, panels.items(), strict=True):
        ax.bar(labels, values, color=colour)
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set(title=title, xlabel=f"{passive} eigenmode")
    axes[0].legend(handles=[Patch(color="tab:blue", label="up-down symmetric"),
                            Patch(color="tab:orange", label="up-down antisymmetric")])  # fmt: skip
    return fig


def plot_passive_currents(results: Results, geometry: Geometry, passive: str = "IVC") -> Figure:
    """The fitted eigenmode amplitudes of `passive`, and the net current in every passive structure."""
    time_ms = results["TIME"] * 1e3
    ok = results["GLOBAL.STATUS"] == 0
    names = list(results["FIT.DOF_NAMES"])
    values = np.where(ok[:, None], results["FIT.DOF_VALUES"], np.nan)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    ax = axes[0]
    modes = [name for name in names if geometry.dof_passive.get(name) == passive and name in geometry.time_constant]
    colours = plt.get_cmap("viridis")(np.linspace(0.0, 0.9, len(modes)))
    for name, colour in zip(modes, colours, strict=True):
        symmetric = up_down_parity(geometry.passives[passive], geometry.dof_current[name]) > 0
        ax.plot(time_ms, values[:, names.index(name)], "-" if symmetric else "--", color=colour,
                label=name.removeprefix(f"{passive}_EIG_"))  # fmt: skip
    ax.set(title=f"Fitted {passive} eigenmode amplitudes (dashed: antisymmetric)", ylabel="amplitude [A]")
    ax.legend(ncol=8, fontsize="small", title="mode")

    ax = axes[1]
    net: dict[str, np.ndarray] = {}
    for name, structure in geometry.dof_passive.items():
        net[structure] = net.get(structure, 0.0) + values[:, names.index(name)] * geometry.dof_current[name].sum()
    for structure, current in net.items():
        ax.plot(time_ms, current * 1e-3, label=structure)
    ax.set(title="Net toroidal current in each passive structure", ylabel="I [kA]", xlabel="t [ms]")
    ax.legend(ncol=4, fontsize="small")
    for ax in axes:
        _shade_failures(ax, time_ms, ok)
    return fig


def plot_signals(signals: Signals, preshot: Preshot) -> Figure:
    """The inputs of `rtgsfit()`: `meas` grouped by signal type, the Rogowski constraints that sensor replacement
    makes from `meas`, and the PF coil currents."""
    names = np.asarray(preshot["MEAS:NAMES"]).astype(str)
    time_ms = signals.time * 1e3
    groups = {
        "PSI_FLOOP_": ("Flux loops", "psi [Wb]", 1.0),
        "B_BPPROBE_": ("Bp probes", "B [T]", 1.0),
        "I_ROG_": ("Rogowski coils", "I [kA]", 1e-3),
        "V_FLOOP_": ("Flux loop voltages", "V [V]", 1.0),
    }
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, constrained_layout=True)
    for ax, (prefix, (title, label, scale)) in zip(axes.flat, groups.items(), strict=False):
        columns = np.flatnonzero(np.char.startswith(names, prefix))
        lines = ax.plot(time_ms, signals.meas[:, columns] * scale, lw=0.8)
        ax.set(title=f"meas: {title} ({len(columns)})", ylabel=label)
        if len(columns) <= 10:
            ax.legend(lines, [_short_name(name) for name in names[columns]], ncol=2, fontsize="small")

    ax = axes[1, 1]
    constraint_names = np.asarray(preshot["CONSTRAINT:NAMES"]).astype(str)
    rows = np.flatnonzero(np.char.startswith(constraint_names, "I_ROG_"))
    lines = ax.plot(time_ms, signals.meas @ preshot["CONSTRAINT:REPLACEMENT"][rows].T * 1e-3, lw=0.8)
    ax.legend(lines, [_short_name(name) for name in constraint_names[rows]], ncol=2, fontsize="small")
    ax.set(title=f"Rogowski constraints ({len(rows)})\n= CONSTRAINT:REPLACEMENT meas", ylabel="I [kA]")

    ax = axes[1, 2]
    lines = ax.plot(time_ms, signals.coil_current(preshot) * 1e-3, lw=1.0)
    ax.legend(lines, list(preshot["COIL:NAMES"]), ncol=2, fontsize="small")
    ax.set(title="PF coil currents\n= COIL:MATRIX coil_signal", ylabel="I [kA]")
    for ax in axes[1]:
        ax.set_xlabel("t [ms]")
    return fig


def plot_time_traces(results: Results) -> Figure:
    """The main 0D quantities against time; failed calls are shaded."""
    time_ms = results["TIME"] * 1e3
    ok = results["GLOBAL.STATUS"] == 0

    def value(key: str) -> np.ndarray:
        return np.where(ok, results[key], np.nan)

    panels = {
        "Plasma current": ("Ip [kA]", {"GLOBAL.IP": ("Ip", 1e-3)}),
        "Major radius": ("R [m]", {"GLOBAL.MAG_AXIS.R": ("magnetic axis", 1.0),
                                   "GLOBAL.CURRENT_CENT.R": ("current centroid", 1.0),
                                   "BOUNDARY.GEO_AXIS.R": ("geometric axis", 1.0)}),
        "Height": ("Z [m]", {"GLOBAL.MAG_AXIS.Z": ("magnetic axis", 1.0),
                             "GLOBAL.CURRENT_CENT.Z": ("current centroid", 1.0),
                             "BOUNDARY.GEO_AXIS.Z": ("geometric axis", 1.0)}),
        "Poloidal flux": ("psi [Wb]", {"GLOBAL.PSI_MAG_AXIS": ("magnetic axis", 1.0),
                                       "BOUNDARY.PSI": ("boundary", 1.0)}),
        "Boundary shape": ("", {"BOUNDARY.ELONGATION": ("elongation", 1.0),
                                "BOUNDARY.TRIANG": ("triangularity", 1.0),
                                "BOUNDARY.MINOR_RADIUS": ("minor radius [m]", 1.0)}),
        "Fit quality": ("chi_sq", {"CONSTRAINTS.FLUX_LOOP.CHI_SQ": ("flux loops", 1.0),
                                   "CONSTRAINTS.BP_PROBE.CHI_SQ": ("Bp probes", 1.0),
                                   "CONSTRAINTS.ROGOWSKI.CHI_SQ": ("Rogowski coils", 1.0),
                                   "CONSTRAINTS.CHI_SQ_MAG": ("all magnetics", 1.0)}),
    }  # fmt: skip
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True, constrained_layout=True)
    for ax, (title, (label, curves)) in zip(axes.flat, panels.items(), strict=True):
        for key, (name, scale) in curves.items():
            ax.plot(time_ms, value(key) * scale, label=name)
        _shade_failures(ax, time_ms, ok)
        ax.set(title=title, ylabel=label)
        ax.legend(fontsize="small")
    axes[2, 1].set_yscale("log")
    for ax in axes[2]:
        ax.set_xlabel("t [ms]")
    return fig


def _index(results: Results, time: float) -> int:
    return int(np.argmin(np.abs(results["TIME"] - time)))


def plot_equilibria(results: Results, times: Sequence[float]) -> Figure:
    """Flux surfaces, toroidal current density, boundary, magnetic axis and x-points at `times` [s]."""
    r, z = results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"]
    indices = [_index(results, time) for time in times]
    j_max = max(float(np.max(np.abs(results["PROFILES_2D.R_Z.J_PHI"][i]))) for i in indices) * 1e-6

    fig, axes = plt.subplots(1, len(indices), figsize=(3.2 * len(indices) + 1, 6.5), sharey=True,
                             constrained_layout=True, squeeze=False)  # fmt: skip
    for ax, i in zip(axes[0], indices, strict=True):
        image = ax.pcolormesh(r, z, results["PROFILES_2D.R_Z.J_PHI"][i] * 1e-6, shading="nearest", cmap="RdBu_r",
                              vmin=-j_max, vmax=j_max)  # fmt: skip
        ax.contour(r, z, results["PROFILES_2D.R_Z.PSI"][i], levels=25, colors="0.45", linewidths=0.6)
        n = results["BOUNDARY.OUTLINE.N"][i]
        ax.plot(results["BOUNDARY.OUTLINE.R"][i, :n], results["BOUNDARY.OUTLINE.Z"][i, :n], "-", color="tab:green",
                lw=2, label="boundary")  # fmt: skip
        ax.plot(results["GLOBAL.MAG_AXIS.R"][i], results["GLOBAL.MAG_AXIS.Z"][i], "+", color="tab:green", ms=12,
                mew=2, label="magnetic axis")  # fmt: skip
        ax.plot(results["XPT.ALL.R"][i], results["XPT.ALL.Z"][i], "x", color="tab:purple", ms=7,
                label="x-points")  # fmt: skip
        _plot_limiter(ax, results["PRESHOT.LIMITER.R"], results["PRESHOT.LIMITER.Z"])
        ax.set_title(f"t = {results['TIME'][i] * 1e3:.1f} ms, Ip = {results['GLOBAL.IP'][i] * 1e-3:.0f} kA\n"
                     f"{describe_status(int(results['GLOBAL.STATUS'][i]))}", fontsize="medium")  # fmt: skip
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
    axes[0, 0].set_ylabel("Z [m]")
    axes[0, 0].legend(loc="lower left", fontsize="small")
    fig.colorbar(image, ax=axes[0].tolist(), label="j_phi [MA / m²]; grey: psi contours", shrink=0.8)
    return fig


def plot_constraints(results: Results, time: float) -> Figure:
    """Measured and reconstructed values of every magnetic constraint at `time` [s]."""
    i = _index(results, time)
    groups = {
        "FLUX_LOOP": ("Flux loops", "psi [Wb]", 1.0),
        "BP_PROBE": ("Bp probes", "B [mT]", 1e3),
        "ROGOWSKI": ("Rogowski coils", "I [kA]", 1e-3),
    }
    names = {group: [_short_name(name) for name in results[f"CONSTRAINTS.{group}.ALL.NAMES"]] for group in groups}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True,
                             gridspec_kw={"width_ratios": [max(len(n), 12) for n in names.values()]})  # fmt: skip
    for ax, (group, (title, label, scale)) in zip(axes, groups.items(), strict=True):
        node = f"CONSTRAINTS.{group}.ALL"
        x = np.arange(len(names[group]))
        ax.plot(x, results[f"{node}.MEASURED"][i] * scale, "o", mfc="none", ms=7, label="measured")
        ax.plot(x, results[f"{node}.RECONSTRUCT"][i] * scale, "x", ms=7, label="reconstructed")
        ax.set_xticks(x, names[group], rotation=90, fontsize="small")
        ax.set(title=f"{title}\nchi_sq = {results[f'CONSTRAINTS.{group}.CHI_SQ'][i]:.3g}", ylabel=label)
    axes[0].legend(fontsize="small")
    fig.suptitle(f"t = {results['TIME'][i] * 1e3:.1f} ms: magnetic chi_sq = {results['CONSTRAINTS.CHI_SQ_MAG'][i]:.3g}")
    return fig


def plot_timing(results: Results) -> Figure | None:
    """Time taken by each section of `rtgsfit()` over the successful calls, and the total per call."""
    ok = results["GLOBAL.STATUS"] == 0
    duration_us = results["TIMING.DURATION"][ok] * 1e6
    if not np.any(duration_us):
        print("No timing: the library was built without build/timing")
        return None
    sections = list(results["TIMING.SECTIONS"])

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True,
                                      gridspec_kw={"width_ratios": [2, 1]})  # fmt: skip
    left.boxplot(duration_us, showfliers=False)
    left.set_xticks(np.arange(1, len(sections) + 1), sections, rotation=45, ha="right")
    left.set_yscale("log")
    left.set(ylabel="time per call [µs]", title=f"Sections of rtgsfit() over {np.count_nonzero(ok)} successful calls")
    right.plot(results["TIME"][ok] * 1e3, duration_us[:, sections.index("TOTAL")], ".", ms=3)
    right.set(xlabel="t [ms]", ylabel="TOTAL [µs]", title="Total time of each call")
    return fig
