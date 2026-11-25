"""Helpers for `example_02_build_preshot_st40_tepcs.ipynb`: checks, comparisons and plots.

The notebook calls the RT-GSFit build, storage and C-library APIs itself. This module holds only the code
that would distract from it: checking the environment and the MDSplus runs, comparing preshots, collecting
the machine geometry from GSFit, and plotting.
"""

from __future__ import annotations

import importlib
import json
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
from rtgsfit.database_readers.rtgsfit_mdsplus import DEFAULT_SERVER
from rtgsfit.database_writers.rtgsfit_mdsplus import check_standard_utility_csv
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings


# ----------------------------------------------------------------------------------------------
# Environment and MDSplus runs
# ----------------------------------------------------------------------------------------------
def check_environment(server: str = DEFAULT_SERVER, *, write_mdsplus: bool = False) -> None:
    """Check the packages, the MDSplus server and the compiler that this example needs.

    Every problem is reported at once, in a single `RuntimeError`.
    """
    problems: list[str] = []
    print(f"rtgsfit {rtgsfit.__version__} from {Path(rtgsfit.__file__).parent}")
    packages = ["gsfit", "gsfit_rs", "shapely", "mdsthin"]
    if write_mdsplus:
        packages.append("standard_utility")
    for name in packages:
        try:
            module = importlib.import_module(name)
        except ImportError as error:
            problems.append(f"cannot import {name} ({error}); install GSFit and the TE extras (see README.md)")
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


def require_free_run(pulse: int, run_name: str, server: str = DEFAULT_SERVER) -> None:
    """Stop if `run_name` already exists at `pulse`: this example never overwrites a run."""
    runs = existing_runs(pulse, server)
    if run_name in runs:
        raise FileExistsError(f"{dd.TREE_NAME} {pulse} already has a run {run_name}; choose an unused PRIVATE_RUN")
    print(f"{dd.TREE_NAME} {pulse}: {run_name} is free (existing runs: {', '.join(runs) or 'none'})")


# ----------------------------------------------------------------------------------------------
# Comparing preshots
# ----------------------------------------------------------------------------------------------
RECIPE_SECTIONS = ("preshot", "gsfit")
"""The sections of the archived settings (`INFO:SETTINGS`) that determine the preshot."""


def assert_identical(actual: Preshot, expected: Preshot) -> None:
    """Every node of `actual` equals `expected` exactly, e.g. a preshot read back after writing it."""
    assert actual.groups.keys() == expected.groups.keys()
    for group, nodes in expected.groups.items():
        assert actual.groups[group].keys() == nodes.keys(), group
        for name, value in nodes.items():
            np.testing.assert_array_equal(
                np.asarray(actual.groups[group][name]), np.asarray(value), err_msg=f"{group}:{name}"
            )
    print(f"All {sum(len(nodes) for nodes in expected.groups.values())} nodes are identical")


def compare_preshots(built: Preshot, stored: Preshot, rtol: float = 1e-12) -> bool:
    """Check that `built` equals `stored`: floats to round-off, everything else exactly.

    `INFO:SETTINGS` must agree in `RECIPE_SECTIONS`; replay and writer settings do not change the preshot.
    Raises `AssertionError` on any real difference. Returns True if every node is bit-for-bit identical.
    """
    assert built.groups.keys() == stored.groups.keys(), "the preshots have different groups"
    n_exact, round_off = 0, {}
    for group, nodes in stored.groups.items():
        assert built.groups[group].keys() == nodes.keys(), f"{group} has different nodes"
        for name, value in nodes.items():
            key, new = f"{group}:{name}", built.groups[group][name]
            if key == "INFO:SETTINGS":
                archived, current = json.loads(str(value)), json.loads(str(new))
                for section in RECIPE_SECTIONS:
                    assert archived[section] == current[section], f"INFO:SETTINGS {section} differs"
                if archived == current:
                    n_exact += 1
                else:
                    print("INFO:SETTINGS: the recipe sections agree; other settings differ")
                continue
            expected, actual = np.asarray(value), np.asarray(new)
            if expected.dtype.kind == "f":
                scale = float(np.abs(expected).max(initial=0.0))
                np.testing.assert_allclose(actual, expected, rtol=rtol, atol=rtol * scale, err_msg=key)
                if np.array_equal(actual, expected):
                    n_exact += 1
                else:
                    round_off[key] = float(np.abs(actual - expected).max() / scale)
            else:
                np.testing.assert_array_equal(actual, expected, err_msg=key)
                n_exact += 1
    n_nodes = sum(len(nodes) for nodes in stored.groups.values())
    print(f"{n_exact} of {n_nodes} nodes are bit-for-bit identical")
    for key, difference in round_off.items():
        print(f"  {key}: differs by round-off, at most {difference:.1e} of its largest value")
    return n_exact == n_nodes


def print_pcs_inputs(preshot: Preshot) -> None:
    """What the caller of `rtgsfit()` must supply each cycle, in the order the preshot fixes."""
    names = np.asarray(preshot["CONSTRAINT:MAG_SIGNALS"]).astype(str)
    prefixes = {"PSI_FLOOP_": "flux loops", "B_BPPROBE_": "Bp probes", "I_ROG_": "Rogowski coils",
                "V_FLOOP_": "flux loop voltages"}  # fmt: skip
    counts = ", ".join(
        f"{np.count_nonzero(np.char.startswith(names, prefix))} {label}" for prefix, label in prefixes.items()
    )
    print(f"mag_signal: {len(names)} values in CONSTRAINT:MAG_SIGNALS order: {counts}")
    print(f"  first {', '.join(names[:3])}; last {', '.join(names[-3:])}")
    signals = list(preshot["COIL:COIL_SIGNALS"])
    print(f"coil signals: {len(signals)} PCS channels in COIL:COIL_SIGNALS order: {', '.join(signals)}")
    print("coil_current = coil_signal @ COIL:MATRIX, one value per COIL:NAMES:")
    for i_coil, name in enumerate(preshot["COIL:NAMES"]):
        print(f"  {name:<6}{_circuit(preshot, i_coil)}")


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
    """The current of coil `i_coil` from the PCS signals: column `i_coil` of COIL:MATRIX."""
    terms = []
    for coefficient, signal in zip(preshot["COIL:MATRIX"][:, i_coil], preshot["COIL:COIL_SIGNALS"], strict=True):
        if coefficient != 0.0:
            sign = "-" if coefficient < 0 else "+"
            size = "" if abs(coefficient) == 1.0 else f"{abs(coefficient):g} "
            terms.append(f"{sign} {size}{signal}")
    if not terms:
        return "I = 0: no signal"
    return "I = " + " ".join(terms).removeprefix("+ ")


def plot_coil_greens(preshot: Preshot, geometry: Geometry) -> Figure:
    """GRID:G_COIL, the poloidal flux per ampere of every PF coil circuit, with its coils in red."""
    r, z = preshot["GRID:R"], preshot["GRID:Z"]
    shape = _grid_shape(preshot)
    names = list(preshot["COIL:NAMES"])
    fig, axes = _panels(len(names), 4, 5.2)
    for i_coil, (ax, name) in enumerate(zip(axes.flat, names, strict=False)):
        image = ax.contourf(r, z, preshot["GRID:G_COIL"][:, i_coil].reshape(shape) * 1e6, levels=20, cmap="viridis")
        fig.colorbar(image, ax=ax, location="bottom", shrink=0.9, label="psi [µWb / A]", ticks=MaxNLocator(5))
        _plot_limiter(ax, preshot["LIMITER:R"], preshot["LIMITER:Z"])
        members = geometry.coil_members[name]
        for coil, filaments in geometry.coils.items():
            in_circuit = coil in members
            _plot_filaments(ax, filaments, color="tab:red" if in_circuit else "0.75", zorder=3 if in_circuit else 2)
        title = name if members == [name] else f"{name} = " + " + ".join(members)
        ax.set_title(f"{title}\n{_circuit(preshot, i_coil)}", fontsize="medium")
    _set_limits(axes, preshot, list(geometry.coils.values()))
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
