"""The preshot: every quantity RT-GSFit needs that can be computed before the pulse.

Real-time use needs `rtgsfit()` to be as fast as possible, so all the work that does not depend on
the measurements is done here, before the pulse, and compiled into C as constants. That includes
the Greens functions, the LU factors of the Grad-Shafranov operator and the weighted response
matrices of the least-squares problem.

A `Preshot` is a nested dictionary `GROUP -> NODE -> value` with exactly the nodes of
`data_dictionary.PRESHOT_NODES`; it is stored in MDSplus at `PRESHOT.GROUP:NODE`. Database readers
compute the physics-specific parts (Greens functions, weights, ...) and call `assemble_preshot`,
which derives everything else, so the derived quantities are computed in one place only.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any
from typing import Literal

import numpy as np
import numpy.typing as npt
import scipy.linalg

from rtgsfit import data_dictionary as dd
from rtgsfit.constants import MU_0
from rtgsfit.constants import N_INTERP
from rtgsfit.constants import N_PLASMA
from rtgsfit.constants import PLASMA_DOF_NAMES

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int32]


class PreshotError(ValueError):
    """The preshot is inconsistent."""


LEGACY_PRESHOT_NAMES = {
    ("SIZES", "N_FIT_TARGET"): "N_B_VEC",
    ("SIZES", "N_FIT_DOF"): "N_X_VEC",
    ("FIT", "DOF_NAMES"): "X_VEC_NAMES",
    ("FIT", "TARGET_NAMES"): "B_VEC_NAMES",
    ("GRID", "INV_R"): "INV_MU0_R",
}


def migrate_legacy_preshot(groups: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Convert earlier preshot names/basis contracts without changing the original mapping."""
    if "N_FIT_COEF" in groups["SIZES"]:
        groups = {group: dict(nodes) for group, nodes in groups.items()}
        for group, old_name, name in (
            ("SIZES", "N_FIT_COEF", "N_FIT_DOF"),
            ("FIT", "COEF_NAMES", "DOF_NAMES"),
        ):
            if old_name not in groups[group] or name in groups[group]:
                raise PreshotError(f"ambiguous or incomplete coefficient-name preshot at {group}:{old_name}")
            groups[group][name] = groups[group].pop(old_name)
        return groups
    if "N_X_VEC" not in groups["SIZES"]:
        return groups
    groups = {group: dict(nodes) for group, nodes in groups.items()}
    for (group, name), old_name in LEGACY_PRESHOT_NAMES.items():
        groups[group][name] = groups[group].pop(old_name)
    r_flat = np.asarray(groups["GRID"]["R_FLAT"])
    if not np.allclose(groups["GRID"]["INV_R"], 1.0 / (r_flat * MU_0), rtol=1e-12, atol=0):
        raise PreshotError("legacy GRID:INV_MU0_R does not match the supported mu0 convention")
    groups["GRID"]["INV_R"] = 1.0 / r_flat
    groups["FIT"]["DOF_NAMES"] = ["FF_PRIME_MU0" if name == "FF_PRIME" else name for name in groups["FIT"]["DOF_NAMES"]]
    source = groups["INFO"]["SOURCE"]
    if source != "analytic" and not source.startswith(("gsfit ", "freegsnke ")):
        raise PreshotError(f"unknown legacy flux convention for source {source!r}")
    groups["INFO"]["MU0"] = MU_0
    groups["INFO"]["PSI_TO_WB"] = 2.0 * np.pi if source == "analytic" else 1.0
    return groups


class Preshot:
    """The preshot, as `GROUP -> NODE -> value`, mirroring `PRESHOT.GROUP:NODE` in MDSplus."""

    def __init__(self, groups: dict[str, dict[str, Any]]) -> None:
        self.groups = {group: dict(nodes) for group, nodes in groups.items()}
        self.validate()

    # Access -----------------------------------------------------------------------------------
    def __getitem__(self, key: str) -> Any:
        """`preshot["GRID:R"]` or `preshot["SIZES:N_R"]`."""
        group, _, node = key.partition(":")
        return self.groups[group][node]

    @property
    def sizes(self) -> dict[str, int]:
        return {name: int(value) for name, value in self.groups["SIZES"].items()}

    def to_nested(self) -> dict[str, Any]:
        """`{"PRESHOT": {...}}`, the layout of the MDSplus tree."""
        return {"PRESHOT": {group: dict(nodes) for group, nodes in self.groups.items()}}

    @classmethod
    def from_nested(cls, nested: dict[str, Any]) -> Preshot:
        groups = nested["PRESHOT"] if "PRESHOT" in nested else nested
        groups = migrate_legacy_preshot(groups)
        converted: dict[str, dict[str, Any]] = {}
        for node in dd.PRESHOT_NODES:
            _, group, name = node.path
            converted.setdefault(group, {})[name] = _convert(groups[group][name], node)
        return cls(converted)

    def content_hash(self) -> str:
        """Hash of every quantity that is compiled into C; used to cache the compiled library."""
        digest = hashlib.sha256()
        for node in dd.C_NODES:
            _, group, name = node.path
            value = np.ascontiguousarray(self.groups[group][name])
            digest.update(node.key.encode())
            digest.update(str(value.dtype).encode())
            digest.update(value.tobytes())
        return digest.hexdigest()[:16]

    def summary(self) -> str:
        return ", ".join(f"{name}={value}" for name, value in self.sizes.items())

    # Validation -------------------------------------------------------------------------------
    def symbolic_sizes(self) -> dict[str, int]:
        sizes = self.sizes
        sizes["N_INTERP"] = N_INTERP
        sizes["N_R_PLUS_2"] = sizes["N_R"] + 2
        sizes["N_COIL_SIGNAL"] = len(self.groups["COIL"]["SIGNALS"])
        return sizes

    def validate(self) -> None:
        """Check every quantity against the data dictionary and the relations between sizes."""
        for node in dd.PRESHOT_NODES:
            _, group, name = node.path
            if group not in self.groups or name not in self.groups[group]:
                raise PreshotError(f"missing {node.mds_path}")
        expected_keys = {node.path[1:] for node in dd.PRESHOT_NODES}
        for group, nodes in self.groups.items():
            for name in nodes:
                if (group, name) not in expected_keys:
                    raise PreshotError(f"PRESHOT.{group}:{name} is not in the data dictionary")

        sizes = self.symbolic_sizes()
        _check_size_relations(sizes)
        for node in dd.PRESHOT_NODES:
            _, group, name = node.path
            value = self.groups[group][name]
            if node.dtype == "str":
                if node.shape:
                    expected_length = sizes[node.shape[0]]
                    if not isinstance(value, list) or len(value) != expected_length:
                        raise PreshotError(f"{node.mds_path} must be a list of {expected_length} strings")
                elif not isinstance(value, str):
                    raise PreshotError(f"{node.mds_path} must be a string")
                continue
            array = np.asarray(value)
            expected_shape = tuple(sizes[dim] for dim in node.shape)
            if array.shape != expected_shape:
                raise PreshotError(f"{node.mds_path} has shape {array.shape}, expected {expected_shape} = {node.shape}")
            if array.dtype != np.dtype(node.dtype):
                raise PreshotError(f"{node.mds_path} has dtype {array.dtype}, expected {node.dtype}")
            if node.dtype == "float64" and not np.all(np.isfinite(array)):
                raise PreshotError(f"{node.mds_path} is not finite")

        for name in ("MU0", "PSI_TO_WB"):
            if float(self["INFO:" + name]) <= 0:
                raise PreshotError(f"INFO:{name} must be positive")
        grid_r, grid_z = self.groups["GRID"]["R"], self.groups["GRID"]["Z"]
        d_r, d_z = float(self.groups["GRID"]["D_R"]), float(self.groups["GRID"]["D_Z"])
        if not (np.allclose(np.diff(grid_r), d_r, rtol=1e-9) and np.allclose(np.diff(grid_z), d_z, rtol=1e-9)):
            raise PreshotError("the grid must be uniform with spacing D_R, D_Z")
        if grid_r[0] <= 0.0:
            raise PreshotError("the grid must have R > 0")
        for name in ("MASK_LIMITER",):
            if not np.all(np.isin(self.groups["GRID"][name], (0, 1))):
                raise PreshotError(f"PRESHOT.GRID:{name} must be 0 or 1")
        if not np.all(np.isin(self.groups["INITIAL"]["MASK"], (0, 1))):
            raise PreshotError("PRESHOT.INITIAL:MASK must be 0 or 1")
        index = self.groups["LIMITER"]["GRID_INDEX"]
        if index.min() < 0 or index.max() >= sizes["N_GRID"]:
            raise PreshotError("PRESHOT.LIMITER:GRID_INDEX is outside the grid")
        permutation = self.groups["POISSON"]["PERMUTATION"]
        if not np.array_equal(np.sort(permutation), np.arange(sizes["N_GRID"])):
            raise PreshotError("PRESHOT.POISSON:PERMUTATION is not a permutation")
        fit_target_names = self.groups["FIT"]["TARGET_NAMES"]
        if fit_target_names[: sizes["N_CONSTRAINT"]] != self.groups["CONSTRAINT"]["NAMES"]:
            raise PreshotError("the first N_CONSTRAINT PRESHOT.FIT:TARGET_NAMES must be PRESHOT.CONSTRAINT:NAMES")
        if np.any(self.groups["FIT"]["WEIGHT"][: sizes["N_CONSTRAINT"]] == 0.0):
            raise PreshotError("constraint weights must be non-zero (fit_prediction is divided by them)")


def _convert(value: Any, node: dd.Node) -> Any:
    """Convert a value read from a database into the type given by the data dictionary."""
    if node.dtype == "str":
        if node.shape:
            return [str(item).strip() for item in np.atleast_1d(value)]
        return str(value)
    array = np.asarray(value, dtype=node.dtype)
    return array.reshape(()) if not node.shape else array


def _check_size_relations(sizes: dict[str, int]) -> None:
    """The same relations as the `_Static_assert`s in `c/src/preshot.h`."""
    relations = {
        "N_GRID == N_R * N_Z": sizes["N_GRID"] == sizes["N_R"] * sizes["N_Z"],
        "N_EDGE == 2 N_R + 2 N_Z - 4": sizes["N_EDGE"] == 2 * sizes["N_R"] + 2 * sizes["N_Z"] - 4,
        "N_R, N_Z >= 4": sizes["N_R"] >= 4 and sizes["N_Z"] >= 4,
        "N_CONSTRAINT == N_FLUX_LOOP + N_BP_PROBE + N_ROGOWSKI": (
            sizes["N_CONSTRAINT"] == sizes["N_FLUX_LOOP"] + sizes["N_BP_PROBE"] + sizes["N_ROGOWSKI"]
        ),
        "N_FIT_TARGET == N_CONSTRAINT + N_REG": sizes["N_FIT_TARGET"] == sizes["N_CONSTRAINT"] + sizes["N_REG"],
        f"N_PLASMA == {N_PLASMA}": sizes["N_PLASMA"] == N_PLASMA,
        "N_FIT_DOF == N_PLASMA + N_PASSIVE": sizes["N_FIT_DOF"] == sizes["N_PLASMA"] + sizes["N_PASSIVE"],
        "N_FIT_TARGET >= N_FIT_DOF": sizes["N_FIT_TARGET"] >= sizes["N_FIT_DOF"],
        "N_LIMITER, N_COIL, N_MEAS, N_XPT_MAX >= 1": min(
            sizes["N_LIMITER"], sizes["N_COIL"], sizes["N_MEAS"], sizes["N_XPT_MAX"]
        )
        >= 1,
    }
    for relation, holds in relations.items():
        if not holds:
            raise PreshotError(f"size relation violated: {relation}")


# ----------------------------------------------------------------------------------------------
# Geometry helpers shared by the database readers
# ----------------------------------------------------------------------------------------------
def rectangle_points(r: FloatArray, z: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Points on the boundary of the rectangle spanned by `r` and `z`, without repeats.

    This is the "edge order" used throughout RT-GSFit: start at the bottom-left corner, go up the
    left side, right along the top, down the right side and left along the bottom. There are
    `2 len(r) + 2 len(z) - 4` points.
    """
    n_r, n_z = len(r), len(z)
    edge_r = np.concatenate(
        ([r[0]], np.full(n_z - 2, r[0]), [r[0]], r[1:-1], [r[-1]], np.full(n_z - 2, r[-1]), [r[-1]], r[-2:0:-1])
    )
    edge_z = np.concatenate(
        ([z[0]], z[1:-1], [z[-1]], np.full(n_r - 2, z[-1]), [z[-1]], z[-2:0:-1], [z[0]], np.full(n_r - 2, z[0]))
    )
    return edge_r, edge_z


def edge_grid_index(n_r: int, n_z: int) -> npt.NDArray[np.int64]:
    """Flat grid index (`i_z * n_r + i_r`) of each edge point, in edge order."""
    i_r, i_z = rectangle_points(np.arange(n_r, dtype=np.float64), np.arange(n_z, dtype=np.float64))
    return i_z.astype(np.int64) * n_r + i_r.astype(np.int64)


def limiter_interpolation(
    grid_r: FloatArray, grid_z: FloatArray, limiter_r: FloatArray, limiter_z: FloatArray
) -> tuple[IntArray, FloatArray]:
    """Bilinear interpolation stencil of each limiter point.

    Returns `(index, weight)`, both `[N_LIMITER, 4]`, ordered (i_r, i_z), (i_r + 1, i_z),
    (i_r, i_z + 1), (i_r + 1, i_z + 1), so that `limiter_psi = sum(weight * psi[index], axis=1)`.
    """
    if np.any((limiter_r < grid_r[0]) | (limiter_r > grid_r[-1]) | (limiter_z < grid_z[0]) | (limiter_z > grid_z[-1])):
        raise PreshotError("every limiter point must be inside the grid")
    n_r = len(grid_r)
    i_r = np.clip(np.searchsorted(grid_r, limiter_r) - 1, 0, n_r - 2)
    i_z = np.clip(np.searchsorted(grid_z, limiter_z) - 1, 0, len(grid_z) - 2)
    r_0, r_1 = grid_r[i_r], grid_r[i_r + 1]
    z_0, z_1 = grid_z[i_z], grid_z[i_z + 1]
    area = (r_1 - r_0) * (z_1 - z_0)
    index = np.stack([n_r * i_z + i_r, n_r * i_z + i_r + 1, n_r * (i_z + 1) + i_r, n_r * (i_z + 1) + i_r + 1], axis=1)
    weight = (
        np.stack(
            [
                (r_1 - limiter_r) * (z_1 - limiter_z),
                (limiter_r - r_0) * (z_1 - limiter_z),
                (r_1 - limiter_r) * (limiter_z - z_0),
                (limiter_r - r_0) * (limiter_z - z_0),
            ],
            axis=1,
        )
        / area[:, np.newaxis]
    )
    return index.astype(np.int32), weight


PoissonOperator = Literal["grad_shafranov", "laplacian"]


def poisson_matrix(grid_r: FloatArray, grid_z: FloatArray, operator: PoissonOperator) -> FloatArray:
    """Dense matrix of the discrete Grad-Shafranov operator, scaled by `d_z**2`.

    Interior rows are eq. (45) of Moret et al., Fusion Eng. Des. 91 (2015) 1:

        psi[i_z+1, i_r] + psi[i_z-1, i_r] + a psi[i_z, i_r+1] + b psi[i_z, i_r-1] - (2 + a + b) psi[i_z, i_r]

    with `a = (d_z/d_r)**2 R / (R + d_r/2)` and `b = (d_z/d_r)**2 R / (R - d_r/2)` for the
    Grad-Shafranov operator Delta*, or `a = b = (d_z/d_r)**2` for the Laplacian (the large aspect
    ratio limit used by the analytic test). Edge rows are identity rows: the edge values are set by
    the free-boundary condition. Returned in Fortran order so that it can be factorised in place.
    """
    n_r, n_z = len(grid_r), len(grid_z)
    n_grid = n_r * n_z
    d_r, d_z = grid_r[1] - grid_r[0], grid_z[1] - grid_z[0]
    if operator == "grad_shafranov":
        a = (d_z / d_r) ** 2 * grid_r / (grid_r + d_r / 2)
        b = (d_z / d_r) ** 2 * grid_r / (grid_r - d_r / 2)
    elif operator == "laplacian":
        a = np.full(n_r, (d_z / d_r) ** 2)
        b = np.full(n_r, (d_z / d_r) ** 2)
    else:
        raise ValueError(f"unknown Poisson operator {operator!r}")
    c = 2.0 + a + b

    matrix = np.zeros((n_grid, n_grid), dtype=np.float64, order="F")
    i_z, i_r = np.divmod(np.arange(n_grid), n_r)
    interior = (i_z > 0) & (i_z < n_z - 1) & (i_r > 0) & (i_r < n_r - 1)
    rows = np.flatnonzero(interior)
    edge = np.flatnonzero(~interior)
    matrix[edge, edge] = 1.0
    matrix[rows, rows] = -c[i_r[rows]]
    matrix[rows, rows + 1] = a[i_r[rows]]
    matrix[rows, rows - 1] = b[i_r[rows]]
    matrix[rows, rows + n_r] = 1.0
    matrix[rows, rows - n_r] = 1.0
    return matrix


def poisson_lu_bands(matrix: FloatArray, n_r: int) -> tuple[FloatArray, FloatArray, IntArray]:
    """Banded LU factors of `matrix` with partial pivoting, as used by `c/src/poisson_solver.c`.

    Returns `(lower_band [N_GRID, N_R], upper_band [N_GRID, N_R + 2], permutation [N_GRID])` with

    * `lower_band[i, k] = L[i, i - N_R + k]` (the unit diagonal of L is not stored);
    * `upper_band[i, k] = U[i, i + k]`, except `upper_band[i, 0] = 1 / U[i, i]`;
    * `(matrix[permutation])` = L U.

    `matrix` is overwritten. Raises PreshotError if L or U do not fit in the bands.
    """
    n_grid = matrix.shape[0]
    lu, pivots = scipy.linalg.lu_factor(matrix, overwrite_a=True, check_finite=False)
    permutation = np.arange(n_grid)
    for i, pivot in enumerate(pivots):
        permutation[i], permutation[pivot] = permutation[pivot], permutation[i]

    rows = np.arange(n_grid)[:, np.newaxis]
    lower_columns = rows - n_r + np.arange(n_r)[np.newaxis, :]
    upper_columns = rows + np.arange(n_r + 2)[np.newaxis, :]
    lower_band = np.where(lower_columns >= 0, lu[rows, np.clip(lower_columns, 0, n_grid - 1)], 0.0)
    upper_band = np.where(upper_columns < n_grid, lu[rows, np.clip(upper_columns, 0, n_grid - 1)], 0.0)

    block = 512
    for start in range(0, n_grid, block):
        stop = min(start + block, n_grid)
        part = lu[start:stop]
        outside_lower = np.tril(part, start - n_r - 1)
        outside_upper = np.triu(part, start + n_r + 2)
        if np.any(outside_lower != 0.0) or np.any(outside_upper != 0.0):
            raise PreshotError("the LU factors of the Poisson matrix do not fit in the bands")

    if np.any(upper_band[:, 0] == 0.0):
        raise PreshotError("the Poisson matrix is singular")
    upper_band[:, 0] = 1.0 / upper_band[:, 0]
    return lower_band, upper_band, permutation.astype(np.int32)


@dataclass(frozen=True)
class SensorReplacement:
    """Replace the faulty sensor `replaced` by `sum(coefficients * replacements)`."""

    replaced: str
    replacements: tuple[str, ...]
    coefficients: tuple[float, ...]


def replacement_matrix(
    constraint_names: Sequence[str], replacements: Sequence[SensorReplacement]
) -> tuple[FloatArray, list[str]]:
    """Build `CONSTRAINT:REPLACEMENT` and `MEAS:NAMES`.

    meas holds the constraints that are not replaced, in order, followed by every sensor used in a
    replacement that is not already a constraint. Each constraint row is either a unit row or the
    replacement coefficients.
    """
    by_name = {item.replaced: item for item in replacements}
    unknown = set(by_name) - set(constraint_names)
    if unknown:
        raise PreshotError(f"cannot replace sensors that are not constraints: {sorted(unknown)}")
    meas_names = [name for name in constraint_names if name not in by_name]
    for item in replacements:
        for name in item.replacements:
            if name in by_name:
                raise PreshotError(f"{name} is replaced, so it cannot be used as a replacement")
            if name not in meas_names:
                meas_names.append(name)

    matrix = np.zeros((len(constraint_names), len(meas_names)))
    for i, name in enumerate(constraint_names):
        if name in by_name:
            item = by_name[name]
            if len(item.replacements) != len(item.coefficients):
                raise PreshotError(f"replacement of {name}: names and coefficients differ in length")
            for replacement, coefficient in zip(item.replacements, item.coefficients, strict=True):
                matrix[i, meas_names.index(replacement)] += coefficient
        else:
            matrix[i, meas_names.index(name)] = 1.0
    return matrix, meas_names


# ----------------------------------------------------------------------------------------------
# Assembly
# ----------------------------------------------------------------------------------------------
@dataclass
class PreshotInputs:
    """The machine- and diagnostic-specific inputs from which `assemble_preshot` builds a Preshot.

    Constraints are ordered flux loops, then Bp probes, then Rogowski coils. Greens functions are
    unweighted: `g_constraint_grid[i, k]` is constraint i per ampere at grid point k.
    """

    mu_0: float
    grid_r: FloatArray
    grid_z: FloatArray
    poisson_operator: PoissonOperator
    poisson_rhs_scale: FloatArray  # [N_GRID] (2 pi / psi_to_weber) mu0 R d_z^2; rhs = -scale * j_phi
    mask_limiter: npt.NDArray[np.bool_]  # [N_GRID]
    limiter_r: FloatArray
    limiter_z: FloatArray
    g_edge_edge: FloatArray  # [N_EDGE, N_EDGE]
    coil_names: list[str]
    g_grid_coil: FloatArray  # [N_GRID, N_COIL]
    passive_names: list[str]
    g_grid_passive: FloatArray  # [N_GRID, N_PASSIVE]
    constraint_names: list[str]
    n_flux_loop: int
    n_bp_probe: int
    n_rogowski: int
    g_constraint_grid: FloatArray  # [N_CONSTRAINT, N_GRID]
    g_constraint_coil: FloatArray  # [N_CONSTRAINT, N_COIL]
    g_constraint_passive: FloatArray  # [N_CONSTRAINT, N_PASSIVE]
    constraint_weight: FloatArray  # [N_CONSTRAINT]
    regularisation_names: list[str]
    regularisation: FloatArray  # [N_REG, N_PASSIVE]
    regularisation_weight: FloatArray  # [N_REG]
    sensor_replacements: list[SensorReplacement]
    coil_signals: list[str]
    coil_matrix: FloatArray  # [N_COIL, N_COIL_SIGNAL]
    initial_psi_norm: FloatArray  # [N_GRID]
    initial_mask: npt.NDArray[np.int32]  # [N_GRID]
    xpt_fraction: float
    tolerance: float
    ip_min: float
    n_xpt_max: int
    settings_json: str
    source: str
    psi_to_weber: float


def assemble_preshot(inputs: PreshotInputs) -> Preshot:
    """Derive every preshot quantity from `inputs` and return the validated Preshot."""
    grid_r = np.asarray(inputs.grid_r, dtype=np.float64)
    grid_z = np.asarray(inputs.grid_z, dtype=np.float64)
    n_r, n_z = len(grid_r), len(grid_z)
    n_grid = n_r * n_z
    d_r = float(np.mean(np.diff(grid_r)))
    d_z = float(np.mean(np.diff(grid_z)))
    r_flat, z_flat = (array.ravel() for array in np.meshgrid(grid_r, grid_z))
    edge_r, edge_z = rectangle_points(grid_r, grid_z)

    limiter_index, limiter_weight = limiter_interpolation(grid_r, grid_z, inputs.limiter_r, inputs.limiter_z)
    lower_band, upper_band, permutation = poisson_lu_bands(poisson_matrix(grid_r, grid_z, inputs.poisson_operator), n_r)

    n_constraint = len(inputs.constraint_names)
    n_passive = len(inputs.passive_names)
    n_reg = len(inputs.regularisation_names)
    if n_constraint != inputs.n_flux_loop + inputs.n_bp_probe + inputs.n_rogowski:
        raise PreshotError("constraint names do not match the number of flux loops, Bp probes and Rogowski coils")
    weight = np.concatenate([inputs.constraint_weight, inputs.regularisation_weight]).astype(np.float64)
    n_fit_target = n_constraint + n_reg

    # Rows of fit_target are weighted; the regularisation rows only act on the passive dofs.
    fit_matrix_grid = np.zeros((n_grid, n_fit_target))
    fit_matrix_grid[:, :n_constraint] = (
        (inputs.g_constraint_grid * inputs.constraint_weight[:, np.newaxis]).T * d_r * d_z
    )
    g_fit_target_passive = np.vstack([inputs.g_constraint_passive, inputs.regularisation]).reshape(
        n_fit_target, n_passive
    )
    fit_matrix_passive = (g_fit_target_passive * weight[:, np.newaxis]).T

    replacement, meas_names = replacement_matrix(inputs.constraint_names, inputs.sensor_replacements)
    fit_coefficient_names = [*PLASMA_DOF_NAMES, *inputs.passive_names]
    sizes = {
        "N_R": n_r, "N_Z": n_z, "N_GRID": n_grid, "N_EDGE": len(edge_r), "N_LIMITER": len(inputs.limiter_r),
        "N_COIL": len(inputs.coil_names), "N_MEAS": len(meas_names), "N_FLUX_LOOP": inputs.n_flux_loop,
        "N_BP_PROBE": inputs.n_bp_probe, "N_ROGOWSKI": inputs.n_rogowski, "N_CONSTRAINT": n_constraint,
        "N_REG": n_reg, "N_FIT_TARGET": n_fit_target, "N_PLASMA": N_PLASMA, "N_PASSIVE": n_passive,
        "N_FIT_DOF": N_PLASMA + n_passive, "N_XPT_MAX": inputs.n_xpt_max,
    }  # fmt: skip
    f64 = np.float64
    groups: dict[str, dict[str, Any]] = {
        "SIZES": {name: np.int32(sizes[name]) for name in dd.SIZE_NAMES},
        "GRID": {
            "R": grid_r, "Z": grid_z, "D_R": f64(d_r), "D_Z": f64(d_z), "R_FLAT": r_flat, "Z_FLAT": z_flat,
            "INV_R": 1.0 / r_flat,
            "MASK_LIMITER": np.asarray(inputs.mask_limiter).astype(np.int32),
            "G_COIL": np.asarray(inputs.g_grid_coil, dtype=f64).reshape(n_grid, -1),
            "G_PASSIVE": np.asarray(inputs.g_grid_passive, dtype=f64).reshape(n_grid, n_passive),
        },
        "EDGE": {
            "R": edge_r, "Z": edge_z, "INV_MU0_R": 1.0 / (edge_r * inputs.mu_0),
            "G_EDGE": np.asarray(inputs.g_edge_edge, dtype=f64),
        },
        "LIMITER": {
            "R": np.asarray(inputs.limiter_r, dtype=f64), "Z": np.asarray(inputs.limiter_z, dtype=f64),
            "GRID_INDEX": limiter_index, "GRID_WEIGHT": limiter_weight,
        },
        "POISSON": {
            "RHS_SCALE": np.asarray(inputs.poisson_rhs_scale, dtype=f64), "LOWER_BAND": lower_band,
            "UPPER_BAND": upper_band, "PERMUTATION": permutation,
        },
        "CONSTRAINT": {
            "NAMES": list(inputs.constraint_names),
            "G_COIL": np.asarray(inputs.g_constraint_coil, dtype=f64).reshape(n_constraint, -1),
            "REPLACEMENT": replacement,
        },
        "FIT": {
            "TARGET_NAMES": [*inputs.constraint_names, *inputs.regularisation_names],
            "DOF_NAMES": fit_coefficient_names,
            "WEIGHT": weight, "A_GRID": fit_matrix_grid, "A_PASSIVE": fit_matrix_passive,
        },
        "SETTINGS": {
            "XPT_FRACTION": f64(inputs.xpt_fraction), "TOLERANCE": f64(inputs.tolerance), "IP_MIN": f64(inputs.ip_min),
        },
        "COIL": {
            "NAMES": list(inputs.coil_names), "SIGNALS": list(inputs.coil_signals),
            "MATRIX": np.asarray(inputs.coil_matrix, dtype=f64).reshape(len(inputs.coil_names), -1),
        },
        "MEAS": {"NAMES": meas_names},
        "INITIAL": {
            "PSI_NORM": np.asarray(inputs.initial_psi_norm, dtype=f64).ravel(),
            "MASK": np.asarray(inputs.initial_mask).astype(np.int32).ravel(),
        },
        "INFO": {"SETTINGS": inputs.settings_json, "SOURCE": inputs.source,
                 "MU0": f64(inputs.mu_0), "PSI_TO_WB": f64(inputs.psi_to_weber)},
    }  # fmt: skip
    return Preshot(groups)


def cosine_bump_psi_norm(grid_r: FloatArray, grid_z: FloatArray, r_0: float, z_0: float, radius: float) -> FloatArray:
    """An initial psi_norm: `(1 - cos(pi rho / radius)) / 2` within `radius` of (r_0, z_0), and 1 outside."""
    r_2d, z_2d = np.meshgrid(grid_r, grid_z)
    rho = np.hypot(r_2d - r_0, z_2d - z_0)
    psi_norm = np.where(rho <= radius, 0.5 * (1.0 - np.cos(np.pi * rho / radius)), 1.0)
    return psi_norm.ravel()
