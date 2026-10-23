"""The RT-GSFit data dictionary: one definition of every stored quantity.

Every quantity has one name and one place, and the same name is used everywhere:

* in C, preshot quantities are called `GROUP_NODE` (e.g. `GRID_G_COIL`) and sizes `N_...`;
* in MDSplus and in local files, the path is `PRESHOT.GROUP:NODE` (e.g. `PRESHOT.GRID:G_COIL`);
* in Python, the nested dictionary key is `data["PRESHOT"]["GRID"]["G_COIL"]`;
* the standard_utility `RTGSFIT.csv`, which creates the MDSplus tree, is generated from this module
  by `python -m rtgsfit write-csv`.

Results follow the IMAS-style layout of GSFit (`database_writers/tokamak_energy_mdsplus_new`) so that
analysis scripts written for GSFit also work for RT-GSFit: `GLOBAL:IP`, `GLOBAL.MAG_AXIS:R`,
`PROFILES_2D.R_Z:PSI`, `CONSTRAINTS.FLUX_LOOP.ALL:MEASURED`, ...

Shapes are written with symbolic sizes. `N_TIME` is the number of time slices and is always the first
dimension of a SIGNAL. Matrices are `[rows, columns]`; C stores them row-major. MATLAB reads MDSplus
arrays with the dimensions reversed, so MATLAB sees the transpose.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from dataclasses import field
from typing import Any
from typing import Literal

from rtgsfit.constants import MDSPLUS_NODE_NAME_MAX

NodeType = Literal["NUMERIC", "SIGNAL", "TEXT"]
DataType = Literal["float64", "int32", "str"]

TREE_NAME = "RTGSFIT"
TIME_DIM = "<>:TIME"


@dataclass(frozen=True)
class Node:
    """One stored quantity.

    `path` is the tuple of MDSplus node names, e.g. `("PRESHOT", "GRID", "R")`; the last entry is the
    data node and the others are structures. `c_name` is set for quantities compiled into C.
    """

    path: tuple[str, ...]
    node_type: NodeType
    unit: str
    description: str
    shape: tuple[str, ...] = ()
    dtype: DataType = "float64"
    c_name: str | None = None
    dims: tuple[str, ...] = field(default=())

    @property
    def mds_path(self) -> str:
        """MDSplus path relative to the run node, e.g. `PRESHOT.GRID:R` or `:TIME`."""
        return ".".join(self.path[:-1]) + ":" + self.path[-1]

    @property
    def key(self) -> str:
        """Flat key used in local files, e.g. `PRESHOT.GRID.R`."""
        return ".".join(self.path)


# ----------------------------------------------------------------------------------------------
# Sizes. These become `enum` constants in the generated `preshot_sizes.h` (issue #97).
# ----------------------------------------------------------------------------------------------
SIZE_DESCRIPTIONS: dict[str, str] = {
    "N_R": "Number of grid points in R",
    "N_Z": "Number of grid points in Z",
    "N_GRID": "Number of grid points, N_R * N_Z",
    "N_EDGE": "Number of grid points on the edge of the grid, 2 N_R + 2 N_Z - 4",
    "N_LIMITER": "Number of limiter points",
    "N_COIL": "Number of PF coil circuits",
    "N_MEAS": "Number of PCS measurements, the length of meas",
    "N_FLUX_LOOP": "Number of flux loop constraints",
    "N_BP_PROBE": "Number of Bp probe constraints",
    "N_ROGOWSKI": "Number of Rogowski coil constraints",
    "N_CONSTRAINT": "Number of constraints, N_FLUX_LOOP + N_BP_PROBE + N_ROGOWSKI",
    "N_REG": "Number of regularisation rows of the least-squares problem",
    "N_B_VEC": "Number of rows of the least-squares problem, N_CONSTRAINT + N_REG",
    "N_PLASMA": "Number of plasma degrees of freedom",
    "N_PASSIVE": "Number of passive-structure degrees of freedom",
    "N_X_VEC": "Number of unknowns of the least-squares problem, N_PLASMA + N_PASSIVE",
    "N_XPT_MAX": "Maximum number of x-points that are stored",
}
SIZE_NAMES: tuple[str, ...] = tuple(SIZE_DESCRIPTIONS)
"""Order of the sizes; the same order as `rtgsfit_get_sizes()` in C."""


def _preshot_nodes() -> list[Node]:
    def p(
        group: str,
        node: str,
        node_type: NodeType,
        unit: str,
        description: str,
        shape: tuple[str, ...] = (),
        dtype: DataType = "float64",
        in_c: bool = True,
    ) -> Node:
        return Node(
            ("PRESHOT", group, node),
            node_type,
            unit,
            description,
            shape,
            dtype,
            c_name=f"{group}_{node}" if in_c else None,
        )

    nodes = [
        Node(("PRESHOT", "SIZES", name), "NUMERIC", "dimensionless", description, (), "int32", c_name=name)
        for name, description in SIZE_DESCRIPTIONS.items()
    ]
    nodes += [
        # GRID
        p("GRID", "R", "NUMERIC", "metre", "Major radius of the grid columns", ("N_R",)),
        p("GRID", "Z", "NUMERIC", "metre", "Height of the grid rows", ("N_Z",)),
        p("GRID", "D_R", "NUMERIC", "metre", "Grid spacing in R"),
        p("GRID", "D_Z", "NUMERIC", "metre", "Grid spacing in Z"),
        p("GRID", "R_FLAT", "NUMERIC", "metre", "R of every grid point; index i_z * N_R + i_r", ("N_GRID",)),
        p("GRID", "Z_FLAT", "NUMERIC", "metre", "Z of every grid point; index i_z * N_R + i_r", ("N_GRID",)),
        p("GRID", "INV_MU0_R", "NUMERIC", "1 / henry", "1 / (mu0 R) at every grid point", ("N_GRID",)),
        p(
            "GRID",
            "MASK_LIMITER",
            "NUMERIC",
            "dimensionless",
            "1 where plasma current is allowed (inside the limiter); 0 elsewhere",
            ("N_GRID",),
            "int32",
        ),
        p(
            "GRID",
            "G_COIL",
            "NUMERIC",
            "weber / ampere",
            "Greens function: psi on the grid per unit coil current",
            ("N_GRID", "N_COIL"),
        ),
        p(
            "GRID",
            "G_PASSIVE",
            "NUMERIC",
            "weber / ampere",
            "Greens function: psi on the grid per unit passive degree of freedom",
            ("N_GRID", "N_PASSIVE"),
        ),
        # EDGE
        p(
            "EDGE",
            "R",
            "NUMERIC",
            "metre",
            "R of the grid edge points: from the bottom-left corner up, then clockwise in (R, Z)",
            ("N_EDGE",),
        ),
        p("EDGE", "Z", "NUMERIC", "metre", "Z of the grid edge points", ("N_EDGE",)),
        p("EDGE", "INV_MU0_R", "NUMERIC", "1 / henry", "1 / (mu0 R) at the edge points", ("N_EDGE",)),
        p(
            "EDGE",
            "G_EDGE",
            "NUMERIC",
            "weber / ampere",
            "Greens function between edge points, used for the von Hagenow free-boundary condition",
            ("N_EDGE", "N_EDGE"),
        ),
        # LIMITER
        p("LIMITER", "R", "NUMERIC", "metre", "R of the limiter points", ("N_LIMITER",)),
        p("LIMITER", "Z", "NUMERIC", "metre", "Z of the limiter points", ("N_LIMITER",)),
        p(
            "LIMITER",
            "GRID_INDEX",
            "NUMERIC",
            "dimensionless",
            "Grid indices of the bilinear interpolation stencil of each limiter point",
            ("N_LIMITER", "N_INTERP"),
            "int32",
        ),
        p(
            "LIMITER",
            "GRID_WEIGHT",
            "NUMERIC",
            "dimensionless",
            "Bilinear interpolation weights of each limiter point",
            ("N_LIMITER", "N_INTERP"),
        ),
        # POISSON
        p(
            "POISSON",
            "RHS_SCALE",
            "NUMERIC",
            "henry * metre ** 2",
            "Factor that turns j_phi into the right-hand side of the discrete Grad-Shafranov equation",
            ("N_GRID",),
        ),
        p(
            "POISSON",
            "LOWER_BAND",
            "NUMERIC",
            "dimensionless",
            "Lower band of the LU factors of the discrete Grad-Shafranov operator",
            ("N_GRID", "N_R"),
        ),
        p(
            "POISSON",
            "UPPER_BAND",
            "NUMERIC",
            "dimensionless",
            "Upper band of the LU factors; the diagonal (column 0) is stored inverted",
            ("N_GRID", "N_R_PLUS_2"),
        ),
        p(
            "POISSON",
            "PERMUTATION",
            "NUMERIC",
            "dimensionless",
            "Row permutation of the LU factorisation",
            ("N_GRID",),
            "int32",
        ),
        # CONSTRAINT
        p(
            "CONSTRAINT",
            "NAMES",
            "TEXT",
            "dimensionless",
            "Constraint names: flux loops, then Bp probes, then Rogowski coils",
            ("N_CONSTRAINT",),
            "str",
            in_c=False,
        ),
        p(
            "CONSTRAINT",
            "G_COIL",
            "NUMERIC",
            "dimensionless",
            "Greens function: constraint value per unit coil current",
            ("N_CONSTRAINT", "N_COIL"),
        ),
        p(
            "CONSTRAINT",
            "REPLACEMENT",
            "NUMERIC",
            "dimensionless",
            "Maps meas onto the constraints; replaces faulty sensors by combinations of others",
            ("N_CONSTRAINT", "N_MEAS"),
        ),
        # FIT
        p(
            "FIT",
            "B_VEC_NAMES",
            "TEXT",
            "dimensionless",
            "Names of the rows of b_vec: constraints, then regularisations",
            ("N_B_VEC",),
            "str",
            in_c=False,
        ),
        p(
            "FIT",
            "X_VEC_NAMES",
            "TEXT",
            "dimensionless",
            "Names of the entries of x_vec: plasma, then passive degrees of freedom",
            ("N_X_VEC",),
            "str",
            in_c=False,
        ),
        p("FIT", "WEIGHT", "NUMERIC", "dimensionless", "Weight of each row of b_vec; 1 / uncertainty", ("N_B_VEC",)),
        p(
            "FIT",
            "A_GRID",
            "NUMERIC",
            "dimensionless",
            "Weighted response of b_vec to the current in each grid cell, including d_r d_z",
            ("N_GRID", "N_B_VEC"),
        ),
        p(
            "FIT",
            "A_PASSIVE",
            "NUMERIC",
            "dimensionless",
            "Weighted response of b_vec to each passive degree of freedom; rows of the transposed a_matrix",
            ("N_PASSIVE", "N_B_VEC"),
        ),
        # SETTINGS
        p(
            "SETTINGS",
            "XPT_FRACTION",
            "NUMERIC",
            "dimensionless",
            "boundary_psi = XPT_FRACTION xpt_psi[0] + (1 - XPT_FRACTION) axis_psi when diverted",
        ),
        p(
            "SETTINGS",
            "TOLERANCE",
            "NUMERIC",
            "dimensionless",
            "Tolerance of the null finder and of boundary_psi == axis_psi",
        ),
        p("SETTINGS", "IP_MIN", "NUMERIC", "ampere", "Smallest ip for which the current centroid is computed"),
        # Not used by C: needed to build the inputs to rtgsfit() and to start the first iteration.
        p("COIL", "NAMES", "TEXT", "dimensionless", "PF coil circuit names", ("N_COIL",), "str", in_c=False),
        p(
            "COIL",
            "SIGNALS",
            "TEXT",
            "dimensionless",
            "PCS signals that the coil currents are made from",
            ("N_COIL_SIGNAL",),
            "str",
            in_c=False,
        ),
        p(
            "COIL",
            "MATRIX",
            "NUMERIC",
            "dimensionless",
            "coil_current = MATRIX coil_signal",
            ("N_COIL", "N_COIL_SIGNAL"),
            in_c=False,
        ),
        p("MEAS", "NAMES", "TEXT", "dimensionless", "PCS signal names of meas", ("N_MEAS",), "str", in_c=False),
        p("INITIAL", "PSI_NORM", "NUMERIC", "dimensionless", "psi_norm at the first call", ("N_GRID",), in_c=False),
        p(
            "INITIAL",
            "MASK",
            "NUMERIC",
            "dimensionless",
            "Plasma mask at the first call",
            ("N_GRID",),
            "int32",
            in_c=False,
        ),
        p("INFO", "SETTINGS", "TEXT", "dimensionless", "XML settings that produced the preshot", (), "str", in_c=False),
        p(
            "INFO",
            "SOURCE",
            "TEXT",
            "dimensionless",
            "Database reader and code version that produced the preshot",
            (),
            "str",
            in_c=False,
        ),
    ]
    return nodes


def _results_nodes() -> list[Node]:
    def s(
        path: tuple[str, ...],
        unit: str,
        description: str,
        shape: tuple[str, ...] = (),
        dtype: DataType = "float64",
        dims: tuple[str, ...] = (),
    ) -> Node:
        return Node(path, "SIGNAL", unit, description, ("N_TIME", *shape), dtype, dims=(TIME_DIM, *dims))

    def text(path: tuple[str, ...], description: str, shape: tuple[str, ...] = ()) -> Node:
        return Node(path, "TEXT", "dimensionless", description, shape, "str")

    grid_dims = ("<>.PROFILES_2D.R_Z.Z", "<>.PROFILES_2D.R_Z.R")
    nodes = [
        Node(("TIME",), "NUMERIC", "second", "Time of each call to rtgsfit()", ("N_TIME",)),
        s(("GLOBAL", "IP"), "ampere", "Plasma current"),
        s(("GLOBAL", "PSI_MAG_AXIS"), "weber", "Poloidal flux at the magnetic axis"),
        s(("GLOBAL", "XPT_DIVERTED"), "dimensionless", "Flag for x-point diverted; 1=diverted; 0=limited", (), "int32"),
        s(("GLOBAL", "STATUS"), "dimensionless", "rtgsfit() status bits; 0 is success; see the manual", (), "int32"),
        s(("GLOBAL", "MAG_AXIS", "R"), "metre", "Major radius of the magnetic axis"),
        s(("GLOBAL", "MAG_AXIS", "Z"), "metre", "Height of the magnetic axis"),
        s(("GLOBAL", "CURRENT_CENT", "R"), "metre", "Major radius of the current centroid, integral(j_phi R dS) / Ip"),
        s(("GLOBAL", "CURRENT_CENT", "Z"), "metre", "Height of the current centroid, integral(j_phi Z dS) / Ip"),
        s(("BOUNDARY", "PSI"), "weber", "Poloidal flux at the plasma boundary"),
        s(
            ("BOUNDARY", "OUTLINE", "R"),
            "metre",
            "Major radius of the boundary outline; NaN padded; post-processed",
            ("N_OUTLINE",),
        ),
        s(
            ("BOUNDARY", "OUTLINE", "Z"),
            "metre",
            "Height of the boundary outline; NaN padded; post-processed",
            ("N_OUTLINE",),
        ),
        s(("BOUNDARY", "OUTLINE", "N"), "dimensionless", "Number of points in the boundary outline", (), "int32"),
        s(("BOUNDARY", "GEO_AXIS", "R"), "metre", "(R_max + R_min) / 2 of the boundary outline"),
        s(("BOUNDARY", "GEO_AXIS", "Z"), "metre", "(Z_max + Z_min) / 2 of the boundary outline"),
        s(("BOUNDARY", "MINOR_RADIUS"), "metre", "(R_max - R_min) / 2 of the boundary outline"),
        s(("BOUNDARY", "ELONGATION"), "dimensionless", "(Z_max - Z_min) / (2 MINOR_RADIUS)"),
        s(("BOUNDARY", "TRIANG"), "dimensionless", "Mean of the upper and lower triangularity"),
        s(("BOUNDARY", "TRIANG_U"), "dimensionless", "(GEO_AXIS:R - R(Z_max)) / MINOR_RADIUS"),
        s(("BOUNDARY", "TRIANG_L"), "dimensionless", "(GEO_AXIS:R - R(Z_min)) / MINOR_RADIUS"),
        s(("CONSTRAINTS", "CHI_SQ_MAG"), "dimensionless", "Weighted squared residual of the magnetic constraints"),
        s(("CONSTRAINTS", "CHI_SQ_TOTAL"), "dimensionless", "Weighted squared residual including regularisation"),
    ]
    for group, unit, size in (
        ("FLUX_LOOP", "weber", "N_FLUX_LOOP"),
        ("BP_PROBE", "tesla", "N_BP_PROBE"),
        ("ROGOWSKI", "ampere", "N_ROGOWSKI"),
    ):
        base = ("CONSTRAINTS", group, "ALL")
        nodes += [
            text((*base, "NAMES"), "Sensor names", (size,)),
            s((*base, "MEASURED"), unit, "Measured value, after sensor replacement", (size,)),
            s((*base, "RECONSTRUCT"), unit, "Value calculated from the reconstructed equilibrium", (size,)),
            s((*base, "WEIGHT"), "dimensionless", "Weight given to the measurement; 1 / uncertainty", (size,)),
        ]
    nodes += [
        text(("CONSTRAINTS", "PF_CURRENT", "ALL", "NAMES"), "PF coil circuit names", ("N_COIL",)),
        s(("CONSTRAINTS", "PF_CURRENT", "ALL", "MEASURED"), "ampere", "PF coil circuit currents", ("N_COIL",)),
        Node(("PROFILES_2D", "R_Z", "R"), "NUMERIC", "metre", "Major radius of the grid columns", ("N_R",)),
        Node(("PROFILES_2D", "R_Z", "Z"), "NUMERIC", "metre", "Height of the grid rows", ("N_Z",)),
        s(("PROFILES_2D", "R_Z", "PSI"), "weber", "Poloidal flux", ("N_Z", "N_R"), dims=grid_dims),
        s(
            ("PROFILES_2D", "R_Z", "PSI_NORM"),
            "dimensionless",
            "Normalised poloidal flux; 0 on axis, 1 on the boundary and outside the plasma",
            ("N_Z", "N_R"),
            dims=grid_dims,
        ),
        s(
            ("PROFILES_2D", "R_Z", "MASK"),
            "dimensionless",
            "1 inside the plasma; 0 outside",
            ("N_Z", "N_R"),
            "int32",
            dims=grid_dims,
        ),
        s(
            ("PROFILES_2D", "R_Z", "J_PHI"),
            "ampere / metre ** 2",
            "Toroidal plasma current density",
            ("N_Z", "N_R"),
            dims=grid_dims,
        ),
        s(("XPOINTS", "ALL", "N"), "dimensionless", "Number of x-points found", (), "int32"),
        s(("XPOINTS", "ALL", "R"), "metre", "Major radius of the x-points, sorted by psi; NaN padded", ("N_XPT_MAX",)),
        s(("XPOINTS", "ALL", "Z"), "metre", "Height of the x-points; NaN padded", ("N_XPT_MAX",)),
        s(("XPOINTS", "ALL", "PSI"), "weber", "Poloidal flux at the x-points; NaN padded", ("N_XPT_MAX",)),
        s(
            ("FIT", "X_VEC"),
            "dimensionless",
            "Solution of the least-squares problem a_matrix x_vec = b_vec",
            ("N_X_VEC",),
        ),
        text(("FIT", "X_VEC_NAMES"), "Names of the entries of x_vec", ("N_X_VEC",)),
        s(("FIT", "LAPACK_INFO"), "dimensionless", "INFO returned by LAPACK dgelss; 0 is success", (), "int32"),
        text(("INPUT", "SETTINGS"), "XML settings of the replay"),
        text(("INPUT", "PRESHOT_SRC"), "Where the preshot was read from"),
        text(("INPUT", "CODE_VERSION"), "RT-GSFit version"),
        s(("INPUT", "MEAS"), "dimensionless", "meas: the PCS measurements passed to rtgsfit()", ("N_MEAS",)),
        s(
            ("INPUT", "COIL_SIGNAL"),
            "ampere",
            "PCS coil signals; coil_current = PRESHOT.COIL:MATRIX coil_signal",
            ("N_COIL_SIGNAL",),
        ),
        text(("TIMING", "SECTIONS"), "Names of the timed sections of rtgsfit()", ("N_TIMING",)),
        s(
            ("TIMING", "DURATION"),
            "second",
            "Thread CPU time of each section; zero unless built with TIMING=1",
            ("N_TIMING",),
        ),
    ]
    return nodes


STRUCTURE_DESCRIPTIONS: dict[tuple[str, ...], str] = {
    ("PRESHOT",): "Everything computed before the pulse; compiled into C for real-time use",
    ("PRESHOT", "SIZES"): "Array sizes; enum constants in C",
    ("PRESHOT", "GRID"): "The computational grid",
    ("PRESHOT", "EDGE"): "The edge of the grid, for the free-boundary condition",
    ("PRESHOT", "LIMITER"): "The limiter",
    ("PRESHOT", "POISSON"): "The discrete Grad-Shafranov operator",
    ("PRESHOT", "CONSTRAINT"): "From the PCS measurements to the constraints",
    ("PRESHOT", "FIT"): "The weighted least-squares problem a_matrix x_vec = b_vec",
    ("PRESHOT", "SETTINGS"): "Scalar algorithm settings",
    ("PRESHOT", "COIL"): "PF coil circuits",
    ("PRESHOT", "MEAS"): "PCS measurements",
    ("PRESHOT", "INITIAL"): "State before the first call",
    ("PRESHOT", "INFO"): "Provenance",
    ("GLOBAL",): "0D parameters of the equilibrium",
    ("GLOBAL", "MAG_AXIS"): "Magnetic axis",
    ("GLOBAL", "CURRENT_CENT"): "Current centroid",
    ("BOUNDARY",): "Plasma boundary",
    ("BOUNDARY", "OUTLINE"): "Contour psi = boundary_psi around the magnetic axis",
    ("BOUNDARY", "GEO_AXIS"): "Geometric axis",
    ("CONSTRAINTS",): "Measurements used to constrain the equilibrium and their reconstructed values",
    ("CONSTRAINTS", "FLUX_LOOP"): "Flux loops",
    ("CONSTRAINTS", "BP_PROBE"): "Poloidal field probes",
    ("CONSTRAINTS", "ROGOWSKI"): "Rogowski coils",
    ("CONSTRAINTS", "PF_CURRENT"): "PF coil circuit currents",
    ("PROFILES_2D",): "2D profiles",
    ("PROFILES_2D", "R_Z"): "On the rectangular (R, Z) grid",
    ("XPOINTS",): "X-points",
    ("XPOINTS", "ALL"): "All x-points that were found",
    ("FIT",): "The least-squares fit",
    ("INPUT",): "Inputs and settings",
    ("TIMING",): "Execution time of rtgsfit()",
}

PRESHOT_NODES: tuple[Node, ...] = tuple(_preshot_nodes())
RESULTS_NODES: tuple[Node, ...] = tuple(_results_nodes())
ALL_NODES: tuple[Node, ...] = RESULTS_NODES + PRESHOT_NODES
NODES_BY_KEY: dict[str, Node] = {node.key: node for node in ALL_NODES}
C_NODES: tuple[Node, ...] = tuple(node for node in PRESHOT_NODES if node.c_name is not None)


def preshot_node(group: str, name: str) -> Node:
    """Return the node `PRESHOT.group:name`."""
    return NODES_BY_KEY[f"PRESHOT.{group}.{name}"]


def check_names() -> None:
    """Check the MDSplus constraints on the node names; raises ValueError."""
    if len(NODES_BY_KEY) != len(ALL_NODES):
        raise ValueError("duplicate node in the data dictionary")
    for node in ALL_NODES:
        for name in node.path:
            if len(name) > MDSPLUS_NODE_NAME_MAX or not name.replace("_", "").isalnum() or not name.isupper():
                raise ValueError(f"invalid MDSplus node name {name!r} in {node.mds_path}")
        if len(node.path) > 6:
            raise ValueError(f"{node.mds_path} is too deep for the CSV format")


# ----------------------------------------------------------------------------------------------
# Nested dictionaries
# ----------------------------------------------------------------------------------------------
def nested_set(data: dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    for key in path[:-1]:
        data = data.setdefault(key, {})
    data[path[-1]] = value


def nested_get(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    for key in path:
        data = data[key]
    return data


def flatten(data: dict[str, Any], prefix: tuple[str, ...] = ()) -> dict[str, Any]:
    """`{"A": {"B": 1}}` -> `{"A.B": 1}`."""
    flat: dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, dict):
            flat.update(flatten(value, (*prefix, key)))
        else:
            flat[".".join((*prefix, key))] = value
    return flat


def unflatten(flat: dict[str, Any]) -> dict[str, Any]:
    data: dict[str, Any] = {}
    for key, value in flat.items():
        nested_set(data, tuple(key.split(".")), value)
    return data


# ----------------------------------------------------------------------------------------------
# documentation/data_dictionary.md
# ----------------------------------------------------------------------------------------------
def generate_markdown() -> str:
    """Return the text of `documentation/data_dictionary.md` (`python -m rtgsfit write-data-dictionary`)."""
    check_names()
    lines = [
        "# RT-GSFit data dictionary",
        "",
        "<!-- Generated by `python -m rtgsfit write-data-dictionary`; do not edit by hand. -->",
        "",
        "Every quantity that RT-GSFit stores, in MDSplus (tree `RTGSFIT`, below the run node) and in",
        "local files. The C name is the name of the compiled-in preshot constant. `N_TIME` is the number",
        "of time slices. See [manual.md](manual.md) for the definitions of the symbols.",
        "",
        "## Sizes",
        "",
        "| Size | Description |",
        "|---|---|",
        *[f"| `{name}` | {description} |" for name, description in SIZE_DESCRIPTIONS.items()],
        "| `N_COIL_SIGNAL` | Number of PCS coil signals (Python only) |",
        "| `N_INTERP` | 4, the points of a bilinear interpolation stencil |",
        "| `N_R_PLUS_2` | N_R + 2 |",
        "| `N_OUTLINE` | Number of points stored for the boundary outline (Python only) |",
        "| `N_TIMING` | Number of timed sections of `rtgsfit()` |",
    ]
    sections = [
        ("Preshot", [node for node in PRESHOT_NODES if node.path[1] != "SIZES"]),
        ("Results", list(RESULTS_NODES)),
    ]
    for title, nodes in sections:
        lines += ["", f"## {title}"]
        group: tuple[str, ...] | None = None
        for node in nodes:
            parents = node.path[:-1]
            if parents != group:
                group = parents
                heading = ".".join(parents) if parents else "(run node)"
                description = STRUCTURE_DESCRIPTIONS.get(parents, "")
                lines += ["", f"### `{heading}`", ""]
                if description:
                    lines += [description + ".", ""]
                lines += ["| Node | C name | Shape | Type | Unit | Description |", "|---|---|---|---|---|---|"]
            shape = "[" + ", ".join(node.shape) + "]" if node.shape else "scalar"
            c_name = f"`{node.c_name}`" if node.c_name else ""
            lines.append(
                f"| `{node.mds_path}` | {c_name} | {shape} | {node.dtype} | {node.unit} | {node.description} |"
            )
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------------
# The standard_utility CSV that defines the MDSplus tree
# ----------------------------------------------------------------------------------------------
N_CSV_NODE_COLUMNS = 5


def generate_csv(run_info: str = "Default run") -> str:
    """Return the text of `standard_utility/csv_files/RTGSFIT.csv`, in the format of `GSFIT.csv`."""
    check_names()
    width = N_CSV_NODE_COLUMNS + 2 + 3 * 4
    rows: list[list[str]] = []

    def row(*cells: str) -> list[str]:
        return [*cells, *[""] * (width - len(cells))]

    rows += [
        row("Root node"),
        row("Root info"),
        row("Run name", "RUN01"),
        row("Script info"),
        row("Script name"),
        row("Tree", TREE_NAME),
        row("Run info", run_info),
    ]
    rows.append(row(*[""] * 6, "Data", "", "", "Time", "", "", "Dimension 1", "", "", "Dimension 2", "", ""))
    rows.append(["Node"] * N_CSV_NODE_COLUMNS + ["Data node", "Node type"] + ["Data", "Unit", "Description"] * 4)
    rows.append(row())

    written_structures: set[tuple[str, ...]] = set()
    previous_top: str | None = None
    for node in ALL_NODES:
        top = node.path[0] if len(node.path) > 1 else ""
        if previous_top is not None and top != previous_top:
            rows.append(row())
        previous_top = top
        for depth in range(1, len(node.path)):
            structure = node.path[:depth]
            if structure not in written_structures:
                written_structures.add(structure)
                cells = [*structure, *[""] * (N_CSV_NODE_COLUMNS - len(structure))]
                rows.append(row(*cells, "", "STRUCTURE", "", "", STRUCTURE_DESCRIPTIONS.get(structure, "")))
        parents = node.path[:-1]
        cells = [*parents, *[""] * (N_CSV_NODE_COLUMNS - len(parents)), node.path[-1], node.node_type]
        data = "data" + "".join(f"['{name}']" for name in node.path)
        cells += [data, node.unit, node.description]
        for dim in node.dims:
            cells += [dim, "", ""]
        rows.append(row(*cells))

    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\n").writerows(rows)
    return buffer.getvalue()
