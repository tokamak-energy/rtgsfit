"""Collect the outputs of a replay into the results layout of the data dictionary.

`Results.nested` is `{"TIME": ..., "GLOBAL": {...}, ..., "PRESHOT": {...}}`, exactly the layout of the
RTGSFIT MDSplus tree (see `data_dictionary.RESULTS_NODES`), so every writer stores the same thing.
`Results.flat` is the same keyed by `Node.key`, e.g. `"GLOBAL.IP"`.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from typing import Any

import numpy as np

import rtgsfit
from rtgsfit import data_dictionary as dd
from rtgsfit.boundary import BoundaryShape
from rtgsfit.boundary import boundary_shape
from rtgsfit.boundary import outline_size
from rtgsfit.database_readers import Signals
from rtgsfit.library import TIMING_SECTIONS
from rtgsfit.library import SliceResult
from rtgsfit.preshot import Preshot


@dataclass
class Results:
    nested: dict[str, Any] = field(default_factory=dict)

    @property
    def flat(self) -> dict[str, Any]:
        return {".".join(path): value for path, value in _walk(self.nested)}

    def __getitem__(self, key: str) -> Any:
        """`results["GLOBAL.IP"]` or `results["PROFILES_2D.R_Z.PSI"]`."""
        return dd.nested_get(self.nested, tuple(key.split(".")))

    @property
    def n_time(self) -> int:
        return len(self.nested["TIME"])


def _walk(data: dict[str, Any], prefix: tuple[str, ...] = ()) -> list[tuple[tuple[str, ...], Any]]:
    items = []
    for key, value in data.items():
        if isinstance(value, dict):
            items += _walk(value, (*prefix, key))
        else:
            items.append(((*prefix, key), value))
    return items


def collect_results(
    preshot: Preshot,
    signals: Signals,
    slices: list[SliceResult],
    settings_xml: str,
    preshot_source: str,
) -> Results:
    """Map the per-slice outputs of `rtgsfit()` onto the results nodes.

    `slices[i]` is the result at `signals.time[i]`; a replay that stopped early has fewer slices.
    """
    sizes = preshot.sizes
    n_time = len(slices)
    n_r, n_z = sizes["N_R"], sizes["N_Z"]
    n_loop, n_probe = sizes["N_FLUX_LOOP"], sizes["N_BP_PROBE"]
    weight = preshot["FIT:WEIGHT"][: sizes["N_CONSTRAINT"]]
    constraint_names = preshot["CONSTRAINT:NAMES"]

    if n_time == 0:
        raise ValueError("a replay must have at least one time slice")

    def stack(name: str, dtype: type = np.float64) -> np.ndarray:
        return np.array([getattr(item, name) for item in slices], dtype=dtype)

    def grid(name: str, dtype: type = np.float64) -> np.ndarray:
        return stack(name, dtype).reshape(n_time, n_z, n_r)

    xpt_n = stack("xpt_n", np.int32)
    xpt_valid = np.arange(sizes["N_XPT_MAX"])[np.newaxis, :] < xpt_n[:, np.newaxis]

    def xpt(name: str) -> np.ndarray:
        return np.where(xpt_valid, stack(name), np.nan)

    measured, reconstructed = stack("constraint_measured"), stack("constraint_reconstructed")
    groups = {"FLUX_LOOP": slice(0, n_loop), "BP_PROBE": slice(n_loop, n_loop + n_probe),
              "ROGOWSKI": slice(n_loop + n_probe, sizes["N_CONSTRAINT"])}  # fmt: skip
    constraints: dict[str, Any] = {
        "CHI_SQ_MAG": stack("chi_sq_magnetic"),
        "CHI_SQ_TOTAL": stack("chi_sq_total"),
    }
    for group, rows in groups.items():
        constraints[group] = {
            "ALL": {
                "NAMES": list(constraint_names[rows]),
                "MEASURED": measured[:, rows],
                "RECONSTRUCT": reconstructed[:, rows],
                "WEIGHT": np.tile(weight[rows], (n_time, 1)),
            }
        }
    coil_signal = signals.coil_signal[:n_time]
    constraints["PF_CURRENT"] = {
        "ALL": {"NAMES": list(preshot["COIL:NAMES"]), "MEASURED": coil_signal @ preshot["COIL:MATRIX"].T}
    }

    r, z = np.asarray(preshot["GRID:R"]), np.asarray(preshot["GRID:Z"])
    shapes = [
        boundary_shape(r, z, item.psi.reshape(n_z, n_r), item.mask.reshape(n_z, n_r), item.boundary_psi,
                       item.axis_psi, item.axis_r, item.axis_z)
        if item.status == 0 else BoundaryShape.empty(outline_size(n_r, n_z))
        for item in slices
    ]  # fmt: skip

    def shape(name: str, dtype: type = np.float64) -> np.ndarray:
        return np.array([getattr(item, name) for item in shapes], dtype=dtype)

    timing = np.zeros((n_time, len(TIMING_SECTIONS)))
    for i, item in enumerate(slices):
        if item.timing is not None:
            timing[i] = item.timing

    nested: dict[str, Any] = {
        "TIME": np.asarray(signals.time[:n_time], dtype=np.float64),
        "GLOBAL": {
            "IP": stack("ip"),
            "PSI_MAG_AXIS": stack("axis_psi"),
            "XPT_DIVERTED": stack("xpt_diverted", np.int32),
            "STATUS": stack("status", np.int32),
            "MAG_AXIS": {"R": stack("axis_r"), "Z": stack("axis_z")},
            "CURRENT_CENT": {"R": stack("current_centroid_r"), "Z": stack("current_centroid_z")},
        },
        "BOUNDARY": {
            "PSI": stack("boundary_psi"),
            "OUTLINE": {"R": shape("outline_r"), "Z": shape("outline_z"), "N": shape("outline_n", np.int32)},
            "GEO_AXIS": {"R": shape("r_geo"), "Z": shape("z_geo")},
            "MINOR_RADIUS": shape("minor_radius"),
            "ELONGATION": shape("elongation"),
            "TRIANG": shape("triangularity"),
            "TRIANG_U": shape("triangularity_upper"),
            "TRIANG_L": shape("triangularity_lower"),
        },
        "CONSTRAINTS": constraints,
        "PROFILES_2D": {
            "R_Z": {
                "R": r,
                "Z": z,
                "PSI": grid("psi"),
                "PSI_NORM": grid("psi_norm"),
                "MASK": grid("mask", np.int32),
                "J_PHI": grid("j_phi"),
            }
        },
        "XPOINTS": {"ALL": {"N": xpt_n, "R": xpt("xpt_r"), "Z": xpt("xpt_z"), "PSI": xpt("xpt_psi")}},
        "FIT": {
            "X_VEC": stack("x_vec"),
            "X_VEC_NAMES": list(preshot["FIT:X_VEC_NAMES"]),
            "LAPACK_INFO": stack("lapack_info", np.int32),
        },
        "INPUT": {
            "SETTINGS": settings_xml,
            "PRESHOT_SRC": preshot_source,
            "CODE_VERSION": rtgsfit.__version__,
            "MEAS": np.asarray(signals.meas[:n_time], dtype=np.float64),
            "COIL_SIGNAL": np.asarray(coil_signal, dtype=np.float64),
        },
        "TIMING": {"SECTIONS": list(TIMING_SECTIONS), "DURATION": timing},
    }
    nested.update(preshot.to_nested())
    results = Results(nested)
    validate_results(results, preshot)
    return results


def validate_results(results: Results, preshot: Preshot) -> None:
    """Check that the results have exactly the nodes, shapes and dtypes of the data dictionary."""
    sizes = preshot.symbolic_sizes()
    sizes["N_TIME"] = results.n_time
    sizes["N_TIMING"] = len(TIMING_SECTIONS)
    sizes["N_OUTLINE"] = outline_size(sizes["N_R"], sizes["N_Z"])
    flat = results.flat
    expected = {node.key for node in dd.ALL_NODES}
    if set(flat) != expected:
        raise ValueError(f"results nodes differ from the data dictionary: {sorted(set(flat) ^ expected)}")
    for node in dd.RESULTS_NODES:
        value = flat[node.key]
        if node.dtype == "str":
            continue
        array = np.asarray(value)
        shape = tuple(sizes[dim] for dim in node.shape)
        if array.shape != shape or array.dtype != np.dtype(node.dtype):
            raise ValueError(f"{node.mds_path}: {array.shape} {array.dtype}, expected {shape} {node.dtype}")
