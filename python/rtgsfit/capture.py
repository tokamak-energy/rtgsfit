"""Portable native captures shared by replay and the PCS logger; no solver calls on readback."""

from __future__ import annotations

import ctypes
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from netCDF4 import Dataset

import rtgsfit
from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers import Signals
from rtgsfit.database_writers.c_source import preshot_checksum as solver_preshot_checksum
from rtgsfit.library import TIMING_SECTIONS
from rtgsfit.library import SliceResult
from rtgsfit.library import _load_library
from rtgsfit.library import build_library
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results
from rtgsfit.results import collect_results
from rtgsfit.settings import parse_settings

SCHEMA = "rtgsfit-capture-2"
Spec = tuple[str, tuple[str, ...]]
OUTPUTS: dict[str, Spec] = {
    "status": ("int32", ()),
    "psi": ("float64", ("N_Z", "N_R")),
    "psi_norm": ("float64", ("N_Z", "N_R")),
    "mask": ("int32", ("N_Z", "N_R")),
    "j_phi": ("float64", ("N_Z", "N_R")),
    "fit_dof": ("float64", ("N_FIT_DOF",)),
    "constraint_measured": ("float64", ("N_CONSTRAINT",)),
    "constraint_reconstructed": ("float64", ("N_CONSTRAINT",)),
    "xpt_r": ("float64", ("N_XPT_MAX",)),
    "xpt_z": ("float64", ("N_XPT_MAX",)),
    "xpt_psi": ("float64", ("N_XPT_MAX",)),
    "xpt_n": ("int32", ()),
    "xpt_diverted": ("int32", ()),
    "ip": ("float64", ()),
    "axis_psi": ("float64", ()),
    "axis_r": ("float64", ()),
    "axis_z": ("float64", ()),
    "boundary_psi": ("float64", ()),
    "current_centroid_r": ("float64", ()),
    "current_centroid_z": ("float64", ()),
    "chi_sq_magnetic": ("float64", ()),
    "chi_sq_total": ("float64", ()),
    "lapack_info": ("int32", ()),
}
INPUTS: dict[str, Spec] = {
    "time": ("float64", ("N_TIME",)),
    "meas": ("float64", ("N_TIME", "N_MEAS")),
    "coil_signal": ("float64", ("N_TIME", "N_COIL_SIGNAL")),
    "coil_current": ("float64", ("N_TIME", "N_COIL")),
}


class _CDimension(ctypes.Structure):
    _fields_ = [("name", ctypes.c_char_p), ("size", ctypes.c_size_t)]


class _CVariable(ctypes.Structure):
    _fields_ = [
        ("group", ctypes.c_char_p),
        ("name", ctypes.c_char_p),
        ("type", ctypes.c_int),
        ("n_dimensions", ctypes.c_size_t),
        ("dimensions", ctypes.POINTER(ctypes.c_char_p)),
        ("values", ctypes.c_void_p),
        ("n_values", ctypes.c_size_t),
    ]


def _full_preshot_checksum(preshot: Preshot) -> int:
    checksum = 0
    for node in dd.PRESHOT_NODES:
        _, group, name = node.path
        checksum = zlib.crc32(node.key.encode() + b"\0", checksum)
        checksum = zlib.crc32(node.dtype.encode() + b"\0", checksum)
        value = preshot[f"{group}:{name}"]
        if node.dtype == "str":
            for item in np.asarray(value, dtype=object).reshape(-1):
                encoded = str(item).encode("utf-8")
                checksum = zlib.crc32(len(encoded).to_bytes(8, "little"), checksum)
                checksum = zlib.crc32(encoded, checksum)
        else:
            dtype = np.dtype(node.dtype).newbyteorder("<")
            data = np.ascontiguousarray(value, dtype=dtype).tobytes()
            checksum = zlib.crc32(data, checksum)
    return checksum


@dataclass
class Capture:
    """Inputs and outputs in native C units, with the exact preshot used for those calls."""

    preshot: Preshot
    signals: Signals
    slices: list[SliceResult]
    settings_json: str
    preshot_source: str
    code_version: str = rtgsfit.__version__

    def validate(self) -> None:
        self.preshot.validate()
        self.signals.validate(self.preshot)
        parse_settings(self.settings_json, name="captured settings")
        if not self.slices or len(self.slices) != len(self.signals.time):
            raise ValueError("capture needs one output slice per input time, and at least one slice")
        if not np.isfinite(self.signals.time).all():
            raise ValueError("capture times must be finite")
        for name in ("time", "meas", "coil_signal"):
            if np.asarray(getattr(self.signals, name)).dtype.kind != "f":
                raise ValueError(f"capture input {name} must contain floating-point values")
        if not self.preshot_source or not self.code_version:
            raise ValueError("capture needs a preshot source and code version")
        sizes = self.preshot.symbolic_sizes()
        timed = self.slices[0].timing is not None
        for i, item in enumerate(self.slices):
            for name, (dtype, dims) in OUTPUTS.items():
                value = np.asarray(getattr(item, name))
                shape = tuple(sizes[dim] for dim in dims)
                if len(shape) == 2:
                    shape = (int(np.prod(shape)),)
                if value.shape != shape:
                    raise ValueError(f"capture slice {i} {name}: shape {value.shape}, expected {shape}")
                if dtype == "int32":
                    if value.dtype.kind not in "iu" or np.any(value < -(2**31)) or np.any(value >= 2**31):
                        raise ValueError(f"capture slice {i} {name} must contain int32 values")
                elif value.dtype.kind != "f":
                    raise ValueError(f"capture slice {i} {name} must contain floating-point values")
            if (item.timing is not None) != timed:
                raise ValueError("timing must be present for every capture slice or absent for every slice")
            if item.timing is not None and item.timing.shape != (len(TIMING_SECTIONS),):
                raise ValueError(f"capture slice {i} timing has the wrong shape")

    def to_results(self) -> Results:
        """Use the same offline processing as replay, without rebuilding or rerunning C."""
        self.validate()
        return collect_results(
            self.preshot,
            self.signals,
            self.slices,
            self.settings_json,
            self.preshot_source,
            code_version=self.code_version,
        )


def _read_variable(group: Any, name: str, spec: Spec, sizes: dict[str, int]) -> Any:
    dtype, dims = spec
    if name not in group.variables:
        raise ValueError(f"capture is missing {group.path}/{name}")
    if not set(dims) <= sizes.keys():
        raise ValueError(f"capture is missing dimensions for {group.path}/{name}: {dims}")
    variable = group.variables[name]
    shape = tuple(sizes[dim] for dim in dims)
    dtype_valid = (
        variable.dtype is str if dtype == "str" else np.dtype(variable.dtype).newbyteorder("=") == np.dtype(dtype)
    )
    if variable.dimensions != dims or variable.shape != shape or not dtype_valid:
        raise ValueError(f"capture {group.path}/{name} must have dimensions {dims}, shape {shape}, dtype {dtype}")
    value = variable[...]
    if dtype == "str":
        return str(value) if not dims else np.asarray(value).tolist()
    return np.asarray(value, dtype=dtype)


def write_capture(path: str | Path, capture: Capture, library_path: str | Path | None = None) -> Path:
    """Write a NetCDF-4 capture with the shared C writer used by replay and PCS."""
    capture.validate()
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    sizes = capture.preshot.symbolic_sizes()
    sizes.update(N_TIME=len(capture.slices), N_TIMING=len(TIMING_SECTIONS))
    timed = capture.slices[0].timing is not None
    if library_path is None:
        library_path = build_library(capture.preshot)
    library = _load_library(Path(library_path))
    library.cdll.rtgsfit_get_preshot_checksum.restype = ctypes.c_uint32
    library.cdll.rtgsfit_get_preshot_checksum.argtypes = []
    if library.cdll.rtgsfit_get_preshot_checksum() != solver_preshot_checksum(capture.preshot):
        raise ValueError("capture preshot does not match the compiled solver")
    function = library.cdll.rtgsfit_capture_write
    function.restype = ctypes.c_int32
    function.argtypes = [
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_int32,
        ctypes.POINTER(_CDimension),
        ctypes.c_size_t,
        ctypes.POINTER(_CVariable),
        ctypes.c_size_t,
    ]
    library.cdll.rtgsfit_capture_last_error.restype = ctypes.c_char_p

    keepalive: list[Any] = []
    dimension_names = [name.encode() for name in sizes]
    keepalive.extend(dimension_names)
    dimensions = (_CDimension * len(sizes))(
        *(_CDimension(name, size) for name, size in zip(dimension_names, sizes.values(), strict=True))
    )

    variables: list[tuple[str, str, str, tuple[str, ...], Any]] = []
    for node in dd.PRESHOT_NODES:
        if node.c_name is not None:
            continue
        _, group, name = node.path
        variables.append((f"PRESHOT/{group}", name, node.dtype, node.shape, capture.preshot[f"{group}:{name}"]))
    for name, (_, dims) in INPUTS.items():
        value = (
            capture.signals.coil_current(capture.preshot) if name == "coil_current" else getattr(capture.signals, name)
        )
        variables.append(("INPUT", name, "float64", dims, value))
    for name, (dtype, dims) in OUTPUTS.items():
        shape = (sizes["N_TIME"], *(sizes[dim] for dim in dims))
        value = np.asarray([getattr(item, name) for item in capture.slices], dtype=dtype).reshape(shape)
        variables.append(("OUTPUT", name, dtype, ("N_TIME", *dims), value))
    if timed:
        variables.extend(
            [
                (
                    "OUTPUT",
                    "timing",
                    "float64",
                    ("N_TIME", "N_TIMING"),
                    np.asarray([item.timing for item in capture.slices], dtype=np.float64),
                ),
                ("OUTPUT", "timing_sections", "str", ("N_TIMING",), TIMING_SECTIONS),
            ]
        )

    c_variables = (_CVariable * len(variables))()
    for i, (group, name, dtype, dims, value) in enumerate(variables):
        group_bytes, name_bytes = group.encode(), name.encode()
        dimension_bytes = tuple(dim.encode() for dim in dims)
        keepalive.extend((group_bytes, name_bytes, dimension_bytes))
        c_dimensions = (ctypes.c_char_p * len(dimension_bytes))(*dimension_bytes)
        keepalive.append(c_dimensions)
        if dtype == "str":
            strings = [str(item).encode("utf-8") for item in np.asarray(value, dtype=object).reshape(-1)]
            c_values = (ctypes.c_char_p * len(strings))(*strings)
            keepalive.extend((strings, c_values))
            c_type = 3
            data_pointer = ctypes.cast(c_values, ctypes.c_void_p)
            n_values = len(strings)
        else:
            array = np.ascontiguousarray(value, dtype=np.float64 if dtype == "float64" else np.int32)
            keepalive.append(array)
            c_type = 1 if dtype == "float64" else 2
            data_pointer = ctypes.c_void_p(array.ctypes.data)
            n_values = array.size
        c_variables[i] = _CVariable(group_bytes, name_bytes, c_type, len(dims), c_dimensions, data_pointer, n_values)
    keepalive.append(c_variables)
    encoded_path = str(path).encode()
    metadata = (
        capture.settings_json.encode(),
        capture.preshot_source.encode(),
        capture.code_version.encode(),
    )
    if (
        function(
            encoded_path,
            *metadata,
            int(timed),
            dimensions,
            len(sizes),
            c_variables,
            len(variables),
        )
        != 0
    ):
        error = library.cdll.rtgsfit_capture_last_error()
        raise OSError(error.decode("utf-8", errors="replace") if error else "C capture writer failed")
    return path


def read_capture(path: str | Path) -> Capture:
    """Read only the versioned v2 capture contract; legacy PCS files need an explicit converter."""
    with Dataset(Path(path).expanduser(), "r") as dataset:
        dataset.set_auto_mask(False)
        if getattr(dataset, "schema", None) != SCHEMA:
            raise ValueError(f"unsupported capture schema; expected {SCHEMA!r} (legacy PCS files are not v2 captures)")
        required = {
            "flux_convention",
            "settings_json",
            "preshot_source",
            "code_version",
            "solver_checksum",
            "preshot_checksum",
            "timing_enabled",
        }
        if not required <= set(dataset.ncattrs()):
            raise ValueError(f"capture is missing attributes: {sorted(required - set(dataset.ncattrs()))}")
        if dataset.flux_convention != "native" or dataset.timing_enabled not in (0, 1):
            raise ValueError("capture needs native flux units and timing_enabled of 0 or 1")
        if set(dataset.groups) != {"PRESHOT", "INPUT", "OUTPUT"} or dataset.variables:
            raise ValueError("capture must contain PRESHOT, INPUT and OUTPUT groups")
        sizes = {name: len(dim) for name, dim in dataset.dimensions.items()}
        if sizes.get("N_TIME", 0) < 1 or sizes.get("N_TIMING") != len(TIMING_SECTIONS):
            raise ValueError("capture has invalid time or timing dimensions")
        expected_groups: dict[str, set[str]] = {}
        for node in dd.PRESHOT_NODES:
            expected_groups.setdefault(node.path[1], set()).add(node.path[2])
        if set(dataset["PRESHOT"].groups) != set(expected_groups) or dataset["PRESHOT"].variables:
            raise ValueError("capture PRESHOT groups do not match the data dictionary")
        for group, names in expected_groups.items():
            if set(dataset[f"PRESHOT/{group}"].variables) != names or dataset[f"PRESHOT/{group}"].groups:
                raise ValueError(f"capture PRESHOT/{group} nodes do not match the data dictionary")
        nested: dict[str, Any] = {}
        for node in dd.PRESHOT_NODES:
            _, group, name = node.path
            if group not in dataset["PRESHOT"].groups:
                raise ValueError(f"capture is missing PRESHOT/{group}")
            value = _read_variable(dataset[f"PRESHOT/{group}"], name, (node.dtype, node.shape), sizes)
            dd.nested_set(nested, node.path, value)
        preshot = Preshot.from_nested(nested)
        expected_sizes = preshot.symbolic_sizes() | {"N_TIME": sizes["N_TIME"], "N_TIMING": len(TIMING_SECTIONS)}
        if sizes != expected_sizes:
            raise ValueError("capture dimensions do not match its preshot")
        if int(dataset.solver_checksum) != solver_preshot_checksum(preshot):
            raise ValueError("capture preshot does not match the solver checksum")
        if int(dataset.preshot_checksum) != _full_preshot_checksum(preshot):
            raise ValueError("capture preshot checksum does not match the embedded preshot")
        values = {name: _read_variable(dataset["INPUT"], name, spec, sizes) for name, spec in INPUTS.items()}
        signals = Signals(time=values["time"], meas=values["meas"], coil_signal=values["coil_signal"])
        if not np.allclose(values["coil_current"], signals.coil_current(preshot), rtol=1e-12, atol=0.0, equal_nan=True):
            raise ValueError("capture coil_current does not match COIL:MATRIX @ coil_signal")
        outputs = {
            name: _read_variable(dataset["OUTPUT"], name, (dtype, ("N_TIME", *dims)), sizes)
            for name, (dtype, dims) in OUTPUTS.items()
        }
        timed = bool(dataset.timing_enabled)
        expected_outputs = set(OUTPUTS) | ({"timing", "timing_sections"} if timed else set())
        if set(dataset["OUTPUT"].variables) != expected_outputs or set(dataset["INPUT"].variables) != set(INPUTS):
            raise ValueError("capture has unexpected input/output variables")
        timing = None
        if timed:
            sections = _read_variable(dataset["OUTPUT"], "timing_sections", ("str", ("N_TIMING",)), sizes)
            if sections != list(TIMING_SECTIONS):
                raise ValueError("capture timing sections do not match this schema")
            timing = _read_variable(dataset["OUTPUT"], "timing", ("float64", ("N_TIME", "N_TIMING")), sizes)
        slices = []
        for i in range(sizes["N_TIME"]):
            item = {
                name: value[i].item() if not OUTPUTS[name][1] else value[i].reshape(-1).copy()
                for name, value in outputs.items()
            }
            slices.append(SliceResult(**item, timing=None if timing is None else timing[i].copy()))
        capture = Capture(preshot, signals, slices, dataset.settings_json, dataset.preshot_source, dataset.code_version)
    capture.validate()
    return capture
