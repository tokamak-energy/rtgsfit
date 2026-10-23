"""Build librtgsfit.so for a preshot and call it through ctypes.

The preshot is compiled into the library, so each preshot needs its own build. Builds are cached
in `$RTGSFIT_CACHE_DIR` (default `~/.cache/rtgsfit`) under a key made from the preshot content, the
C sources and the build options, so replaying the same preshot again does not recompile.
"""

from __future__ import annotations

import ctypes
import fcntl
import hashlib
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from rtgsfit import data_dictionary as dd
from rtgsfit.database_writers.c_source import write_c_source
from rtgsfit.preshot import Preshot

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int32]

TIMING_SECTIONS = (
    "CONSTRAINTS", "BASIS", "A_MATRIX", "LEAST_SQUARES", "CURRENT", "RESIDUAL", "POISSON",
    "EXTERNAL_PSI", "NULLS", "BOUNDARY", "FLOOD_FILL", "NORMALISE", "TOTAL",
)  # fmt: skip
"""Order of `enum rtgsfit_timer` in c/src/rtgsfit.h."""

STATUS_FLAGS = {
    1: "AXIS_EQUALS_BOUNDARY",
    2: "NO_AXIS",
    4: "TOO_MANY_XPTS",
    8: "TOO_MANY_OPTS",
    16: "AXIS_OUTSIDE_CORE",
    32: "BOUNDARY_ABOVE_AXIS",
    64: "ARGUMENT_SIZE",
    128: "LAPACK",
}
"""`enum rtgsfit_status` in c/src/rtgsfit.h."""


def describe_status(status: int) -> str:
    if status == 0:
        return "OK"
    return " | ".join(name for bit, name in STATUS_FLAGS.items() if status & bit)


def c_directory() -> Path:
    """The `c/` directory holding the Makefile and sources (`$RTGSFIT_C_DIR` overrides)."""
    if "RTGSFIT_C_DIR" in os.environ:
        return Path(os.environ["RTGSFIT_C_DIR"])
    path = Path(__file__).resolve().parents[2] / "c"
    if not (path / "Makefile").is_file():
        raise FileNotFoundError(f"cannot find the C sources at {path}; set RTGSFIT_C_DIR")
    return path


def cache_directory() -> Path:
    return Path(os.environ.get("RTGSFIT_CACHE_DIR", Path.home() / ".cache" / "rtgsfit"))


def _source_hash(c_dir: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted([c_dir / "Makefile", *(c_dir / "src").glob("*.[ch]")]):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()[:16]


@dataclass(frozen=True)
class BuildOptions:
    timing: bool = False
    arch_flags: str = "-march=native"
    debug: bool = False

    @classmethod
    def from_settings(cls, settings: dict[str, Any]) -> BuildOptions:
        build = settings.get("build", {})
        return cls(
            timing=bool(build.get("timing", False)),
            arch_flags=str(build.get("arch_flags", "-march=native")),
            debug=bool(build.get("debug", False)),
        )


DEFAULT_BUILD_OPTIONS = BuildOptions()


def _resolved_arch(arch_flags: str) -> str:
    """The target that `-march=native` resolves to on this machine, so that hosts with different
    CPUs sharing a cache directory never load each other's binaries."""
    if "native" not in arch_flags:
        return ""
    compiler = os.environ.get("CC", "gcc")
    result = subprocess.run(
        [compiler, *arch_flags.split(), "-Q", "--help=target"], capture_output=True, text=True, check=False
    )
    return hashlib.sha256(result.stdout.encode()).hexdigest() if result.returncode == 0 else os.uname().nodename


def build_library(preshot: Preshot, options: BuildOptions = DEFAULT_BUILD_OPTIONS, verbose: bool = False) -> Path:
    """Compile librtgsfit.so for `preshot` (or reuse a cached build) and return its path."""
    c_dir = c_directory()
    identity = f"{preshot.content_hash()} {_source_hash(c_dir)} {options} {_resolved_arch(options.arch_flags)}"
    key = hashlib.sha256(identity.encode()).hexdigest()[:20]
    build_root = cache_directory() / key
    library = build_root / "build" / "librtgsfit.so"
    build_root.mkdir(parents=True, exist_ok=True)
    with open(build_root / ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not library.is_file():
            write_c_source(preshot, build_root / "preshot")
            command = [
                "make", "-C", str(c_dir), f"PRESHOT_DIR={build_root / 'preshot'}", f"BUILD_DIR={build_root / 'build'}",
                f"TIMING={int(options.timing)}", f"DEBUG={int(options.debug)}", f"ARCH_FLAGS={options.arch_flags}",
                f"PYTHON={sys.executable}",
            ]  # fmt: skip
            result = subprocess.run(command, capture_output=True, text=True)
            if verbose or result.returncode != 0:
                print(result.stdout, result.stderr, sep="\n")
            if result.returncode != 0:
                raise RuntimeError(f"building librtgsfit failed: {' '.join(command)}")
    return library


@dataclass
class SliceResult:
    """The outputs of one call to `rtgsfit()`, named as in c/src/rtgsfit.h."""

    status: int
    psi: FloatArray
    psi_norm: FloatArray
    mask: IntArray
    j_phi: FloatArray
    x_vec: FloatArray
    constraint_measured: FloatArray
    constraint_reconstructed: FloatArray
    xpt_r: FloatArray
    xpt_z: FloatArray
    xpt_psi: FloatArray
    xpt_n: int
    xpt_diverted: int
    ip: float
    axis_psi: float
    axis_r: float
    axis_z: float
    boundary_psi: float
    current_centroid_r: float
    current_centroid_z: float
    chi_sq_magnetic: float
    chi_sq_total: float
    lapack_info: int
    timing: FloatArray | None  # [s], ordered as TIMING_SECTIONS


_DOUBLE_P = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_INT32_P = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
_I32 = ctypes.c_int32


class Rtgsfit:
    """One compiled RT-GSFit library plus the iteration state (`psi_norm` and `mask`).

    ```python
    rtgsfit = Rtgsfit.build(preshot)
    result = rtgsfit(meas, coil_current)
    ```
    The library keeps a static workspace, so only one `Rtgsfit` per library should be used at a time.
    """

    def __init__(self, library_path: str | Path, preshot: Preshot, timing: bool = False) -> None:
        self.library_path = Path(library_path)
        self.preshot = preshot
        self.timing = timing
        self._lib = ctypes.CDLL(str(self.library_path), mode=ctypes.RTLD_LOCAL)
        self._lib.rtgsfit.restype = _I32
        self._lib.rtgsfit.argtypes = (
            [_DOUBLE_P, _I32] * 2 + [_DOUBLE_P, _I32, _INT32_P, _I32] + [_DOUBLE_P, _I32] * 8
            + [ctypes.POINTER(_I32)] * 2 + [ctypes.POINTER(ctypes.c_double)] * 9 + [ctypes.POINTER(_I32)]
        )  # fmt: skip
        self._lib.rtgsfit_get_sizes.restype = _I32
        self._lib.rtgsfit_get_sizes.argtypes = [_INT32_P, _I32]
        self._lib.rtgsfit_timing_get.restype = _I32
        self._lib.rtgsfit_timing_get.argtypes = [_DOUBLE_P, _I32]
        self._lib.rtgsfit_timing_reset.restype = None
        self._lib.rtgsfit_timing_reset.argtypes = []
        self._check_sizes()
        self.sizes = preshot.sizes
        self.reset()

    @classmethod
    def build(cls, preshot: Preshot, options: BuildOptions = DEFAULT_BUILD_OPTIONS) -> Rtgsfit:
        return cls(build_library(preshot, options), preshot, timing=options.timing)

    def _check_sizes(self) -> None:
        sizes = np.zeros(len(dd.SIZE_NAMES), dtype=np.int32)
        if self._lib.rtgsfit_get_sizes(sizes, len(sizes)) != 0:
            raise RuntimeError("rtgsfit_get_sizes failed: the library and the data dictionary disagree")
        compiled = dict(zip(dd.SIZE_NAMES, sizes.tolist(), strict=True))
        if compiled != self.preshot.sizes:
            raise RuntimeError(f"{self.library_path} was built for sizes {compiled}, not {self.preshot.sizes}")

    def reset(self) -> None:
        """Restart from the initial conditions in the preshot."""
        self.psi_norm = np.array(self.preshot["INITIAL:PSI_NORM"], dtype=np.float64)
        self.mask = np.array(self.preshot["INITIAL:MASK"], dtype=np.int32)

    def __call__(self, meas: FloatArray, coil_current: FloatArray) -> SliceResult:
        sizes = self.sizes
        meas = np.ascontiguousarray(meas, dtype=np.float64)
        coil_current = np.ascontiguousarray(coil_current, dtype=np.float64)
        if meas.shape != (sizes["N_MEAS"],) or coil_current.shape != (sizes["N_COIL"],):
            raise ValueError(
                f"rtgsfit() needs meas of shape ({sizes['N_MEAS']},) and coil_current of shape "
                f"({sizes['N_COIL']},); got {meas.shape} and {coil_current.shape}"
            )
        arrays = {
            "psi": sizes["N_GRID"], "j_phi": sizes["N_GRID"], "x_vec": sizes["N_X_VEC"],
            "constraint_measured": sizes["N_CONSTRAINT"], "constraint_reconstructed": sizes["N_CONSTRAINT"],
            "xpt_r": sizes["N_XPT_MAX"], "xpt_z": sizes["N_XPT_MAX"], "xpt_psi": sizes["N_XPT_MAX"],
        }  # fmt: skip
        out = {name: np.zeros(size) for name, size in arrays.items()}
        int_scalars = {name: _I32(0) for name in ("xpt_n", "xpt_diverted", "lapack_info")}
        scalars = {
            name: ctypes.c_double(0.0)
            for name in (
                "ip", "axis_psi", "axis_r", "axis_z", "boundary_psi", "current_centroid_r",
                "current_centroid_z", "chi_sq_magnetic", "chi_sq_total",
            )
        }  # fmt: skip
        if self.timing:
            self._lib.rtgsfit_timing_reset()
        array_args: list[Any] = []
        for name in arrays:
            array_args += [out[name], len(out[name])]
        status = self._lib.rtgsfit(
            meas, len(meas), coil_current, len(coil_current),
            self.psi_norm, len(self.psi_norm), self.mask, len(self.mask),
            *array_args,
            ctypes.byref(int_scalars["xpt_n"]), ctypes.byref(int_scalars["xpt_diverted"]),
            *(ctypes.byref(value) for value in scalars.values()),
            ctypes.byref(int_scalars["lapack_info"]),
        )  # fmt: skip
        timing = None
        if self.timing:
            timing = np.zeros(len(TIMING_SECTIONS))
            self._lib.rtgsfit_timing_get(timing, len(timing))
            timing *= 1e-6
        return SliceResult(
            status=int(status),
            psi_norm=self.psi_norm.copy(),
            mask=self.mask.copy(),
            **out,
            **{name: int(value.value) for name, value in int_scalars.items()},
            **{name: float(value.value) for name, value in scalars.items()},
            timing=timing,
        )
