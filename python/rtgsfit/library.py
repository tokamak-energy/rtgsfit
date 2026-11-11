"""Build librtgsfit.so for a preshot and call it through ctypes.

The preshot is compiled into the library. `build_library` writes it to `c/build/preshot.bin` and runs
`make`, which generates `c/src/preshot_sizes.h` and `c/src/preshot_data.c` with the C generator
`c/tools/write_preshot.c` and compiles `c/build/librtgsfit.so`. Make recompiles only what changed, so
building the same preshot again is quick. There is one build, of the last preshot built, so build from
one process at a time.
"""

from __future__ import annotations

import _ctypes
import ctypes
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from rtgsfit import data_dictionary as dd
from rtgsfit.database_writers.c_source import preshot_checksum
from rtgsfit.database_writers.c_source import write_preshot_file
from rtgsfit.preshot import Preshot

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int32]

TIMING_SECTIONS = (
    "CONSTRAINTS", "BASIS", "FIT_MATRIX", "LEAST_SQUARES", "CURRENT", "RESIDUAL", "POISSON",
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
    256: "FIT_INVALID_ARGUMENT",
    512: "FIT_NO_CONVERGENCE",
    1024: "FIT_NONFINITE_INPUT",
    2048: "FIT_RESOURCE_FAILURE",
}
"""`enum rtgsfit_status` in c/src/rtgsfit.h."""


def describe_status(status: int) -> str:
    if status == 0:
        return "OK"
    names = [name for bit, name in STATUS_FLAGS.items() if status & bit]
    unknown = status & ~sum(STATUS_FLAGS)
    if unknown:
        names.append(f"UNKNOWN_BITS({unknown:#x})")
    return " | ".join(names)


def c_directory() -> Path:
    """The `c/` directory holding the Makefile and sources (`$RTGSFIT_C_DIR` overrides)."""
    if "RTGSFIT_C_DIR" in os.environ:
        return Path(os.environ["RTGSFIT_C_DIR"])
    path = Path(__file__).resolve().parents[2] / "c"
    if not (path / "Makefile").is_file():
        raise FileNotFoundError(f"cannot find the C sources at {path}; set RTGSFIT_C_DIR")
    return path


def cache_directory() -> Path:
    """Generated outputs live beside the C sources, independent of the working directory."""
    if "RTGSFIT_CACHE_DIR" in os.environ:
        return Path(os.environ["RTGSFIT_CACHE_DIR"])
    return c_directory().resolve() / "cache"


def build_directory() -> Path:
    """`c/build`, which holds the library, its objects and the preshot file (`make clean` removes it)."""
    return c_directory().resolve() / "build"


PRESHOT_FILE = "preshot.bin"
"""The preshot file that `build_library` writes into `build_directory()` for `make`."""
LIBRARY = "librtgsfit.so"


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


def run_make(arguments: list[str], verbose: bool = False) -> None:
    """Run `make -C c <arguments>`; raise with make's error output if it fails."""
    command = ["make", "-C", str(c_directory()), *arguments]
    result = subprocess.run(command, capture_output=True, text=True)
    if verbose:
        print(result.stdout, result.stderr, sep="\n")
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed:\n{result.stderr[-4000:]}")


def build_library(preshot: Preshot, options: BuildOptions = DEFAULT_BUILD_OPTIONS, verbose: bool = False) -> Path:
    """Compile `c/build/librtgsfit.so` for `preshot` and return its path."""
    preshot_file = write_preshot_file(preshot, build_directory() / PRESHOT_FILE)
    run_make(
        [
            f"PRESHOT_FILE={preshot_file}", f"TIMING={int(options.timing)}", f"DEBUG={int(options.debug)}",
            f"ARCH_FLAGS={options.arch_flags}",
        ],
        verbose=verbose,
    )  # fmt: skip
    return build_directory() / LIBRARY


@dataclass
class _LoadedLibrary:
    file_id: tuple[int, ...]
    cdll: ctypes.CDLL
    unloaded: bool = False


_LOADED: dict[Path, _LoadedLibrary] = {}


def _load_library(path: Path) -> _LoadedLibrary:
    """Load `path`, or reuse it if it is loaded and unchanged.

    The dynamic loader returns the library already loaded from a path even if the file has been rebuilt
    since, so a rebuilt library is first unloaded. `Rtgsfit` objects that still use it then raise.
    """
    path = path.resolve()
    status = path.stat()
    file_id = (status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns)
    loaded = _LOADED.get(path)
    if loaded is not None and loaded.file_id == file_id:
        return loaded
    if loaded is not None:
        loaded.unloaded = True
        _ctypes.dlclose(loaded.cdll._handle)
    _LOADED[path] = _LoadedLibrary(file_id, ctypes.CDLL(str(path), mode=ctypes.RTLD_LOCAL))
    return _LOADED[path]


@dataclass
class SliceResult:
    """The outputs of one call to `rtgsfit()`, named as in c/src/rtgsfit.h."""

    status: int
    psi: FloatArray
    psi_norm: FloatArray
    mask: IntArray
    j_phi: FloatArray
    fit_dof: FloatArray
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
    The library keeps a static workspace, so use one `Rtgsfit` at a time. Creating an `Rtgsfit` for a
    rebuilt library (another preshot or other options) unloads the old library; earlier `Rtgsfit`
    objects then raise.
    """

    def __init__(self, library_path: str | Path, preshot: Preshot, timing: bool = False) -> None:
        self.library_path = Path(library_path)
        self.preshot = preshot
        self.timing = timing
        self._loaded = _load_library(self.library_path)
        self._lib = self._loaded.cdll
        self._lib.rtgsfit.restype = _I32
        self._lib.rtgsfit.argtypes = (
            [_DOUBLE_P, _I32] * 2 + [_DOUBLE_P, _I32, _INT32_P, _I32] + [_DOUBLE_P, _I32] * 8
            + [ctypes.POINTER(_I32)] * 2 + [ctypes.POINTER(ctypes.c_double)] * 9 + [ctypes.POINTER(_I32)]
        )  # fmt: skip
        self._lib.rtgsfit_get_sizes.restype = _I32
        self._lib.rtgsfit_get_sizes.argtypes = [_INT32_P, _I32]
        self._lib.rtgsfit_get_preshot_checksum.restype = ctypes.c_uint32
        self._lib.rtgsfit_get_preshot_checksum.argtypes = []
        self._lib.rtgsfit_timing_get.restype = _I32
        self._lib.rtgsfit_timing_get.argtypes = [_DOUBLE_P, _I32]
        self._lib.rtgsfit_timing_reset.restype = None
        self._lib.rtgsfit_timing_reset.argtypes = []
        self._check_preshot()
        self.sizes = preshot.sizes
        self.reset()

    @classmethod
    def build(cls, preshot: Preshot, options: BuildOptions = DEFAULT_BUILD_OPTIONS) -> Rtgsfit:
        return cls(build_library(preshot, options), preshot, timing=options.timing)

    def _check_preshot(self) -> None:
        """Check that the library was compiled from `self.preshot`: first the sizes, then every value."""
        sizes = np.zeros(len(dd.SIZE_NAMES), dtype=np.int32)
        if self._lib.rtgsfit_get_sizes(sizes, len(sizes)) != 0:
            raise RuntimeError("rtgsfit_get_sizes failed: the library and the data dictionary disagree")
        compiled = dict(zip(dd.SIZE_NAMES, sizes.tolist(), strict=True))
        if compiled != self.preshot.sizes:
            raise RuntimeError(f"{self.library_path} was built for sizes {compiled}, not {self.preshot.sizes}")
        checksum = int(self._lib.rtgsfit_get_preshot_checksum())
        expected = preshot_checksum(self.preshot)
        if checksum != expected:
            raise RuntimeError(
                f"{self.library_path} was built from another preshot (checksum {checksum:08X}, expected "
                f"{expected:08X}); build it with Rtgsfit.build(preshot)"
            )

    def reset(self) -> None:
        """Restart from the initial conditions in the preshot."""
        self.psi_norm = np.array(self.preshot["INITIAL:PSI_NORM"], dtype=np.float64)
        self.mask = np.array(self.preshot["INITIAL:MASK"], dtype=np.int32)

    def __call__(self, meas: FloatArray, coil_current: FloatArray) -> SliceResult:
        if self._loaded.unloaded:
            raise RuntimeError(
                f"{self.library_path} was rebuilt and reloaded after this Rtgsfit was created; "
                "create a new one with Rtgsfit.build(preshot)"
            )
        sizes = self.sizes
        meas = np.ascontiguousarray(meas, dtype=np.float64)
        coil_current = np.ascontiguousarray(coil_current, dtype=np.float64)
        if meas.shape != (sizes["N_MEAS"],) or coil_current.shape != (sizes["N_COIL"],):
            raise ValueError(
                f"rtgsfit() needs meas of shape ({sizes['N_MEAS']},) and coil_current of shape "
                f"({sizes['N_COIL']},); got {meas.shape} and {coil_current.shape}"
            )
        arrays = {
            "psi": sizes["N_GRID"], "j_phi": sizes["N_GRID"], "fit_dof": sizes["N_FIT_DOF"],
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
