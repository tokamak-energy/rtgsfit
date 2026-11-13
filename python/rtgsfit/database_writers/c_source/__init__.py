"""Compile the preshot into C: `c/src/preshot_sizes.h` and `c/src/preshot_data.c`.

The C files are written by the C generator `c/tools/write_preshot.c`, which `make` runs. It reads the
preshot either from MDSplus (`make SHOT=<pulse> RUN_NAME=<run>`, which is what the PCS uses, so no
Python is needed there) or from a preshot file written by `write_preshot_file`, which is how a
preshot that is only in Python (analytic, FreeGSNKE, ...) gets into C.

The preshot file is binary: `RTGSFIT PRESHOT 1` and `SOURCE <INFO:SOURCE>` lines, then for every node of
`data_dictionary.C_NODES`, in order, a line `<GROUP>:<NAME> <float64|int32> <count>` followed by the
values, raw and little-endian, and finally `END`.
"""

from __future__ import annotations

import zlib
from collections.abc import Iterator
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from rtgsfit import data_dictionary as dd
from rtgsfit.database_writers import DatabaseWriter
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

if TYPE_CHECKING:
    from rtgsfit.results import Results

PRESHOT_FILE_VERSION = 1
_MAX_SOURCE_BYTES = 400


def _c_values(preshot: Preshot) -> Iterator[tuple[dd.Node, np.ndarray]]:
    """Every value compiled into C, in the order of `data_dictionary.C_NODES`, as little-endian arrays."""
    for node in dd.C_NODES:
        _, group, name = node.path
        dtype = np.dtype(node.dtype).newbyteorder("<")
        yield node, np.ascontiguousarray(preshot[f"{group}:{name}"], dtype=dtype)


def preshot_checksum(preshot: Preshot) -> int:
    """The CRC-32 of the values compiled into C; `rtgsfit_get_preshot_checksum()` returns the same."""
    checksum = 0
    for _, values in _c_values(preshot):
        checksum = zlib.crc32(values.tobytes(), checksum)
    return checksum


def write_preshot_file(preshot: Preshot, path: str | Path) -> Path:
    """Write the values compiled into C to `path`, for `make PRESHOT_FILE=<path>`."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    source = " ".join(str(preshot["INFO:SOURCE"]).split()).encode("utf-8")[:_MAX_SOURCE_BYTES]
    source = source.decode("utf-8", errors="ignore").encode("utf-8")
    with open(path, "wb") as file:
        file.write(f"RTGSFIT PRESHOT {PRESHOT_FILE_VERSION}\n".encode())
        file.write(b"SOURCE " + source + b"\n")
        for node, values in _c_values(preshot):
            file.write(f"{node.path[1]}:{node.path[2]} {node.dtype} {values.size}\n".encode())
            file.write(values.tobytes())
        file.write(b"END\n")
    return path


def write_c_source(preshot: Preshot, verbose: bool = False) -> Path:
    """Generate `preshot_sizes.h` and `preshot_data.c` in `c/src` (without compiling); returns `c/src`.

    The files are rewritten only when their content changes, so `make` recompiles only what changed.
    """
    from rtgsfit.library import PRESHOT_FILE
    from rtgsfit.library import build_directory
    from rtgsfit.library import c_directory
    from rtgsfit.library import run_make

    preshot_file = write_preshot_file(preshot, build_directory() / PRESHOT_FILE)
    run_make(["preshot", f"PRESHOT_FILE={preshot_file}"], verbose=verbose)
    return c_directory() / "src"


class CSourceWriter(DatabaseWriter):
    """Makes the preshot the active one in `c/src`, ready for `make` (see `write_c_source`)."""

    def write_preshot(self, preshot: Preshot, settings: Settings) -> None:
        print(f"Wrote the C preshot to {write_c_source(preshot)}")

    def write_results(self, results: Results, settings: Settings) -> None:
        raise NotImplementedError("the c_source writer only writes the preshot")
