"""Write the preshot or the results to a local `.npz` file, keyed by data dictionary node.

Keys are `Node.key`, e.g. `GLOBAL.IP` or `PRESHOT.GRID.R`, so a local file has exactly the nodes of
the MDSplus tree. Strings are stored as numpy unicode arrays, so the file can be read with
`numpy.load(path)` without pickling; `rtgsfit.database_readers.local` converts them back.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import numpy as np

from rtgsfit import data_dictionary as dd
from rtgsfit.database_writers import DatabaseWriter
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

if TYPE_CHECKING:
    from rtgsfit.results import Results


def _to_numpy(value: Any) -> np.ndarray:
    if isinstance(value, (list, tuple)):
        return np.array([str(item) for item in value], dtype=np.str_)
    if isinstance(value, str):
        return np.array(value, dtype=np.str_)
    return np.asarray(value)


def write_npz(path: str | Path, nested: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    flat = dd.flatten(nested)
    np.savez_compressed(path, **{key: _to_numpy(value) for key, value in flat.items()})
    return path if path.suffix == ".npz" else path.with_name(path.name + ".npz")


def _output_path(settings: Settings, default_name: str) -> Path:
    local = settings.get("writers", {}).get("local", {})
    return Path(local.get("path", default_name))


class LocalWriter(DatabaseWriter):
    def write_preshot(self, preshot: Preshot, settings: Settings) -> None:
        local = settings.get("writers", {}).get("local", {})
        path = Path(local.get("preshot_path", "rtgsfit_preshot.npz"))
        written = write_npz(path, preshot.to_nested())
        print(f"Wrote the preshot to {written}")

    def write_results(self, results: Results, settings: Settings) -> None:
        written = write_npz(_output_path(settings, "rtgsfit_results.npz"), results.nested)
        print(f"Wrote the results to {written}")
