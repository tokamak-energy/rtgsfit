"""Read a preshot or results file written by the `local` writer.

Settings for the preshot reader: `preshot/local/path`.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import numpy as np

from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers import PreshotReader
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

if TYPE_CHECKING:
    from rtgsfit.results import Results

_STR_NODES = {node.key for node in dd.ALL_NODES if node.dtype == "str"}
_LIST_NODES = {node.key for node in dd.ALL_NODES if node.dtype == "str" and node.shape}


def read_npz(path: str | Path) -> dict[str, Any]:
    """The nested dictionary stored by `rtgsfit.database_writers.local.write_npz`."""
    flat: dict[str, Any] = {}
    with np.load(Path(path), allow_pickle=False) as data:
        for key in data.files:
            value = data[key]
            if key in _LIST_NODES:
                flat[key] = [str(item) for item in value]
            elif key in _STR_NODES:
                flat[key] = str(value)
            else:
                flat[key] = value
    return dd.unflatten(flat)


def read_results(path: str | Path) -> Results:
    from rtgsfit.results import Results

    return Results(read_npz(path))


def read_preshot(path: str | Path) -> Preshot:
    return Preshot.from_nested(read_npz(path))


class LocalPreshotReader(PreshotReader):
    def read_preshot(self, settings: Settings) -> Preshot:
        return read_preshot(settings["preshot"]["local"]["path"])
