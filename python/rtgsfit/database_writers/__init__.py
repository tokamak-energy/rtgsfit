"""Database writers: where the preshot and the results go.

Writers are selected by name in the settings (`writers/names`), and each has its own settings
section, e.g. `writers/local/path`:

| writer            | preshot | results | needs                      |
|-------------------|---------|---------|----------------------------|
| `local`           | yes     | yes     | nothing (writes a `.npz`)  |
| `c_source`        | yes     |         | a C compiler (writes c/src)|
| `rtgsfit_mdsplus` | yes     | yes     | the TE network             |

Every writer stores the same nodes, those of `rtgsfit.data_dictionary`, so a local file and the
MDSplus tree hold the same names and values.
"""

from __future__ import annotations

import importlib
from abc import ABC
from abc import abstractmethod
from typing import TYPE_CHECKING

from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

if TYPE_CHECKING:
    from rtgsfit.results import Results


class DatabaseWriter(ABC):
    @abstractmethod
    def write_preshot(self, preshot: Preshot, settings: Settings) -> None:
        """Store the preshot alone (the `PRESHOT` subtree)."""

    @abstractmethod
    def write_results(self, results: Results, settings: Settings) -> None:
        """Store the results of a replay, including a copy of the preshot."""


_WRITERS = {
    "local": "rtgsfit.database_writers.local:LocalWriter",
    "c_source": "rtgsfit.database_writers.c_source:CSourceWriter",
    "rtgsfit_mdsplus": "rtgsfit.database_writers.rtgsfit_mdsplus:RtgsfitMdsplusWriter",
}


def get_database_writer(name: str) -> DatabaseWriter:
    if name not in _WRITERS:
        raise ValueError(f"unknown database writer {name!r}; available: {', '.join(_WRITERS)}")
    module_name, class_name = _WRITERS[name].split(":")
    writer = getattr(importlib.import_module(module_name), class_name)()
    assert isinstance(writer, DatabaseWriter)
    return writer
