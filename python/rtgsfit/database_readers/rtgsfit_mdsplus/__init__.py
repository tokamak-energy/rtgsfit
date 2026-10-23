"""Read the preshot from the RTGSFIT MDSplus tree (TE network only).

Settings, `preshot/rtgsfit_mdsplus/`: `pulse` (e.g. 52_014_922), `run_name` (e.g. `RUN21`) and
optionally `server` (default `smaug`). The preshot is written there by
`python -m rtgsfit preshot <settings>` with the `rtgsfit_mdsplus` writer.

`mdsthin` is imported lazily because it is only needed here.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers import PreshotReader
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

DEFAULT_SERVER = "smaug"


def _decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode().strip()
    if isinstance(value, np.ndarray) and value.dtype.kind in "SO":
        return np.array([_decode(item) for item in value.ravel()]).reshape(value.shape)
    return value


def read_mdsplus_nodes(
    pulse: int, run_name: str, nodes: tuple[dd.Node, ...], server: str = DEFAULT_SERVER
) -> dict[str, Any]:
    """Read `nodes` of `\\RTGSFIT::TOP.<run_name>` into a nested dictionary."""
    import mdsthin

    nested: dict[str, Any] = {}
    with mdsthin.Connection(server) as connection:
        connection.openTree(dd.TREE_NAME, pulse)
        for node in nodes:
            path = f"\\{dd.TREE_NAME}::TOP.{run_name}.{node.mds_path}"
            dd.nested_set(nested, node.path, _decode(connection.get(path).data()))
    return nested


class RtgsfitMdsplusReader(PreshotReader):
    def read_preshot(self, settings: Settings) -> Preshot:
        section = settings["preshot"]["rtgsfit_mdsplus"]
        nested = read_mdsplus_nodes(
            int(section["pulse"]), str(section["run_name"]), dd.PRESHOT_NODES, section.get("server", DEFAULT_SERVER)
        )
        return Preshot.from_nested(nested)
