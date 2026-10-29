"""Write the preshot or the results to the RTGSFIT MDSplus tree (TE network only).

The tree layout is created from `standard_utility/csv_files/RTGSFIT.csv`, which is generated from
`rtgsfit.data_dictionary` by `python -m rtgsfit write-csv`, so the nodes written here are exactly the
nodes of a local `.npz` file.

Settings, `writers/rtgsfit_mdsplus/`:

* `pulse`: pulse number to write to (for a preshot, e.g. 52_014_922);
* `run_name`: e.g. `RUN01`;
* `run_description`: free text;
* `link_best`: point `BEST` at this run;
* `enabled`: must be `true`. Writing to the shared tree is never done by accident.

`standard_utility` is imported lazily because it is specific to Tokamak Energy.
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


def _section(settings: Settings) -> dict[str, Any]:
    section = settings.get("writers", {}).get("rtgsfit_mdsplus", {})
    if not section.get("enabled", False):
        raise PermissionError(
            "writing to MDSplus is disabled; set writers/rtgsfit_mdsplus/enabled to true in the settings"
        )
    for key in ("pulse", "run_name"):
        if key not in section:
            raise ValueError(f"writers/rtgsfit_mdsplus/{key} must be set")
    return section


def check_standard_utility_csv(standard_utility: Any) -> None:
    """Refuse to write if the installed RTGSFIT.csv differs from the data dictionary: the tree would
    be created with the wrong nodes."""
    path = Path(standard_utility.__file__).parent / "csv_files" / "RTGSFIT.csv"
    if not path.is_file() or path.read_text() != dd.generate_csv():
        raise RuntimeError(
            f"{path} does not match rtgsfit.data_dictionary. Update it with `python -m rtgsfit write-csv`, "
            "commit it to standard_utility, and reinstall standard_utility in this environment."
        )


def _to_mdsplus(data: dict[str, Any]) -> dict[str, Any]:
    """Text lists become numpy string arrays, as GSFit writes them."""
    converted: dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, dict):
            converted[key] = _to_mdsplus(value)
        elif isinstance(value, (list, tuple)):
            converted[key] = np.array([str(item) for item in value])
        else:
            converted[key] = value
    return converted


def write_to_mdsplus(data: dict[str, Any], settings: Settings) -> None:
    import standard_utility as util

    section = _section(settings)
    check_standard_utility_csv(util)
    pulse = int(section["pulse"])
    run_name = str(section["run_name"])
    run_description = str(section.get("run_description", ""))
    util.create_script_nodes(
        script_name=dd.TREE_NAME,
        pulseNo_write=pulse,
        run_name=run_name,
        run_info=run_description,
        link_best=bool(section.get("link_best", False)),
        use_mdsthin=True,
    )
    util.write_script_data(
        script_name=dd.TREE_NAME,
        pulseNo_write=pulse,
        data_to_write=_to_mdsplus(data),
        run_name=run_name,
        run_description=run_description,
        force_write=True,
        use_mdsthin=True,
    )
    print(f"Wrote to MDSplus {dd.TREE_NAME} pulse {pulse} {run_name}")


class RtgsfitMdsplusWriter(DatabaseWriter):
    def write_preshot(self, preshot: Preshot, settings: Settings) -> None:
        write_to_mdsplus(preshot.to_nested(), settings)

    def write_results(self, results: Results, settings: Settings) -> None:
        write_to_mdsplus(results.nested, settings)
