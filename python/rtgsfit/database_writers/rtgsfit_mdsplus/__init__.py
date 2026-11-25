"""Write the preshot or the results to the RTGSFIT MDSplus tree (TE network only).

The tree layout is created from `standard_utility/csv_files/RTGSFIT.csv`, which is generated from
`rtgsfit.data_dictionary` by `python -m rtgsfit write-csv`, so the physics data match a local `.npz`
file. The common `CODE_VERSION` metadata records the code environment and the UTC start of the data write.

Settings, `writers/rtgsfit_mdsplus/`:

* `pulse`: pulse number to write to (for a preshot, e.g. 52_014_922);
* `run_name`: e.g. `RUN01`;
* `run_description`: free text;
* `link_best`: point `BEST` at this run;
* `enabled`: must be `true`. Writing to the shared tree is never done by accident.

`standard_utility` is imported lazily because it is specific to Tokamak Energy.
"""

from __future__ import annotations

import platform
import socket
import subprocess
from datetime import UTC
from datetime import datetime
from getpass import getuser
from importlib.metadata import distributions
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import numpy as np

from rtgsfit import __version__
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
            converted[key] = np.array([str(item) for item in value], dtype=np.str_)
        else:
            converted[key] = value
    return converted


def _code_version() -> dict[str, str]:
    git_id = ""
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            git_id = subprocess.run(
                ["git", "-C", str(parent), "rev-parse", "--short", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            break
    library = sorted(
        f"{distribution.metadata['Name']}=={distribution.version}"
        for distribution in distributions()
        if distribution.metadata.get("Name")
    )
    return {
        "COMPUTER": socket.gethostname(),
        "DATETIME": datetime.now(UTC).isoformat(),
        "GIT_ID": git_id,
        "LIBRARY": ", ".join(library),
        "MATLAB": "N/A",
        "PYTHON": platform.python_version(),
        "USER": getuser(),
        "VERSION": __version__,
    }


def _write_quality(pulse: int, run_name: str, quality: dict[str, Any]) -> None:
    import standard_utility as util

    util.create_write_nodes_from_csv_file(
        pulse,
        "DATA_Q",
        tree=dd.TREE_NAME,
        base_node=run_name,
        run_name="",
        selected_nodes="select_from_data",
        operation_type="write",
        data_to_write=quality,
        force_write=True,
        use_mdsthin=True,
    )


def write_quality_only(quality: dict[str, Any], settings: Settings) -> None:
    """Write a failed-input report without replacing an existing run or linking BEST."""
    import mdsthin
    import standard_utility as util
    from mdsthin.exceptions import TreeFOPENR
    from mdsthin.exceptions import TreeNNF

    section = _section(settings)
    pulse, run_name = int(section["pulse"]), str(section["run_name"])
    if not run_name.isidentifier() or len(run_name) > 12 or run_name.upper() == "BEST":
        raise ValueError(f"invalid MDSplus run name {run_name!r}")
    check_standard_utility_csv(util)
    with mdsthin.Connection("smaug", timeout=30) as connection:
        try:
            connection.openTree(dd.TREE_NAME, pulse)
        except TreeFOPENR:
            pass
        else:
            try:
                try:
                    connection.get(f'getnci("\\\\{dd.TREE_NAME}::TOP.{run_name}", "NID_NUMBER")')
                except TreeNNF:
                    pass
                else:
                    raise FileExistsError(f"{dd.TREE_NAME} pulse {pulse} {run_name} already exists")
            finally:
                connection.closeTree(dd.TREE_NAME, pulse)
    util.create_script_nodes(
        script_name=dd.TREE_NAME,
        pulseNo_write=pulse,
        run_name=run_name,
        run_info=str(section.get("run_description", "")),
        link_best=False,
        use_mdsthin=True,
    )
    util.write_script_data(
        script_name=dd.TREE_NAME,
        pulseNo_write=pulse,
        data_to_write={"CODE_VERSION": _code_version()},
        run_name=run_name,
        run_description=str(section.get("run_description", "")),
        force_write=True,
        use_mdsthin=True,
    )
    _write_quality(pulse, run_name, quality)


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
        link_best=False,
        use_mdsthin=True,
    )
    data_to_write = _to_mdsplus(data)
    quality = data_to_write.pop("DATA_Q", None)
    code_version = data_to_write.setdefault("CODE_VERSION", {})
    for key, value in _code_version().items():
        if key == "DATETIME":
            code_version[key] = value
        else:
            code_version.setdefault(key, value)
    util.write_script_data(
        script_name=dd.TREE_NAME,
        pulseNo_write=pulse,
        data_to_write=data_to_write,
        run_name=run_name,
        run_description=run_description,
        force_write=True,
        use_mdsthin=True,
    )
    if quality is not None:
        _write_quality(pulse, run_name, quality)
    if section.get("link_best", False):
        util.link_BEST_direct(pulse, f"{dd.TREE_NAME}.{run_name}", use_mdsthin=True)
    print(f"Wrote to MDSplus {dd.TREE_NAME} pulse {pulse} {run_name}")


class RtgsfitMdsplusWriter(DatabaseWriter):
    def write_preshot(self, preshot: Preshot, settings: Settings) -> None:
        write_to_mdsplus(preshot.to_nested(), settings)

    def write_results(self, results: Results, settings: Settings) -> None:
        write_to_mdsplus(results.nested, settings)
