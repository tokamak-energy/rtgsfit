"""The data dictionary: the single definition of every name stored locally, in MDSplus and in C."""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

import pytest
from rtgsfit import data_dictionary as dd
from rtgsfit.constants import MDSPLUS_NODE_NAME_MAX
from rtgsfit.database_writers.rtgsfit_mdsplus import check_standard_utility_csv
from rtgsfit.library import c_directory
from rtgsfit.results import Results


def test_names_are_valid_mdsplus_names() -> None:
    dd.check_names()
    for node in dd.ALL_NODES:
        assert all(len(name) <= MDSPLUS_NODE_NAME_MAX for name in node.path), node.mds_path


def test_every_c_node_is_declared_in_preshot_h() -> None:
    header = (c_directory() / "src" / "preshot.h").read_text()
    declared = set(re.findall(r"extern const \w+ (\w+)\s*[\[;]", header))
    sizes_enum = set(dd.SIZE_NAMES)
    expected = {node.c_name for node in dd.C_NODES if node.c_name not in sizes_enum}
    assert declared == expected


def test_documentation_is_up_to_date() -> None:
    path = c_directory().parent / "documentation" / "data_dictionary.md"
    assert path.read_text() == dd.generate_markdown(), "run: python -m rtgsfit write-data-dictionary " + str(path)


def test_every_shape_size_is_documented() -> None:
    documented = set(re.findall(r"^\| `(N_\w+)` \|", dd.generate_markdown(), flags=re.MULTILINE))
    used = {size for node in dd.ALL_NODES for size in node.shape} - {"N_TIME"}
    assert used <= documented


def test_results_have_exactly_the_dictionary_nodes(analytic_results: Results) -> None:
    """`collect_results` validates this too; this guards against that check being removed."""
    assert set(analytic_results.flat) == {node.key for node in dd.ALL_NODES}


def test_flatten_round_trip() -> None:
    nested = {"A": {"B": 1, "C": {"D": [2, 3]}}, "E": "text"}
    flat = dd.flatten(nested)
    assert flat == {"A.B": 1, "A.C.D": [2, 3], "E": "text"}
    assert dd.unflatten(flat) == nested


def _csv_rows() -> list[list[str]]:
    return list(csv.reader(io.StringIO(dd.generate_csv())))


def test_csv_has_every_node_once() -> None:
    rows = _csv_rows()
    references = [cell for row in rows for cell in row if cell.startswith("data[")]
    expected = ["data" + "".join(f"['{name}']" for name in node.path) for node in dd.ALL_NODES]
    assert sorted(references) == sorted(expected)
    assert len({len(row) for row in rows}) == 1


def test_csv_loads_with_standard_utility(tmp_path: Path) -> None:
    functions = pytest.importorskip("standard_utility.functions.functions_node_information")
    path = tmp_path / "RTGSFIT.csv"
    path.write_text(dd.generate_csv())
    functions.load_node_information_from_csv_file(str(path))


def test_standard_utility_csv_is_up_to_date() -> None:
    """The rtgsfit_mdsplus writer creates the tree from the installed CSV, so it must match."""
    standard_utility = pytest.importorskip("standard_utility")
    check_standard_utility_csv(standard_utility)
