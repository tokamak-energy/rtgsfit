"""Local files and generated C: round trips and validation."""

from __future__ import annotations

import copy
import re
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers.local import read_preshot
from rtgsfit.database_readers.local import read_results
from rtgsfit.database_writers.c_source import data_source
from rtgsfit.database_writers.c_source import sizes_header
from rtgsfit.database_writers.c_source import write_c_source
from rtgsfit.database_writers.local import write_npz
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import PreshotError
from rtgsfit.results import Results


def _assert_same(actual: object, expected: object, key: str) -> None:
    if isinstance(expected, (str, list)):
        assert actual == expected, key
    else:
        assert_array_equal(actual, expected, err_msg=key)
        assert np.asarray(actual).dtype == np.asarray(expected).dtype, key


def test_results_round_trip(analytic_results: Results, tmp_path: Path) -> None:
    path = write_npz(tmp_path / "results.npz", analytic_results.nested)
    again = read_results(path).flat
    expected = analytic_results.flat
    assert set(again) == set(expected)
    for key, value in expected.items():
        _assert_same(again[key], value, key)


def test_preshot_round_trip(analytic_preshot: Preshot, tmp_path: Path) -> None:
    path = write_npz(tmp_path / "preshot.npz", analytic_preshot.to_nested())
    again = read_preshot(path)
    assert again.content_hash() == analytic_preshot.content_hash()
    assert again["FIT:B_VEC_NAMES"] == analytic_preshot["FIT:B_VEC_NAMES"]


def _modified(preshot: Preshot, key: str, value: object) -> dict:
    groups = copy.deepcopy(preshot.groups)
    group, _, name = key.partition(":")
    groups[group][name] = value
    return groups


@pytest.mark.parametrize(
    ("key", "change"),
    [
        ("GRID:R", lambda value: value[:-1]),
        ("GRID:G_COIL", lambda value: value.astype(np.float32)),
        ("FIT:WEIGHT", lambda value: np.where(np.arange(len(value)) == 0, np.nan, value)),
        ("GRID:MASK_LIMITER", lambda value: value * 2),
        ("POISSON:PERMUTATION", lambda value: np.zeros_like(value)),
        ("LIMITER:GRID_INDEX", lambda value: value - 1_000_000),
        ("CONSTRAINT:NAMES", lambda value: value[::-1]),
    ],
)
def test_preshot_validation_rejects(analytic_preshot: Preshot, key: str, change) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(PreshotError):
        Preshot(_modified(analytic_preshot, key, change(analytic_preshot[key])))


def test_preshot_validation_rejects_unknown_node(analytic_preshot: Preshot) -> None:
    with pytest.raises(PreshotError):
        Preshot(_modified(analytic_preshot, "GRID:EXTRA", np.zeros(1)))


def test_c_sizes_are_enum_constants(analytic_preshot: Preshot) -> None:
    """Issue #97: array sizes must be integer constant expressions, so they are enum members."""
    header = sizes_header(analytic_preshot)
    enum_body = header[header.index("enum {") : header.index("};")]
    values = dict(re.findall(r"(\w+)\s*=\s*(\d+),", enum_body))
    assert {name: int(value) for name, value in values.items()} == analytic_preshot.sizes


def test_c_arrays_have_the_dictionary_lengths(analytic_preshot: Preshot) -> None:
    source = data_source(analytic_preshot)
    lengths = {name: int(length) for name, length in re.findall(r"const \w+ (\w+)\[(\d+)\] = \{", source)}
    sizes = analytic_preshot.symbolic_sizes()
    for node in dd.C_NODES:
        if node.shape:
            assert lengths[node.c_name] == int(np.prod([sizes[dim] for dim in node.shape])), node.c_name


def test_c_source_is_only_rewritten_when_it_changes(analytic_preshot: Preshot, tmp_path: Path) -> None:
    assert len(write_c_source(analytic_preshot, tmp_path)) == 2
    assert write_c_source(analytic_preshot, tmp_path) == []
