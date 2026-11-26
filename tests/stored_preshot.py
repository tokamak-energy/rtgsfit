"""Checks shared by the ST40 mode tests."""

from __future__ import annotations

import json

import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.replay import build_preshot
from rtgsfit.replay import read_preshot
from rtgsfit.settings import Settings


def assert_mdsplus_preshot_is_the_gsfit_build(settings: Settings) -> None:
    """The PRESHOT in MDSplus (`preshot/reader`) equals a fresh build from GSFit (`preshot/builder`).

    Floating-point nodes are compared to round-off rather than by `content_hash`: the Poisson LU
    bands come from threaded LAPACK, whose rounding depends on the BLAS thread count. The archived
    `INFO:SETTINGS` must agree in the sections that make the preshot; replay and writer settings may
    change without publishing a new preshot.
    """
    stored, built = read_preshot(settings), build_preshot(settings)
    assert built.groups.keys() == stored.groups.keys()
    for group, nodes in stored.groups.items():
        assert built.groups[group].keys() == nodes.keys(), group
        for name, value in nodes.items():
            if f"{group}:{name}" == "INFO:SETTINGS":
                archived, current = json.loads(str(value)), json.loads(str(built.groups[group][name]))
                for section in ("preshot", "gsfit"):
                    assert archived[section] == current[section], f"INFO:SETTINGS {section}"
                continue
            expected, actual = np.asarray(value), np.asarray(built.groups[group][name])
            if expected.dtype.kind == "f":
                scale = float(np.abs(expected).max(initial=0.0))
                assert_allclose(actual, expected, rtol=1e-12, atol=1e-12 * scale, err_msg=f"{group}:{name}")
            else:
                assert_array_equal(actual, expected, err_msg=f"{group}:{name}")
