"""Shared fixtures and options.

* Tests marked `te_network` need the Tokamak Energy network (MDSplus, GSFit's ST40 readers) and are
  skipped unless `--te-network` is given.
* Tests marked `freegsnke` are skipped when GSFit or FreeGSNKE cannot be imported.
* Everything else runs anywhere with numpy, scipy, contourpy, a C compiler, OpenBLAS and LAPACKE.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
from rtgsfit.preshot import Preshot
from rtgsfit.replay import read_preshot
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

DATA_DIR = Path(__file__).parent / "data"


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--te-network", action="store_true", help="run the tests that need the Tokamak Energy network")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    have_freegsnke = all(importlib.util.find_spec(name) is not None for name in ("gsfit", "freegsnke"))
    for item in items:
        if "te_network" in item.keywords and not config.getoption("--te-network"):
            item.add_marker(pytest.mark.skip(reason="needs the TE network; run with --te-network"))
        if "freegsnke" in item.keywords and not have_freegsnke:
            item.add_marker(pytest.mark.skip(reason="needs gsfit and freegsnke"))


@pytest.fixture(scope="session")
def analytic_preshot() -> Preshot:
    """Building it takes ~15 s (a dense LU of the Poisson matrix), so build it once; it does not
    depend on the plasma current."""
    return read_preshot(load_settings("analytic"))


@pytest.fixture(scope="session")
def analytic_results(analytic_preshot: Preshot) -> Results:
    """The analytic replay: 10 calls of `rtgsfit()` with the same measurements."""
    return replay(load_settings("analytic"), preshot=analytic_preshot, verbose=False)


@pytest.fixture(scope="session")
def analytic_zero_results(analytic_preshot: Preshot) -> Results:
    """The analytic replay with no plasma current: every call must fail cleanly with NO_AXIS."""
    settings = load_settings("analytic").copy_with({"analytic/plasma/ip": 0.0})
    return replay(settings, preshot=analytic_preshot, verbose=False)
