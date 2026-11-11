"""Shared fixtures and options.

* Tests marked `te_network` need the Tokamak Energy network (MDSplus, GSFit's ST40 readers) and are
  skipped unless `--te-network` is given.
* Tests marked `freegsnke` are skipped when GSFit or FreeGSNKE cannot be imported.
* Tests marked `timing` (the ST40 timing benchmark) are skipped unless `--timing` is given.
* Tests marked `figures` draw diagnostic figures (see `plotting.py`) and are skipped unless
  `--figures` is given. Each run writes them to a new `tests/figures/<UTC time>/`.
* Everything else runs anywhere with numpy, scipy, contourpy, a C compiler, OpenBLAS and LAPACKE.
"""

from __future__ import annotations

import importlib.util
import itertools
import platform
import shlex
import socket
import subprocess
from datetime import UTC
from datetime import datetime
from pathlib import Path

import pytest
from rtgsfit.preshot import Preshot
from rtgsfit.replay import read_preshot
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

DATA_DIR = Path(__file__).parent / "data"
FIGURES_DIR = Path(__file__).parent / "figures"
FIGURE_DIR_KEY = pytest.StashKey[Path]()


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--te-network", action="store_true", help="run the tests that need the Tokamak Energy network")
    parser.addoption(
        "--timing", action="store_true", help="run the timing benchmark; writes a new CSV to c/cache/timing/"
    )
    parser.addoption(
        "--figures", action="store_true", help="draw diagnostic figures into a new tests/figures/<UTC time>/"
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    have_freegsnke = all(importlib.util.find_spec(name) is not None for name in ("gsfit", "freegsnke"))
    for item in items:
        # Markers, not `item.keywords`: keywords include parametrize ids such as `[freegsnke]`.
        if item.get_closest_marker("te_network") is not None and not config.getoption("--te-network"):
            item.add_marker(pytest.mark.skip(reason="needs the TE network; run with --te-network"))
        if item.get_closest_marker("freegsnke") is not None and not have_freegsnke:
            item.add_marker(pytest.mark.skip(reason="needs gsfit and freegsnke"))
        if item.get_closest_marker("timing") is not None and not config.getoption("--timing"):
            item.add_marker(pytest.mark.skip(reason="timing benchmark; run with --timing"))
        if item.get_closest_marker("figures") is not None and not config.getoption("--figures"):
            item.add_marker(pytest.mark.skip(reason="diagnostic figures; run with --figures"))


def pytest_terminal_summary(terminalreporter: pytest.TerminalReporter, config: pytest.Config) -> None:
    if FIGURE_DIR_KEY in config.stash:
        terminalreporter.write_line(f"Figures written to {config.stash[FIGURE_DIR_KEY]}")


@pytest.fixture(scope="session")
def git_commit() -> str:
    """Short hash of HEAD, with `-dirty` if the working tree has changes."""
    root = Path(__file__).resolve().parents[1]
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, check=True
        ).stdout.strip()
        changes = subprocess.run(
            ["git", "status", "--porcelain"], cwd=root, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return f"{commit}-dirty" if changes else commit


@pytest.fixture(scope="session")
def figure_dir(request: pytest.FixtureRequest, git_commit: str) -> Path:
    """`tests/figures/<UTC time>/`, new for each pytest run, so that earlier figures are never overwritten.

    `run.txt` in it records when, where and from which commit the figures were made.
    """
    now = datetime.now(UTC)
    stamp = now.strftime("%Y%m%dT%H%M%SZ")
    for n in itertools.count():
        path = FIGURES_DIR / (stamp if n == 0 else f"{stamp}_{n}")
        try:
            path.mkdir(parents=True)
        except FileExistsError:
            continue
        break
    path.joinpath("run.txt").write_text(
        f"time: {now:%Y-%m-%dT%H:%M:%SZ}\n"
        f"commit: {git_commit}\n"
        f"host: {socket.gethostname()}\n"
        f"python: {platform.python_version()}\n"
        f"command: pytest {shlex.join(request.config.invocation_params.args)}\n"
    )
    request.config.stash[FIGURE_DIR_KEY] = path
    return path


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
