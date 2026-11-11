"""Timing benchmark of `rtgsfit()` on ST40: opt-in, run with `pytest --te-network --timing`.

As in v1 (#100), each vs-GSFit case is converged at a fixed time and `rtgsfit()` is then called
`N_CALLS` times from that converged state. The PCS path is timed: the MDSplus preshot and the
PCS-equivalent signals. The thread CPU time of each section (`TIMING_SECTIONS`) is summarised and
printed. Each pytest run writes a new `c/cache/timing/rtgsfit_timing_<UTC time>.csv` (not
committed), with rows tagged by time, commit and host so that runs can be compared.
"""

from __future__ import annotations

import csv
import socket
from datetime import UTC
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit.library import TIMING_SECTIONS
from rtgsfit.library import BuildOptions
from rtgsfit.library import Rtgsfit
from rtgsfit.library import cache_directory
from rtgsfit.replay import read_preshot
from rtgsfit.replay import read_signals
from rtgsfit.settings import load_settings

pytestmark = [pytest.mark.te_network, pytest.mark.timing]

CASES = [(pulse, time) for pulse in (13343, 13345, 13346) for time in (0.030, 0.100)]  # as test_st40.py
N_WARM_UP = 10
N_CALLS = 1_000
STATISTICS = ("mean", "median", "std", "min", "max", "p95", "p99")
COLUMNS = ["timestamp", "commit", "hostname", "pulse", "time_s", "n_calls", "section"]
COLUMNS += [f"{statistic}_us" for statistic in STATISTICS]


def _statistics(column: np.ndarray) -> dict[str, float]:
    return {
        "mean": column.mean(), "median": np.median(column), "std": column.std(), "min": column.min(),
        "max": column.max(), "p95": np.percentile(column, 95), "p99": np.percentile(column, 99),
    }  # fmt: skip


@pytest.fixture(scope="module")
def timing_csv() -> Path:
    """A new CSV for this run, with the header written; `x` mode never overwrites an earlier run."""
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    path = cache_directory() / "timing" / f"rtgsfit_timing_{stamp}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="") as file:
        csv.DictWriter(file, fieldnames=COLUMNS).writeheader()
    return path


@pytest.mark.parametrize(("pulse", "time"), CASES, ids=[f"{pulse}@{time}" for pulse, time in CASES])
def test_timing(pulse: int, time: float, timing_csv: Path, git_commit: str, capsys: pytest.CaptureFixture[str]) -> None:
    settings = load_settings("st40").copy_with({
        "build/timing": True,
        "replay/pulse": pulse,
        "replay/time_start": time,
        "replay/time_end": time,
        "replay/repeat": 1,
    })  # fmt: skip
    preshot = read_preshot(settings)
    signals = read_signals(settings, preshot)
    meas, coil_current = signals.meas[0], signals.coil_current(preshot)[0]
    rtgsfit = Rtgsfit.build(preshot, BuildOptions.from_settings(settings))
    for _ in range(N_WARM_UP):
        rtgsfit(meas, coil_current)
    psi_norm, mask = rtgsfit.psi_norm.copy(), rtgsfit.mask.copy()

    timing = np.empty((N_CALLS, len(TIMING_SECTIONS)))
    status = np.empty(N_CALLS, dtype=np.int64)
    for i in range(N_CALLS):
        # Restart every call from the converged state, so that each call does the same work.
        rtgsfit.psi_norm[:] = psi_norm
        rtgsfit.mask[:] = mask
        result = rtgsfit(meas, coil_current)
        status[i] = result.status
        assert result.timing is not None
        timing[i] = result.timing * 1e6
    assert_array_equal(status, 0)

    # The sections are timed back to back, so they add up to (almost all of) TOTAL.
    total = TIMING_SECTIONS.index("TOTAL")
    parts = np.delete(timing, total, axis=1).sum(axis=1)
    assert np.all(timing >= 0.0)
    assert np.all(parts <= timing[:, total] + 1e-3)
    assert np.all(parts >= 0.9 * timing[:, total])

    common = {
        "timestamp": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit": git_commit,
        "hostname": socket.gethostname(),
        "pulse": pulse,
        "time_s": time,
        "n_calls": N_CALLS,
    }
    rows: list[dict[str, object]] = []
    lines = [f"\nrtgsfit() timing, pulse {pulse} at {time * 1e3:g} ms, {N_CALLS} calls [us]"]
    lines.append(f"{'section':<14}" + "".join(f"{statistic:>9}" for statistic in STATISTICS))
    for section, column in zip(TIMING_SECTIONS, timing.T, strict=True):
        statistics = _statistics(column)
        rows.append(common | {"section": section} | {f"{k}_us": f"{v:.3f}" for k, v in statistics.items()})
        lines.append(f"{section:<14}" + "".join(f"{value:9.1f}" for value in statistics.values()))
    with timing_csv.open("a", newline="") as file:
        csv.DictWriter(file, fieldnames=COLUMNS).writerows(rows)
    with capsys.disabled():
        print("\n".join(lines) + f"\nWritten to {timing_csv}")
