"""Replay RT-GSFit: run the real-time code on stored or synthetic data, exactly as the PCS would.

```python
from rtgsfit import load_settings
from rtgsfit.replay import run_replay

results = run_replay(load_settings("analytic"))
results["GLOBAL.IP"]
```

A replay is

1. read the preshot with the `preshot/reader`;
2. compile it into `librtgsfit.so` (cached; see `rtgsfit.library`);
3. read the signals with the `signals/reader`;
4. call `rtgsfit()` once per time slice, passing `psi_norm` and `mask` from one call to the next;
5. collect the results in the data-dictionary layout and pass them to every `writers/names` writer.
"""

from __future__ import annotations

import json
import time as wall_clock
from pathlib import Path

import numpy as np

from rtgsfit.capture import Capture
from rtgsfit.capture import write_capture
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import get_preshot_reader
from rtgsfit.database_readers import get_signals_reader
from rtgsfit.database_writers import get_database_writer
from rtgsfit.library import TIMING_SECTIONS
from rtgsfit.library import BuildOptions
from rtgsfit.library import Rtgsfit
from rtgsfit.library import SliceResult
from rtgsfit.library import cache_directory
from rtgsfit.library import describe_status
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results
from rtgsfit.settings import Settings


def preshot_location(settings: Settings, section: str = "preshot", key: str = "reader") -> str:
    """Where the preshot comes from, e.g. `rtgsfit_mdsplus {"pulse": 52000230, "run_name": "TEST03"}`."""
    name = str(settings[section][key])
    options = settings[section].get(name, {})
    return f"{name} {json.dumps(options, default=str)}" if options else name


def read_preshot(settings: Settings) -> Preshot:
    """The preshot used by a replay: `preshot/reader`."""
    return get_preshot_reader(str(settings["preshot"]["reader"])).read_preshot(settings)


def build_preshot(settings: Settings) -> Preshot:
    """A new preshot: `preshot/builder`, or `preshot/reader` if there is no builder."""
    preshot_settings = settings["preshot"]
    name = str(preshot_settings.get("builder", preshot_settings["reader"]))
    return get_preshot_reader(name).read_preshot(settings)


def read_signals(settings: Settings, preshot: Preshot) -> Signals:
    signals = get_signals_reader(str(settings["signals"]["reader"])).read_signals(settings, preshot)
    n_repeat = int(settings["replay"].get("repeat", 1))
    if n_repeat > 1:
        signals = signals.repeat(n_repeat)
    signals.validate(preshot)
    return signals


def _stop_reason(result: SliceResult, settings: Settings) -> str | None:
    replay = settings["replay"]
    if replay.get("stop_on_status", False) and result.status != 0:
        return f"status {describe_status(result.status)}"
    stop_below_ip = float(replay.get("stop_below_ip", 0.0))
    if stop_below_ip > 0.0 and abs(result.ip) < stop_below_ip:
        return f"|ip| = {abs(result.ip):.4g} A < {stop_below_ip:.4g} A"
    return None


def print_timing(slices: list[SliceResult]) -> None:
    """Per-section timing statistics over the successful slices."""
    timed = np.array([item.timing for item in slices if item.timing is not None and item.status == 0]) * 1e6
    if len(timed) == 0:
        return
    print(f"\nTiming of {len(timed)} successful slices [us]")
    print(f"{'section':<14}{'mean':>10}{'median':>10}{'p99':>10}{'max':>10}")
    for i, section in enumerate(TIMING_SECTIONS):
        column = timed[:, i]
        print(f"{section:<14}{column.mean():10.1f}{np.median(column):10.1f}"
              f"{np.percentile(column, 99):10.1f}{column.max():10.1f}")  # fmt: skip


def replay(
    settings: Settings,
    preshot: Preshot | None = None,
    signals: Signals | None = None,
    verbose: bool = True,
) -> Results:
    """Replay without result writers; optionally write native `capture/path` before processing.

    `preshot` and `signals` override the readers in `settings`.
    """
    if preshot is None:
        preshot = read_preshot(settings)
        source = preshot_location(settings)
    else:
        source = "supplied by the caller"
    if signals is None:
        signals = read_signals(settings, preshot)
    else:
        signals.validate(preshot)

    rtgsfit = Rtgsfit.build(preshot, BuildOptions.from_settings(settings))
    if verbose:
        print(f"Using C library: {rtgsfit.library_path}")
    coil_current = signals.coil_current(preshot)

    slices: list[SliceResult] = []
    start = wall_clock.perf_counter()
    for i, time in enumerate(signals.time):
        result = rtgsfit(signals.meas[i], coil_current[i])
        slices.append(result)
        if verbose:
            print(f"t = {time * 1e3:8.2f} ms  ip = {result.ip:11.4e} A  chi_sq = {result.chi_sq_magnetic:9.3e}  "
                  f"status = {describe_status(result.status)}")  # fmt: skip
        reason = _stop_reason(result, settings)
        if reason is not None:
            if verbose:
                print(f"Stopping the replay at t = {time:.4f} s: {reason}")
            break
    if verbose:
        elapsed = wall_clock.perf_counter() - start
        print(f"{len(slices)} slices in {elapsed:.3f} s")
        print_timing(slices)

    n_time = len(slices)
    logged = Signals(signals.time[:n_time], signals.meas[:n_time], signals.coil_signal[:n_time])
    capture = Capture(preshot, logged, slices, settings.json_text, source)
    if "capture" in settings:
        path = Path(settings["capture"]["path"]).expanduser()
        if not path.is_absolute():
            path = cache_directory() / path
        written = write_capture(path, capture, library_path=rtgsfit.library_path)
        if verbose:
            print(f"Wrote native capture to {written}")
    return capture.to_results()


def write_results(results: Results, settings: Settings) -> None:
    for name in settings.get("writers", {}).get("names", []):
        get_database_writer(name).write_results(results, settings)


def write_preshot(preshot: Preshot, settings: Settings) -> None:
    for name in settings.get("writers", {}).get("names", []):
        get_database_writer(name).write_preshot(preshot, settings)


def run_replay(settings: Settings, verbose: bool = True) -> Results:
    """Replay as described by `settings` and write the results with every `writers/names` writer."""
    results = replay(settings, verbose=verbose)
    write_results(results, settings)
    return results
