"""RT-GSFit replays in PPAC, ST40's post-pulse analysis chain: rebuild the preshot, replay the pulse, write MDSplus.

This replays a pulse after the event, so it can be run on old pulses. It is not the import of what the PCS
itself computed during the pulse. PPAC imports the package named in its `ppac_scripts.py` entry and runs

```python
controller = rtgsfit.RtgsfitReplay(
    pulseNo=13343, run_name="RUN01_POST", run_description="...", pulseNo_write=13343, link_run_to_best=True
)
controller.run()
```

`run()` is the whole pipeline for one pulse, with no state from earlier pulses:

1. build a new preshot with GSFit for this pulse (`preshot/reader = gsfit`, `gsfit/pulse`), so the
   machine description and the MAG `BEST` sensors are this pulse's;
2. read this pulse's MAG `BEST` and PSU signals by their PCS names (`st40_mdsplus`);
3. compile the preshot into a private copy of the C library (`isolated_build`) and call `rtgsfit()`
   every millisecond over the PCS's window, exactly as the PCS would;
4. write the results, which contain the preshot, settings and input signals, to the `RTGSFIT` tree.

Every other setting is that of the `st40` configuration. A setting can be changed before `run()`:

```python
controller = RtgsfitReplay(pulseNo=13343, run_name="RUN99", write_to_mds=False)
controller.settings = controller.settings.copy_with({"replay/time_end": 0.2})
results = controller.run()
```
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from rtgsfit.library import isolated_build
from rtgsfit.replay import replay
from rtgsfit.replay import write_results
from rtgsfit.results import Results
from rtgsfit.settings import Settings
from rtgsfit.settings import load_settings

POST_SHOT_SETTINGS: dict[str, Any] = {
    "preshot/reader": "gsfit",
    "signals/reader": "st40_mdsplus",
    "signals/st40_mdsplus/mag_run": "BEST",
    "replay/time_start": 0.0055,
    "replay/time_end": 0.6,
    "replay/time_step": 0.001,
    "replay/repeat": 1,
    "replay/stop_on_status": False,
    "replay/stop_below_ip": 0.0,
}
"""Changes to the base settings for every pulse.

The PCS runs `rtgsfit()` every millisecond from 5.5 ms to 600 ms, starting from the initial conditions
each time it is enabled (`RT_Diagnostics.c` enables it for PCS clock times 1.1055 s to 1.7 s, and the
PCS clock starts at pulse time -1.1 s), so the replay does the same. Do not start earlier: from 0 ms,
pulse 14922 locks onto a one-cell 'plasma' with status 0 and never recovers. Status 0 is not proof of a
plasma: on 14922 a 5 kA, 0.001 m^2 filament also has status 0 from 519.5 ms. A failed slice leaves
`psi_norm` and `mask` unchanged, and the `st40_mdsplus` reader drops times that a signal did not
record.
"""


class RtgsfitReplay:
    """Post-shot RT-GSFit replay of one ST40 pulse, with the constructor PPAC calls.

    The run is `RUN01_POST` of the `RTGSFIT` tree unless `run_name` says otherwise, and is linked to
    `BEST` by default. `pulseNo_write` defaults to `pulseNo`.
    `write_to_mds=False` runs everything but the MDSplus write.
    `settings` is the base configuration, a packaged name or a JSON path or `Settings`, and `overrides`
    are applied last, as slash-separated paths.
    """

    def __init__(
        self,
        pulseNo: int,
        run_name: str = "RUN01_POST",
        run_description: str = "RT-GSFit post-shot replay",
        pulseNo_write: int | None = None,
        link_run_to_best: bool = True,
        write_to_mds: bool = True,
        settings: str | Path | Settings = "st40",
        overrides: dict[str, Any] | None = None,
        verbose: bool = True,
    ) -> None:
        self.pulseNo = int(pulseNo)
        self.pulseNo_write = self.pulseNo if pulseNo_write is None else int(pulseNo_write)
        self.run_name = run_name
        self.run_description = run_description
        self.link_run_to_best = bool(link_run_to_best)
        self.write_to_mds = bool(write_to_mds)
        self.verbose = verbose
        base = settings if isinstance(settings, Settings) else load_settings(settings)
        pulse_settings = {
            **POST_SHOT_SETTINGS,
            "replay/pulse": self.pulseNo,
            "gsfit/pulse": self.pulseNo,
            "writers/names": ["rtgsfit_mdsplus"] if self.write_to_mds else [],
            "writers/rtgsfit_mdsplus": {
                "enabled": self.write_to_mds,
                "pulse": self.pulseNo_write,
                "run_name": run_name,
                "run_description": run_description,
                "link_best": self.link_run_to_best,
            },
        }
        self.settings = base.copy_with({**pulse_settings, **(overrides or {})})
        self.results: Results | None = None

    def run(self) -> Results:
        """Build the preshot, replay the pulse and write the results; return them."""
        if self.verbose:
            print(f"RT-GSFit post-shot replay of pulse {self.pulseNo}")
        with isolated_build():
            results = replay(self.settings, verbose=self.verbose)
        write_results(results, self.settings)
        self.results = results
        if self.verbose:
            status = np.asarray(results["GLOBAL.STATUS"])
            time = np.asarray(results["TIME"])
            print(
                f"Pulse {self.pulseNo}: {np.count_nonzero(status == 0)} of {len(status)} slices from "
                f"{time[0]:.4f} s to {time[-1]:.4f} s have status 0"
            )
        return results
