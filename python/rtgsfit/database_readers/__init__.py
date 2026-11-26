"""Database readers: where the preshot and the measured signals come from.

A *preshot reader* returns the `Preshot`; a *signals reader* returns the `Signals` passed to
`rtgsfit()` at each time slice. Readers are selected by name in the settings:

```json
{
    "preshot": {"reader": "analytic"},
    "signals": {"reader": "analytic"}
}
```

| reader            | preshot | signals | needs                                      |
|-------------------|---------|---------|--------------------------------------------|
| `analytic`        | yes     | yes     | nothing                                    |
| `local`           | yes     |         | a file written by the `local` writer       |
| `gsfit`           | yes     |         | GSFit (and MDSplus for ST40)               |
| `freegsnke`       | yes     | yes     | GSFit and FreeGSNKE (no TE network)        |
| `rtgsfit_mdsplus` | yes     |         | the TE network                             |
| `st40_magrun16`   |         | yes     | processed MAG RUN16 and PSU2COIL (TE network) |
| `st40_tepcs`      |         | yes     | TEPCS counts and per-pulse PSS (TE network) |

ST40 is replayed in one of two ways, each a packaged settings file named after its signals reader:
`st40_magrun16` feeds RT-GSFit what GSFit's post-shot reconstructions (GSFIT RUN04) fit, and
`st40_tepcs` reconstructs what the PCS sees in an experiment.

Readers are imported lazily, so that optional dependencies are only needed when used.
"""

from __future__ import annotations

import importlib
from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from dataclasses import field

import numpy as np
import numpy.typing as npt

from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

FloatArray = npt.NDArray[np.float64]


class InputDataError(ValueError):
    """The recorded data of a pulse cannot support a replay (e.g. a missing or inconsistent signal)."""


@dataclass
class Signals:
    """The inputs to `rtgsfit()` at every time slice, before the coil matrix is applied."""

    time: FloatArray  # [N_TIME]
    mag_signal: FloatArray  # [N_TIME, N_MAG_SIGNAL], ordered as PRESHOT.CONSTRAINT:MAG_SIGNALS
    coil_signal: FloatArray  # [N_TIME, N_COIL_SIGNAL], ordered as PRESHOT.COIL:COIL_SIGNALS
    input_quality: dict[str, str] = field(default_factory=dict)

    def validate(self, preshot: Preshot) -> None:
        n_time = len(self.time)
        n_mag_signal = preshot.sizes["N_MAG_SIGNAL"]
        n_coil_signal = len(preshot["COIL:COIL_SIGNALS"])
        if self.mag_signal.shape != (n_time, n_mag_signal):
            raise ValueError(f"mag_signal has shape {self.mag_signal.shape}, expected ({n_time}, {n_mag_signal})")
        if self.coil_signal.shape != (n_time, n_coil_signal):
            raise ValueError(f"coil_signal has shape {self.coil_signal.shape}, expected ({n_time}, {n_coil_signal})")
        if n_time > 1 and np.any(np.diff(self.time) < 0.0):
            raise ValueError("time must be non-decreasing")

    def repeat(self, n_repeat: int) -> Signals:
        """Every time slice `n_repeat` times: more iterations of `rtgsfit()` per time."""
        return Signals(
            time=np.repeat(self.time, n_repeat),
            mag_signal=np.repeat(self.mag_signal, n_repeat, axis=0),
            coil_signal=np.repeat(self.coil_signal, n_repeat, axis=0),
            input_quality=self.input_quality.copy(),
        )

    def coil_current(self, preshot: Preshot) -> FloatArray:
        """[N_TIME, N_COIL]: coil_current = coil_signal @ PRESHOT.COIL:MATRIX."""
        return self.coil_signal @ preshot["COIL:MATRIX"]


class PreshotReader(ABC):
    @abstractmethod
    def read_preshot(self, settings: Settings) -> Preshot:
        """Return the preshot described by `settings`."""


class SignalsReader(ABC):
    @abstractmethod
    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        """Return the signals of the replay described by `settings["replay"]`."""


_PRESHOT_READERS = {
    "analytic": "rtgsfit.database_readers.analytic:AnalyticReader",
    "local": "rtgsfit.database_readers.local:LocalPreshotReader",
    "gsfit": "rtgsfit.database_readers.gsfit:GsfitReader",
    "freegsnke": "rtgsfit.database_readers.freegsnke:FreegsnkeReader",
    "rtgsfit_mdsplus": "rtgsfit.database_readers.rtgsfit_mdsplus:RtgsfitMdsplusReader",
}
_SIGNALS_READERS = {
    "analytic": "rtgsfit.database_readers.analytic:AnalyticReader",
    "freegsnke": "rtgsfit.database_readers.freegsnke:FreegsnkeReader",
    "st40_magrun16": "rtgsfit.database_readers.st40_magrun16:St40Magrun16Reader",
    "st40_tepcs": "rtgsfit.database_readers.st40_tepcs:St40TepcsReader",
}


def _load(registry: dict[str, str], name: str, kind: str) -> object:
    if name not in registry:
        raise ValueError(f"unknown {kind} reader {name!r}; available: {', '.join(registry)}")
    module_name, class_name = registry[name].split(":")
    return getattr(importlib.import_module(module_name), class_name)()


def get_preshot_reader(name: str) -> PreshotReader:
    reader = _load(_PRESHOT_READERS, name, "preshot")
    assert isinstance(reader, PreshotReader)
    return reader


def get_signals_reader(name: str) -> SignalsReader:
    reader = _load(_SIGNALS_READERS, name, "signals")
    assert isinstance(reader, SignalsReader)
    return reader
