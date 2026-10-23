"""Database readers: where the preshot and the measured signals come from.

A *preshot reader* returns the `Preshot`; a *signals reader* returns the `Signals` passed to
`rtgsfit()` at each time slice. Readers are selected by name in the settings:

```xml
<preshot><reader type="str">analytic</reader></preshot>
<signals><reader type="str">analytic</reader></signals>
```

| reader            | preshot | signals | needs                                      |
|-------------------|---------|---------|--------------------------------------------|
| `analytic`        | yes     | yes     | nothing                                    |
| `local`           | yes     |         | a file written by the `local` writer       |
| `gsfit`           | yes     | yes     | GSFit (and MDSplus for ST40)               |
| `freegsnke`       | yes     | yes     | GSFit and FreeGSNKE (no TE network)        |
| `rtgsfit_mdsplus` | yes     |         | the TE network                             |
| `st40_mdsplus`    |         | yes     | the TE network                             |

Readers are imported lazily, so that optional dependencies are only needed when used.
"""

from __future__ import annotations

import importlib
from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings

FloatArray = npt.NDArray[np.float64]


@dataclass
class Signals:
    """The inputs to `rtgsfit()` at every time slice, before the coil matrix is applied."""

    time: FloatArray  # [N_TIME]
    meas: FloatArray  # [N_TIME, N_MEAS], ordered as PRESHOT.MEAS:NAMES
    coil_signal: FloatArray  # [N_TIME, N_COIL_SIGNAL], ordered as PRESHOT.COIL:SIGNALS

    def validate(self, preshot: Preshot) -> None:
        n_time = len(self.time)
        n_meas = preshot.sizes["N_MEAS"]
        n_coil_signal = len(preshot["COIL:SIGNALS"])
        if self.meas.shape != (n_time, n_meas):
            raise ValueError(f"meas has shape {self.meas.shape}, expected ({n_time}, {n_meas})")
        if self.coil_signal.shape != (n_time, n_coil_signal):
            raise ValueError(f"coil_signal has shape {self.coil_signal.shape}, expected ({n_time}, {n_coil_signal})")
        if n_time > 1 and np.any(np.diff(self.time) < 0.0):
            raise ValueError("time must be non-decreasing")

    def repeat(self, n_repeat: int) -> Signals:
        """Every time slice `n_repeat` times: more iterations of `rtgsfit()` per time."""
        return Signals(
            time=np.repeat(self.time, n_repeat),
            meas=np.repeat(self.meas, n_repeat, axis=0),
            coil_signal=np.repeat(self.coil_signal, n_repeat, axis=0),
        )

    def coil_current(self, preshot: Preshot) -> FloatArray:
        """[N_TIME, N_COIL]: coil_current = PRESHOT.COIL:MATRIX coil_signal."""
        return self.coil_signal @ preshot["COIL:MATRIX"].T


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
    "gsfit": "rtgsfit.database_readers.gsfit:GsfitReader",
    "freegsnke": "rtgsfit.database_readers.freegsnke:FreegsnkeReader",
    "st40_mdsplus": "rtgsfit.database_readers.st40_mdsplus:St40MdsplusReader",
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
