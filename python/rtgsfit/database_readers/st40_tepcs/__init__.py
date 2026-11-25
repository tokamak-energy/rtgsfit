"""The ST40 PCS inputs, reconstructed from the TEPCS counts and the pulse's PSS calibration (TE network).

This is what the PCS itself does in an experiment: the counts of every analogue channel are calibrated
with `PSS.AI`, the magnetic sensors have their offset and linear drift subtracted and the PSS
compensation matrix applied, as in `RT_Diagnostics.c`. With the PCS preshot (PSU coil signals and its
sensor replacements) this is the closest replay to the real-time reconstruction.

Like the PCS, every cycle reads each channel's last sample before the cycle time, counting samples from
the start of the record. The stored times are not used to choose it because they can drift from the
sample count: by 0.6 s in pulse 14211, ACQ2106_033's stored times are 34 µs behind it.

Settings: `signals/st40_tepcs/server` (default `smaug`). Replay times must lie on the PCS cycle grid.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from typing import Any

import numpy as np

from rtgsfit.database_readers import InputDataError
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings
from rtgsfit.settings import replay_times

if TYPE_CHECKING:
    import mdsthin

LOG = logging.getLogger(__name__)
PCS_DT = 5e-5
PCS_T0 = -1.1
PCS_AVERAGE = 600
PCS_FIRST = 600
PCS_SECOND = 10000


def channel(index: int) -> str:
    if not 0 <= index < 288:
        raise InputDataError(f"invalid PSS analogue channel index {index}")
    if index < 160:
        return f"\\TEPCS::TOP.ACQ2106_032:CH{index + 1:03d}"
    return f"\\TEPCS::TOP.ACQ2106_033:CH{index - 159:03d}"


def pcs_signals(connection: mdsthin.Connection, settings: Settings, preshot: Preshot) -> Signals:
    """Reconstruct inputs using this pulse's PSS and the current PCS correction algorithm."""
    from mdsthin.exceptions import MdsException
    from mdsthin.exceptions import TreeNNF

    pulse = int(settings["replay"]["pulse"])
    time = replay_times(settings)
    cycles_float = (time - PCS_T0) / PCS_DT
    cycles = np.rint(cycles_float).astype(int)
    if not len(cycles) or not np.allclose(cycles_float, cycles, rtol=0, atol=1e-6) or np.any(cycles < 1):
        raise InputDataError("replay times do not lie on the PCS cycle grid after its first cycle")

    try:
        connection.openTree("PSS", pulse)
        try:

            def get(name: str) -> Any:
                return connection.get(f"\\PSS::TOP.{name}").data()

            names = [str(name).strip() for name in np.ravel(get("AI.LIST"))]
            offsets = np.ravel(get("AI.OFFSETCOUNTS")).astype(float)
            gains = np.ravel(get("AI.COUNT2MKS")).astype(float)
            if not np.isclose(float(get("TIMES.T_OFFSET")), PCS_T0, rtol=0, atol=1e-7):
                raise InputDataError("PSS clock offset differs from the reconstruction")
            if not np.isclose(float(get("TIMES.DT_CYCLE")), PCS_DT, rtol=0, atol=1e-9):
                raise InputDataError("PSS cycle time differs from the reconstruction")
            nonempty_names = [name for name in names if name]
            if (
                len(set(nonempty_names)) != len(nonempty_names)
                or len(offsets) != len(names)
                or len(gains) != len(names)
                or not np.isfinite(offsets).all()
                or not np.isfinite(gains).all()
            ):
                raise InputDataError("PSS AI names and calibration arrays are inconsistent or non-finite")
            optional = {}
            missing = []
            for key in ("MAG.SUB_NAMES", "MAG.COMP_MATRIX", "MAG.INDEX_COMP"):
                try:
                    optional[key] = get(key)
                except TreeNNF:
                    missing.append(key)
            if len(missing) == 3:
                LOG.warning("Pulse %s has no PSS MAG metadata; using AI calibration only", pulse)
                subset, comp_index, matrix = [], np.array([], dtype=int), np.empty((0, 0))
                correction = "unknown: PSS MAG metadata absent; calibration only"
            elif missing:
                raise InputDataError(f"pulse {pulse} has incomplete PSS MAG metadata: {missing}")
            else:
                subset = [str(name).strip() for name in np.ravel(optional["MAG.SUB_NAMES"])]
                comp_index = np.ravel(optional["MAG.INDEX_COMP"]).astype(int) - 1
                matrix = np.asarray(optional["MAG.COMP_MATRIX"], dtype=float)
                if (
                    len(set(subset)) != len(subset)
                    or matrix.shape != (len(comp_index), len(subset))
                    or np.any(comp_index < 0)
                    or np.any(comp_index >= len(names))
                    or not np.isfinite(matrix).all()
                ):
                    raise InputDataError("PSS MAG metadata has inconsistent names, indices or shape")
                correction = "PSS MAG compensation and current PCS offset/drift algorithm"
        finally:
            connection.closeTree("PSS", pulse)

        requested = [str(name) for name in preshot["MEAS:NAMES"]] + [str(name) for name in preshot["COIL:SIGNALS"]]
        unknown = sorted(set(requested) - set(names))
        if unknown:
            raise InputDataError(f"pulse {pulse} PSS AI.LIST lacks required signals: {unknown}")
        connection.openTree("TEPCS", pulse)
        try:
            cache = {}

            def mks(index: int) -> np.ndarray:
                if index not in cache:
                    path = channel(index)
                    try:
                        counts = np.asarray(connection.get(f"raw_of({path})").data(), dtype=float)
                        times = np.asarray(connection.get(f"dim_of({path})").data(), dtype=float)
                    except MdsException as exc:
                        raise InputDataError(
                            f"pulse {pulse} TEPCS signal {names[index]} ({path}) read failed: "
                            f"{type(exc).__name__}: {exc}"
                        ) from exc
                    if counts.ndim != 1 or times.ndim != 1 or len(counts) != len(times) or len(times) < 2:
                        raise InputDataError(f"{path} has inconsistent counts and time")
                    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
                        raise InputDataError(f"{path} has an invalid time base")
                    per_cycle = round(PCS_DT * (len(times) - 1) / (times[-1] - times[0]))
                    if per_cycle < 1 or abs(times[0] - PCS_T0) > PCS_DT / per_cycle / 2:
                        raise InputDataError(f"{path} does not start on the PCS cycle grid")
                    if np.diff(times).max() > 1.5 * PCS_DT / per_cycle:
                        raise InputDataError(f"{path} has gaps in the recorded time base")
                    n_cycle = max(int(cycles.max()), PCS_SECOND) + 1
                    if len(counts) < per_cycle * (n_cycle - 1):
                        raise InputDataError(f"{path} is too short for the requested cycles and drift calibration")
                    # Cycle 0 has no earlier sample.
                    samples = per_cycle * np.arange(1, n_cycle) - 1
                    cache[index] = np.concatenate([[np.nan], (counts[samples] + offsets[index]) * gains[index]])
                return cache[index]

            output = {}
            for name in dict.fromkeys(requested):
                x = mks(names.index(name))
                y = x[cycles].copy()
                if name in subset:
                    k = subset.index(name)
                    avg1 = x[PCS_FIRST - PCS_AVERAGE + 1 : PCS_FIRST + 1].mean()
                    avg2 = x[PCS_SECOND - PCS_AVERAGE + 1 : PCS_SECOND + 1].mean()
                    on = cycles >= PCS_SECOND
                    y[on] -= avg2 + (avg2 - avg1) / 0.47 * PCS_DT * (cycles[on] - PCS_SECOND)
                    for j, index in enumerate(comp_index):
                        y += matrix[j, k] * mks(int(index))[cycles]
                if not np.isfinite(y).all():
                    raise InputDataError(f"pulse {pulse} PCS input {name} contains non-finite values")
                output[name] = y
        finally:
            connection.closeTree("TEPCS", pulse)
    except InputDataError:
        raise
    except (MdsException, ValueError) as exc:
        raise InputDataError(f"pulse {pulse} PSS/TEPCS read failed: {type(exc).__name__}: {exc}") from exc

    n_meas = len(preshot["MEAS:NAMES"])
    return Signals(
        time=time,
        meas=np.column_stack([output[name] for name in requested[:n_meas]]),
        coil_signal=np.column_stack([output[name] for name in requested[n_meas:]]),
        input_quality={"source": "TEPCS counts and per-pulse PSS calibration", "correction": correction},
    )


class St40TepcsReader(SignalsReader):
    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        import mdsthin

        section = settings.get("signals", {}).get("st40_tepcs", {})
        with mdsthin.Connection(section.get("server", "smaug"), timeout=30) as connection:
            return pcs_signals(connection, settings, preshot)
