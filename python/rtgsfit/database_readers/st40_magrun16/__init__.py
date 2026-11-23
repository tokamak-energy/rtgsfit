"""The ST40 measurements that GSFit's post-shot reconstructions (GSFIT RUN04) fit (TE network).

GSFit's `st40_mdsplus_with_digital_filter` reader fits the processed magnetic signals of MAG RUN16
(offset and drift corrected and digitally filtered) and the PF coil currents of PSU2COIL. This reader
gives RT-GSFit the same inputs, so that the two codes differ only in their algorithms:

* `PRESHOT.MEAS:NAMES` are PCS sensor names: `PSI_FLOOP_007` is `\\MAG::TOP.RUN16.FLOOP.L007:PSI_PROCESS`,
  `B_BPPROBE_101` is `BPPROBE.P101:B_PROCESS` and `I_ROG_INIVC000` is `ROG.INIVC000:I_PROCESS`;
* `PRESHOT.COIL:SIGNALS` are `I_<COIL>_PSU2COIL`: the current of the PSU2COIL channel that carries every
  PF coil of `gsfit/coils/groups/<COIL>`.

Settings: `signals/st40_magrun16/psu2coil_run` (default `BEST`, as GSFIT RUN04 used) and `server`
(default `smaug`). The signals are always those of `replay/pulse`.

Every signal is linearly interpolated onto the replay times. Replay times outside the time recorded by
every signal are dropped, so a generous `replay/time_end` cannot replay made-up data.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from rtgsfit.database_readers import InputDataError
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings
from rtgsfit.settings import replay_times

if TYPE_CHECKING:
    import mdsthin

MAG_ROOT = "\\MAG::TOP.RUN16"
PSU2COIL_SUFFIX = "_PSU2COIL"


def mag_node(pcs_name: str) -> str:
    """The processed MAG RUN16 node of a PCS sensor name: `PSI_FLOOP_007` -> `FLOOP.L007:PSI_PROCESS`."""
    parts = pcs_name.split("_")
    prefix = {"FLOOP": "L", "BPPROBE": "P", "ROG": ""}.get(parts[1]) if len(parts) == 3 else None
    if prefix is None:
        raise ValueError(f"{pcs_name!r} is not a flux loop, Bp probe or Rogowski coil name such as PSI_FLOOP_007")
    # MDSplus node names have at most 12 characters: GSFit's PSI_PROCESSED is the node PSI_PROCESS.
    return f"{parts[1]}.{prefix}{parts[2]}:{parts[0]}_PROCESS"


def psu2coil_channels(signals: list[str], groups: dict[str, list[str]], channel_coils: np.ndarray) -> list[int]:
    """The PSU2COIL channel (row of `PF.ALL:COILS`) of every `I_<COIL>_PSU2COIL` signal.

    Every PF coil of a group must be in the same channel, and every PSU2COIL coil must be in a group,
    so that no current is lost. Several groups can share a channel (MCT and MCB before pulse 12551).
    """
    channel_of = {}
    for channel, row in enumerate(np.asarray(channel_coils).reshape(len(channel_coils), -1)):
        for coil in (str(name).strip() for name in row):
            if coil:
                channel_of[coil] = channel
    unmodelled = sorted(set(channel_of) - {coil for members in groups.values() for coil in members})
    if unmodelled:
        raise InputDataError(f"PSU2COIL drives PF coils that gsfit/coils/groups does not have: {unmodelled}")
    channels = []
    for signal in signals:
        group = signal.removeprefix("I_").removesuffix(PSU2COIL_SUFFIX)
        if not signal.startswith("I_") or not signal.endswith(PSU2COIL_SUFFIX) or group not in groups:
            raise ValueError(f"{signal!r} is not I_<COIL>{PSU2COIL_SUFFIX} for a coil of gsfit/coils/groups")
        missing = [coil for coil in groups[group] if coil not in channel_of]
        if missing:
            raise InputDataError(f"PSU2COIL has no current for the PF coils {missing} of {group}")
        group_channels = {channel_of[coil] for coil in groups[group]}
        if len(group_channels) != 1:
            raise InputDataError(f"the PF coils of {group} are in different PSU2COIL channels")
        channels.append(group_channels.pop())
    return channels


def _get(connection: mdsthin.Connection, path: str, pulse: int) -> np.ndarray:
    from mdsthin.exceptions import MdsException

    try:
        return np.asarray(connection.get(path).data())
    except MdsException as exc:
        raise InputDataError(f"pulse {pulse} {path}: {type(exc).__name__}") from exc


def _check_time(time: np.ndarray, n_samples: int, path: str) -> None:
    if time.ndim != 1 or len(time) < 2 or len(time) != n_samples:
        raise InputDataError(f"{path} has inconsistent values and time")
    if not np.isfinite(time).all() or np.any(np.diff(time) <= 0.0):
        raise InputDataError(f"{path} has an invalid time base")


def magrun16_signals(connection: mdsthin.Connection, settings: Settings, preshot: Preshot) -> Signals:
    pulse = int(settings["replay"]["pulse"])
    psu2coil_run = str(settings.get("signals", {}).get("st40_magrun16", {}).get("psu2coil_run", "BEST"))
    time = replay_times(settings)
    recorded = np.ones(len(time), dtype=bool)

    def interpolate(signal_time: np.ndarray, values: np.ndarray) -> np.ndarray:
        nonlocal recorded
        recorded &= (time >= signal_time[0]) & (time <= signal_time[-1])
        return np.interp(time, signal_time, values)

    connection.openTree("ST40", pulse)
    try:
        meas = []
        for name in preshot["MEAS:NAMES"]:
            path = f"{MAG_ROOT}.{mag_node(str(name))}"
            values = _get(connection, f"data({path})", pulse).astype(np.float64)
            signal_time = _get(connection, f"dim_of({path})", pulse).astype(np.float64)
            _check_time(signal_time, len(values), path)
            meas.append(interpolate(signal_time, values))
    finally:
        connection.closeTree("ST40", pulse)

    root = f"\\PSU2COIL::TOP.{psu2coil_run}"
    connection.openTree("PSU2COIL", pulse)
    try:
        channel_coils = _get(connection, f"{root}.PF.ALL:COILS", pulse)
        channels = psu2coil_channels(
            [str(name) for name in preshot["COIL:SIGNALS"]], settings["gsfit"]["coils"]["groups"], channel_coils
        )
        coil_time = _get(connection, f"{root}:TIME", pulse).astype(np.float64)
        currents = _get(connection, f"{root}.PF.ALL:I", pulse).astype(np.float64)
        if currents.shape != (len(coil_time), len(channel_coils)):
            raise InputDataError(f"pulse {pulse} {root}.PF.ALL:I does not match its TIME and PF.ALL:COILS")
        _check_time(coil_time, len(currents), f"{root}:TIME")
        coil_signal = [interpolate(coil_time, currents[:, channel]) for channel in channels]
        if psu2coil_run == "BEST":
            psu2coil_run = f"BEST ({str(_get(connection, f'{root}:BEST_RUN', pulse)).strip()})"
    finally:
        connection.closeTree("PSU2COIL", pulse)

    if not recorded.any():
        raise InputDataError(f"pulse {pulse} has no replay times recorded by every MAG RUN16 and PSU2COIL signal")
    signals = Signals(
        time=time[recorded],
        meas=np.column_stack(meas)[recorded],
        coil_signal=np.column_stack(coil_signal)[recorded],
        input_quality={"source": f"MAG RUN16 processed signals and PSU2COIL {psu2coil_run}"},
    )
    if not np.isfinite(signals.meas).all() or not np.isfinite(signals.coil_signal).all():
        raise InputDataError(f"pulse {pulse} has non-finite MAG RUN16 or PSU2COIL values at the replay times")
    return signals


class St40Magrun16Reader(SignalsReader):
    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        import mdsthin

        server = settings.get("signals", {}).get("st40_magrun16", {}).get("server", "smaug")
        with mdsthin.Connection(server) as connection:
            return magrun16_signals(connection, settings, preshot)
