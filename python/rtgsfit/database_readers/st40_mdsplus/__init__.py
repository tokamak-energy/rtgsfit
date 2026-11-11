"""Read the ST40 measurements (MAG) and PF power supply currents (PSU) from MDSplus (TE network only).

Settings:

* `replay/pulse`: the ST40 pulse;
* `signals/st40_mdsplus/mag_run` (default `BEST`) and `server` (default `smaug`).

The signal names come from the preshot: `PRESHOT.MEAS:NAMES` are PCS names such as
`PSI_FLOOP_007`, `B_BPPROBE_101` or `I_ROG_INIVC000`, and `PRESHOT.COIL:SIGNALS` are PCS names such as
`I_SOL_PSU` or `I_ROG_MCWIRE`. Every signal is linearly interpolated onto the replay times. Replay
times outside the time range recorded by every signal are dropped, rather than given the first or last
recorded value, so a generous `replay/time_end` cannot replay made-up data.

All reads use `mdsthin`; neither the native MDSplus bindings nor `st40_database` are needed.
"""

from __future__ import annotations

import numpy as np

from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.preshot import Preshot
from rtgsfit.settings import Settings
from rtgsfit.settings import replay_times

DEFAULT_SERVER = "smaug"


def mag_node(pcs_name: str) -> str:
    """The MAG node of a PCS sensor name: `PSI_FLOOP_007` -> `FLOOP.L007.PSI`."""
    parts = pcs_name.split("_")
    if len(parts) != 3:
        raise ValueError(f"{pcs_name!r} is not a PCS sensor name of the form QUANTITY_SENSOR_ID")
    prefix = {"FLOOP": "L", "BPPROBE": "P", "ROG": ""}.get(parts[1])
    if prefix is None:
        raise ValueError(f"{pcs_name!r} is not a flux loop, Bp probe or Rogowski coil")
    return f"{parts[1]}.{prefix}{parts[2]}.{parts[0]}"


def psu_node(pcs_name: str) -> str:
    """The PSU node of a PCS coil signal name: `I_SOL_PSU` -> `\\ST40::TOP.PSU.SOL:I`."""
    parts = pcs_name.split("_")
    if len(parts) != 3 or parts[0] != "I" or parts[2] != "PSU":
        raise ValueError(f"{pcs_name!r} is not a PSU signal name of the form I_<PSU>_PSU")
    return f"\\ST40::TOP.PSU.{parts[1]}:I"


class St40MdsplusReader(SignalsReader):
    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        import mdsthin

        section = settings.get("signals", {}).get("st40_mdsplus", {})
        server = section.get("server", DEFAULT_SERVER)
        pulse = int(settings["replay"]["pulse"])
        time = replay_times(settings)

        with mdsthin.Connection(server) as connection:
            # MAG links refer to RAW_MAG in the ST40 tree. Its shared time base is the
            # dimension of BPPROBE.P101:B (the same reference as st40_database's MAG template).
            connection.openTree("ST40", pulse)
            root = f"\\MAG::TOP.{section.get('mag_run', 'BEST')}"
            mag_time = np.asarray(connection.get(f"dim_of({root}.BPPROBE.P101:B)").data(), dtype=np.float64)
            recorded = (time >= mag_time[0]) & (time <= mag_time[-1])

            def from_mag(pcs_name: str) -> np.ndarray:
                structure, _, quantity = mag_node(pcs_name).rpartition(".")
                path = f"{root}.{structure}:{quantity}"
                values = np.asarray(connection.get(path).data(), dtype=np.float64)
                return np.interp(time, mag_time, values)

            meas = np.column_stack([from_mag(name) for name in preshot["MEAS:NAMES"]])
            coil_signal = np.zeros((len(time), len(preshot["COIL:SIGNALS"])))
            for i, name in enumerate(preshot["COIL:SIGNALS"]):
                if "_ROG_" in name:
                    coil_signal[:, i] = from_mag(name)
                    continue
                path = psu_node(name)
                value = np.asarray(connection.get(path).data(), dtype=np.float64)
                value_time = np.asarray(connection.get(f"dim_of({path})").data(), dtype=np.float64)
                recorded &= (time >= value_time[0]) & (time <= value_time[-1])
                coil_signal[:, i] = np.interp(time, value_time, value)

        if not recorded.any():
            raise ValueError(
                f"pulse {pulse} has no replay time from {time[0]} s to {time[-1]} s recorded by every MAG "
                "and PSU signal"
            )
        return Signals(time=time[recorded], meas=meas[recorded], coil_signal=coil_signal[recorded])
