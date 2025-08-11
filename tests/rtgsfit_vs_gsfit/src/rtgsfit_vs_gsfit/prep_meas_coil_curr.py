"""
Contains `prep_meas()` and `prep_coil_curr()` which read measurement data from MAG
and PSU2COIL on MDS+ to supply inputs for `replay_rtgsfit()`.
"""

import numpy as np
import mdsthin

import st40_database

def convert_sens_name_to_mag(sens_name: str) -> str:
    """
    Convert a PCS sensor name to the format used on MDS+
    for MAG.
    """
    parts = sens_name.split("_")
    if "FLOOP" in sens_name:
        append_character = "L"
    elif "BPPROBE" in sens_name:
        append_character = "P"
    elif "ROG" in sens_name:
        append_character = ""
    else:
        raise ValueError(f"Sensor name {sens_name} does not match expected patterns.")
    return f"{parts[1]}.{append_character}{parts[2]}.{parts[0]}"

def convert_coil_name_to_psu2coil(coil_name: str) -> str:
    """
    Convert a coil name from the COIL_NAMES list to the format PSU2COIL expects.
    """
    return f"PF.{coil_name}.I"

def convert_coil_signal_to_psu(coil_signal: str) -> str:
    """
    Convert a coil signal name to the PSU format expected on MDS+.
    """
    if "ROG" in coil_signal:
        return convert_sens_name_to_mag(coil_signal)
    else:
        parts = coil_signal.split("_")
        return parts[1]


def prep_meas_pcs(cfg: dict) -> np.ndarray:
    """
    Prepare the meas_pcs input array for RTGSFIT by reading data from MAG on MDS+.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:SENS_NAMES").data()

    mag = st40_database.GetData(cfg["pulse_num"], "MAG#BEST")
    time_array_mag = mag.get("TIME")

    sens_names_mag = [convert_sens_name_to_mag(sens_name) for sens_name in sens_names]
    meas_pcs = np.zeros(len(sens_names), dtype=np.float64)
    for i, sens_name_mag in enumerate(sens_names_mag):
        meas_mag = mag.get(sens_name_mag)
        meas_pcs[i] = np.interp(cfg["time"], time_array_mag, meas_mag)

    return meas_pcs

def prep_coil_curr(cfg: dict) -> np.ndarray:

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        coil_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:COIL_NAMES").data()
        coil_matrix = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:COIL_MATRIX").data()
        coil_signals = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:COIL_SIGNALS").data()

    mag = st40_database.GetData(cfg["pulse_num"], "MAG#BEST")
    time_array_mag = mag.get("TIME")

    if cfg["use_psu2coil"]:
        coil_names_psu2coil = [convert_coil_name_to_psu2coil(coil_name) for coil_name in coil_names]
        psu2coil = st40_database.GetData(cfg["pulse_num"], f'PSU2COIL#{cfg["psu2coil_run_name"]}')
        coil_curr = np.zeros(len(coil_names), dtype=np.float64)
        for i, coil_name_psu2coil in enumerate(coil_names_psu2coil):
            coil_curr_psu2coil = psu2coil.get(coil_name_psu2coil)
            coil_curr[i] = np.interp(cfg["time"], time_array_mag, coil_curr_psu2coil)
    else:
        # Use coil_matrix to compute coil currents from the PSU signals
        psu_currents = np.zeros(len(coil_signals), dtype=np.float64)
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("ST40", cfg["pulse_num"])
            for i, coil_signal in enumerate(coil_signals):
                coil_signal_psu = convert_coil_signal_to_psu(coil_signal)
                if "ROG" not in coil_signal:
                    psu_current = conn.get(f"\\ST40::TOP.PSU.{coil_signal_psu}:I").data()
                else:
                    # One of the signals is a Rogowski coil which we read from MAG instead of PSU
                    psu_current = mag.get(coil_signal_psu)
                psu_currents[i] = np.interp(cfg["time"], time_array_mag, psu_current)
            coil_curr = psu_currents @ coil_matrix.T

    return coil_curr

if __name__ == "__main__":
    from rtgsfit_vs_gsfit import config_loader

    cfg = config_loader.load_and_prepare_config()
    meas_pcs = prep_meas_pcs(cfg)
    coil_curr_with_psu2coil = prep_coil_curr(cfg)
    # Allow more columns to be printed for better visibility
    np.set_printoptions(threshold=np.inf, linewidth=200)
    print("Prepared meas_pcs:")
    print(meas_pcs)
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        coil_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:COIL_NAMES").data()
    print("Coil names:", coil_names)
    print("Prepared coil_curr with PSU2COIL:", coil_curr_with_psu2coil)

    cfg = config_loader.load_and_prepare_config()
    cfg["use_psu2coil"] = False
    _ = prep_meas_pcs(cfg)
    coil_curr_no_psu2coil = prep_coil_curr(cfg)
    print("Prepared coil_curr wout PSU2COIL:", coil_curr_no_psu2coil)
