import numpy as np
import mdsthin

import st40_database

from rtgsfit_vs_gsfit import cnst

mag = st40_database.GetData(cnst.PULSE_NUM, "MAG#BEST")
TIME_EXPERIMENTAL = mag.get('TIME')

def convert_fl_list_format(inpsi_list):
    outpsi_list = []
    for item in inpsi_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"FLOOP.L{number_part}.PSI"
        outpsi_list.append(new_format)
    return outpsi_list

def convert_bp_list_format(inbp_list):
    outbp_list = []
    for item in inbp_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"BPPROBE.P{number_part}.B"
        outbp_list.append(new_format)
    return outbp_list

def convert_rog_list_format(inrog_list):
    outrog_list = []
    for name in inrog_list:
        name = name.strip()
        if name.startswith("I_ROG_"):
            body = name[len("I_ROG_"):]
            outrog_list.append(f"ROG.{body}.I")
        else:
            raise ValueError(f"Unexpected format: {name}")
    return outrog_list

def convert_coil_list_format(incoil_list):
    return [f"PF.{name.strip()}.I" for name in incoil_list]

def prep_sens_meas(time: float,
                   sens_names: list) -> np.ndarray:
    """
    Prepare the sensor measurements for RTGSFIT.
    """

    sens_meas = np.zeros(len(sens_names), dtype=np.float64)
    for i, sens_name in enumerate(sens_names):
        sens_experiment = mag.get(sens_name)
        sens_interp = np.interp(time, TIME_EXPERIMENTAL, sens_experiment)
        sens_meas[i] = sens_interp

    return sens_meas

def prep_meas(time: float) -> np.ndarray:
    """
    Prepare the measurements for RTGSFIT.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        n_meas = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_MEAS").data()
        n_flux_loops = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_F_LOOPS").data()
        n_bp_probes = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_BP_PROBES").data()
        n_rogowski_coils = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_ROG_COILS").data()
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:SENS_NAMES").data()

    fl_range = range(n_flux_loops)
    bp_range = range(n_flux_loops, n_flux_loops + n_bp_probes)
    rogowski_range = range(n_flux_loops + n_bp_probes,
                           n_flux_loops + n_bp_probes + n_rogowski_coils)
    fl_names = convert_fl_list_format(sens_names[fl_range])
    bp_probe_names = convert_bp_list_format(sens_names[bp_range])
    rogowski_coil_names = convert_rog_list_format(sens_names[rogowski_range])

    meas = np.zeros(n_meas, dtype=np.float64)
    meas[fl_range] = prep_sens_meas(time, fl_names)
    meas[bp_range] = prep_sens_meas(time, bp_probe_names)
    meas[rogowski_range] = prep_sens_meas(time, rogowski_coil_names)

    print("Flux loop measurements:", meas[fl_range])
    print("BP probe measurements:", meas[bp_range])
    print("Rogowski coil measurements:", meas[rogowski_range])
    return meas

def prep_coil_curr(time: float) -> np.ndarray:
    """
    Prepare the coil currents for RTGSFIT.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        coil_names = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:COIL_NAMES").data()

    coil_names = convert_coil_list_format(coil_names)
    psu2coil = st40_database.GetData(cnst.PULSE_NUM, f'PSU2COIL#{cnst.PSU2COIL_RUN_NAME}')

    coil_curr = np.zeros(len(coil_names), dtype=np.float64)
    for i, coil_name in enumerate(coil_names):
        coil_experiment = psu2coil.get(coil_name)
        coil_interp = np.interp(time, TIME_EXPERIMENTAL, coil_experiment)
        coil_curr[i] = coil_interp

    print("Coil names:", coil_names)
    print("Coil currents:", coil_curr)
    return coil_curr