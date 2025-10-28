import os
import sys

import numpy as np
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
import mdsthin
from netCDF4 import Dataset
import numpy as np
import os


def write_data_to_mdsplus(
    pulseNo: int,
    run_name: str = "RUN01",
    run_description: str = "Standard run with default settings",
    settings_path: str = "default",
    write_to_mds: bool = True,
    pulseNo_write: int | None = None,
    pulse_num_preshot: int = 99_000_230,
    run_name_preshot: str = "RUN08",
) -> None:
    """
    Write RT-GSFit results to MDSplus

    :param pulseNo: pulse number
    :param run_name: run_name to save to MDSplus
    :param run_description: help string for MDSplus Tree
    :param settings_path: location where code inputs are stored
    :param write_to_mds: flag to turn on / off writing to MDSplus
    :param pulseNo_write: pulse number in which data is to be written, if different from the current pulse

    :return: None
    """
    data_file_name = f"/home/pcs.user/ops_double_threaded/st40pcs_dtacq/results/rtgsfit_results_{pulseNo}.nc"
    data_file_name = f"/home/filip.janky/ops/pcs/model/ST40PCS/results/rtgsfit_results_{pulseNo}.nc"
    data_file_name = f"/home/alex.prokopyszyn/Data/filip_data/rtgsfit_results_{pulseNo}.nc"

    # If `pulseNo_write` is not specified, we will write to `pulseNo`
    if pulseNo_write is None:
        pulseNo_write = pulseNo

    # Load data from netCDF file
    with Dataset(data_file_name, "r") as nc:
        print("Variables in netCDF file new:")
        for var in nc.variables:
            print(var)
        time = np.array(nc.variables["time"][:])
        flux_norm = np.array(nc.variables["flux_norm"][:])
        mask = np.array(nc.variables["mask"][:])
        chi_sq_err = np.array(nc.variables["chi_sq_err"][:])
        lcfs_n = np.array(nc.variables["lcfs_n"][:])
        flux_total = np.array(nc.variables["flux_total"][:])
        lcfs_r = np.array(nc.variables["lcfs_r"][:])
        lcfs_z = np.array(nc.variables["lcfs_z"][:])
        coef = np.array(nc.variables["coef"][:])
        flux_boundary = np.array(nc.variables["flux_boundary"][:])
        plasma_current = np.array(nc.variables["plasma_current"][:])
        lcfs_err_code = np.array(nc.variables["lcfs_err_code"][:])
        lapack_dgelss_info = np.array(nc.variables["lapack_dgelss_info"][:])
        meas_model = np.array(nc.variables["meas_model"][:])
        sensors_IN = np.array(nc.variables["sensors_IN"][:])
        I_PF_IN = np.array(nc.variables["I_PF_IN"][:])
        r = np.array(nc.variables["r"][:])
        z = np.array(nc.variables["z"][:])
        
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", pulse_num_preshot)
        coil_names = conn.get(f"\\RTGSFIT::TOP.{run_name_preshot}.PRESHOT:COIL_NAMES").data()
        meas_names = conn.get(f"\\RTGSFIT::TOP.{run_name_preshot}.PRESHOT:MEAS_NAMES").data()
        sens_rep_mat = conn.get(f"\\RTGSFIT::TOP.{run_name_preshot}.PRESHOT:SENS_REP_MAT").data()
        n_sens = int(np.sqrt(len(sens_rep_mat)))
        sens_rep_mat = sens_rep_mat.reshape((n_sens, n_sens))
        coef_names = conn.get(f"\\RTGSFIT::TOP.{run_name_preshot}.PRESHOT:COEF_NAMES").data()
        weight = conn.get(f"\\RTGSFIT::TOP.{run_name_preshot}.PRESHOT:WEIGHT").data()
    
    flux_loop_indices = []
    bp_probe_indices = []
    rog_coil_indices = []
    for i, meas_name in enumerate(meas_names):
        if meas_name.startswith("L") and meas_name[1:].isdigit():
            flux_loop_indices.append(i)
        elif meas_name.startswith("P") and meas_name[1:].isdigit():
            bp_probe_indices.append(i)
        elif meas_name.startswith("IVC"):
            continue
        elif meas_name.startswith("OVC"):
            continue
        else:
            rog_coil_indices.append(i)
    ivc_indices = []
    case_indices = []
    psr_indices = []
    for i, coef_name in enumerate(coef_names):
        if coef_name.startswith("eig_"):
            ivc_indices.append(i)
        elif coef_name == "OVC":
            ovc_index = i
        elif "CASE" in coef_name:
            case_indices.append(i)
        elif "PSR" in coef_name:
            psr_indices.append(i)

    # Create a nested dictionary to store data
    results = NestedDict()

    ## Store results
    results["TIME"] = time
    results["TWO_D"]["MASK"] = mask
    results["GLOBAL"]["CHIT"] = chi_sq_err
    results["P_BOUNDARY"]["NBND"] = lcfs_n
    results["TWO_D"]["PSI"] = flux_total
    results["P_BOUNDARY"]["RBND"] = lcfs_r
    results["P_BOUNDARY"]["ZBND"] = lcfs_z
    results["GLOBAL"]["PSI_B"] = flux_boundary
    results["GLOBAL"]["IP"] = plasma_current
    results["GLOBAL"]["LCFS_ERR"] = lcfs_err_code
    results["GLOBAL"]["DGELSS_INFO"] = lapack_dgelss_info
    results["TWO_D"]["RGRID"] = r
    results["TWO_D"]["ZGRID"] = z
        
    # Constraint and coil names
    results["CONSTRAINTS"]["COIL"]["NAME"] = coil_names
    results["CONSTRAINTS"]["FLOOP"]["NAME"] = meas_names[flux_loop_indices]
    results["CONSTRAINTS"]["BPPROBE"]["NAME"] = meas_names[bp_probe_indices]
    results["CONSTRAINTS"]["ROGOWSKI"]["NAME"] = meas_names[rog_coil_indices]

    # Constraint and coil mvalues and weights
    meas_no_reg = sensors_IN @ sens_rep_mat.T
    results["CONSTRAINTS"]["FLOOP"]["MVALUE"] = meas_no_reg[:, flux_loop_indices]
    results["CONSTRAINTS"]["FLOOP"]["WEIGHT"] = weight[flux_loop_indices]
    results["CONSTRAINTS"]["BPPROBE"]["MVALUE"] = meas_no_reg[:, bp_probe_indices]
    results["CONSTRAINTS"]["BPPROBE"]["WEIGHT"] = weight[bp_probe_indices]
    results["CONSTRAINTS"]["ROGOWSKI"]["MVALUE"] = meas_no_reg[:, rog_coil_indices]
    results["CONSTRAINTS"]["ROGOWSKI"]["WEIGHT"] = weight[rog_coil_indices]
    results["CONSTRAINTS"]["COIL"]["MVALUE"] = I_PF_IN

    # Constraints CVALUEs
    weight_matrix = np.tile(weight, (len(time), 1))
    meas_model_no_weight = meas_model / weight_matrix
    results["CONSTRAINTS"]["FLOOP"]["CVALUE"] = meas_model_no_weight[:, flux_loop_indices]
    results["CONSTRAINTS"]["BPPROBE"]["CVALUE"] = meas_model_no_weight[:, bp_probe_indices]
    results["CONSTRAINTS"]["ROGOWSKI"]["CVALUE"] = meas_model_no_weight[:, rog_coil_indices]

    # Passive dof values
    for i, ivc_idx in enumerate(ivc_indices):
        results["PASSIVES"]["IVC"]["DOF"][f"EIG_{i+1:02d}"] = coef[:, ivc_idx]
    for case_idx in case_indices:
        case_name = coef_names[case_idx]
        results["PASSIVES"][case_name]["DOF"]["CONSTANT_J"] = coef[:, case_idx]
    for psr_idx in psr_indices:
        psr_name = coef_names[psr_idx]
        results["PASSIVES"][psr_name]["DOF"]["CONSTANT_J"] = coef[:, psr_idx]
    results["PASSIVES"]["OVC"]["DOF"]["CONSTANT_J"] = coef[:, ovc_index]

    # Plasma dof values
    pls0_idx = np.where(coef_names == "pls0")[0][0]
    results["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["P_PRIME_DOF"] = coef[:, pls0_idx]
    pls1_idx = np.where(coef_names == "pls1")[0][0]
    results["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["FF_PRIM_DOF"] = coef[:, pls1_idx]
    pls2_idx = np.where(coef_names == "pls2")[0][0]
    results["GLOBAL"]["ASYM_Z_DOF"] = coef[:, pls2_idx]

    # Calculate PSI_A using
    # flux_norm = (psi_a - flux_total) / (psi_a - psi_b)
    # Rearranging gives:
    # psi_a = (flux_norm * psi_b - flux_total) / (flux_norm - 1)
    flux_norm_flat = flux_norm.reshape(flux_norm.shape[0], -1)
    flux_total_flat = flux_total.reshape(flux_total.shape[0], -1)
    min_indices = np.argmin(flux_norm_flat, axis=1)
    flux_norm_argmin = flux_norm_flat[np.arange(flux_norm_flat.shape[0]), min_indices]
    flux_total_argmin = flux_total_flat[np.arange(flux_total_flat.shape[0]), min_indices]
    psi_a = (flux_norm_argmin * flux_boundary - flux_total_argmin) \
        / (flux_norm_argmin - 1 + (flux_norm_argmin == 1))
    # Note that we added (flux_norm_argmin == 1) to the denominator to avoid division by zero
    results["GLOBAL"]["PSI_A"] = psi_a

    util.create_script_nodes(
        script_name="RTGSFIT",
        pulseNo_write=pulseNo_write,
        pulseNo_cal=None,
        run_name=run_name,
        run_info=run_description,
        workflows=None,
        link_best=False,
    )
    util.write_script_data(
        script_name="RTGSFIT",
        pulseNo_write=pulseNo_write,
        data_to_write=results.to_dictionary(),
        pulseNo_cal=None,
        run_name=run_name,
        run_description=run_description,
    )


if __name__ == "__main__":
    """
    Main function to run the `write_data_to_mdsplus` from the command line, with arguments

    Usage: `python write_data_to_mdsplus.py <pulseNo> [<run_name>]`
    Note: `run_name` argument is optional
    """
    args = sys.argv
    pulseNo = int(args[1])
    print(args[0])
    print(args[1])
    print(len(args))

    username = os.getlogin()
    # only PCS user should write to the real 5 digit pulse
    if username == "pcs.user":
        pulseNo_write = None
    elif username == "filip.janky":
        # Write to Filip's million pulse range
        pulseNo_write = pulseNo + 30_000_000
    elif username == "alex.prokopyszyn":
        # Write to Alex's million pulse range
        pulseNo_write = pulseNo + 52_000_000

    if len(args) > 2:
        run_name = args[2]
        write_data_to_mdsplus(pulseNo, run_name, pulseNo_write=pulseNo_write)
    else:
        write_data_to_mdsplus(pulseNo, pulseNo_write=pulseNo_write)
