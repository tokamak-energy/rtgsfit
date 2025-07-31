"""
This module contains the function to run GSFIT and save the results to a dictionary.
"""

import os
import typing

from gsfit import Gsfit
import st40_database
import numpy as np

from rtgsfit_vs_gsfit import cnst

def replay_gsfit():
    """
    Run GSFIT, then save the results to a dictionary.
    """

    from gsfit import Gsfit
    
    gsfit_controller = Gsfit(
        pulseNo=cnst.PULSE_NUM,
        run_name=cnst.RUN_NAME,
        run_description=cnst.RUN_DESCRIPTION,
        write_to_mds=True,
        pulseNo_write=cnst.PULSE_NUM_WRITE,
        settings_path=cnst.settings_path,
    )

    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["method"] = "user_defined"
    gsfit_controller.settings["GSFIT_code_settings.json"]["timeslices"]["user_defined"] = [cnst.TIME]
    gsfit_controller.settings["GSFIT_code_settings.json"]["database_writer"]["method"] = "tokamak_energy_mdsplus"
    elmag_run_name = gsfit_controller.settings["GSFIT_code_settings.json"]["database_reader"]["st40_mdsplus"]["workflow"]["elmag"]["run_name"]
    ivc_n_dof = gsfit_controller.settings["passive_dof_regularisation.json"]["IVC"]["n_dof"]

    gsfit_controller.run()

    plasma = gsfit_controller.plasma

    flux_loops = gsfit_controller.flux_loops
    flux_loops_to_include = flux_loops.get_vec_bool(["*", "fit_settings", "include"])
    flux_loop_names = np.array(flux_loops.keys())
    flux_loop_names = flux_loop_names[flux_loops_to_include]
    flux_loop_r = np.zeros(len(flux_loop_names))
    flux_loop_z = np.zeros(len(flux_loop_names))
    for i, name in enumerate(flux_loop_names):
        flux_loop_r[i] = flux_loops.get_array1([name, "geometry", "r"])[0]
        flux_loop_z[i] = flux_loops.get_array1([name, "geometry", "z"])[0]
    
    bp_probes = gsfit_controller.bp_probes
    bp_probes_to_include = bp_probes.get_vec_bool(["*", "fit_settings", "include"])
    bp_probe_names = np.array(bp_probes.keys())
    bp_probe_names = bp_probe_names[bp_probes_to_include]
    bp_probe_r = np.zeros(len(bp_probe_names))
    bp_probe_z = np.zeros(len(bp_probe_names))
    bp_probe_angle = np.zeros(len(bp_probe_names))
    for i, name in enumerate(bp_probe_names):
        bp_probe_r[i] = bp_probes.get_array1([name, "geometry", "r"])[0]
        bp_probe_z[i] = bp_probes.get_array1([name, "geometry", "z"])[0]
        bp_probe_angle[i] = bp_probes.get_array1([name, "geometry", "angle_pol"])[0]

    rogowski_coils = gsfit_controller.rogowski_coils
    rogowski_coils_to_include = rogowski_coils.get_vec_bool(["*", "fit_settings", "include"])
    rogowski_coil_names = np.array(rogowski_coils.keys())
    rogowski_coil_names = rogowski_coil_names[rogowski_coils_to_include]

    passives = gsfit_controller.passives
    ivc_r = passives.get_array1(["IVC", "geometry", "r"])
    ivc_z = passives.get_array1(["IVC", "geometry", "z"])
    ivc_dr = passives.get_array1(["IVC", "geometry", "d_r"])
    ivc_dz = passives.get_array1(["IVC", "geometry", "d_z"])
    ivc_area = passives.get_array1(["IVC", "geometry", "area"])
    ivc_angle1 = passives.get_array1(["IVC", "geometry", "angle1"])
    ivc_angle2 = passives.get_array1(["IVC", "geometry", "angle2"])

    ovc_r = passives.get_array1(["OVC", "geometry", "r"])
    ovc_z = passives.get_array1(["OVC", "geometry", "z"])

    elmag = st40_database.GetData(11012050, f"ELMAG#{elmag_run_name}")
    coil_names = typing.cast(list[str], elmag.get("COILS.COIL_NAMES"))
    coils = gsfit_controller.coils

    output_dict = {}
    output_dict["grid"] = {}
    output_dict["grid"]["r"] = plasma.get_array1(["grid", "r"])
    output_dict["grid"]["z"] = plasma.get_array1(["grid", "z"])
    output_dict["two_d"] = {}
    output_dict["two_d"]["psi"] = plasma.get_array3(["two_d", "psi"])[0, :, :]
    output_dict["two_d"]["psi_n"] = plasma.get_array3(["two_d", "psi_n"])[0, :, :]
    output_dict["p_boundary"] = {}
    output_dict["p_boundary"]["rbnd"] = plasma.get_array2(["p_boundary", "rbnd"])[0, :]
    output_dict["p_boundary"]["zbnd"] = plasma.get_array2(["p_boundary", "zbnd"])[0, :]
    output_dict["p_boundary"]["nbnd"] = plasma.get_vec_usize(["p_boundary", "nbnd"])[0]
    output_dict["flux_loops"] = {}
    output_dict["flux_loops"]["names"] = flux_loop_names
    output_dict["flux_loops"]["r"] = flux_loop_r
    output_dict["flux_loops"]["z"] = flux_loop_z
    output_dict["bp_probes"] = {}
    output_dict["bp_probes"]["names"] = bp_probe_names
    output_dict["bp_probes"]["r"] = bp_probe_r
    output_dict["bp_probes"]["z"] = bp_probe_z
    output_dict["bp_probes"]["angle"] = bp_probe_angle
    output_dict["rogowski_coils"] = {}
    for name in rogowski_coil_names:
        output_dict["rogowski_coils"][name] = {}
        output_dict["rogowski_coils"][name]["r"] = rogowski_coils.get_array1([name, "geometry", "r"])
        output_dict["rogowski_coils"][name]["z"] = rogowski_coils.get_array1([name, "geometry", "z"]) 
    output_dict["IVC"] = {}
    output_dict["IVC"]["r"] = ivc_r
    output_dict["IVC"]["z"] = ivc_z
    output_dict["IVC"]["d_r"] = ivc_dr
    output_dict["IVC"]["d_z"] = ivc_dz
    output_dict["IVC"]["area"] = ivc_area
    output_dict["IVC"]["angle1"] = ivc_angle1
    output_dict["IVC"]["angle2"] = ivc_angle2
    output_dict["IVC"]["dof"] = {}
    for i_ivc_dof in range(ivc_n_dof):
        output_dict["IVC"]["dof"][f"eig_{i_ivc_dof+1:02d}"] = {}
        output_dict["IVC"]["dof"][f"eig_{i_ivc_dof+1:02d}"]["current_distribution"] = \
            passives.get_array1(["IVC", "dof", f"eig_{i_ivc_dof+1:02d}", "current_distribution"])
        output_dict["IVC"]["dof"][f"eig_{i_ivc_dof+1:02d}"]["calculated"] = \
            passives.get_array1(["IVC", "dof", f"eig_{i_ivc_dof+1:02d}", "calculated"])[0]
    output_dict["OVC"] = {}
    output_dict["OVC"]["r"] = ovc_r
    output_dict["OVC"]["z"] = ovc_z
    output_dict["OVC"]["current_distribution"] = passives.get_array1(["OVC", "dof", "*", "current_distribution"])
    output_dict["OVC"]["calculated"] = passives.get_array1(["OVC", "dof", "*", "calculated"])[0]
        
    output_dict["coils"] = {}
    for name in coil_names:
        output_dict["coils"][name] = {}
        output_dict["coils"][name]["r"] = coils.get_array1(["pf", name, "geometry", "r"])
        output_dict["coils"][name]["z"] = coils.get_array1(["pf", name, "geometry", "z"])
        output_dict["coils"][name]["d_r"] = coils.get_array1(["pf", name, "geometry", "d_r"])
        output_dict["coils"][name]["d_z"] = coils.get_array1(["pf", name, "geometry", "d_z"])
        coil_curr_experimental = coils.get_array1(["pf", name, "i", "measured_experimental"])
        coil_time_experimental = coils.get_array1(["pf", name, "i", "time_experimental"])
        coil_curr_interp = np.interp(cnst.TIME, coil_time_experimental, coil_curr_experimental)
        output_dict["coils"][name]["i"] = coil_curr_interp

    output_file = os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy')
    os.makedirs(cnst.DATA_DIR, exist_ok=True)
    with open(output_file, 'wb') as f:
        np.save(f, output_dict, allow_pickle=True)

if __name__ == "__main__":
    replay_gsfit()
