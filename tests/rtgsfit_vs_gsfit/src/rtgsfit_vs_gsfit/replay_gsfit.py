"""
This module contains the function to run GSFIT and save the results to a dictionary.
"""

import os

from gsfit import Gsfit
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

    gsfit_controller.run()

    plasma = gsfit_controller.plasma
    flux_loops = gsfit_controller.flux_loops
    flux_loops_to_include = flux_loops.get_vec_bool(["*", "fit_settings", "include"])
    flux_loop_names = np.array(flux_loops.keys())
    flux_loop_names = flux_loop_names[flux_loops_to_include]
    flux_loop_r = np.zeros(len(flux_loop_names))
    flux_loop_z = np.zeros(len(flux_loop_names))
    for i, name in enumerate(flux_loop_names):
        flux_loop_r[i] = flux_loops.get_array1([name, "geometry", "r"])
        flux_loop_z[i] = flux_loops.get_array1([name, "geometry", "z"])

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
    # Save the output dictionary to a file
    output_file = os.path.join(cnst.DATA_DIR, 'gsfit_output_dict.npy')
    os.makedirs(cnst.DATA_DIR, exist_ok=True)
    with open(output_file, 'wb') as f:
        np.save(f, output_dict, allow_pickle=True)

if __name__ == "__main__":
    replay_gsfit()
