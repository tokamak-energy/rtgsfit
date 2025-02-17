import re
import sys

import numpy as np
from MDSplus import Connection


def const_to_file(const_dict, file_read="constants.h", file_write="constants.c"):
    with open(file_read, "r") as file:
        file_contents = file.read()

    # Define a function to perform variable replacements
    def replace_variables(match):
        assert (match.lastindex == 4) or (match.lastindex == 5)

        var_name = match.group(4).lower()
        var_name_in_file = " ".join([match.group(xx) for xx in range(2, 5)])
        if match.lastindex == 5:
            var_name_in_file = var_name_in_file + match.group(match.lastindex)

        if var_name in const_dict:
            replacement = const_dict[var_name]
            # Check the data type of the replacement value
            if isinstance(replacement, bool):
                return f"{var_name_in_file} = {int(replacement)};"
            elif isinstance(replacement, (int, float)):
                return f"{var_name_in_file} = {replacement};"
            elif isinstance(replacement, (np.int32, float)):
                return f"{var_name_in_file} = {replacement};"
            elif isinstance(replacement, np.ndarray) and (match.lastindex == 5):
                data_list = map(str, replacement.flatten().tolist())
                data_list = [x if ((ii + 1) % 20) else x + "\n" for ii, x in enumerate(data_list)]
                array_values = ", ".join(data_list)
                return f"{var_name_in_file} = {{{array_values}}};"
            else:
                raise Exception(f"{type(replacement)} is not catered for")
        else:
            raise Exception(f"{var_name} is not in dictionary")

    pattern = r"\b(\w+)\s*\b(\w+)\s*\b(\w+)\s*\b(\w+)(\[\s*\d*\s*\])?\s*;"

    # Replace variables in the C code
    new_c_code = re.sub(pattern, replace_variables, file_contents)

    with open(file_write, "w") as file:
        file.write(new_c_code)


if __name__ == "__main__":
    args = sys.argv[1:]
    for arg in args:
        key, value = arg.split("=")
        if key.upper() == "SHOT":
            shot = int(value)
        elif key.upper() == "RUN_NAME":
            run_name = value
        else:
            raise Exception(f"Unknown argument: {key}")
    conn = Connection("smaug")
    conn.openTree("RTGSFIT", shot)

    # Read from MDSplus and add to dictionary. Note, this is not qute a 1:1 mapping of variable names
    # because MDSplus restricts node names to 12 characters.
    data_dictionary = {}
    data_dictionary["dr"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:DR").data()
    data_dictionary["dz"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:DZ").data()
    data_dictionary["frac"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:FRAC").data()
    data_dictionary["inv_r_ltrb_mu0"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:INV_R_L_MU0").data()
    data_dictionary["inv_r_mu0"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:INV_R_MU0").data()
    data_dictionary["limit_idx"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:LIMIT_IDX").data()
    data_dictionary["limit_weight"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:LIMIT_W").data()
    data_dictionary["lower_band"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:LOWER_BAND").data()
    data_dictionary["mask_lim"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:MASK_LIM").data()
    data_dictionary["n_bp_probes"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_BP_PROBES").data()
    data_dictionary["n_coef"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_COEF").data()
    data_dictionary["n_coil"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_COIL").data()
    data_dictionary["n_flux_loops"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_F_LOOPS").data()
    data_dictionary["n_grid"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_GRID").data()
    data_dictionary["n_intrp"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_INTRP").data()
    data_dictionary["n_lcfs_max"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_LCFS_MAX").data()
    data_dictionary["n_limit"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_LIMIT").data()
    data_dictionary["n_ltrb"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_LTRB").data()
    data_dictionary["n_meas"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_MEAS").data()
    data_dictionary["n_pls"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_PLS").data()
    data_dictionary["n_r"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_R").data()
    data_dictionary["n_rogowski_coils"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_ROG_COILS").data()
    data_dictionary["n_vess"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_VESS").data()
    data_dictionary["n_xpt_max"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_XPT_MAX").data()
    data_dictionary["n_z"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:N_Z").data()
    data_dictionary["perm_idx"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:PERM_IDX").data()
    data_dictionary["r_grid"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:R_GRID").data()
    data_dictionary["r_mu0_dz2"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:R_MU0_DZ2").data()
    data_dictionary["r_vec"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:R_VEC").data()
    data_dictionary["thresh"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:THRESH").data()
    data_dictionary["upper_band"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:UPPER_BAND").data()
    data_dictionary["weight"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:WEIGHT").data()
    data_dictionary["z_grid"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:Z_GRID").data()
    data_dictionary["z_vec"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:Z_VEC").data()
    data_dictionary["g_coef_meas_weight"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.COEF_MEAS_W").data()
    data_dictionary["g_grid_coil"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.GRID_COIL").data()
    data_dictionary["g_grid_meas_weight"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.GRID_MEAS_W").data()
    data_dictionary["g_grid_vessel"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.GRID_VESSEL").data()
    data_dictionary["g_ltrb"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.LTRB").data()
    data_dictionary["g_meas_coil"] = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT:GREENS.MEAS_COIL").data()

    const_to_file(data_dictionary)
