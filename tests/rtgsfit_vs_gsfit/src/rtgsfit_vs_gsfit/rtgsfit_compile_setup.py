"""
Module for interfacing with MDS+ and building RTGSFIT.

This module provides utility functions to:
- `rtgsfit_mds_nodeclear`: Delete existing MDS+ node for RTGSFIT.
- `initialise_rtgsfit_node`: Generate and populate the MDS+ node with the required data.
- `compile_rtgsfit`: Compile the RTGSFIT code using its `Makefile`.
"""
import os
import subprocess

import MDSplus
import numpy as np

from gsfit import Gsfit

def rtgsfit_mds_nodeclear(cfg: dict):
    """
    Delete the existing RTGSFIT MDSplus node and ensure it is removed.
    """

    # Delete the existing RTGSFIT MDSplus node (if it exists)
    def delete_node_recursive(tree, node):
        for child in node.getNodeWild('*'):
            delete_node_recursive(tree, child)
        print(f"Deleting node {node.getPath()}")
        tree.deleteNode(node.getPath())
    tree = MDSplus.Tree("RTGSFIT", cfg["pulse_num_write"], "EDIT")
    try:
        node = tree.getNode(f":{cfg['run_name']}")
        delete_node_recursive(tree, node)
        print(f"Deleted node {cfg['run_name']} and its children.")
    except Exception as e:
        print(f"Failed to delete node {cfg['run_name']}: {e}")
    tree.write()
    tree.close()

    # Re-open tree to check if node exists
    def node_exists(tree, path):
        try:
            tree.getNode(path)
            return True
        except Exception:
            return False
    tree_check = MDSplus.Tree("RTGSFIT", cfg["pulse_num_write"])
    if not node_exists(tree_check, f":{cfg['run_name']}"):
        deleted = True
    else:
        deleted = False
    assert deleted, f"Node {cfg['run_name']} was not deleted successfully."

def initialise_rtgsfit_node(cfg: dict):
    """
    Initialise the RTGSFIT MDSplus node by running GSFIT with the specified settings.
    This function sets up the GSFIT controller, modifies the settings, and runs the analysis
    to generate the necessary data for RTGSFIT.
    """

    # Construct the GSFit object; using the "st40_setup_for_rtgsfit" settings
    gsfit_controller = Gsfit(
        pulseNo=cfg['pulse_num'],
        run_name=cfg['run_name'],
        run_description=cfg['run_description'],
        write_to_mds=True,
        pulseNo_write=cfg['pulse_num_write'],
        settings_path=cfg['settings_path'],
    )

    # Change the analysis_name, so that GSFit writes into RT-GSFit MDSplus tree
    gsfit_controller.analysis_name = "RTGSFIT"

    # Add a list of signals to be read using PCS formatting
    gsfit_controller.results["PRESHOT"]["COIL_SIGNALS"] = np.array(cfg["coil_signals"])
    # coil_matrix = np.array(
    #     [
    #         # BVL_PSU, BVUB_PSU, BVUT_PSU, DIV_PSU, MCVC_PSU, PSH_PSU, ROG_MCWIRE, SOL_PSU
    #         [0.0,      0.0,      0.0,      0.0,     0.0,      0.0,     0.0,        1.0],  # SOL coil
    #         [0.0,      0.0,      0.0,      0.0,     0.0,      0.0,     1.0,        0.0],  # MCT coil
    #         [0.0,      0.0,      0.0,      0.0,     1.0,      0.0,     1.0,        0.0],  # MCB coil
    #         [0.0,      0.0,      0.0,      1.0,     0.0,      0.0,     0.0,        0.0],  # DIV coil
    #         [1.0,      0.0,      0.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVL coil
    #         [0.0,      0.0,      1.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVUT coil
    #         [0.0,      1.0,      0.0,      0.0,     0.0,      0.0,     0.0,        0.0],  # BVUB coil
    #         [0.0,      0.0,      0.0,      0.0,     0.0,      1.0,     0.0,        0.0],  # PSH coil
    #     ]
    # )
    coil_matrix = np.array(cfg["coil_matrix"])
    gsfit_controller.results["PRESHOT"]["COIL_MATRIX"] = coil_matrix

    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_r"] = cfg["n_r"]
    gsfit_controller.settings["GSFIT_code_settings.json"]["grid"]["n_z"] = cfg["n_z"]

    # Run
    gsfit_controller.run()
    cfg["rtgsfit_node_initialised"] = True

    # Create IVC dictionary and save it, need to save this to MDS+ instead
    ivc_dict = {}
    passives = gsfit_controller.passives
    ivc_dict["r"] = passives.get_array1(["IVC", "geometry", "r"])
    ivc_dict["z"] = passives.get_array1(["IVC", "geometry", "z"])
    ivc_dict["dr"] = passives.get_array1(["IVC", "geometry", "d_r"])
    ivc_dict["dz"] = passives.get_array1(["IVC", "geometry", "d_z"])
    n_eigs = gsfit_controller.settings["passive_dof_regularisation.json"]["IVC"]["n_dof"]
    n_segs = len(passives.get_array1(["IVC", "dof", f"eig_01", "current_distribution"]))
    ivc_dict["current_distributions"] = np.zeros((n_eigs, n_segs))
    for eig_num in range(n_eigs):
        ivc_dict["current_distributions"][eig_num, :] = \
            passives.get_array1(["IVC", "dof", f"eig_{eig_num + 1:02d}", "current_distribution"])
    np.save(cfg["ivc_dict_path"], ivc_dict, allow_pickle=True)

    # Create coef_names list and save it
    passives = gsfit_controller.passives
    coef_names = ["pls0", "pls1", "pls2"]
    passive_names = passives.keys()
    for passive_name in passive_names:
        dof_names = passives.keys([passive_name, "dof"])
        for dof_name in dof_names:
            if dof_name == "constant_current_density":
                coef_names.append(passive_name)
            elif dof_name.startswith("eig_"):
                coef_names.append(dof_name)
            else:
                raise ValueError(f"Unknown DoF name: {dof_name}")
    # Save coef_names to a file
    with open(cfg["coef_names_path"], "w") as f:
        f.writelines(name + "\n" for name in coef_names)

    # create meas_names list and save it
    meas_names = []
    flux_loops = gsfit_controller.flux_loops
    bp_probes = gsfit_controller.bp_probes
    rogowski_coils = gsfit_controller.rogowski_coils
    # Get the "included" sensors
    flux_loops_to_include = flux_loops.get_vec_bool(["*", "fit_settings", "include"])
    bp_probes_to_include = bp_probes.get_vec_bool(["*", "fit_settings", "include"])
    rogowski_coils_to_include = rogowski_coils.get_vec_bool(["*", "fit_settings", "include"])
    meas_names = []
    for i_flux_loop, floop_name in enumerate(flux_loops.keys()):
        if flux_loops_to_include[i_flux_loop]:
            meas_names.append(floop_name)
    for i_bp_probe, bp_name in enumerate(bp_probes.keys()):
        if bp_probes_to_include[i_bp_probe]:
            meas_names.append(bp_name)
    rogowski_coils_names_rtgsfit_order = [
        "INIVC000",
        "BVLT",
        "BVLB",
        "GASBFLT",
        "GASBFLB",
        "HFSPSRT",
        "HFSPSRB",
        "DIVPSRT",
        "DIVPSRB",
    ]
    for rog_name in rogowski_coils_names_rtgsfit_order:
        meas_names.append(rog_name)
    for passive_name in passive_names:
        passive_regularisation_local = passives.get_array2([passive_name, "regularisations"])
        [n_reg_local, _] = passive_regularisation_local.shape
        for i_reg in range(n_reg_local):
            meas_names.append(f"{passive_name}_reg_{i_reg}")
    # save meas_names to file
    with open(cfg["meas_names_path"], "w") as f:
        f.writelines(name + "\n" for name in meas_names)

def compile_rtgsfit(cfg: dict):
    """
    Clean and then compile RTGSFIT.
    """

    # Check the RTGSFIT node is initialised
    if not cfg["rtgsfit_node_initialised"]:
        raise RuntimeError("RTGSFIT node is not initialised. Please run `initialise_rtgsfit_node()` first.")

    # Clean the RTGSFIT source directory
    os.chdir(cfg['rtgsfit_src_path'])
    subprocess.run(["make", "clean"], cwd=cfg['rtgsfit_src_path'], check=True)
    # Check if the object files are removed
    object_files = [f for f in os.listdir(cfg['rtgsfit_src_path']) if f.endswith('.o')]
    assert not object_files, "Object files were not removed during cleaning."
    # Check constants.c was removed
    constants_c_path = os.path.join(cfg['rtgsfit_src_path'], 'constants.c')
    assert not os.path.exists(constants_c_path), f"constants.c was not removed: {constants_c_path}"

    # Compile RTGSFIT
    os.chdir(cfg['rtgsfit_src_path'])
    subprocess.run(
        f"make SHOT={cfg['pulse_num_write']} RUN_NAME={cfg['run_name']} DEBUG=1",
        cwd=cfg['rtgsfit_src_path'],
        check=True,
        shell=True
    )
    cfg["rtgsfit_compiled"] = True

if __name__ == "__main__":

    from rtgsfit_vs_gsfit import config_loader

    cfg = config_loader.load_and_prepare_config()
    rtgsfit_mds_nodeclear(cfg)
    initialise_rtgsfit_node(cfg)
    compile_rtgsfit(cfg)