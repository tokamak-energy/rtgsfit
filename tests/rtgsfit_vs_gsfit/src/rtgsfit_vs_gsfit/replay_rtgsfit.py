import os

import ctypes
import mdsthin
import numpy as np

from rtgsfit_vs_gsfit import prep_meas_coil_curr

def initial_flux_norm(r_vec, z_vec, 
                      r_axis, z_axis,
                      rho_boundary):
    """
    Create concentric circles for the initial normalised flux centred at (r_axis, z_axis)
    with a radius of rho_boundary.
    The flux is normalised such that it is 0 at the boundary and 1 at the centre.

    Parameters:
    r_vec : np.ndarray
        1D Array of radial coordinates.
    z_vec : np.ndarray
        1D Array of vertical coordinates.
    r_axis : float
        Radial coordinate of the axis of the plasma.
    z_axis : float
        Vertical coordinate of the axis of the plasma.
    rho_boundary : float
        Radius of the boundary of the plasma in metres.

    Returns:
    np.ndarray
        Flattened 2D array of normalised flux values.
    """
    
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    rho = np.sqrt((r_grid - r_axis)**2 + (z_grid - z_axis)**2)
    flux_norm = (1 - np.cos(np.pi * rho / rho_boundary)) / 2 * (rho <= rho_boundary) \
              + (rho > rho_boundary)

    return flux_norm.flatten()

def replay_rtgsfit(cfg: dict):
    """
    Run the RTGSFIT code for a single snapshot over multiple iterations.
    """

    if not cfg["rtgsfit_node_initialised"]:
        raise RuntimeError("RTGSFIT node is not initialised. Please run `initialise_rtgsfit_node()` first.")
    if not cfg["rtgsfit_compiled"]:
        raise RuntimeError("RTGSFIT code is not compiled. Please run `compile_rtgsfit()` first.")

    # Ensure the RTGSFIT library is loaded
    # and the function signatures are set correctly
    librtgsfit_path = os.path.join(cfg["rtgsfit_path"], 'lib', 'librtgsfit.so')
    rtgsfit_lib = ctypes.CDLL(librtgsfit_path)
    rtgsfit_lib.rtgsfit.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # meas_pcs
        ctypes.POINTER(ctypes.c_double),  # coil_curr
        ctypes.POINTER(ctypes.c_double),  # flux_norm
        ctypes.POINTER(ctypes.c_int32),   # mask
        ctypes.POINTER(ctypes.c_double),  # flux_total
        ctypes.POINTER(ctypes.c_double),  # chi_sq_err
        ctypes.POINTER(ctypes.c_double),  # lcfs_r
        ctypes.POINTER(ctypes.c_double),  # lcfs_z
        ctypes.POINTER(ctypes.c_int32),   # lcfs_n
        ctypes.POINTER(ctypes.c_double),  # coef
        ctypes.POINTER(ctypes.c_double),  # flux_boundary
        ctypes.POINTER(ctypes.c_double),  # plasma_current
        ctypes.POINTER(ctypes.c_int32),   # lcfs_err_code
        ctypes.POINTER(ctypes.c_int32),   # lapack_dgelss_info
        ctypes.POINTER(ctypes.c_double),  # meas_model
        ctypes.c_int32,                   # n_meas_model
    ]
    rtgsfit_lib.rtgsfit.restype = None

    meas_pcs = prep_meas_coil_curr.prep_meas_pcs(cfg)
    coil_curr = prep_meas_coil_curr.prep_coil_curr(cfg)
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:Z_VEC").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_COEF").data()
        n_lcfs_max = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_LCFS_MAX").data()
        n_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_COIL").data()
        n_meas = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_MEAS").data()

    n_r = len(r_vec)
    n_z = len(z_vec)
    n_grid = n_r * n_z
    n_meas_pcs = len(meas_pcs)
    flux_norm = initial_flux_norm(r_vec, z_vec,
                                  cfg["r_axis0"], cfg["z_axis0"],
                                  cfg["rho_boundary0"])
    mask = np.ones(n_grid, dtype=np.int32)
    mask = mask.flatten().astype(np.int32)
    flux_total = np.zeros(n_grid, dtype=np.float64)
    coef = np.zeros(n_coef, dtype=np.float64)
    chi_sq_err = np.array([0.0], dtype=np.float64)
    lcfs_r = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_z = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_n = np.array([0], dtype=np.int32)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)
    lcfs_err_code = np.array([0], dtype=np.int32)
    lapack_dgelss_info = np.array([0], dtype=np.int64)
    meas_model = np.zeros(n_meas, dtype=np.float64)
    n_meas_model = np.array([n_meas], dtype=np.int32)

    output_dict = {"meas_pcs" : np.zeros((cfg["n_iters"] + 1, n_meas_pcs), dtype=np.float64),
                   "coil_curr" : np.zeros((cfg["n_iters"] + 1, n_coil), dtype=np.float64),
                   "flux_norm" : np.zeros((cfg["n_iters"] + 1, n_grid), dtype=np.float64),
                   "mask" : np.zeros((cfg["n_iters"] + 1, n_grid), dtype=np.int32),
                   "flux_total" : np.zeros((cfg["n_iters"] + 1, n_grid), dtype=np.float64),
                   "chi_sq_err" : np.zeros((cfg["n_iters"] + 1), dtype=np.float64),
                   "lcfs_r" : np.zeros((cfg["n_iters"] + 1, n_lcfs_max), dtype=np.float64),
                   "lcfs_z" : np.zeros((cfg["n_iters"] + 1, n_lcfs_max), dtype=np.float64),
                   "lcfs_n" : np.zeros((cfg["n_iters"] + 1), dtype=np.int32),
                   "coef" : np.zeros((cfg["n_iters"] + 1, n_coef), dtype=np.float64),
                   "flux_boundary" : np.zeros((cfg["n_iters"] + 1, 1), dtype=np.float64),
                   "plasma_current" : np.zeros((cfg["n_iters"] + 1, 1), dtype=np.float64),
                   "lcfs_err_code" : np.zeros((cfg["n_iters"] + 1, 1), dtype=np.int32),
                   "lapack_dgelss_info" : np.zeros((cfg["n_iters"] + 1), dtype=np.int32),
                   "meas_model" : np.zeros((cfg["n_iters"] + 1, n_meas), dtype=np.float64),
                   "n_meas_model": np.zeros((cfg["n_iters"] + 1), dtype=np.int32),
                   }
    output_dict["meas_pcs"][0, :] = meas_pcs
    output_dict["coil_curr"][0, :] = coil_curr
    output_dict["flux_norm"][0, :] = flux_norm
    output_dict["mask"][0, :] = mask
    output_dict["flux_total"][0, :] = flux_total
    output_dict["chi_sq_err"][0] = chi_sq_err[0]
    output_dict["lcfs_r"][0, :] = lcfs_r
    output_dict["lcfs_z"][0, :] = lcfs_z
    output_dict["lcfs_n"][0] = lcfs_n[0]
    output_dict["coef"][0, :] = coef
    output_dict["flux_boundary"][0, :] = flux_boundary
    output_dict["plasma_current"][0, :] = plasma_current
    output_dict["lcfs_err_code"][0] = lcfs_err_code[0]
    output_dict["lapack_dgelss_info"][0] = lapack_dgelss_info[0]
    output_dict["meas_model"][0, :] = meas_model
    output_dict["n_meas_model"][0] = n_meas_model[0]

    for i_iter in range(cfg["n_iters"]):

        meas_pcs_copy = meas_pcs.copy()
        coil_curr_copy = coil_curr.copy()

        # Call the rtgsfit function
        rtgsfit_lib.rtgsfit(
            meas_pcs_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            coil_curr_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            flux_norm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            flux_total.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            chi_sq_err.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_n.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            flux_boundary.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            plasma_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_err_code.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            lapack_dgelss_info.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            meas_model.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(n_meas)
        )

        print("Iteration:", i_iter + 1)

        output_dict["meas_pcs"][i_iter + 1, :] = meas_pcs
        output_dict["coil_curr"][i_iter + 1, :] = coil_curr
        output_dict["flux_norm"][i_iter + 1, :] = flux_norm
        output_dict["mask"][i_iter + 1, :] = mask
        output_dict["flux_total"][i_iter + 1, :] = flux_total
        output_dict["chi_sq_err"][i_iter + 1] = chi_sq_err[0]
        output_dict["lcfs_r"][i_iter + 1, :] = lcfs_r
        output_dict["lcfs_z"][i_iter + 1, :] = lcfs_z
        output_dict["lcfs_n"][i_iter + 1] = lcfs_n[0]
        output_dict["coef"][i_iter + 1, :] = coef
        output_dict["flux_boundary"][i_iter + 1, :] = flux_boundary
        output_dict["plasma_current"][i_iter + 1, :] = plasma_current
        output_dict["lcfs_err_code"][i_iter + 1] = lcfs_err_code[0]
        output_dict["lapack_dgelss_info"][i_iter + 1] = lapack_dgelss_info[0]
        output_dict["meas_model"][i_iter + 1, :] = meas_model
        output_dict["n_meas_model"][i_iter + 1] = n_meas_model[0]        

    # Save output_dict to a file
    np.save(cfg["rtgsfit_output_dict_path"], output_dict, allow_pickle=True)
    cfg["rtgsfit_replayed"] = True

if __name__ == "__main__":

    from rtgsfit_vs_gsfit import config_loader, rtgsfit_compile_setup

    cfg = config_loader.load_and_prepare_config()

    rtgsfit_compile_setup.rtgsfit_mds_nodeclear(cfg)
    rtgsfit_compile_setup.initialise_rtgsfit_node(cfg)
    rtgsfit_compile_setup.compile_rtgsfit(cfg)
    # cfg["rtgsfit_node_initialised"] = True
    # cfg["rtgsfit_compiled"] = True

    replay_rtgsfit(cfg)
    print("RTGSFIT replay completed successfully.")


