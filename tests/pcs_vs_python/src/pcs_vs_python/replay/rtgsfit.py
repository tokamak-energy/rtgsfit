import os

import ctypes
import mdsthin
import numpy as np
from pcs_vs_python.replay.input_signals import prep_coil_curr_2d, prep_meas_pcs_2d
import standard_utility

from pcs_vs_python.replay import input_signals, out_data_rtgsfit

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

    # Ensure the RTGSFIT library is loaded
    # and the function signatures are set correctly
    librtgsfit_path = os.path.join(cfg["rtgsfit_path"], 'lib', 'librtgsfit.so')
    rtgsfit_lib = ctypes.CDLL(librtgsfit_path)
    rtgsfit_lib.rtgsfit.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # meas_pcs
        ctypes.POINTER(ctypes.c_double),  # coil_curr
        ctypes.POINTER(ctypes.c_double),  # flux_norm
        ctypes.POINTER(ctypes.c_int),     # mask
        ctypes.POINTER(ctypes.c_double),  # flux_total
        ctypes.POINTER(ctypes.c_double),  # error
        ctypes.POINTER(ctypes.c_double),  # lcfs_r
        ctypes.POINTER(ctypes.c_double),  # lcfs_z
        ctypes.POINTER(ctypes.c_int),     # lcfs_n
        ctypes.POINTER(ctypes.c_double),  # coef
        ctypes.POINTER(ctypes.c_double),  # flux_boundary
        ctypes.POINTER(ctypes.c_double)   # plasma_current
    ]
    rtgsfit_lib.rtgsfit.restype = ctypes.c_int

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_preshot"])
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:Z_VEC").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_COEF").data()
        n_lcfs_max = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_LCFS_MAX").data()
        n_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_COIL").data()
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:SENS_NAMES").data()

    n_r = len(r_vec)
    n_z = len(z_vec)
    n_grid = n_r * n_z
    n_meas_pcs = len(sens_names)

    meas_pcs_2d = input_signals.prep_meas_pcs_2d(cfg)
    coil_curr_2d = input_signals.prep_coil_curr_2d(cfg)

    meas_pcs = np.zeros(n_meas_pcs, dtype=np.float64)
    coil_curr = np.zeros(n_coil, dtype=np.float64)
    flux_norm = initial_flux_norm(r_vec, z_vec,
                                  cfg["r_axis0"], cfg["z_axis0"],
                                  cfg["rho_boundary0"])
    mask = np.ones(n_grid, dtype=np.int32)
    flux_total = np.zeros(n_grid, dtype=np.float64)
    coef = np.zeros(n_coef, dtype=np.float64)
    error = np.array([0.0], dtype=np.float64)
    lcfs_r = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_z = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_n = np.array([0], dtype=np.int32)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)

    time_array_replay = np.linspace(cfg["t_min"], cfg["t_max"], cfg["n_t"],
                                    dtype=np.float64)

    out = out_data_rtgsfit.OutputDataRTGSFITClass(cfg)

    for iteration, time in enumerate(time_array_replay):

        meas_pcs = meas_pcs_2d[iteration, :]
        coil_curr = coil_curr_2d[iteration, :]
        out.update_mvalues(iteration, meas_pcs, coil_curr)
        print("iteration:", iteration)

        # Call the rtgsfit function
        result = rtgsfit_lib.rtgsfit(
            meas_pcs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            coil_curr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            flux_norm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            flux_total.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            error.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_n.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            flux_boundary.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            plasma_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )

        out.update_cvalues(iteration, flux_norm, mask, flux_total, error, lcfs_r, lcfs_z, lcfs_n,
                           coef, flux_boundary, plasma_current)
        print("Iteration:", iteration + 1)
        print("Time (ms):", time * 1e3)
        print("Result:", result)
        print("Plasma current:", plasma_current)

    print("Writing to MDSplus")
    print("Creating script nodes")
    # Write to MDSplus
    standard_utility.create_script_nodes(
        script_name="RTGSFIT",
        pulseNo_write=cfg["pulse_num_replay"],
        pulseNo_cal=None,
        run_name=cfg["run_name_replay"],
        run_info=cfg["run_description_replay"],
        workflows=None,
        link_best=False,
    )
    print("Writing data to MDSplus")
    standard_utility.write_script_data(
        script_name="RTGSFIT",
        pulseNo_write=cfg["pulse_num_replay"],
        data_to_write=out.data_dict,
        pulseNo_cal=None,
        run_name=cfg["run_name_replay"],
        run_description=cfg["run_description_replay"],
    )
    print("Finished writing to MDSplus")

if __name__ == "__main__":

    from pcs_vs_python import config_loader

    cfg = config_loader.load_and_prepare_config()
    replay_rtgsfit(cfg)
    print("RTGSFIT replay completed successfully.")


