import os

import ctypes
import mdsthin
import numpy as np

from rtgsfit_vs_gsfit import cnst, rtgsfit_signalprep

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

def replay_rtgsfit(time: float):
    """
    Run the RTGSFIT code for a single snapshot over multiple iterations.
    """

    # Ensure the RTGSFIT library is loaded
    # and the function signatures are set correctly
    librtgsfit_path = os.path.join(cnst.RTGSFIT_PATH, 'lib', 'librtgsfit.so')
    rtgsfit_lib = ctypes.CDLL(librtgsfit_path)
    rtgsfit_lib.rtgsfit.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # meas
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

    meas = rtgsfit_signalprep.prep_meas(time)
    coil_curr = rtgsfit_signalprep.prep_coil_curr(time)
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cnst.PULSE_NUM_WRITE)
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:Z_VEC").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COEF").data()
        n_lcfs_max = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_LCFS_MAX").data()
        n_meas = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_MEAS").data()
        n_coil = conn.get(f"\\RTGSFIT::TOP.{cnst.RUN_NAME}.PRESHOT:N_COIL").data()
    n_r = len(r_vec)
    n_z = len(z_vec)
    n_grid = n_r * n_z
    flux_norm = initial_flux_norm(r_vec, z_vec,
                                  cnst.R_AXIS0, cnst.Z_AXIS0,
                                  cnst.RHO_BOUNDARY0)
    mask = np.ones(n_grid, dtype=np.int32)
    flux_total = np.zeros(n_grid, dtype=np.float64)
    coef = np.zeros(n_coef, dtype=np.float64)
    error = np.array([0.0], dtype=np.float64)
    lcfs_r = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_z = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_n = np.array([0], dtype=np.int32)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)

    output_dict = {"meas" : np.zeros((cnst.N_ITERS + 1, n_meas), dtype=np.float64),
                   "coil_curr" : np.zeros((cnst.N_ITERS + 1, n_coil), dtype=np.float64),
                   "flux_norm" : np.zeros((cnst.N_ITERS + 1, n_grid), dtype=np.float64),
                   "mask" : np.zeros((cnst.N_ITERS + 1, n_grid), dtype=np.int32),
                   "flux_total" : np.zeros((cnst.N_ITERS + 1, n_grid), dtype=np.float64),
                   "error" : np.zeros((cnst.N_ITERS + 1, 1), dtype=np.float64),
                   "lcfs_r" : np.zeros((cnst.N_ITERS + 1, n_lcfs_max), dtype=np.float64),
                   "lcfs_z" : np.zeros((cnst.N_ITERS + 1, n_lcfs_max), dtype=np.float64),
                   "lcfs_n" : np.zeros((cnst.N_ITERS + 1, 1), dtype=np.int32),
                   "coef" : np.zeros((cnst.N_ITERS + 1, n_coef), dtype=np.float64),
                   "flux_boundary" : np.zeros((cnst.N_ITERS + 1, 1), dtype=np.float64),
                   "plasma_current" : np.zeros((cnst.N_ITERS + 1, 1), dtype=np.float64)}
    output_dict["meas"][0, :] = meas
    output_dict["coil_curr"][0, :] = coil_curr
    output_dict["flux_norm"][0, :] = flux_norm
    output_dict["mask"][0, :] = mask
    output_dict["flux_total"][0, :] = flux_total
    output_dict["error"][0, :] = error
    output_dict["lcfs_r"][0, :] = lcfs_r
    output_dict["lcfs_z"][0, :] = lcfs_z
    output_dict["lcfs_n"][0, :] = lcfs_n
    output_dict["coef"][0, :] = coef
    output_dict["flux_boundary"][0, :] = flux_boundary
    output_dict["plasma_current"][0, :] = plasma_current

    for i_iter in range(cnst.N_ITERS):

        meas_copy = meas.copy()
        coil_curr_copy = coil_curr.copy()

        # Call the rtgsfit function
        result = rtgsfit_lib.rtgsfit(
            meas_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            coil_curr_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
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

        print("Iteration:", i_iter + 1)
        print("Result:", result)

        output_dict["meas"][i_iter + 1, :] = meas
        output_dict["coil_curr"][i_iter + 1, :] = coil_curr
        output_dict["flux_norm"][i_iter + 1, :] = flux_norm
        output_dict["mask"][i_iter + 1, :] = mask
        output_dict["flux_total"][i_iter + 1, :] = flux_total
        output_dict["error"][i_iter + 1, :] = error
        output_dict["lcfs_r"][i_iter + 1, :] = lcfs_r
        output_dict["lcfs_z"][i_iter + 1, :] = lcfs_z
        output_dict["lcfs_n"][i_iter + 1, :] = lcfs_n
        output_dict["coef"][i_iter + 1, :] = coef
        output_dict["flux_boundary"][i_iter + 1, :] = flux_boundary
        output_dict["plasma_current"][i_iter + 1, :] = plasma_current

    # Save output_dict to a file
    output_file = os.path.join(cnst.DATA_DIR, 'rtgsfit_output_dict.npy')
    np.save(output_file, output_dict, allow_pickle=True)

if __name__ == "__main__":
    replay_rtgsfit(cnst.TIME)