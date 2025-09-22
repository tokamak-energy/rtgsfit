import os
import re

import ctypes
import matplotlib.pyplot as plt
import numpy as np

from rtgsfit_verify_analytic import cnst, measurements, read_constants_c

def initial_flux_norm(r_vec, z_vec):
    """
    Create the initial array for the normalised flux.
    """
    
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    opt_r = cnst.ANALYTIC_RO + 0.4 * cnst.ANALYTIC_RHOB
    opt_z = cnst.ANALYTIC_ZO + cnst.ANALYTIC_RHOB
    flux_norm_radius = 0.25
    rho = np.sqrt((r_grid - opt_r)**2 + (z_grid - opt_z)**2)
    flux_norm = (1 - np.cos(np.pi * rho / flux_norm_radius)) / 2 * (rho <= flux_norm_radius) \
              + (rho > flux_norm_radius)
    
    return flux_norm.flatten()

def replay_rtgsfit():

    librtgsfit_path = os.path.join(cnst.RTGSFIT_PATH, 'lib', 'librtgsfit.so')
    constants_c_path = os.path.join(cnst.REPO_PATH, 'data', 'constants.c')

    rtgsfit_lib = ctypes.CDLL(librtgsfit_path)

    # Define the argument types for the rtgsfit function
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
        ctypes.POINTER(ctypes.c_double),   # plasma_current
        ctypes.POINTER(ctypes.c_int32)   # lcfs_err_code
    ]
    # Define the return type for the rtgsfit function
    rtgsfit_lib.rtgsfit.restype = ctypes.c_int

    constants_c_dict = read_constants_c.constants_c_dict(constants_c_path)
    n_r = constants_c_dict["n_r"]
    n_z = constants_c_dict["n_z"]
    n_meas = constants_c_dict["n_meas"]
    n_coil = constants_c_dict["n_coil"]
    n_lcfs_max = constants_c_dict["n_lcfs_max"]
    n_coef = constants_c_dict["n_coef"]
    r_vec = constants_c_dict["r_vec"]
    z_vec = constants_c_dict["z_vec"]
    n_grid = n_r * n_z

    # Read in the flux loop coordinates and the BP probe coordinates
    data = np.loadtxt(os.path.join(cnst.DATA_DIR, 'flux_loops.txt'), dtype=np.float64, skiprows=1)
    fl_r_coords = data[:, 0]
    fl_z_coords = data[:, 1]
    fl_coords = np.zeros((len(fl_r_coords), 2), dtype=np.float64)
    fl_coords[:, 0] = fl_r_coords
    fl_coords[:, 1] = fl_z_coords
    data = np.loadtxt(os.path.join(cnst.DATA_DIR, 'bp_probes.txt'), dtype=np.float64, skiprows=1)
    bp_r_coords = data[:, 0]
    bp_z_coords = data[:, 1]
    bp_alpha_coords = data[:, 2]
    bp_probe_coords = np.zeros((len(bp_r_coords), 3), dtype=np.float64)
    bp_probe_coords[:, 0] = bp_r_coords
    bp_probe_coords[:, 1] = bp_z_coords
    bp_probe_coords[:, 2] = bp_alpha_coords

    meas = measurements.generate_measurements(
        fl_coords,
        bp_probe_coords,
        cnst.ANALYTIC_RO,
        cnst.ANALYTIC_ZO,
        cnst.ANALYTIC_PLASMA_CURRENT
    )
    coil_curr = np.zeros(n_coil, dtype=np.float64)
    flux_norm = initial_flux_norm(r_vec, z_vec)
    mask = np.ones(n_grid, dtype=np.int32)
    flux_total = np.zeros(n_grid, dtype=np.float64)
    coef = np.zeros(n_coef, dtype=np.float64)
    error = np.array([0.1], dtype=np.float64) # PROKOPYSZYN: Change to 0.0
    lcfs_r = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_z = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_n = np.array([0], dtype=np.int32)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)
    lcfs_err_code = np.array([0], dtype=np.int32)

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

        meas = measurements.generate_measurements(
            fl_coords,
            bp_probe_coords,
            cnst.ANALYTIC_RO,
            cnst.ANALYTIC_ZO,
            cnst.ANALYTIC_PLASMA_CURRENT
        )
        coil_curr = np.zeros(n_coil, dtype=np.float64)

        # Call the rtgsfit function
        result = rtgsfit_lib.rtgsfit(
            meas.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
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
            plasma_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lcfs_err_code.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
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
    output_file = os.path.join(cnst.DATA_DIR, 'output_dict.npy')
    np.save(output_file, output_dict, allow_pickle=True)

    # # Save just the last iteration data to a separate file
    # # for regression testing
    # reference_output_dict = {}

    # for key, value in output_dict.items():
    #     # Get the last iteration's 1D array
    #     last_iter_value = value[-1]
    #     print(key)
    #     print(np.shape(last_iter_value))
    #     reference_output_dict[key] = last_iter_value.copy()

    # reference_file = os.path.join(cnst.DATA_DIR, 'reference_output_dict.npy')
    # np.save(reference_file, reference_output_dict, allow_pickle=True)
