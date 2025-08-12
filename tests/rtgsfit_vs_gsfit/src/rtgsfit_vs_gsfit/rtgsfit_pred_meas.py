"""
This module contains calc_pred_meas() which computes RTGSFIT's
least-squares solution corresponding meas values.
"""

import mdsthin
import numpy as np

def calc_pred_meas(cfg: dict,
                   flux_norm: np.ndarray,
                   coef: np.ndarray,
                   coil_curr: np.ndarray) -> np.ndarray:
    """
    Calculate the predicted measurements from RTGSFIT output.

    Parameters:
        flux_norm (np.ndarray): Normalized flux values from RTGSFIT. With shape (N_GRID, )
        coef (np.ndarray): Coefficients for the linear combination. With shape (N_COEF, )
        coil_curr (np.ndarray): Coil currents from RTGSFIT. With shape (N_COIL, )

    Returns:
        np.ndarray: Predicted measurements. With shape (N_MEAS, )
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        n_grid = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_GRID").data()
        n_meas = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_MEAS").data()
        n_coef = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_COEF").data()
        n_pls = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_PLS").data()
        n_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_COIL").data()
        n_r = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_R").data()
        n_z = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:N_Z").data()
        g_grid_meas_weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT.GREENS:GRID_MEAS_W").data()
        g_coef_meas_weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT.GREENS:COEF_MEAS_W").data()
        g_meas_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT.GREENS:MEAS_COIL").data()
        weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:WEIGHT").data()
        r_flat = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:R_GRID").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name"]}.PRESHOT:Z_VEC").data()

    g_grid_meas_weight = g_grid_meas_weight.reshape((n_grid, n_meas))
    g_coef_meas_weight = g_coef_meas_weight.reshape((n_coef, n_meas))
    g_meas_coil = g_meas_coil.reshape((n_meas, n_coil)).T

    # Calculate g_pls_grid
    g_pls_grid = np.zeros((n_pls, n_grid), dtype=np.float64)
    g_pls_grid[0, :] = (1 - flux_norm) * r_flat
    g_pls_grid[1, :] = (1 - flux_norm)  / (r_flat * cfg["mu_0"])
    dflux_norm_dz = np.gradient(flux_norm.reshape((n_z, n_r)), z_vec, axis=0)
    g_pls_grid[2, :] = dflux_norm_dz.flatten()

    # Calculate first n_pls rows of g_coef_meas_weight
    g_coef_meas_weight[:n_pls, :] = g_pls_grid @ g_grid_meas_weight

    # Now calculate the predicted measurements from the plasma current
    pred_meas = coef @ g_coef_meas_weight
    pred_meas = pred_meas / weight

    # Now need to calculate the contribution from the PF coils
    pred_meas += coil_curr @ g_meas_coil

    return pred_meas
