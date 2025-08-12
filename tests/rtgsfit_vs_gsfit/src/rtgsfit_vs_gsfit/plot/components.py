"""
Provides reusable subplot helper functions that take an `ax` argument
and draw individual plot elements such as contours or annotations.
"""

import os
import typing

import matplotlib.pyplot as plt
import mdsthin
import numpy as np

def psi_contours(iteration: int,
                 ax: plt.Axes,
                 cfg: dict,
                 levels: int = 30):
    """
    Draw the ψ [Wb] contours on the given axis.
    """
    
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        n_r = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:N_R").data()
        n_z = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:N_Z").data()
        r_vec = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:R_VEC").data()
        z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:Z_VEC").data()
    rtgsfit_output_dict =np.load(cfg["rtgsfit_output_dict_path"],
                                 allow_pickle=True).item()
    psi_rtgsfit = rtgsfit_output_dict["flux_total"][iteration].reshape(n_z, n_r)
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        psi_gsfit = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.TWO_D:PSI").data()[0, :, :]

    ax.contour(r_vec, z_vec, psi_rtgsfit,
               levels=levels,
               colors="tab:blue")
    ax.contour(r_vec, z_vec, psi_gsfit,
               levels=levels,
               colors="tab:orange")
    ax.set_aspect('equal')
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("ψ [Wb]")

def limiter_scatter(ax: plt.Axes,
                    cfg: dict):
    """
    Draw the limiter points on the given axis.
    """
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        limiter_r = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:LIMIT_R").data()
        limiter_z = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:LIMIT_Z").data()
    
    ax.scatter(limiter_r, limiter_z,
               color="tab:green")
    
def flux_loop_annotate(ax: plt.Axes,
                       cfg: dict):
    """
    Annotate the flux loop positions on the given axis.
    """

    from st40_database import GetData

    mag =  GetData(cfg["pulse_num"], f"MAG#BEST", is_fail_quiet=False)
    names_long = typing.cast(list[str], mag.get("FLOOP.ALL.NAMES"))
    names_short = np.char.replace(names_long, "FLOOP_", "L")

    r_vec = typing.cast(np.typing.NDArray[np.float64], mag.get("FLOOP.ALL.R"))
    z_vec = typing.cast(np.typing.NDArray[np.float64], mag.get("FLOOP.ALL.Z"))

    for j, name in enumerate(names_short):
        ax.annotate(name, xy=(r_vec[j], z_vec[j]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, color='black', ha='center')
        
def bp_probe_annotate(ax: plt.Axes,
                      cfg: dict):
    """
    Annotate the Bp probe positions on the given axis.
    """

    from st40_database import GetData

    mag = GetData(cfg["pulse_num"], f"MAG#BEST", is_fail_quiet=False)
    names_long = typing.cast(list[str], mag.get("BPPROBE.ALL.NAMES"))
    names_short = np.char.replace(names_long, "BPPROBE_", "P")

    r_vec = typing.cast(np.typing.NDArray[np.float64], mag.get("BPPROBE.ALL.R"))
    z_vec = typing.cast(np.typing.NDArray[np.float64], mag.get("BPPROBE.ALL.Z"))

    for j, name in enumerate(names_short):
        ax.annotate(name, xy=(r_vec[j], z_vec[j]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, color='black', ha='center')
        
def flux_loop_line(iteration: int,
                   ax: plt.Axes,
                   cfg: dict):
    """
    Plot measured and predicted flux loop values against flux loop names.

    The plot includes:
      • Measured flux loop values.
      • Predicted values from RTGSFIT.
      • Predicted values from GSFIT.
    """

    from rtgsfit_vs_gsfit import rtgsfit_pred_meas

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:SENS_NAMES").data()
    flux_loop_indices = np.where(np.char.startswith(sens_names, "PSI_FLOOP_"))[0]
    long_names = sens_names[flux_loop_indices]
    short_names = np.char.replace(long_names, "PSI_FLOOP_", "L")

    rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                  allow_pickle=True).item()
    flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
    coef = rtgsfit_output_dict["coef"][iteration, :]
    coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
    pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        fl_include = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:INCLUDE").data() == 1
        gsfit_fl_meas = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:MVALUE").data()[0, fl_include]
        gsfit_fl_pred = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.FLOOP:CVALUE").data()[0, fl_include]


    ax.plot(short_names,
            pred_meas_rtgsfit[flux_loop_indices],
            label="RTGSFIT Predicted",
            color='tab:blue')
    ax.plot(short_names,
            gsfit_fl_pred,
            label="GSFIT Predicted",
            color='tab:orange')
    ax.plot(short_names,
            gsfit_fl_meas,
            label="Measured",
            color='tab:green')
    ax.set_ylabel("Flux Loop Value [Wb]")
    ax.tick_params(axis='x', rotation=90)
    ax.legend()

def bp_probe_line(iteration: int,
                  ax: plt.Axes,
                  cfg: dict):
    """
    Plot measured and predicted Bp probe values against Bp probe names.

    The plot includes:
      • Measured Bp probe values.
      • Predicted values from RTGSFIT.
      • Predicted values from GSFIT.
    """

    from rtgsfit_vs_gsfit import rtgsfit_pred_meas

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", cfg["pulse_num_write"])
        sens_names = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name']}.PRESHOT:SENS_NAMES").data()
    bp_probe_indices = np.where(np.char.startswith(sens_names, "B_BPPROBE_"))[0]
    long_names = sens_names[bp_probe_indices]
    short_names = np.char.replace(long_names, "B_BPPROBE_", "P")

    rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                  allow_pickle=True).item()
    flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
    coef = rtgsfit_output_dict["coef"][iteration, :]
    coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
    pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        bp_include = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:INCLUDE").data() == 1
        gsfit_bp_meas = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:MVALUE").data()[0, bp_include]
        gsfit_bp_pred = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.BPPROBE:CVALUE").data()[0, bp_include]

    ax.plot(short_names,
            pred_meas_rtgsfit[bp_probe_indices],
            label="RTGSFIT Predicted",
            color='tab:blue')
    ax.plot(short_names,
            gsfit_bp_pred,
            label="GSFIT Predicted",
            color='tab:orange')
    ax.plot(short_names,
            gsfit_bp_meas,
            label="GSFIT Measured",
            color='tab:green')
    ax.set_ylabel("Bp Probe Value [T]")
    ax.tick_params(axis='x', rotation=90)
    ax.legend()
