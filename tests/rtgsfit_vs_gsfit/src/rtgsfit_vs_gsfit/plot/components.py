"""
Provides reusable subplot helper functions that take an `ax` argument
and draw individual plot elements such as contours or annotations.

These functions are then called by the examples.py module.
"""

import os
import typing

import matplotlib.pyplot as plt
import mdsthin
import numpy as np

def psi_contours(iteration: int,
                 ax: plt.Axes,
                 cfg: dict,
                 levels: int = 30,
                 plot_rtgsfit = True,
                 plot_gsfit = True):
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
    rtgsfit_lcfs_n = rtgsfit_output_dict["lcfs_n"][iteration, 0]
    rtgsfit_lcfs_r = rtgsfit_output_dict["lcfs_r"][iteration, :rtgsfit_lcfs_n]
    rtgsfit_lcfs_z = rtgsfit_output_dict["lcfs_z"][iteration, :rtgsfit_lcfs_n]
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        psi_gsfit = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.TWO_D:PSI").data()[0, :, :]
        gsfit_lcfs_n = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.P_BOUNDARY:NBND").data()[0]
        gsfit_lcfs_r = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.P_BOUNDARY:RBND").data()[0, :gsfit_lcfs_n]
        gsfit_lcfs_z = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.P_BOUNDARY:ZBND").data()[0, :gsfit_lcfs_n]

    dot_size = 5
    if plot_rtgsfit:
        ax.contour(r_vec, z_vec, psi_rtgsfit,
                   levels=levels,
                   colors="tab:blue")
        ax.scatter(rtgsfit_lcfs_r, rtgsfit_lcfs_z,
                   label='RTGSFIT LCFS',
                   color='tab:blue',
                   s=dot_size)
    if plot_gsfit:
        ax.contour(r_vec, z_vec, psi_gsfit,
                levels=levels,
                colors="tab:orange")
        ax.scatter(gsfit_lcfs_r, gsfit_lcfs_z,
                   label='GSFIT LCFS',
                   color='tab:orange',
                   s=dot_size)
    ax.set_aspect('equal')
    ax.set_xlabel("R [m]")
    ax.set_ylabel("z [m]")

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

def ivc_j_filled_polygon(ivc_j: np.ndarray, ivc_dict: dict, ax: plt.Axes, cfg: dict):
    """
    Plot filled polygons representing the current density in the IVC.
    """
    if cfg["j_vrange"] is None:
        cfg["j_vrange"] = [ivc_j.min(), ivc_j.max()]
    else:
        cfg["j_vrange"][0] = min(cfg["j_vrange"][0], ivc_j.min())
        cfg["j_vrange"][1] = max(cfg["j_vrange"][1], ivc_j.max())

    for idx in range(len(ivc_dict["r"])):
        # if iteration == 0: continue
        r = ivc_dict["r"][idx]
        z = ivc_dict["z"][idx]
        dr = ivc_dict["dr"][idx]
        dz = ivc_dict["dz"][idx]
        filament_r = [r - 0.5 * dr, r + 0.5 * dr, r + 0.5 * dr, r - 0.5 * dr]
        filament_z = [z - 0.5 * dz, z - 0.5 * dz, z + 0.5 * dz, z + 0.5 * dz]
        color = plt.cm.viridis((ivc_j[idx] - cfg["j_vrange"][0]) /
                               (cfg["j_vrange"][1] - cfg["j_vrange"][0]))
        ax.fill(filament_r, filament_z, color=color, edgecolor='none', linewidth=0.5)
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.viridis,
        norm=plt.Normalize(vmin=cfg["j_vrange"][0], vmax=cfg["j_vrange"][1])
    )
    ax.set_aspect("equal")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Current Density (A/m²)')
    cbar.ax.tick_params(labelsize=8)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()

def ivc_j_rtgsfit(iteration, ax: plt.Axes, cfg: dict):
    """
    Plot the IVC current density from RTGSFIT.
    """

    from rtgsfit_vs_gsfit import rtgsfit_pred_meas

    ivc_dict = np.load(cfg["ivc_dict_path"], allow_pickle=True).item()
    areas = ivc_dict["dr"] * ivc_dict["dz"]

    with open(cfg["coef_names_path"], "r") as f:
        coef_names = [line.strip() for line in f]
    ivc_indices = [i for i, name in enumerate(coef_names) if name.startswith("eig_")]

    rtgsfit_output_dict =np.load(cfg["rtgsfit_output_dict_path"],
                                 allow_pickle=True).item()
    coef = rtgsfit_output_dict["coef"][iteration, :]
    coef_ivc = coef[ivc_indices]

    eigenvector_currents = ivc_dict["current_distributions"] * coef_ivc[:, np.newaxis]
    total_j = np.sum(eigenvector_currents, axis=0) / areas

    ivc_j_filled_polygon(total_j, ivc_dict, ax, cfg)

    ax.set_title("RTGSFIT \n"
                 f"IVC Current = {np.sum(total_j * areas):.2e} A" "\n"
                 f"|IVC Current| = {np.sum(np.abs(total_j) * areas):.2e} A")

def ivc_j_gsfit(iteration, ax: plt.Axes, cfg: dict):
    """
    Plot the IVC current density from GSFIT.
    """

    ivc_dict = np.load(cfg["ivc_dict_path"], allow_pickle=True).item()
    areas = ivc_dict["dr"] * ivc_dict["dz"]
    n_eigs = ivc_dict["current_distributions"].shape[0]

    eigenvalues = np.zeros(n_eigs)
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        for eig_num in range(1, n_eigs + 1):
            eigenvalues[eig_num - 1] = \
                conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.PASSIVES.IVC.DOF:EIG_{eig_num:02d}")[0]
    
    eigenvector_currents = ivc_dict["current_distributions"] * eigenvalues[:, np.newaxis]
    total_j = np.sum(eigenvector_currents, axis=0) / areas

    ivc_j_filled_polygon(total_j, ivc_dict, ax, cfg)

    ax.set_title("GSFIT \n"
                 f"IVC Current = {np.sum(total_j * areas):.2e} A" "\n"
                 f"|IVC Current| = {np.sum(np.abs(total_j) * areas):.2e} A")
    
def get_passive_support_ring_dict(cfg) -> dict:
    """
    Get the passive support ring geometry and area information.
    """

    from shapely.geometry import Polygon

    passive_support_ring_names = ["DIVPSRB", "DIVPSRT", "HFSPSRB", "HFSPSRT"]
    passive_support_ring_dict = {}
    for name in passive_support_ring_names:
        passive_support_ring_dict[name] = {}
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("MAG", cfg["pulse_num"])
        for name in passive_support_ring_names:
            passive_support_ring_dict[name]["r_path"] = \
                conn.get(f"\\MAG::TOP.BEST.ROG.{name}.R_PATH").data()[:4]
            passive_support_ring_dict[name]["z_path"] = \
                conn.get(f"\\MAG::TOP.BEST.ROG.{name}.Z_PATH").data()[:4]
            polygon = []
            for r, z in zip(passive_support_ring_dict[name]["r_path"],
                            passive_support_ring_dict[name]["z_path"]):
                polygon.append((r, z))
            passive_support_ring_dict[name]["area"] = Polygon(polygon).area

    return passive_support_ring_dict

def passive_j_filled_polygon(psr_dict: dict,
                             ax: plt.Axes,
                             cfg: dict):
    """
    Plot the passive support ring current density.
    Using an input passive support ring dictionary.
    """
    max_j = -np.inf
    min_j = np.inf
    for name in psr_dict.keys():
        min_j = min(min_j, psr_dict[name]["pred_j"])
        max_j = max(max_j, psr_dict[name]["pred_j"])

    if cfg["j_vrange"] is None:
        cfg["j_vrange"] = [min_j, max_j]
    else:
        cfg["j_vrange"][0] = min(cfg["j_vrange"][0], min_j)
        cfg["j_vrange"][1] = max(cfg["j_vrange"][1], max_j)

    title = ax.get_title()
    for name in psr_dict.keys():
        color = plt.cm.viridis((psr_dict[name]["pred_j"] - cfg["j_vrange"][0]) /
                               (cfg["j_vrange"][1] - cfg["j_vrange"][0]))
        ax.fill(psr_dict[name]["r_path"], psr_dict[name]["z_path"],
                color=color, edgecolor='none', linewidth=0.5)
        title += "\n" f"{name} current = " \
                 f"{psr_dict[name]['current']:.2e} A, {psr_dict[name]['pred_j']:.2e} A/m²"
    title += "\n"
    ax.set_title(title)

def passive_j_rtgsfit(iteration: int, ax: plt.Axes, cfg: dict):
    """
    Plot the passive support ring current density from RTGSFIT.
    """

    from rtgsfit_vs_gsfit import rtgsfit_pred_meas

    passive_support_ring_dict = get_passive_support_ring_dict(cfg)

    rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                  allow_pickle=True).item()
    flux_norm = rtgsfit_output_dict["flux_norm"][iteration, :]
    coef = rtgsfit_output_dict["coef"][iteration, :]
    coil_curr = rtgsfit_output_dict["coil_curr"][iteration, :]
    pred_meas_rtgsfit = rtgsfit_pred_meas.calc_pred_meas(cfg, flux_norm, coef, coil_curr)

    with open(cfg["meas_names_path"], "r") as f:
        meas_names = [line.strip() for line in f]

    for psr_name in passive_support_ring_dict.keys():
        for i, meas_name in enumerate(meas_names):
            if meas_name == psr_name:
                passive_support_ring_dict[psr_name]["pred_j"] = \
                    pred_meas_rtgsfit[i] / passive_support_ring_dict[psr_name]["area"]
                passive_support_ring_dict[psr_name]["current"] = \
                    pred_meas_rtgsfit[i]

    passive_j_filled_polygon(passive_support_ring_dict, ax, cfg)

def passive_j_gsfit(iteration: int, ax: plt.Axes, cfg: dict):
    """
    Plot the passive support ring current density from GSFIT.
    """

    passive_support_ring_dict = get_passive_support_ring_dict(cfg)

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        rog_names = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.ROG:NAME").data()
        pred_currents = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.CONSTRAINTS.ROG:CVALUE").data()[0]

    for psr_name in passive_support_ring_dict.keys():
        for i, rog_name in enumerate(rog_names):
            if psr_name in rog_name:
                passive_support_ring_dict[psr_name]["current"] = pred_currents[i]
                passive_support_ring_dict[psr_name]["pred_j"] = \
                    pred_currents[i] / passive_support_ring_dict[psr_name]["area"]

    passive_j_filled_polygon(passive_support_ring_dict, ax, cfg)

def ovc_area(cfg: dict) -> np.ndarray:
    """
    Get the area of the OVC filaments.
    RTGSFIT and GSFIT assume a uniform current in the OVC, 
    so we only need the area to determine the total current.
    """

    with mdsthin.Connection('smaug') as conn:
        conn.openTree("ELMAG", cfg["pulse_num"])
        vessel_d_r = conn.get(f"\\ELMAG::TOP.BEST.VESSEL:DR").data()
        vessel_d_z = conn.get(f"\\ELMAG::TOP.BEST.VESSEL:DZ").data()
        fils2passive = conn.get(f"\\ELMAG::TOP.BEST.VESSEL:FILS2PASSIVE").data()
        passive_names = conn.get(f"\\ELMAG::TOP.BEST.VESSEL:PASSIVE_NAME").data()
        for i_passive, passive_name in enumerate(passive_names):
            if passive_name == "OVC":
                indices = fils2passive[:, i_passive].astype(int) == True
                d_r = vessel_d_r[indices]
                d_z = vessel_d_z[indices]
                areas = d_r * d_z

    return np.sum(areas)

def ovc_current_rtgsfit(iteration: int, cfg: dict) -> float:
    """
    Get the Total OVC current from RTGSFIT.
    """
    rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                  allow_pickle=True).item()
    coef = rtgsfit_output_dict["coef"][iteration, :]
    with open(cfg["coef_names_path"], "r") as f:
        coef_names = [line.strip() for line in f]

    for i, name in enumerate(coef_names):
        if name == "OVC":
            return coef[i] * ovc_area(cfg)

def ovc_current_gsfit(cfg: dict) -> float:
    """
    Get the Total OVC current from GSFIT.
    """
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("GSFIT", cfg["pulse_num_write"])
        ovc_current = conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.PASSIVES.OVC.DOF:CONSTANT_J").data()[0]

    return ovc_current * ovc_area(cfg)

def ivc_eigenvalue_bar_chart(ax: plt.Axes,
                             iteration: int,
                             cfg: dict):
    """
    Create a bar chart comparing the IVC eigenvalues from RTGSFIT and GSFIT.
    """

    import pandas as pd

    def ivc_rtgsfit_row(iteration: int,
                        cfg: dict) -> np.ndarray:
        rtgsfit_output_dict = np.load(cfg["rtgsfit_output_dict_path"],
                                    allow_pickle=True).item()
        coef = rtgsfit_output_dict["coef"][iteration, :]
        with open(cfg["coef_names_path"], "r") as f:
            coef_names = [line.strip() for line in f]
        ivc_indices = [i for i, name in enumerate(coef_names) if name.startswith("eig_")]
        return coef[ivc_indices]
    
    def ivc_gsfit_row(cfg: dict) -> np.ndarray:
        ivc_dict = np.load(cfg["ivc_dict_path"], allow_pickle=True).item()
        n_eigs = ivc_dict["current_distributions"].shape[0]
        eigenvalues = np.zeros(n_eigs)
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("GSFIT", cfg["pulse_num_write"])
            for eig_num in range(1, n_eigs + 1):
                eigenvalues[eig_num - 1] = \
                    conn.get(f"\\GSFIT::TOP.{cfg['run_name']}.PASSIVES.IVC.DOF:EIG_{eig_num:02d}")[0]
        return eigenvalues
    
    rtgsfit_vals = ivc_rtgsfit_row(iteration, cfg)
    gsfit_vals = ivc_gsfit_row(cfg)

    n_eigs = len(rtgsfit_vals)
    cols = [f"eig_{i:02d}" for i in range(1, n_eigs + 1)]

    df = pd.DataFrame([rtgsfit_vals, gsfit_vals], 
                      index=["RTGSFIT", "GSFIT"],
                      columns=cols).T
        
    df.plot(kind="bar", ax=ax)
    ax.set_ylabel("Eigenvalue")
    ax.set_xlabel("Eigenvector Index")
    ax.set_title("IVC Eigenvalues: RTGSFIT vs GSFIT")
