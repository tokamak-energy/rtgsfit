"""
Contains full example plotting routines that assemble the subplot
components in components.py into complete figures.
"""

import os
import time

import matplotlib.pyplot as plt
import mdsthin
import numpy as np

def psi_fluxloop_bp(cfg: dict,
                    iterations = None):
    """
    Create a 2×2 subplot layout showing ψ [Wb] contours, flux loop ψ profiles,
    and Bp probe measurements.

    Layout:
      • Top-left & bottom-left (merged into one large axis):
          Contour plot of poloidal flux ψ [Wb] from both GSFIT and RTGSFIT,
          including limiter points and annotated positions of flux loops and Bp probes.

      • Top-right axis:
          Line plot of ψ [Wb] along the flux loop positions.

      • Bottom-right axis:
          Line plot of poloidal magnetic field Bp [T] along the Bp probe positions.
    """

    import matplotlib.gridspec as gridspec
    from pcs_vs_python.plot.components import BpProbes, FluxLoops, Limiter, Psi

    this_plot_dir = os.path.join(cfg["plot_dir_this_run"], "psi_fluxloop_bp")
    os.makedirs(this_plot_dir, exist_ok=True)

    iterations  = [0, 1, 2, 10, 100, 200]

    psi_class = Psi(cfg)
    limiter_class = Limiter(cfg)
    flux_loops_class = FluxLoops(cfg)
    bp_probe_class = BpProbes(cfg)

    for iteration in iterations:

        fig = plt.figure(figsize=(12, 7))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

        ax_contour = fig.add_subplot(gs[:, 0])
        ax_flux_loop = fig.add_subplot(gs[0, 1])
        ax_bp_probe = fig.add_subplot(gs[1, 1])

        psi_class.contour(iteration, ax_contour)
        limiter_class.scatter(ax_contour)
        flux_loops_class.line(iteration, ax_flux_loop)
        bp_probe_class.line(iteration, ax_bp_probe)

        fig.suptitle(f"{cfg['pulse_num']},  run_name_replay: {cfg['run_name_replay']}"
                     f", run_name_preshot: {cfg['run_name_preshot']}" "\n"
                     f"time = {cfg['time'][iteration]:.1e} s, Iteration {iteration:02d}")

        filename = f"psi_fluxloop_bp_{cfg['pulse_num']}_{cfg['run_name_replay']}_{iteration:02d}.png"
        fig.savefig(os.path.join(this_plot_dir, filename),
                    bbox_inches='tight',
                    dpi=300)
        plt.close("all")

def just_psi(cfg: dict,
             iterations = None):
    """
    Create a single plot showing ψ [Wb] contours.

    Layout:
      • Single axis:
          Contour plot of poloidal flux ψ [Wb] from both GSFIT and RTGSFIT,
          including limiter points and annotated positions of flux loops and Bp probes.
    """

    from pcs_vs_python.plot.components import Psi

    this_plot_dir = os.path.join(cfg["plot_dir_this_run"], "just_psi")
    os.makedirs(this_plot_dir, exist_ok=True)

    # iterations  = [0, 1, 2, 10, 100, 200]
    iterations = np.arange(0, 68000, 1000)

    psi_class = Psi(cfg, plasma_boundary=False)

    for iteration in iterations:

        fig, ax = plt.subplots()

        psi_class.contour(iteration, ax)

        fig.suptitle(f"{cfg['pulse_num']},  run_name_replay: {cfg['run_name_replay']}" "\n"
                     f"Iteration {iteration:02d}")

        filename = f"just_psi_{cfg['pulse_num']}_{cfg['run_name_replay']}_{iteration:02d}.png"
        fig.savefig(os.path.join(this_plot_dir, filename),
                    bbox_inches='tight',
                    dpi=300)
        plt.close("all")

if __name__ == "__main__":

    from pcs_vs_python import config_loader

    # start_time = time.time()
    # cfg = config_loader.load_and_prepare_config()
    # psi_fluxloop_bp(cfg)
    # end_time = time.time()
    # print(f"Plotting completed in {end_time - start_time:.2e} seconds.")

    start_time = time.time()
    cfg = config_loader.load_and_prepare_config()
    cfg["pulse_num_replay"] = 30_013_343
    cfg["run_name_replay"] = "RUN01"
    cfg["plot_dir_this_run"] = os.path.join(cfg["plot_dir"], f'{cfg["pulse_num"]}_{cfg["run_name_replay"]}')
    just_psi(cfg)
    end_time = time.time()
    print(f"Plotting completed in {end_time - start_time:.2e} seconds.")
