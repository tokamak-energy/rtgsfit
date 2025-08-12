"""
Contains full example plotting routines that assemble the subplot
components in components.py into complete figures.
"""

import os
import time

import matplotlib.pyplot as plt
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
    from rtgsfit_vs_gsfit.plot.components import bp_probe_line, flux_loop_annotate, \
    flux_loop_line, limiter_scatter, psi_contours

    this_plot_dir = os.path.join(cfg["plots_this_run_dir"], "psi_fluxloop_bp")
    os.makedirs(this_plot_dir, exist_ok=True)

    if iterations is None:
        iterations = np.arange(cfg["n_iters"], dtype=int)

    for iteration in iterations:

        fig = plt.figure(figsize=(12, 7))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

        ax_contour = fig.add_subplot(gs[:, 0])
        ax_flux_loop = fig.add_subplot(gs[0, 1])
        ax_bp_probe = fig.add_subplot(gs[1, 1])

        psi_contours(iteration, ax_contour, cfg)
        limiter_scatter(ax_contour, cfg)
        flux_loop_annotate(ax_contour, cfg)
        flux_loop_line(iteration, ax_flux_loop, cfg)
        bp_probe_line(iteration, ax_bp_probe, cfg)

        fig.suptitle(f"{cfg['pulse_num']},  {cfg['run_name']}" "\n"
                     f"time = {cfg['time']:.1e} s, Iteration {iteration:02d}")

        filename = f"psi_fluxloop_bp_{cfg['pulse_num']}_{cfg['run_name']}_{iteration:02d}.png"
        fig.savefig(os.path.join(this_plot_dir, filename),
                    bbox_inches='tight',
                    dpi=300)
        plt.close("all")
        end_time = time.time()

if __name__ == "__main__":

    from rtgsfit_vs_gsfit import config_loader, rtgsfit_compile_setup

    cfg = config_loader.load_and_prepare_config()

    iterations = [0, 1, 8]
    start_time = time.time()
    psi_fluxloop_bp(cfg, iterations=iterations)
    end_time = time.time()
    print(f"psi_fluxloop_bp took {end_time - start_time:.2e} seconds")

