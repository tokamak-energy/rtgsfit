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

class Psi:
    def __init__(self,
                 cfg: dict):
        self.cfg = cfg

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_replay"])
            self.psi_py = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.TWO_D:PSI").data()
            self.r_vec = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.TWO_D:RGRID").data()
            self.z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.TWO_D:ZGRID").data()
            self.rbnd_py = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.P_BOUNDARY:RBND").data()
            self.zbnd_py = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.P_BOUNDARY:ZBND").data()
            self.nbnd_py = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.P_BOUNDARY:NBND").data()

        self.dot_size = 5

    def contour(self,
                iteration: int,
                ax: plt.Axes,
                levels: int = 30,
                plot_rtgsfit = True,
                plot_gsfit = True):
        """
        Draw the ψ [Wb] contours on the given axis.
        """

        if plot_rtgsfit:
            ax.contour(self.r_vec, self.z_vec,
                       self.psi_py[iteration, :, :],
                       levels=levels,
                       colors="tab:blue")
            ax.scatter(self.rbnd_py[iteration, :self.nbnd_py[iteration]],
                       self.zbnd_py[iteration, :self.nbnd_py[iteration]],
                       label='RTGSFIT LCFS',
                       color='tab:blue',
                       s=self.dot_size)
        ax.set_aspect('equal')
        ax.set_xlabel("R [m]")
        ax.set_ylabel("z [m]")

class Limiter:
    def __init__(self,
                 cfg: dict):
        self.cfg = cfg
        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_preshot"])
            self.limiter_r = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_preshot']}.PRESHOT:LIMIT_R").data()
            self.limiter_z = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_preshot']}.PRESHOT:LIMIT_Z").data()

    def scatter(self,
                ax: plt.Axes):
        """
        Draw the limiter points on the given axis.
        """
        ax.scatter(self.limiter_r, self.limiter_z,
                   color="tab:green")

class FluxLoops:
    def __init__(self,
                 cfg: dict):
        self.cfg = cfg

        with mdsthin.Connection('smaug') as conn:
            #  \RTGSFIT::TOP.TEST01.CONSTRAINTS.FLOOP:CVALUE
            conn.openTree("RTGSFIT", cfg["pulse_num_replay"])
            self.pred_py = \
                conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.FLOOP:CVALUE").data()
            self.meas_py = \
                conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.FLOOP:MVALUE").data()
            self.name = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.FLOOP:NAME").data()

    def line(self,
             iteration: int,
             ax: plt.Axes):
        """
        Draw the flux loop values on the given axis.
        """

        ax.plot(self.name,
                self.pred_py[iteration],
                label="Predicted",
                color='tab:blue')
        ax.plot(self.name,
                self.meas_py[iteration],
                label="Measured",
                color='tab:green')
        ax.set_ylabel("Flux Loop Value [Wb]")
        ax.tick_params(axis='x', rotation=90)
        ax.legend()

class BpProbes:
    def __init__(self,
                 cfg: dict):
        self.cfg = cfg

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_replay"])
            self.pred_py = \
                conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.BPPROBE:CVALUE").data()
            self.meas_py = \
                conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.BPPROBE:MVALUE").data()
            self.name = conn.get(f"\\RTGSFIT::TOP.{cfg['run_name_replay']}.CONSTRAINTS.BPPROBE:NAME").data()

    def line(self,
             iteration: int,
             ax: plt.Axes):
        """
        Draw the Bp probe values on the given axis.
        """

        ax.plot(self.name,
                self.pred_py[iteration],
                label="Predicted",
                color='tab:blue')
        ax.plot(self.name,
                self.meas_py[iteration],
                label="Measured",
                color='tab:green')
        ax.set_ylabel("Bp Probe Value [T]")
        ax.tick_params(axis='x', rotation=90)
        ax.legend()