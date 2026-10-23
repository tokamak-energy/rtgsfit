"""Analytic solution of the large aspect ratio Grad-Shafranov equation used to verify RT-GSFit.

A circular plasma of radius `rho_b` centred on (r_0, z_0), carrying `ip`, in the large aspect ratio
limit where the Grad-Shafranov equation becomes  Laplacian(psi) = -mu0 R j_phi  with R ~ r_0. With
j_phi proportional to (1 - psi_norm) R (the p' basis function), the solution is a Bessel function
inside the plasma and a vacuum (logarithmic) solution outside:

    psi = psi_b + (psi_0 - psi_b) J0(a rho / rho_b)                       rho <= rho_b
    psi = psi_b - (psi_0 - psi_b) a J1(a) ln(rho / rho_b)                 rho >  rho_b

where `a` is the first zero of J0, and

    psi_0 = (mu0 / 2 pi) ip r_0 (ln(d / rho_b) + 1 / (a J1(a))),  psi_b = (mu0 / 2 pi) ip r_0 ln(d / rho_b)

with d = |(r_0, z_0)|. See the manual for the derivation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.special

from rtgsfit.constants import MU_0

FloatArray = npt.NDArray[np.float64]

A_STAR = float(scipy.special.jn_zeros(0, 1)[0])
"""First zero of J0."""
J1_A_STAR = float(scipy.special.j1(A_STAR))


@dataclass(frozen=True)
class AnalyticPlasma:
    r_0: float
    z_0: float
    rho_b: float
    ip: float

    @property
    def mu0_bar(self) -> float:
        return MU_0 / (2.0 * np.pi)

    @property
    def boundary_psi(self) -> float:
        return self.mu0_bar * self.ip * self.r_0 * np.log(np.hypot(self.r_0, self.z_0) / self.rho_b)

    @property
    def axis_psi(self) -> float:
        return self.boundary_psi + self.mu0_bar * self.ip * self.r_0 / (A_STAR * J1_A_STAR)

    def rho(self, r: FloatArray, z: FloatArray) -> FloatArray:
        return np.asarray(np.hypot(np.subtract(r, self.r_0), np.subtract(z, self.z_0)))

    def psi(self, r: FloatArray, z: FloatArray) -> FloatArray:
        rho = self.rho(r, z)
        delta_psi = self.axis_psi - self.boundary_psi
        inside = delta_psi * scipy.special.j0(A_STAR * rho / self.rho_b)
        outside = -delta_psi * A_STAR * J1_A_STAR * np.log(np.where(rho == 0.0, 1.0, rho) / self.rho_b)
        return np.asarray(np.where(rho <= self.rho_b, inside, outside) + self.boundary_psi)

    def psi_norm(self, r: FloatArray, z: FloatArray) -> FloatArray:
        rho = self.rho(r, z)
        psi_norm = (self.axis_psi - self.psi(r, z)) / (self.axis_psi - self.boundary_psi)
        return np.asarray(np.where(rho <= self.rho_b, psi_norm, 1.0))

    def b_theta(self, r: FloatArray, z: FloatArray) -> FloatArray:
        """Poloidal field, positive anticlockwise around the plasma in the (R, Z) plane."""
        rho = self.rho(r, z)
        delta_psi = self.axis_psi - self.boundary_psi
        inside = -delta_psi * A_STAR * scipy.special.j1(A_STAR * rho / self.rho_b) / (self.rho_b * np.asarray(r))
        outside = -delta_psi * A_STAR * J1_A_STAR / (np.where(rho == 0.0, 1.0, rho) * np.asarray(r))
        return np.asarray(np.where(rho <= self.rho_b, inside, outside))

    def j_phi(self, r: FloatArray, z: FloatArray) -> FloatArray:
        a_raw = self.mu0_bar * A_STAR * self.ip * self.r_0 / (J1_A_STAR * self.rho_b**2)
        return np.asarray(a_raw * (1.0 - self.psi_norm(r, z)) / (np.asarray(r) * MU_0))


def b_theta_from_psi(psi_2d: FloatArray, grid_r: FloatArray, grid_z: FloatArray, r_0: float, z_0: float) -> FloatArray:
    """Poloidal field around (r_0, z_0) from a `[N_Z, N_R]` psi, by finite differences."""
    d_psi_d_z, d_psi_d_r = np.gradient(psi_2d, grid_z, grid_r)
    r_2d, z_2d = np.meshgrid(grid_r, grid_z)
    b_r = -d_psi_d_z / r_2d
    b_z = d_psi_d_r / r_2d
    theta = np.arctan2(z_2d - z_0, r_2d - r_0)
    return np.asarray(-b_r * np.sin(theta) + b_z * np.cos(theta))
