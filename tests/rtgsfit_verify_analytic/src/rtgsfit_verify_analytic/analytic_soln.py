"""
This file generates the measurements which RTGSFIT uses to fit the plasma equilibrium.
We generate them using our analytic solution.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.special

from rtgsfit_verify_analytic import cnst

def analytic_psi(r, z, r_o, z_o, rho_b, total_current):
    """
    Computes the analytic solution ψ(r, z) to the Grad-Shafranov equation 
    in cylindrical coordinates, assuming radial symmetry.

    The solution is defined piecewise in terms of the radial coordinate ρ:

        ρ = sqrt((r - r₀)^2 + (z - z₀)^2)

    For ρ ≤ ρ_b:
        ψ(ρ) = Δψ * J₀(a* * ρ / ρ_b)

    For ρ > ρ_b:
        ψ(ρ) = -Δψ * a* * J₁(a*) * ln(ρ / ρ_b)

    where:
        - Δψ = ψ₀ - ψ_b
        - ψ₀ is the value of the magnetic flux at the o-point
        - ψ_b is the value of the magnetic flux at the boundary
        - a* is the first positive zero of J₀(x) (~2.4048255576957727)
        - ρ_b is the boundary radius
        - (r₀, z₀) is the coordinate of the o-point
        - J₀ and J₁ are Bessel functions of the first kind (orders 0 and 1)

    Parameters:
        r (float): Radial coordinate
        z (float): vertical coordinate
        r_o (float): Center of the plasma in the radial direction
        z_o (float): Center of the plasma in the vertical direction
        rho_b (float): Boundary radius
        delta_psi (float): Difference in magnetic flux between the o-point and the boundary

    Returns:
        float: Value of the analytic solution ψ at the point (r, z)
    """

    rho = np.sqrt((r-r_o)**2 + (z-z_o)**2)
    rho_log = rho + (rho == 0)  # Avoid log(0)
    d = np.sqrt(r_o**2 + z_o**2)
    psi_o = cnst.MU0_BAR * total_current * r_o * (np.log(d / rho_b) + 1 / (cnst.A_J1_A))
    psi_b = cnst.MU0_BAR * total_current * r_o * np.log(d / rho_b)
    delta_psi = psi_o - psi_b

    soln = delta_psi * scipy.special.j0(cnst.A_STAR * rho / rho_b) * (rho <= rho_b) \
         - delta_psi * cnst.A_J1_A * np.log(rho_log / rho_b) * (rho > rho_b) \
         + psi_b

    return soln

def analytic_b_theta(r, z, r_o, z_o, rho_b, total_current):
    """
    Computes the analytic solution B_θ(r, z) to the Grad-Shafranov equation 
    in cylindrical coordinates, assuming radial symmetry.

    The solution is defined piecewise in terms of the radial coordinate ρ:

        ρ = sqrt((r - r₀)^2 + (z - z₀)^2)

    For ρ ≤ ρ_b:
        B_θ(ρ) = -Δψ * a* * J₁(a* * ρ / ρ_b) / (ρ_b * r)

    For ρ > ρ_b:
        B_θ(ρ) = -Δψ * a* * J₁(a*) / (ρ * r)

    where:
        - Δψ = ψ₀ - ψ_b
        - ψ₀ is the value of the magnetic flux at the o-point
        - ψ_b is the value of the magnetic flux at the boundary
        - a* is the first positive zero of J₀(x) (~2.4048255576957727)
        - ρ_b is the boundary radius
        - (r₀, z₀) is the coordinate of the o-point
        - J₁ is Bessel function of the first kind (order 1)

    Parameters:
        r (float): Radial coordinate
        z (float): vertical coordinate
        r_o (float): Center of the plasma in the radial direction
        z_o (float): Center of the plasma in the vertical direction
        rho_b (float): Boundary radius
        delta_psi (float): Difference in magnetic flux between the o-point and the boundary

    Returns:
        float: Value of the analytic solution B_θ at the point (r, z)
    """

    rho = np.sqrt((r-r_o)**2 + (z-z_o)**2)
    rho_inv = rho + (rho == 0)  # Avoid division by zero
    d = np.sqrt(r_o**2 + z_o**2)
    psi_o = cnst.MU0_BAR * total_current * r_o * (np.log(d / rho_b) + 1 / (cnst.A_J1_A))
    psi_b = cnst.MU0_BAR * total_current * r_o * np.log(d / rho_b)
    delta_psi = psi_o - psi_b

    b_theta = -delta_psi * cnst.A_STAR * scipy.special.j1(cnst.A_STAR * rho / rho_b) / (rho_b * r) * (rho <= rho_b) \
            - delta_psi * cnst.A_J1_A / (rho_inv * r) * (rho > rho_b)
    
    return b_theta


def add(a, b):
    return a + b

if __name__ == "__main__":
    # Example usage
    r_o = 1.5  # Center of the plasma in the radial direction
    z_o = 0.5  # Center of the plasma in the vertical direction
    rho_b = 0.5  # Boundary radius
    r_min = r_o - 3 * rho_b
    r_max = r_o + 3 * rho_b
    z_min = z_o - 3 * rho_b
    z_max = z_o + 3 * rho_b
    n_r = 1024
    n_z = 1024
    delta_psi = 1.0  # Difference in magnetic flux between the o-point and the boundary

    r_vec= np.linspace(r_min, r_max, n_r)  # Radial coordinates
    z_vec = np.linspace(z_min, z_max, n_z)   # Vertical coordinates
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)

    psi_values = analytic_psi(r_grid, z_grid, r_o, z_o, rho_b, delta_psi)

    # Make a contour plot of psi
    fig, ax = plt.subplots()
    contour = ax.contour(r_grid, z_grid, psi_values, levels=50)
    cbar = fig.colorbar(contour)
    cbar.set_label('Magnetic Flux (ψ)')
    ax.set_xlabel('Radial Coordinate (r)')
    ax.set_ylabel('Vertical Coordinate (z)')
    ax.set_title('Contour Plot of Magnetic Flux (ψ)')
    ax.set_aspect('equal', adjustable='box')
    fig.savefig('plots/analytic_psi_contour.png', dpi=300, bbox_inches='tight')

    # Make line plot along the midplane
    fig, ax = plt.subplots()
    ax.plot(r_vec, psi_values[n_z // 2, :])
    ax.set_xlabel('Radial Coordinate (r)')
    ax.set_ylabel('Magnetic Flux (ψ)')
    ax.set_title(f'Magnetic Flux (ψ) along z = {z_vec[n_z//2]:.2f}')
    fig.savefig('plots/analytic_psi_midplane.png', dpi=300, bbox_inches='tight')