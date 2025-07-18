"""
This module provides functions to calculate the mutual inductance between coaxial circular wires in the large
aspect-ratio limit.
The primary aim it to provide the Green's matrices needed for the RTGSFIT code.
"""

import numpy as np
import matplotlib.pyplot as plt

MU_0 = 4e-7 * np.pi # Permeability of free space in Tm/A

def mutual_inductance_psi(
    r: float,
    z: float,
    r_prime: float,
    z_prime: float) -> float:
    """
    This function computes the mutual inductance between two coaxial circular wires
    located at (R, z) and (R', z') using the formula:

        M = -(μ₀ / (2 * π)) * ln(ρ / d) * R'
          = -2e-7 * ln(ρ / d) * R'

    where

        - we assume |R - R'| << R'
        - ρ = √[(R - R')² + (z - z')²]
        - d = √[R'² + z'²]
        - μ₀ is the permeability of free space.

    Note 1: The mutual inductance --> infinity as the distance between the wires approaches zero.
            When calling this function ensure (r, z) is never equal to (r_prime, z_prime).

    Note 2: Our definition of ψ is that it satisfies the following equations:
        ∂ψ/∂z = - R * B_R
        ∂ψ/∂R = R * B_z
        and we choose the integration constant such that ψ = 0 at (R, z) = (0, 0).

    Parameters:
        r (float): Radial coordinate of the first wire.
        z (float): Axial coordinate of the first wire.
        r_prime (float): Radial coordinate of the second wire.
        z_prime (float): Axial coordinate of the second wire.

    Returns:
        float: The mutual inductance between the two wires.
    """

    rho = np.sqrt((r - r_prime) ** 2 + (z - z_prime) ** 2)
    dist = np.sqrt(r_prime ** 2 + z_prime ** 2)

    return -(MU_0 / (2 * np.pi)) * np.log(rho / dist) * r_prime

def self_inductance_linear_avg(r: float, z: float, delta: float) -> float:
    """
    Average the mutual inductance for the case where (R, z) = (R', z') over a small 
    linear distance Δ.

    The mutual inductance is averaged over a small interval using the formula:

        <M> = ∫ M(R, z, R', z) dR' / Δ
            = ∫ M(R, z, R, z') dz' / Δ
            = -(μ₀ / 2π) [ln(Δ / (2d)) - 1] * R
            = -2e-7 * [ln(Δ / (2d)) - 1] * R

    where:
        - the lower limit of integration is R' = R - Δ/2 (or z' = z - Δ/2 in the second case),
        - the upper limit of integration is R' = R + Δ/2 (or z' = z + Δ/2 in the second case),
        - μ₀ is the permeability of free space,
        - Δ is the averaging interval,
        - d = √(r² + z²)

    Parameters:
        r (float): Radial coordinate of the observation point.
        z (float): Axial coordinate of the observation point.
        delta (float): Distance over which to average the mutual inductance.

    Returns:
        float: The averaged self-inductance in Tm²/A.
    """

    dist = np.sqrt(r ** 2 + z ** 2)

    return -(MU_0 / (2 * np.pi)) * (np.log(0.5 * delta / dist) - 1) * r

def mutual_inductance_br(
    r: np.ndarray,
    z: np.ndarray,
    r_prime: float,
    z_prime: float) -> np.ndarray:
    """
    Calculate the mutual inductance between two circular coaxial wires based on the radial component
    of the magnetic field, given by -(∂M/∂z) / R.

    The expression used is:

        -(∂M/∂z) / R = (μ_0 / (2 * π)) * {sin(θ) / ρ} * R'/R
                     = (μ_0 / (2 * π)) * {(z - z') / [(R - R')² + (z - z')²]} * R'/R

    where:
        sin(θ) = (z - z') / √[(R - R')² + (z - z')²]
        ρ      = √[(R - R')² + (z - z')²]

    Parameters:
        r_vec (np.ndarray): Array of radial coordinates for the grid points.
        z_vec (np.ndarray): Array of axial coordinates for the grid points.
        r_prime (float): Radial coordinate of the source wire.
        z_prime (float): Axial coordinate of the source wire.

    Returns:
        np.ndarray: Array of mutual inductance values computed at each grid point.
    """
    
    return MU_0 / (2 * np.pi) * (z - z_prime) / ((r - r_prime) ** 2 + (z - z_prime) ** 2) * r_prime / r

def mutual_inductance_bz(
    r: np.ndarray,
    z: np.ndarray,
    r_prime: float,
    z_prime: float) -> np.ndarray:
    """
    Calculate the mutual inductance between two circular coaxial wires based on the radial component
    of the magnetic field, given by (∂M/∂R) / R.

    The expression used is:

        +(∂M/∂R) / R = -(μ_0 / (2 * π)) * {cos(θ) / ρ} * R'/R
                     = -(μ_0 / (2 * π)) * {(R - R') / [(R - R')² + (z - z')²]} * R'/R

    where:
        cos(θ) = (R - R') / √[(R - R')² + (z - z')²]
        ρ      = √[(R - R')² + (z - z')²]

    Parameters:
        r_vec (np.ndarray): Array of radial coordinates for the grid points.
        z_vec (np.ndarray): Array of axial coordinates for the grid points.
        r_prime (float): Radial coordinate of the source wire.
        z_prime (float): Axial coordinate of the source wire.

    Returns:
        np.ndarray: Array of mutual inductance values computed at each grid point.
    """
    
    return -MU_0 / (2 * np.pi) * (r - r_prime) / ((r - r_prime) ** 2 + (z - z_prime) ** 2) * r_prime / r

if __name__ == "__main__":

    r_min = 0.5
    r_max = 3.5
    z_min = -1.5
    z_max = 1.5
    n_r, n_z = 32, 32
    r_v = np.linspace(r_min, r_max, n_r + 1) # Cell vertices
    z_v = np.linspace(z_min, z_max, n_z + 1) # Cell vertices
    r_c = 0.5 * (r_v[:-1] + r_v[1:]) # Cell centers
    z_c = 0.5 * (z_v[:-1] + z_v[1:]) # Cell centers
    rho_b = 1
    r_o = 2  # o-point coordinate
    z_o = 0.0  # o-point coordinate

    r_v_grid, z_v_grid = np.meshgrid(r_v, z_v)
    r_c_grid, z_c_grid = np.meshgrid(r_c, z_c)

    psi = np.zeros((n_z + 1, n_r + 1), dtype=np.float64) # Psi a cell vertices
    rho_v_grid = np.sqrt((r_v_grid - r_o) ** 2 + (z_v_grid - z_o) ** 2)
    rho_c_grid = np.sqrt((r_c_grid - r_o) ** 2 + (z_c_grid - z_o) ** 2)

    j_phi = rho_c_grid <= rho_b
    plasma_current = 1e6
    j_phi = plasma_current * j_phi / np.sum(j_phi)

    for i, z in enumerate(z_v):
        for j, r in enumerate(r_v):
            psi[i, j] = np.sum(mutual_inductance_psi(r, z, r_c_grid, z_c_grid, rho_b * 1.5) * j_phi)

    fig, ax = plt.subplots()
    ax.scatter(r_v_grid, z_v_grid, c=psi.flatten())
    ax.set_xlabel('Radial Coordinate (R)')
    ax.set_ylabel('Axial Coordinate (Z)')
    ax.set_title('Mutual Inductance (Psi) Distribution')
    plt.colorbar(ax.collections[0], ax=ax, label='Mutual Inductance (Psi)')
    fig.savefig('mutual_inductance_psi.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    ax.scatter(r_c_grid, z_c_grid, c=j_phi.flatten())
    ax.set_xlabel('Radial Coordinate (R)')
    ax.set_ylabel('Axial Coordinate (Z)')
    ax.set_title('Plasma Current Density Distribution')
    plt.colorbar(ax.collections[0], ax=ax, label='Plasma current (A)')
    fig.savefig('plasma_current_density.png', dpi=300)
    plt.close()

    # Line plot psi along the radial coordinate at z = 0
    psi_at_z0 = psi[n_z // 2, :]  # Take the middle
    print("z_v[n_z // 2] = ", z_v[n_z // 2])
    fig, ax = plt.subplots()
    ax.plot(r_v, psi[n_z//2, :])
    ax.set_xlabel('Radial Coordinate (R)')
    ax.set_ylabel('Mutual Inductance (Psi)')
    ax.set_title('Mutual Inductance (Psi) at Z = 0')
    fig.savefig('mutual_inductance_psi_radial.png', dpi=300)
    plt.close()

    