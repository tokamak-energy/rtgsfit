"""COCOS 13 normalisation checks independent of result-export conversion."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.constants import MU_0
from rtgsfit.database_readers.analytic import analytic_plasma
from rtgsfit.database_readers.analytic import greens
from rtgsfit.database_readers.analytic.solution import AnalyticPlasma
from rtgsfit.preshot import Preshot
from rtgsfit.replay import read_signals
from rtgsfit.settings import load_settings


def test_filament_total_flux_and_field_derivatives() -> None:
    r, z, source_r, source_z = 10.4, 0.3, 10.0, 0.0
    rho = np.hypot(r - source_r, z - source_z)
    expected = -MU_0 * source_r * np.log(rho / np.hypot(source_r, source_z))
    assert greens.psi(r, z, source_r, source_z) == pytest.approx(expected)
    step = 1e-5
    d_r = (greens.psi(r + step, z, source_r, source_z) - greens.psi(r - step, z, source_r, source_z)) / (2 * step)
    d_z = (greens.psi(r, z + step, source_r, source_z) - greens.psi(r, z - step, source_r, source_z)) / (2 * step)
    assert_allclose(greens.b_r(r, z, source_r, source_z), -d_z / (2 * np.pi * r), rtol=1e-8)
    assert_allclose(greens.b_z(r, z, source_r, source_z), d_r / (2 * np.pi * r), rtol=1e-8)


def test_analytic_preshot_and_raw_measurements_use_total_flux(analytic_preshot: Preshot) -> None:
    preshot = analytic_preshot
    assert preshot["INFO:PSI_TO_WB"] == 1.0
    assert_allclose(
        preshot["POISSON:RHS_SCALE"],
        2 * np.pi * MU_0 * preshot["GRID:R_FLAT"] * preshot["GRID:D_Z"] ** 2,
        rtol=1e-14,
    )
    settings = load_settings("analytic")
    signals = read_signals(settings, preshot)
    from rtgsfit.database_readers.analytic import flux_loop_geometry

    r, z = flux_loop_geometry(settings)
    plasma = analytic_plasma(settings)
    assert_allclose(signals.meas[0, : len(r)], plasma.psi(r, z), rtol=1e-14)


def test_reference_total_flux_satisfies_poisson_equation() -> None:
    plasma = analytic_plasma(load_settings("analytic"))
    r, z, step = 10.2, 0.1, 1e-3
    centre = plasma.psi(r, z)
    laplacian = (
        plasma.psi(r + step, z)
        + plasma.psi(r - step, z)
        + plasma.psi(r, z + step)
        + plasma.psi(r, z - step)
        - 4 * centre
    ) / step**2
    assert_allclose(laplacian, -2 * np.pi * MU_0 * r * plasma.j_phi(r, z), rtol=1e-5)


def test_zero_current_reference_has_no_axis_and_no_current() -> None:
    plasma = AnalyticPlasma(r_0=0.5, z_0=0.0, rho_b=0.15, ip=0.0)
    r, z = np.array([0.5, 0.55, 0.8]), np.zeros(3)  # axis, inside, outside
    with np.errstate(all="raise"):
        assert_array_equal(plasma.psi_norm(r, z), [np.nan, np.nan, 1.0])
        assert_array_equal(plasma.j_phi(r, z), 0.0)
