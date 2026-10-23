"""Replay a FreeGSNKE equilibrium (MAST-U-like machine): no TE network, but needs GSFit and FreeGSNKE.

FreeGSNKE solves a diverted equilibrium, GSFit turns its machine description and synthetic
measurements into a preshot and signals, and RT-GSFit must reconstruct the equilibrium.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import load_settings

pytestmark = pytest.mark.freegsnke


@pytest.fixture(scope="module")
def settings():  # type: ignore[no-untyped-def]
    settings = load_settings("freegsnke")
    section = settings["freegsnke"]
    if (
        not Path(section["machine_dir"]).expanduser().is_dir()
        or not Path(section["currents_path"]).expanduser().is_file()
    ):
        pytest.skip("needs the FreeGSNKE repository (machine pickles and example currents)")
    return settings


@pytest.fixture(scope="module")
def results(settings) -> Results:  # type: ignore[no-untyped-def]
    return replay(settings, verbose=False)


@pytest.fixture(scope="module")
def equilibrium(settings, results):  # type: ignore[no-untyped-def]
    from rtgsfit.database_readers.freegsnke import freegsnke_equilibrium

    return freegsnke_equilibrium(settings)


def test_converges_without_errors(results: Results) -> None:
    assert_array_equal(results["GLOBAL.STATUS"][-5:], 0)
    assert_array_equal(results["FIT.LAPACK_INFO"], 0)
    psi = results["PROFILES_2D.R_Z.PSI"]
    assert np.max(np.abs(psi[-1] - psi[-2])) < 1e-6 * np.max(np.abs(psi[-1]))


def test_matches_freegsnke(results: Results, equilibrium) -> None:  # type: ignore[no-untyped-def]
    # FreeGSNKE's psi is per radian
    assert results["GLOBAL.IP"][-1] == pytest.approx(equilibrium.plasmaCurrent(), rel=1e-2)
    assert results["BOUNDARY.PSI"][-1] == pytest.approx(2.0 * np.pi * equilibrium.psi_bndry, rel=1e-2)
    axis_r, axis_z = equilibrium.magneticAxis()[:2]
    d_r = np.diff(results["PROFILES_2D.R_Z.R"])[0]
    d_z = np.diff(results["PROFILES_2D.R_Z.Z"])[0]
    assert abs(results["GLOBAL.MAG_AXIS.R"][-1] - axis_r) < d_r
    assert abs(results["GLOBAL.MAG_AXIS.Z"][-1] - axis_z) < d_z
    assert results["GLOBAL.XPT_DIVERTED"][-1] == 1
