"""Offline-only products, exported conventions and exact diagnostic bookkeeping."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from rtgsfit.library import describe_status
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results


def test_fit_diagnostics_close(analytic_results: Results) -> None:
    results = analytic_results
    assert_allclose(results["FIT.RESIDUAL"], results["FIT.PREDICTION"] - results["FIT.TARGET"])
    assert_allclose(
        results["FIT.CHI_SQ"].sum(axis=1),
        results["CONSTRAINTS.CHI_SQ_TOTAL"],
        rtol=1e-8,
        atol=1e-12,
    )
    magnetic = sum(results[f"CONSTRAINTS.{name}.CHI_SQ"] for name in ("FLUX_LOOP", "BP_PROBE", "ROGOWSKI"))
    assert_allclose(magnetic, results["CONSTRAINTS.CHI_SQ_MAG"], rtol=1e-8, atol=1e-12)
    assert_allclose(
        results["CONSTRAINTS.CHI_SQ_TOTAL"],
        magnetic + results["CONSTRAINTS.CHI_SQ_REG"],
        rtol=1e-8,
        atol=1e-12,
    )
    assert_allclose(
        results["FIT.FF_PRIME"],
        results["FIT.DOF_VALUES"][:, 1] * results["PRESHOT.INFO.MU0"],
    )


def test_weighted_target_is_distinct_from_physical_dof(analytic_results: Results, analytic_preshot: Preshot) -> None:
    results, preshot = analytic_results, analytic_preshot
    sizes = preshot.sizes
    measured = np.concatenate(
        [results[f"CONSTRAINTS.{name}.ALL.MEASURED"] for name in ("FLUX_LOOP", "BP_PROBE", "ROGOWSKI")],
        axis=1,
    )
    coil = results["CONSTRAINTS.PF_CURRENT.ALL.MEASURED"] @ preshot["CONSTRAINT:G_COIL"].T
    n_constraint = sizes["N_CONSTRAINT"]
    assert_allclose(
        results["FIT.TARGET"][:, :n_constraint],
        (measured - coil) * preshot["FIT:WEIGHT"][:n_constraint],
    )
    assert_allclose(results["FIT.TARGET"][:, n_constraint:], 0.0)
    assert sizes["N_FIT_TARGET"] != sizes["N_FIT_DOF"]
    assert results["FIT.DOF_VALUES"].shape == (len(results["TIME"]), sizes["N_FIT_DOF"])
    assert results["FIT.DOF_NAMES"] == preshot["FIT:DOF_NAMES"]
    assert results["FIT.DOF_UNITS"] == preshot["FIT:DOF_UNITS"]
    assert len(results["FIT.DOF_UNITS"]) == sizes["N_FIT_DOF"]
    assert not {"COEFFICIENTS", "COEF_NAMES", "COEF_UNITS"} & results.nested["FIT"].keys()


def test_exported_field_and_geometry(analytic_results: Results) -> None:
    results = analytic_results
    psi = results["PROFILES_2D.R_Z.PSI"]
    r, z = results["PROFILES_2D.R_Z.R"], results["PROFILES_2D.R_Z.Z"]
    d_z, d_r = np.gradient(psi, z, r, axis=(1, 2), edge_order=2)
    assert_allclose(results["PROFILES_2D.R_Z.B_FIELD_R"], -d_z / (2 * np.pi * r), atol=1e-13)
    assert_allclose(results["PROFILES_2D.R_Z.B_FIELD_Z"], d_r / (2 * np.pi * r), atol=1e-13)
    radius = 0.5
    assert results["GLOBAL.AREA"][-1] == pytest.approx(np.pi * radius**2, rel=0.02)
    assert results["GLOBAL.VOLUME"][-1] == pytest.approx(2 * np.pi**2 * 10 * radius**2, rel=0.02)
    assert results["GLOBAL.LENGTH_POL"][-1] == pytest.approx(2 * np.pi * radius, rel=0.02)


def test_failed_slices_do_not_acquire_plausible_derived_values(analytic_zero_results: Results) -> None:
    for key in (
        "FIT.PREDICTION",
        "FIT.RESIDUAL",
        "FIT.FF_PRIME",
        "GLOBAL.AREA",
        "GLOBAL.VOLUME",
        "GLOBAL.LENGTH_POL",
        "PROFILES_2D.R_Z.B_FIELD_R",
        "CONSTRAINTS.CHI_SQ_REG",
    ):
        assert np.all(np.isnan(analytic_zero_results[key])), key


def test_logged_outputs_can_be_processed_without_running_the_solver(
    analytic_preshot: Preshot, monkeypatch: pytest.MonkeyPatch
) -> None:
    from rtgsfit.database_readers import Signals
    from rtgsfit.library import Rtgsfit
    from rtgsfit.replay import read_signals
    from rtgsfit.results import collect_results
    from rtgsfit.settings import load_settings

    settings = load_settings("analytic")
    signals = read_signals(settings, analytic_preshot)
    solver = Rtgsfit.build(analytic_preshot)
    captured = solver(signals.meas[0], signals.coil_current(analytic_preshot)[0])
    native_psi = captured.psi.copy()

    def reject_solver_call(*args: object, **kwargs: object) -> None:
        pytest.fail("offline post-processing called the solver")

    monkeypatch.setattr(Rtgsfit, "__call__", reject_solver_call)
    logged = Signals(signals.time[:1], signals.meas[:1], signals.coil_signal[:1])
    results = collect_results(analytic_preshot, logged, [captured], settings.json_text, "captured log")
    assert_allclose(
        results["PROFILES_2D.R_Z.PSI"][0].ravel(),
        native_psi * analytic_preshot["INFO:PSI_TO_WB"],
    )
    assert_allclose(captured.psi, native_psi)
    assert np.all(np.isfinite(results["PROFILES_2D.R_Z.B_FIELD_R"]))


def test_already_processed_capture_is_rejected(analytic_results: Results, analytic_preshot: Preshot) -> None:
    from rtgsfit.post_processing import post_process

    with pytest.raises(ValueError, match="already"):
        post_process(analytic_results.nested, analytic_preshot)


def test_dof_units_must_match_dof_count(analytic_results: Results, analytic_preshot: Preshot) -> None:
    from copy import deepcopy

    from rtgsfit.results import validate_results

    invalid = Results(deepcopy(analytic_results.nested))
    invalid.nested["FIT"]["DOF_UNITS"] = ["dimensionless"]
    with pytest.raises(ValueError, match="DOF_UNITS"):
        validate_results(invalid, analytic_preshot)


def test_dof_units_must_match_preshot(analytic_results: Results, analytic_preshot: Preshot) -> None:
    from copy import deepcopy

    from rtgsfit.results import validate_results

    invalid = Results(deepcopy(analytic_results.nested))
    invalid.nested["FIT"]["DOF_UNITS"][0] = "ampere"
    with pytest.raises(ValueError, match=r"PRESHOT.FIT:DOF_UNITS"):
        validate_results(invalid, analytic_preshot)


@pytest.mark.parametrize(
    ("dof", "expected"),
    [
        ("constant_current_density", "ampere / metre ** 2"),
        ("eig_01", "ampere"),
    ],
)
def test_gsfit_passive_dof_units(dof: str, expected: str) -> None:
    from rtgsfit.database_readers.gsfit import _passive_dof_unit

    assert _passive_dof_unit(dof) == expected


@pytest.mark.parametrize(
    "info, bits",
    [
        (0, 0),
        (-1, 128 | 256),
        (-12, 128 | 256),
        (-1010, 128 | 2048),
        (-1011, 128 | 2048),
        (1, 128 | 512),
        (27, 128 | 512),
    ],
)
def test_lapack_status_classification(analytic_preshot: Preshot, info: int, bits: int) -> None:
    import ctypes

    from rtgsfit.library import Rtgsfit

    library = Rtgsfit.build(analytic_preshot)
    # A separate CDLL would keep the old preshot loaded after Rtgsfit rebuilds it.
    function = library._lib.rtgsfit_status_from_lapack_info
    function.argtypes = [ctypes.c_int32]
    function.restype = ctypes.c_int32
    assert function(info) == bits


def test_unknown_status_bits_are_visible() -> None:
    assert describe_status(4096) == "UNKNOWN_BITS(0x1000)"
    assert describe_status(128 | 512) == "LAPACK | FIT_NO_CONVERGENCE"
