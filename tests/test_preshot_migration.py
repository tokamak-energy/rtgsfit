"""Legacy preshots are converted explicitly, not silently treated as the new basis convention."""

from __future__ import annotations

import sys
from copy import deepcopy
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from rtgsfit import data_dictionary as dd
from rtgsfit.constants import MU_0
from rtgsfit.database_readers.rtgsfit_mdsplus import RtgsfitMdsplusReader
from rtgsfit.preshot import LEGACY_FIT_GREENS
from rtgsfit.preshot import LEGACY_PRESHOT_NAMES
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import PreshotError
from rtgsfit.settings import load_settings


def transposed_fit_greens(groups: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Store the fit Greens functions in their original A_* [source, N_FIT_TARGET] layout."""
    for (group, name), old_name in LEGACY_FIT_GREENS.items():
        groups[group][old_name] = np.ascontiguousarray(np.asarray(groups[group].pop(name)).T)
    return groups


def legacy_groups(preshot: Preshot) -> dict[str, dict[str, Any]]:
    groups = transposed_fit_greens(deepcopy(preshot.groups))
    # Recreate the old analytic Wb/rad contract, including its unweighted flux-loop rows.
    n_loop = preshot.sizes["N_FLUX_LOOP"]
    groups["POISSON"]["RHS_SCALE"] /= 2.0 * np.pi
    groups["GRID"]["G_COIL"] /= 2.0 * np.pi
    groups["GRID"]["G_PASSIVE"] /= 2.0 * np.pi
    groups["CONSTRAINT"]["G_COIL"][:n_loop] /= 2.0 * np.pi
    groups["FIT"]["WEIGHT"][:n_loop] *= 2.0 * np.pi
    for (group, name), old_name in LEGACY_PRESHOT_NAMES.items():
        groups[group][old_name] = groups[group].pop(name)
    groups["GRID"]["INV_MU0_R"] = 1.0 / (groups["GRID"]["R_FLAT"] * MU_0)
    groups["FIT"]["X_VEC_NAMES"][1] = "FF_PRIME"
    del groups["INFO"]["MU0"]
    del groups["INFO"]["PSI_TO_WB"]
    return groups


def test_legacy_preshot_has_the_same_compiled_contract(analytic_preshot: Preshot) -> None:
    old = legacy_groups(analytic_preshot)
    old["INFO"]["SETTINGS"] = "<parameters><legacy_settings/></parameters>"
    migrated = Preshot.from_nested(old)
    assert migrated.content_hash() != analytic_preshot.content_hash()
    np.testing.assert_allclose(
        migrated["POISSON:RHS_SCALE"] * migrated["INFO:PSI_TO_WB"],
        analytic_preshot["POISSON:RHS_SCALE"],
    )
    assert migrated["INFO:PSI_TO_WB"] == 2 * np.pi
    assert migrated["FIT:DOF_NAMES"][1] == "FF_PRIME_MU0"
    assert "N_X_VEC" in old["SIZES"]  # Original data was not mutated.
    assert migrated["INFO:SETTINGS"] == old["INFO"]["SETTINGS"]


def test_preshot_without_dof_units_requires_rebuild(analytic_preshot: Preshot) -> None:
    old = deepcopy(analytic_preshot.groups)
    del old["FIT"]["DOF_UNITS"]
    with pytest.raises(PreshotError, match="rebuild the preshot"):
        Preshot.from_nested(old)


def test_unknown_legacy_units_are_rejected(analytic_preshot: Preshot) -> None:
    old = legacy_groups(analytic_preshot)
    old["INFO"]["SOURCE"] = "unknown model"
    with pytest.raises(PreshotError, match="unknown legacy flux"):
        Preshot.from_nested(old)


def test_unsupported_legacy_permeability_is_rejected(analytic_preshot: Preshot) -> None:
    old = legacy_groups(analytic_preshot)
    old["GRID"]["INV_MU0_R"] *= 2
    with pytest.raises(PreshotError, match="mu0 convention"):
        Preshot.from_nested(old)


def test_coefficient_named_preshot_is_renamed_without_rescaling(analytic_preshot: Preshot) -> None:
    old = deepcopy(analytic_preshot.groups)
    old["SIZES"]["N_FIT_COEF"] = old["SIZES"].pop("N_FIT_DOF")
    old["FIT"]["COEF_NAMES"] = old["FIT"].pop("DOF_NAMES")
    migrated = Preshot.from_nested(old)
    assert migrated.content_hash() == analytic_preshot.content_hash()
    assert migrated["FIT:DOF_NAMES"] == analytic_preshot["FIT:DOF_NAMES"]
    assert "N_FIT_COEF" in old["SIZES"]
    assert "DOF_NAMES" not in old["FIT"]


def test_ambiguous_coefficient_named_preshot_is_rejected(analytic_preshot: Preshot) -> None:
    old = deepcopy(analytic_preshot.groups)
    old["SIZES"]["N_FIT_COEF"] = old["SIZES"]["N_FIT_DOF"]
    with pytest.raises(PreshotError, match="ambiguous"):
        Preshot.from_nested(old)


def test_transposed_fit_greens_are_migrated(analytic_preshot: Preshot) -> None:
    old = transposed_fit_greens(deepcopy(analytic_preshot.groups))
    migrated = Preshot.from_nested(old)
    assert migrated.content_hash() == analytic_preshot.content_hash()
    assert migrated["FIT:G_GRID"].shape == (analytic_preshot.sizes["N_FIT_TARGET"], analytic_preshot.sizes["N_GRID"])
    assert "A_GRID" in old["FIT"]  # Original data was not mutated.
    assert "G_GRID" not in old["FIT"]


def test_ambiguous_fit_greens_are_rejected(analytic_preshot: Preshot) -> None:
    old = deepcopy(analytic_preshot.groups)
    old["FIT"]["A_GRID"] = old["FIT"]["G_GRID"].T
    with pytest.raises(PreshotError, match="ambiguous"):
        Preshot.from_nested(old)


def test_legacy_thin_reader_requests_only_original_nodes(
    analytic_preshot: Preshot, monkeypatch: pytest.MonkeyPatch
) -> None:
    old = legacy_groups(analytic_preshot)

    class Connection:
        def __init__(self, server: str) -> None:
            pass

        def __enter__(self) -> Connection:
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def openTree(self, tree: str, pulse: int) -> None:
            pass

        def get(self, expression: str) -> SimpleNamespace:
            key = expression.split("RUN21.")[1].replace(":", ".")
            value = dd.nested_get({"PRESHOT": old}, tuple(key.split(".")))
            return SimpleNamespace(data=lambda: value)

    monkeypatch.setitem(sys.modules, "MDSplus", None)
    monkeypatch.setitem(sys.modules, "mdsthin", SimpleNamespace(Connection=Connection))
    settings = load_settings("st40").copy_with({
        "preshot/rtgsfit_mdsplus/run_name": "RUN21",
        "preshot/rtgsfit_mdsplus/legacy_schema": True,
    })  # fmt: skip
    actual = RtgsfitMdsplusReader().read_preshot(settings)
    assert actual.content_hash() == Preshot.from_nested(old).content_hash()
