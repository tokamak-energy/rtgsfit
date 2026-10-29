"""Offline diagnostics from captured RT-GSFit inputs and outputs; never called by the C solver.

The capture uses native flux units. This module converts exported flux to weber, but leaves the
embedded preshot and INPUT arrays in the native convention required to reproduce the C call.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtgsfit.preshot import Preshot


def post_process(nested: dict[str, Any], preshot: Preshot) -> None:
    """Enrich a collected capture in place, once, before validation or writing.

    Use ``collect_results`` for logged ``Signals`` and ``SliceResult`` objects. Repeated application
    to an already processed mapping is rejected: conversion must never silently happen twice.
    """
    fit = nested["FIT"]
    if "TARGET" in fit:
        raise ValueError("capture has already been post-processed")
    sizes = preshot.sizes
    n_time = len(nested["TIME"])
    n_constraint = sizes["N_CONSTRAINT"]
    n_plasma = sizes["N_PLASMA"]
    valid = nested["GLOBAL"]["STATUS"] == 0
    constraints = nested["CONSTRAINTS"]
    family_names = ("FLUX_LOOP", "BP_PROBE", "ROGOWSKI")
    measured = np.concatenate([constraints[name]["ALL"]["MEASURED"] for name in family_names], axis=1)
    reconstructed = np.concatenate([constraints[name]["ALL"]["RECONSTRUCT"] for name in family_names], axis=1)
    coil_current = constraints["PF_CURRENT"]["ALL"]["MEASURED"]
    coil_response = coil_current @ preshot["CONSTRAINT:G_COIL"].T
    weight = np.asarray(preshot["FIT:WEIGHT"])
    target = np.zeros((n_time, sizes["N_FIT_TARGET"]))
    target[:, :n_constraint] = (measured - coil_response) * weight[:n_constraint]
    prediction = np.zeros_like(target)
    prediction[:, :n_constraint] = (reconstructed - coil_response) * weight[:n_constraint]
    prediction[:, n_constraint:] = fit["DOF_VALUES"][:, n_plasma:] @ preshot["FIT:A_PASSIVE"][:, n_constraint:]
    prediction[~valid] = np.nan
    fit.update(
        TARGET_NAMES=list(preshot["FIT:TARGET_NAMES"]),
        TARGET=target,
        PREDICTION=prediction,
        RESIDUAL=prediction - target,
        CHI_SQ=(prediction - target) ** 2,
    )
    fit["FF_PRIME"] = np.where(valid, fit["DOF_VALUES"][:, 1] * float(preshot["INFO:MU0"]), np.nan)
    fit["DOF_UNITS"] = [
        "ampere / metre ** 3",
        "ampere / metre",
        "ampere / metre",
        *(["ampere"] * sizes["N_PASSIVE"]),
    ]
    for name in family_names:
        values = constraints[name]["ALL"]
        residual = values["RECONSTRUCT"] - values["MEASURED"]
        residual[~valid] = np.nan
        values["RESIDUAL"] = residual
        values["RESIDUAL_W"] = residual * values["WEIGHT"]
        values["CHI_SQ"] = values["RESIDUAL_W"] ** 2
        constraints[name]["CHI_SQ"] = np.where(valid, np.sum(values["CHI_SQ"], axis=1), np.nan)
    constraints["CHI_SQ_REG"] = np.where(valid, np.sum(fit["CHI_SQ"][:, n_constraint:], axis=1), np.nan)
    constraints["PF_PASSIVE"] = {
        "ALL": {
            "NAMES": list(preshot["FIT:DOF_NAMES"][n_plasma:]),
            "RECONSTRUCT": np.where(valid[:, None], fit["DOF_VALUES"][:, n_plasma:], np.nan),
        }
    }

    profile = nested["PROFILES_2D"]["R_Z"]
    flux_factor = float(preshot["INFO:PSI_TO_WB"])
    d_psi_d_z, d_psi_d_r = np.gradient(profile["PSI"], profile["Z"], profile["R"], axis=(1, 2), edge_order=2)
    field_factor = flux_factor / (2.0 * np.pi * profile["R"])
    profile["B_FIELD_R"] = np.where(valid[:, None, None], -d_psi_d_z * field_factor, np.nan)
    profile["B_FIELD_Z"] = np.where(valid[:, None, None], d_psi_d_r * field_factor, np.nan)
    profile["PSI"] *= flux_factor
    nested["GLOBAL"]["PSI_MAG_AXIS"] *= flux_factor
    nested["BOUNDARY"]["PSI"] *= flux_factor
    nested["XPOINTS"]["ALL"]["PSI"] *= flux_factor
    flux_loops = constraints["FLUX_LOOP"]["ALL"]
    for name in ("MEASURED", "RECONSTRUCT", "RESIDUAL"):
        flux_loops[name] *= flux_factor
    flux_loops["WEIGHT"] /= flux_factor

    area = np.full(n_time, np.nan)
    volume = np.full(n_time, np.nan)
    perimeter = np.full(n_time, np.nan)
    outline = nested["BOUNDARY"]["OUTLINE"]
    for i_time, n_outline in enumerate(outline["N"]):
        if not valid[i_time] or n_outline < 4:
            continue
        r = outline["R"][i_time, :n_outline]
        z = outline["Z"][i_time, :n_outline]
        r_next, z_next = np.roll(r, -1), np.roll(z, -1)
        cross = r * z_next - r_next * z
        area[i_time] = 0.5 * abs(np.sum(cross))
        volume[i_time] = np.pi / 3.0 * abs(np.sum((r + r_next) * cross))
        perimeter[i_time] = np.sum(np.hypot(r_next - r, z_next - z))
    nested["GLOBAL"].update(AREA=area, VOLUME=volume, LENGTH_POL=perimeter)
    nested["INPUT"]["SCHEMA"] = "rtgsfit-4; DOF naming; exported flux in Wb; capture/preshot in native units"
    nested["INPUT"]["COCOS"] = np.int32(13)
