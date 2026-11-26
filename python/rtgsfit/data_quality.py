"""Conservative replay reports in standard_utility's five-level DATA_Q layout."""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from rtgsfit.results import Results
from rtgsfit.settings import Settings

SCORE_INFORMATION = {
    -2: "Missing data - no data present in MDS, so no data quality score",
    -1: "Invalid data - Should NOT be communicated internally or externally",
    0: "Default flag - Data not yet validated, should not be communicated internally or externally",
    1: (
        "Automatically validated data - Data is deemed trustworthy for quick internal analyses. "
        "Not suitable for external publication or communication"
    ),
    2: "Data is validated in detail - The data can be communicated externally and used in more complicated analysis",
}


def quality_report(
    settings: Settings, checks: dict[str, Any], *, score: int, input_reader: str | None = None
) -> dict[str, Any]:
    if score not in SCORE_INFORMATION:
        raise ValueError(f"unknown DATA_Q score {score}")
    return {
        "SCORE": score,
        "SCORE_INFO": SCORE_INFORMATION[score],
        "RESULTS": json.dumps(checks, sort_keys=True),
        "CONFIG": json.dumps(
            {
                "input_reader": input_reader or settings["signals"]["reader"],
                "pulse": int(settings["replay"]["pulse"]),
                "replay": settings["replay"],
            },
            sort_keys=True,
        ),
    }


def unvalidated_quality(
    settings: Settings, checks: dict[str, Any], *, input_reader: str | None = None
) -> dict[str, Any]:
    """A score-0 report for callers that have not validated their data."""
    return quality_report(settings, checks, score=0, input_reader=input_reader)


def add_replay_quality(
    results: Results, settings: Settings, *, input_reader: str | None = None, input_source: str | None = None
) -> None:
    """Flag failed slices as invalid; a run with any unvalidated slice remains unvalidated."""
    status = np.asarray(results["GLOBAL.STATUS"])
    codes, counts = np.unique(status, return_counts=True)
    score_1d = np.where(status == 0, 0, -1).astype(np.int32)
    if len(status) == 0:
        score = -2
    elif np.any(status == 0):
        score = 0
    else:
        score = -1
    input_quality = results.input_quality.copy()
    if input_source is not None:
        input_quality["source"] = input_source
    quality = quality_report(
        settings,
        {
            "completed": True,
            "input": input_quality,
            "solver_status_counts": {str(code): int(count) for code, count in zip(codes, counts, strict=True)},
            "note": "Solver status is not a physics-validation or plasma-presence test.",
        },
        score=score,
        input_reader=input_reader,
    )
    quality["TIME"] = np.asarray(results["TIME"])
    quality["SCORE_1D"] = score_1d
    results.nested["DATA_Q"] = quality
