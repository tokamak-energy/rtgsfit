"""Unvalidated replay reports in standard_utility's DATA_Q layout."""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from rtgsfit.results import Results
from rtgsfit.settings import Settings


def unvalidated_quality(
    settings: Settings, checks: dict[str, Any], *, input_reader: str | None = None
) -> dict[str, Any]:
    return {
        "SCORE": 0,
        "SCORE_INFO": "Not yet validated; input and solver checks do not establish equilibrium quality",
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


def add_replay_quality(
    results: Results, settings: Settings, *, input_reader: str | None = None, input_source: str | None = None
) -> None:
    status = np.asarray(results["GLOBAL.STATUS"])
    codes, counts = np.unique(status, return_counts=True)
    input_quality = results.input_quality.copy()
    if input_source is not None:
        input_quality["source"] = input_source
    quality = unvalidated_quality(
        settings,
        {
            "completed": True,
            "input": input_quality,
            "solver_status_counts": {str(code): int(count) for code, count in zip(codes, counts, strict=True)},
            "note": "Solver status is not a physics-validation or plasma-presence test.",
        },
        input_reader=input_reader,
    )
    quality["TIME"] = np.asarray(results["TIME"])
    quality["SCORE_1D"] = np.zeros(len(status), dtype=np.int32)
    results.nested["DATA_Q"] = quality
