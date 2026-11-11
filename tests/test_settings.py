"""JSON settings, schema validation, overrides and reproducible archived values."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from rtgsfit.database_readers import _PRESHOT_READERS
from rtgsfit.database_readers import _SIGNALS_READERS
from rtgsfit.database_writers import _WRITERS
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results
from rtgsfit.settings import available_settings
from rtgsfit.settings import load_settings
from rtgsfit.settings import parse_overrides
from rtgsfit.settings import parse_settings
from rtgsfit.settings import replay_times


def test_shipped_settings() -> None:
    assert {"analytic", "freegsnke", "st40"} <= set(available_settings())


@pytest.mark.parametrize("name", available_settings())
def test_shipped_settings_load_and_name_known_readers(name: str) -> None:
    settings = load_settings(name)
    assert settings["preshot"]["reader"] in _PRESHOT_READERS
    assert settings["preshot"].get("builder", settings["preshot"]["reader"]) in _PRESHOT_READERS
    assert settings["signals"]["reader"] in _SIGNALS_READERS
    assert set(settings["writers"]["names"]) <= set(_WRITERS)
    assert len(replay_times(settings)) >= 1


def test_native_json_types() -> None:
    settings = load_settings("st40")
    assert type(settings["replay"]["pulse"]) is int
    assert type(settings.get_path("replay/time_step")) is float
    assert type(settings["build"]["timing"]) is bool
    assert isinstance(settings["gsfit"]["rogowski_names"], list)
    matrix = np.asarray(settings["gsfit"]["coils"]["matrix"])
    assert matrix.shape == (8, 8)
    assert_allclose(matrix @ np.arange(8), [7, 6, 10, 3, 0, 2, 1, 0])


@pytest.mark.parametrize(
    "text",
    [
        "[]",
        "null",
        '{"build": {"timing": true, "timing": false}}',
        '{"build": {"timing": NaN}}',
        '{"replay": {"time_step": Infinity}}',
        '{"replay": {"time_step": 1e999}}',
        '{"build": {"timing": "false"}}',
        '{"unknown_section": {}}',
        '{"build": {"no_such_setting": true}}',
        '{"replay": {"repeat": true}}',
        '{"replay": {"time_step": 0}}',
        '{"replay": {"time_step": 0.1},}',
        "<parameters/>",
    ],
)
def test_invalid_json_or_schema_is_rejected(text: str) -> None:
    with pytest.raises(ValueError):
        parse_settings(text)


def test_schema_errors_name_the_setting() -> None:
    with pytest.raises(ValueError, match="replay/repeat"):
        load_settings("analytic").copy_with({"replay/repeat": 1.5})


def test_missing_setting_raises_key_error() -> None:
    with pytest.raises(KeyError):
        load_settings("analytic").get_path("replay/missing")


def test_overrides_use_json_literals_and_schema_types() -> None:
    settings = load_settings("analytic")
    overrides = parse_overrides(
        settings,
        [
            "replay/pulse=13345",
            "build/debug=true",
            "analytic/passives/r=[7.0, 8.0]",
            "writers/local/path=path with spaces.npz",
            'writers/names=["local"]',
        ],
    )
    assert overrides["replay/pulse"] == 13345
    assert overrides["analytic/passives/r"] == [7.0, 8.0]
    assert overrides["build/debug"] is True
    assert overrides["writers/local/path"] == "path with spaces.npz"
    assert overrides["writers/names"] == ["local"]


@pytest.mark.parametrize(
    "assignment",
    [
        "replay/pulse",
        "replay/pulse=word",
        "build/timing=False",
        "build/timing=0",
        "replay/time_step=NaN",
        "replay/time_step=1e999",
        "build/timing:bool=true",
        "build/no_such_setting=true",
        "build/timing/child=1",
        "/build/timing=true",
    ],
)
def test_invalid_overrides_are_rejected(assignment: str) -> None:
    with pytest.raises(ValueError):
        parse_overrides(load_settings("analytic"), [assignment])


def test_schema_types_allow_optional_settings_without_type_suffixes() -> None:
    settings = load_settings("analytic")
    overrides = parse_overrides(settings, ["preshot/reader=local", "preshot/local/path=shot.npz"])
    changed = settings.copy_with(overrides)
    assert changed["preshot"]["local"]["path"] == "shot.npz"


def test_copy_with_updates_archived_json_independently() -> None:
    settings = load_settings("analytic")
    changed = settings.copy_with(
        {
            "replay/pulse": 13346,
            "analytic/passives/r": [7.0],
            "writers/local/path": "shot.npz",
        }
    )
    assert settings["replay"]["pulse"] == 0
    reparsed = parse_settings(changed.json_text)
    assert reparsed["replay"]["pulse"] == changed["replay"]["pulse"] == 13346
    assert reparsed["analytic"]["passives"]["r"] == [7.0]
    changed["replay"]["pulse"] = 7
    assert json.loads(changed.json_text)["replay"]["pulse"] == 7


def test_file_load_round_trip_and_missing_path(tmp_path: Path) -> None:
    path = tmp_path / "run.json"
    path.write_text(load_settings("analytic").json_text)
    assert load_settings(path) == load_settings("analytic")
    with pytest.raises(FileNotFoundError):
        load_settings(tmp_path / "missing.json")


def test_new_preshot_and_replay_archive_effective_json(analytic_preshot: Preshot, analytic_results: Results) -> None:
    settings = load_settings("analytic")
    assert parse_settings(analytic_preshot["INFO:SETTINGS"]) == settings
    assert parse_settings(analytic_results["INPUT.SETTINGS"]) == settings


def test_missing_sections_are_rejected() -> None:
    with pytest.raises(ValueError, match="required"):
        parse_settings("{}")


def test_schema_hint_cannot_override_runtime_validation() -> None:
    settings = load_settings("analytic")
    values = dict(settings)
    values["$schema"] = "https://invalid.example/untrusted-schema.json"
    values["build"] = {"timing": "true"}
    with pytest.raises(ValueError, match="build/timing"):
        parse_settings(json.dumps(values))


def test_replay_times_include_the_end() -> None:
    settings = load_settings("analytic").copy_with(
        {
            "replay/time_start": 0.0,
            "replay/time_end": 0.1,
            "replay/time_step": 0.025,
        }
    )
    assert_allclose(replay_times(settings), [0.0, 0.025, 0.05, 0.075, 0.1])
    assert_allclose(replay_times(settings.copy_with({"replay/time_end": 0.0})), [0.0])
    with pytest.raises(ValueError):
        replay_times(settings.copy_with({"replay/time_end": -0.1}))
