"""XML settings: parsing, overrides and the settings shipped with the package."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from rtgsfit.database_readers import _PRESHOT_READERS
from rtgsfit.database_readers import _SIGNALS_READERS
from rtgsfit.database_writers import _WRITERS
from rtgsfit.settings import available_settings
from rtgsfit.settings import load_settings
from rtgsfit.settings import parse_overrides
from rtgsfit.settings import parse_settings
from rtgsfit.settings import replay_times

XML = """<?xml version="1.0"?>
<parameters>
  <replay>
    <pulse type="int">13343</pulse>
    <time_start type="float" unit="s">0.0</time_start>
    <time_end type="float" unit="s">0.1</time_end>
    <time_step type="float" unit="s">0.025</time_step>
  </replay>
  <fit>
    <names type="str_array">A B C</names>
    <weights type="float_array">1.0 2.5</weights>
    <index type="int_array">3 4</index>
    <enabled type="bool">yes</enabled>
  </fit>
</parameters>
"""


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


def test_parse_types() -> None:
    settings = parse_settings(XML)
    assert settings["replay"]["pulse"] == 13343
    assert settings.get_path("replay/time_step") == 0.025
    assert settings["fit"]["names"] == ["A", "B", "C"]
    assert_array_equal(settings["fit"]["weights"], [1.0, 2.5])
    assert settings["fit"]["index"].dtype == np.int64
    assert settings["fit"]["enabled"] is True


@pytest.mark.parametrize(
    "xml",
    [
        "<settings><a type='int'>1</a></settings>",
        "<parameters><a>1</a></parameters>",
        "<parameters><a type='complex'>1</a></parameters>",
        "<parameters><a type='int'>1</a><a type='int'>2</a></parameters>",
        "<parameters><a type='bool'>maybe</a></parameters>",
    ],
)
def test_parse_rejects_invalid_xml(xml: str) -> None:
    with pytest.raises(ValueError):
        parse_settings(xml)


def test_missing_setting_raises_key_error() -> None:
    with pytest.raises(KeyError):
        parse_settings(XML).get_path("replay/missing")


def test_overrides_use_the_declared_type() -> None:
    settings = parse_settings(XML)
    overrides = parse_overrides(settings, ["replay/pulse=13345", "fit/weights=3 4", "new/flag:bool=true"])
    assert overrides["replay/pulse"] == 13345
    assert_array_equal(overrides["fit/weights"], [3.0, 4.0])
    assert overrides["new/flag"] is True
    with pytest.raises(ValueError):
        parse_overrides(settings, ["new/flag=true"])
    with pytest.raises(ValueError):
        parse_overrides(settings, ["replay/pulse"])


def test_copy_with_updates_values_and_xml() -> None:
    settings = parse_settings(XML)
    changed = settings.copy_with({"replay/pulse": 13346, "fit/weights": np.array([7.0]), "new/name": "x"})
    assert settings["replay"]["pulse"] == 13343
    reparsed = parse_settings(changed.xml_text)
    assert reparsed["replay"]["pulse"] == changed["replay"]["pulse"] == 13346
    assert_array_equal(reparsed["fit"]["weights"], [7.0])
    assert reparsed["new"]["name"] == "x"


def test_replay_times_include_the_end() -> None:
    assert_allclose(replay_times(parse_settings(XML)), [0.0, 0.025, 0.05, 0.075, 0.1])
    single = parse_settings(XML).copy_with({"replay/time_end": 0.0})
    assert_allclose(replay_times(single), [0.0])
    with pytest.raises(ValueError):
        replay_times(parse_settings(XML).copy_with({"replay/time_step": 0.0}))
