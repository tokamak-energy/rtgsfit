"""Thin-client transport only; these tests never contact or modify the shared database."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit.database_readers.st40_mdsplus import St40MdsplusReader
from rtgsfit.database_writers.rtgsfit_mdsplus import write_to_mdsplus
from rtgsfit.preshot import Preshot
from rtgsfit.settings import load_settings


def test_writer_explicitly_uses_mdsthin(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    package = tmp_path / "standard_utility"
    (package / "csv_files").mkdir(parents=True)
    (package / "csv_files" / "RTGSFIT.csv").write_text(dd.generate_csv())
    calls: list[dict[str, Any]] = []

    def record(**kwargs: Any) -> None:
        calls.append(kwargs)

    monkeypatch.setitem(sys.modules, "MDSplus", None)
    monkeypatch.setitem(
        sys.modules,
        "standard_utility",
        SimpleNamespace(
            __file__=str(package / "__init__.py"),
            create_script_nodes=record,
            write_script_data=record,
        ),
    )
    settings = load_settings("st40").copy_with({"writers/rtgsfit_mdsplus/enabled": True})
    write_to_mdsplus({"TIME": np.array([0.0])}, settings)
    assert len(calls) == 2
    assert all(call["use_mdsthin"] is True for call in calls)


def _fake_st40(
    monkeypatch: pytest.MonkeyPatch, preshot: Preshot, psu_time: list[float], calls: list[str] | None = None
) -> None:
    """MAG is recorded at 0, 10 and 20 ms with values 1, 2 and 3; the PSU at `psu_time` with the
    same values."""

    class Connection:
        def __init__(self, server: str) -> None:
            assert server == "smaug"

        def __enter__(self) -> Connection:
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def openTree(self, tree: str, pulse: int) -> None:
            if calls is not None:
                calls.append(tree)

        def get(self, expression: str) -> SimpleNamespace:
            if calls is not None:
                calls.append(expression)
            if expression.startswith("dim_of("):
                values = np.array([0.0, 0.01, 0.02] if "\\MAG::" in expression else psu_time)
            else:
                values = np.array([1.0, 2.0, 3.0] if "\\MAG::" in expression else [1.0, 2.0, 3.0][: len(psu_time)])
            return SimpleNamespace(data=lambda: values)

    monkeypatch.setitem(sys.modules, "MDSplus", None)
    monkeypatch.setitem(sys.modules, "st40_database", None)
    monkeypatch.setitem(sys.modules, "mdsthin", SimpleNamespace(Connection=Connection))
    monkeypatch.setitem(preshot.groups["MEAS"], "NAMES", ["PSI_FLOOP_007"] * preshot.sizes["N_MEAS"])
    monkeypatch.setitem(preshot.groups["COIL"], "SIGNALS", ["I_SOL_PSU", "I_ROG_INIVC000"])


ST40_TIMES = {"replay/time_start": 0.0, "replay/time_end": 0.02, "replay/time_step": 0.01}


def test_st40_reader_needs_only_mdsthin(
    monkeypatch: pytest.MonkeyPatch,
    analytic_preshot: Preshot,
) -> None:
    calls: list[str] = []
    _fake_st40(monkeypatch, analytic_preshot, [0.0, 0.01, 0.02], calls)
    settings = load_settings("st40").copy_with(ST40_TIMES)
    signals = St40MdsplusReader().read_signals(settings, analytic_preshot)
    assert_array_equal(signals.meas[:, 0], [1.0, 2.0, 3.0])
    assert_array_equal(signals.coil_signal, [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    assert "\\MAG::TOP.BEST.FLOOP.L007:PSI" in calls
    assert "\\MAG::TOP.BEST.ROG.INIVC000:I" in calls


def test_st40_reader_drops_times_a_signal_did_not_record(
    monkeypatch: pytest.MonkeyPatch,
    analytic_preshot: Preshot,
) -> None:
    _fake_st40(monkeypatch, analytic_preshot, [0.0, 0.01])
    settings = load_settings("st40").copy_with(ST40_TIMES)
    signals = St40MdsplusReader().read_signals(settings, analytic_preshot)
    assert_array_equal(signals.time, [0.0, 0.01])
    assert_array_equal(signals.meas[:, 0], [1.0, 2.0])
    assert_array_equal(signals.coil_signal, [[1.0, 1.0], [2.0, 2.0]])


def test_st40_reader_rejects_replay_times_no_signal_covers(
    monkeypatch: pytest.MonkeyPatch,
    analytic_preshot: Preshot,
) -> None:
    _fake_st40(monkeypatch, analytic_preshot, [0.05, 0.06])
    settings = load_settings("st40").copy_with(ST40_TIMES)
    with pytest.raises(ValueError, match="no replay time"):
        St40MdsplusReader().read_signals(settings, analytic_preshot)
