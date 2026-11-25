"""Thin-client transport only; these tests never contact or modify the shared database."""

from __future__ import annotations

import sys
from datetime import UTC
from datetime import datetime
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit.database_writers.rtgsfit_mdsplus import _write_quality
from rtgsfit.database_writers.rtgsfit_mdsplus import write_quality_only
from rtgsfit.database_writers.rtgsfit_mdsplus import write_to_mdsplus
from rtgsfit.settings import load_settings


@pytest.mark.parametrize("existing_metadata", [False, True])
def test_writer_explicitly_uses_mdsthin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, existing_metadata: bool
) -> None:
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
    settings = load_settings("st40_magrun16").copy_with({"writers/rtgsfit_mdsplus/enabled": True})
    original: dict[str, Any] = {"TIME": np.array([0.0]), "PASSIVE": {"NAMES": []}, "COIL": {"NAMES": ["SOL"]}}
    if existing_metadata:
        original["CODE_VERSION"] = {"VERSION": "0.2.0", "DATETIME": "2000-01-01T00:00:00+00:00"}
    before = datetime.now(UTC)
    write_to_mdsplus(original, settings)
    after = datetime.now(UTC)
    assert len(calls) == 2
    assert all(call["use_mdsthin"] is True for call in calls)
    data = calls[1]["data_to_write"]
    assert data["PASSIVE"]["NAMES"].dtype.kind == "U"
    assert data["PASSIVE"]["NAMES"].shape == (0,)
    assert_array_equal(data["COIL"]["NAMES"], ["SOL"])
    code_version = data["CODE_VERSION"]
    assert {
        "COMPUTER",
        "DATETIME",
        "GIT_ID",
        "LIBRARY",
        "MATLAB",
        "PYTHON",
        "USER",
        "VERSION",
    } <= set(code_version)
    assert all(code_version[key] for key in ("COMPUTER", "GIT_ID", "LIBRARY", "PYTHON", "USER", "VERSION"))
    assert code_version["MATLAB"] == "N/A"
    written_at = datetime.fromisoformat(data["CODE_VERSION"]["DATETIME"])
    assert written_at.tzinfo == UTC
    assert before <= written_at <= after
    if existing_metadata:
        assert data["CODE_VERSION"]["VERSION"] == "0.2.0"
        assert original["CODE_VERSION"]["DATETIME"] == "2000-01-01T00:00:00+00:00"
    else:
        assert "CODE_VERSION" not in original


@pytest.mark.parametrize("existing", [False, True])
def test_quality_only_does_not_replace_a_run_or_link_best(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, existing: bool
) -> None:
    from mdsthin.exceptions import TreeNNF

    package = tmp_path / "standard_utility"
    (package / "csv_files").mkdir(parents=True)
    (package / "csv_files" / "RTGSFIT.csv").write_text(dd.generate_csv())
    calls: list[dict[str, Any]] = []

    class Connection:
        def __init__(self, server: str, timeout: int) -> None:
            pass

        def __enter__(self) -> Connection:
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def openTree(self, tree: str, pulse: int) -> None:
            pass

        def closeTree(self, tree: str, pulse: int) -> None:
            pass

        def get(self, expression: str) -> None:
            if not existing:
                raise TreeNNF()

    monkeypatch.setitem(sys.modules, "mdsthin", SimpleNamespace(Connection=Connection))
    monkeypatch.setitem(
        sys.modules,
        "standard_utility",
        SimpleNamespace(
            __file__=str(package / "__init__.py"),
            create_script_nodes=lambda **kwargs: calls.append(kwargs),
            write_script_data=lambda **kwargs: calls.append(kwargs),
            create_write_nodes_from_csv_file=lambda *args, **kwargs: calls.append(kwargs),
        ),
    )
    settings = load_settings("st40_magrun16").copy_with({"writers/rtgsfit_mdsplus/enabled": True})
    if existing:
        with pytest.raises(FileExistsError, match="already exists"):
            write_quality_only({"SCORE": 0}, settings)
        assert calls == []
    else:
        write_quality_only({"SCORE": 0}, settings)
        assert calls[0]["link_best"] is False
        assert set(calls[1]["data_to_write"]["CODE_VERSION"]) == {
            "COMPUTER",
            "DATETIME",
            "GIT_ID",
            "LIBRARY",
            "MATLAB",
            "PYTHON",
            "USER",
            "VERSION",
        }
        assert calls[2]["data_to_write"] == {"SCORE": 0}


def test_quality_writer_uses_standard_data_quality_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("standard_utility")
    module = import_module("standard_utility.read_write.create_write_nodes_from_csv_file")
    calls = []

    class NodeWriting:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            calls.append(kwargs)

        def write_nodes(self, data: dict[str, Any], force_write: bool) -> None:
            calls.append(data)

    monkeypatch.setattr(module.mdswrapper, "NodeWriting", NodeWriting)
    quality = {
        "SCORE": 0,
        "SCORE_INFO": "unvalidated",
        "RESULTS": "{}",
        "CONFIG": "{}",
        "SCORE_1D": np.array([0]),
        "TIME": np.array([0.0055]),
    }
    _write_quality(52_014_928, "DQTEST", quality)
    assert set(calls[0]["selected_nodes"]) == {"SCORE", "SCORE_INFO", "RESULTS", "CONFIG", "SCORE_1D", "TIME"}
    assert calls[1] is quality
