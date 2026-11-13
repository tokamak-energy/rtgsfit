"""The command line interface, run in-process."""

from __future__ import annotations

from pathlib import Path

import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit.__main__ import main
from rtgsfit.database_readers.local import read_results
from rtgsfit.database_writers.local import write_npz
from rtgsfit.preshot import Preshot
from rtgsfit.results import Results


def test_list_settings(capsys: pytest.CaptureFixture[str]) -> None:
    main(["list-settings"])
    assert {"analytic", "freegsnke", "st40"} <= set(capsys.readouterr().out.split())


def test_write_csv(tmp_path: Path) -> None:
    path = tmp_path / "RTGSFIT.csv"
    main(["write-csv", str(path)])
    assert path.read_text() == dd.generate_csv()


def test_replay_from_a_local_preshot(
    analytic_preshot: Preshot, analytic_results: Results, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """`rtgsfit replay analytic` with the preshot read from a file and the results written to one."""
    preshot_path = write_npz(tmp_path / "preshot.npz", analytic_preshot.to_nested())
    results_path = tmp_path / "results.npz"
    main([
        "replay", "analytic", "--quiet",
        "--set", "preshot/reader=local",
        "--set", f"preshot/local/path={preshot_path}",
        "--set", f"writers/local/path={results_path}",
    ])  # fmt: skip
    assert str(results_path) in capsys.readouterr().out
    results = read_results(results_path)
    assert_array_equal(results["PROFILES_2D.R_Z.PSI"], analytic_results["PROFILES_2D.R_Z.PSI"])
    assert "local" in results["INPUT.PRESHOT_SRC"]
    assert str(results_path) in results["INPUT.SETTINGS"]


def test_unknown_setting_is_an_error() -> None:
    with pytest.raises(ValueError):
        main(["replay", "analytic", "--set", "replay/no_such_setting=1"])
