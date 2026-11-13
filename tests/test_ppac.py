"""The PPAC replay class: its constructor contract, its settings and a complete run without the TE network."""

from __future__ import annotations

import os
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from numpy.testing import assert_array_equal
from rtgsfit import data_dictionary as dd
from rtgsfit import ppac
from rtgsfit.library import LIBRARY
from rtgsfit.library import build_directory
from rtgsfit.library import isolated_build
from rtgsfit.preshot import Preshot
from rtgsfit.replay import replay
from rtgsfit.results import Results
from rtgsfit.settings import Settings

ANALYTIC_READERS = {
    "preshot/reader": "analytic",
    "signals/reader": "analytic",
    "replay/time_start": 0.0,
    "replay/time_end": 0.009,
}
"""Replace the ST40 readers so that `run()` replays the analytic case, as `analytic_results` does."""


def test_ppac_constructs_and_runs_it_as_slurm_python_job_does() -> None:
    # neo_ppac/slurm_python_job.py: getattr(import_module(package_name), class_name)(...), then run().
    package_name, class_name = "rtgsfit", "RtgsfitReplay"
    class_ = getattr(import_module(package_name), class_name)
    controller = class_(
        pulseNo=13343, run_name="RUN01_POST", run_description="PPAC run", pulseNo_write=None, link_run_to_best=True
    )
    assert callable(controller.run)
    settings = controller.settings
    assert settings.get_path("replay/pulse") == 13343
    assert settings.get_path("gsfit/pulse") == 13343
    assert settings.get_path("preshot/reader") == "gsfit"
    assert settings.get_path("signals/reader") == "st40_mdsplus"
    assert settings.get_path("signals/st40_mdsplus/mag_run") == "BEST"
    assert settings.get_path("writers/names") == ["rtgsfit_mdsplus"]
    assert settings.get_path("writers/rtgsfit_mdsplus") == {
        "enabled": True, "pulse": 13343, "run_name": "RUN01_POST", "run_description": "PPAC run", "link_best": True,
    }  # fmt: skip


def test_ppac_replay_defaults_to_run01_post_linked_to_best() -> None:
    settings = ppac.RtgsfitReplay(pulseNo=13343).settings
    assert settings.get_path("writers/rtgsfit_mdsplus/run_name") == "RUN01_POST"
    assert settings.get_path("writers/rtgsfit_mdsplus/link_best") is True


def test_ppac_replays_the_pcs_window() -> None:
    settings = ppac.RtgsfitReplay(pulseNo=13343).settings
    assert settings.get_path("replay/time_start") == 0.0055
    assert settings.get_path("replay/time_end") == 0.6
    assert settings.get_path("replay/time_step") == 0.001


def test_ppac_settings_options() -> None:
    controller = ppac.RtgsfitReplay(
        pulseNo=13343,
        pulseNo_write=45_013_343,
        write_to_mds=False,
        overrides={"replay/time_end": 0.2, "build/arch_flags": "-march=x86-64-v3"},
    )
    assert controller.pulseNo_write == 45_013_343
    assert controller.settings.get_path("writers/names") == []
    assert controller.settings.get_path("writers/rtgsfit_mdsplus/enabled") is False
    assert controller.settings.get_path("writers/rtgsfit_mdsplus/pulse") == 45_013_343
    assert controller.settings.get_path("replay/time_end") == 0.2
    assert controller.settings.get_path("build/arch_flags") == "-march=x86-64-v3"


def test_ppac_run_builds_privately_and_writes_to_mdsplus(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, analytic_preshot: Preshot, analytic_results: Results
) -> None:
    package = tmp_path / "standard_utility"
    (package / "csv_files").mkdir(parents=True)
    (package / "csv_files" / "RTGSFIT.csv").write_text(dd.generate_csv())
    calls: list[dict[str, Any]] = []
    monkeypatch.setitem(
        sys.modules,
        "standard_utility",
        SimpleNamespace(
            __file__=str(package / "__init__.py"),
            create_script_nodes=lambda **kwargs: calls.append(kwargs),
            write_script_data=lambda **kwargs: calls.append(kwargs),
        ),
    )
    copies: list[Path] = []

    @contextmanager
    def recorded_isolated_build() -> Iterator[Path]:
        with isolated_build() as copy:
            copies.append(copy)
            yield copy

    def replay_analytic(settings: Settings, verbose: bool) -> Results:
        # The analytic preshot takes ~15 s to build, so reuse the session's.
        return replay(settings, preshot=analytic_preshot, verbose=verbose)

    monkeypatch.setattr(ppac, "isolated_build", recorded_isolated_build)
    monkeypatch.setattr(ppac, "replay", replay_analytic)
    shared_library = build_directory() / LIBRARY
    shared_before = shared_library.stat().st_mtime_ns if shared_library.exists() else None
    environment = {name: os.environ.get(name) for name in ("RTGSFIT_C_DIR", "RTGSFIT_CACHE_DIR")}

    controller = ppac.RtgsfitReplay(
        pulseNo=0, run_name="RUN02_POST", run_description="test", pulseNo_write=45_000_000, link_run_to_best=True,
        settings="analytic", overrides=ANALYTIC_READERS, verbose=False,
    )  # fmt: skip
    results = controller.run()

    assert controller.results is results
    assert_array_equal(results["GLOBAL.IP"], analytic_results["GLOBAL.IP"])
    assert_array_equal(results["PROFILES_2D.R_Z.PSI"], analytic_results["PROFILES_2D.R_Z.PSI"])
    assert len(copies) == 1 and not copies[0].exists()
    assert (shared_library.stat().st_mtime_ns if shared_library.exists() else None) == shared_before
    assert {name: os.environ.get(name) for name in environment} == environment
    create, write = calls
    assert create["pulseNo_write"] == write["pulseNo_write"] == 45_000_000
    assert create["run_name"] == write["run_name"] == "RUN02_POST"
    assert create["link_best"] is True
    assert create["use_mdsthin"] is write["use_mdsthin"] is True
