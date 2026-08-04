"""
Timing benchmark for RTGSFIT.

Runs RTGSFIT N_REPEATS times for each of the 6 test cases, records per-section
CPU time (CLOCK_THREAD_CPUTIME_ID) in microseconds, uploads a summary CSV row
to a GitHub Gist, and prints a table to stdout.

Runs automatically as part of `pytest -s tests` when executed from the
`tests/rtgsfit_vs_gsfit/` sub-repository.
Environment variables:
  GITHUB_TOKEN      Personal-access token with gist scope.  If absent, the
                    upload step is skipped and results are only printed.
                    The Gist is found automatically by description; a new one
                    is created if none exists.
"""

import ctypes
import os
import socket
import subprocess
from datetime import datetime, timezone

import numpy as np
import pytest
import requests

from rtgsfit_vs_gsfit import config_loader, replay_rtgsfit, rtgsfit_compile_setup

T_NTIMERS = 18
N_REPEATS = 1_000
GIST_FILENAME = "rtgsfit_timing.csv"

CSV_HEADER = (
    "timestamp,commit,hostname,pulse,time_s,section,"
    "mean_us,median_us,std_us,min_us,max_us,p95_us,p99_us\n"
)

test_cases = [
    (13_343, 0.030),
    (13_343, 0.100),
    (13_345, 0.030),
    (13_345, 0.100),
    (13_346, 0.030),
    (13_346, 0.100),
]


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


GIST_DESCRIPTION = "RTGSFIT per-section timing benchmarks"


def _find_gist_id(headers: dict) -> str:
    """Return the ID of the existing timing Gist, or '' if not found."""

    page = 1
    while True:
        resp = requests.get(
            "https://api.github.com/gists",
            headers=headers,
            params={"per_page": 100, "page": page},
            timeout=10,
        )
        resp.raise_for_status()
        gists = resp.json()
        if not gists:
            return ""
        for g in gists:
            if g.get("description") == GIST_DESCRIPTION and GIST_FILENAME in g.get("files", {}):
                return g["id"]
        page += 1


def _upload_to_gist(csv_rows: str) -> None:
    """Append csv_rows to the Gist CSV file, or skip if no token."""
    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        print("\n[timing] GITHUB_TOKEN not set — skipping Gist upload.")
        return

    headers = {"Authorization": f"token {token}", "Accept": "application/vnd.github+json"}
    gist_id = _find_gist_id(headers)

    if gist_id:
        # Fetch existing content and append
        resp = requests.get(f"https://api.github.com/gists/{gist_id}", headers=headers, timeout=10)
        resp.raise_for_status()
        existing = resp.json()["files"].get(GIST_FILENAME, {}).get("content", CSV_HEADER)
        if not existing.startswith("timestamp"):
            existing = CSV_HEADER
        new_content = existing + csv_rows
        resp = requests.patch(
            f"https://api.github.com/gists/{gist_id}",
            headers=headers,
            json={"files": {GIST_FILENAME: {"content": new_content}}},
            timeout=10,
        )
        resp.raise_for_status()
        print(f"\n[timing] Gist updated: https://gist.github.com/{gist_id}")
    else:
        # Create a new gist
        resp = requests.post(
            "https://api.github.com/gists",
            headers=headers,
            json={
                "description": GIST_DESCRIPTION,
                "public": False,
                "files": {GIST_FILENAME: {"content": CSV_HEADER + csv_rows}},
            },
            timeout=10,
        )
        resp.raise_for_status()
        new_id = resp.json()["id"]
        print(f"\n[timing] New Gist created: https://gist.github.com/{new_id}")


def _run_timing_benchmark(cfg: dict) -> np.ndarray:
    """
    Run N_REPEATS calls to rtgsfit() starting from the converged state saved
    by the warm-up run.  Returns a (N_REPEATS, T_NTIMERS) float64 array of
    per-section CPU times in microseconds.

    All array sizes are derived from the output_dict — no extra MDS+ query.
    """
    librtgsfit_path = os.path.join(cfg["rtgsfit_path"], "lib", "librtgsfit.so")
    lib = ctypes.CDLL(librtgsfit_path)

    lib.rtgsfit.argtypes = [
        ctypes.POINTER(ctypes.c_double),   # meas_pcs
        ctypes.POINTER(ctypes.c_double),   # coil_curr
        ctypes.POINTER(ctypes.c_double),   # flux_norm
        ctypes.POINTER(ctypes.c_int32),    # mask
        ctypes.POINTER(ctypes.c_double),   # flux_total
        ctypes.POINTER(ctypes.c_double),   # chi_sq_err
        ctypes.POINTER(ctypes.c_double),   # lcfs_r
        ctypes.POINTER(ctypes.c_double),   # lcfs_z
        ctypes.POINTER(ctypes.c_int32),    # lcfs_n
        ctypes.POINTER(ctypes.c_double),   # coef
        ctypes.POINTER(ctypes.c_double),   # flux_boundary
        ctypes.POINTER(ctypes.c_double),   # plasma_current
        ctypes.POINTER(ctypes.c_int32),    # lcfs_err_code
        ctypes.POINTER(ctypes.c_int32),    # lapack_dgelss_info
        ctypes.POINTER(ctypes.c_double),   # meas_model
        ctypes.c_int32,                    # n_meas
        ctypes.POINTER(ctypes.c_double),   # r_mag_axis
        ctypes.POINTER(ctypes.c_double),   # z_mag_axis
        ctypes.POINTER(ctypes.c_double),   # mag_axis_flux
        ctypes.POINTER(ctypes.c_double),   # r_cur_centroid
        ctypes.POINTER(ctypes.c_double),   # z_cur_centroid
        ctypes.POINTER(ctypes.c_double),   # xpt_r
        ctypes.POINTER(ctypes.c_double),   # xpt_z
        ctypes.POINTER(ctypes.c_double),   # xpt_flux
        ctypes.c_int32,                    # n_xpt_max
        ctypes.POINTER(ctypes.c_int32),    # xpt_n
        ctypes.POINTER(ctypes.c_int32),    # xpt_diverted
    ]
    lib.rtgsfit.restype = None
    lib.rtgsfit_timing_reset.argtypes = []
    lib.rtgsfit_timing_reset.restype = None
    lib.rtgsfit_timing_get.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    lib.rtgsfit_timing_get.restype = None

    # Load warm-up output and derive all sizes from it (no extra MDS+ query)
    output_dict = np.load(cfg["rtgsfit_output_dict_path"], allow_pickle=True).item()
    last = cfg["n_iters"]          # last row: output after final iteration
    n_grid = output_dict["flux_norm"].shape[1]
    n_coef = output_dict["coef"].shape[1]
    n_lcfs_max = output_dict["lcfs_r"].shape[1]
    n_meas = output_dict["meas_model"].shape[1]
    n_xpt_max = output_dict["xpt_r"].shape[1]

    # Fixed inputs (same every call — raw sensor readings at this time slice)
    meas_pcs_fixed = output_dict["meas_pcs"][0].copy().astype(np.float64)
    coil_curr_fixed = output_dict["coil_curr"][0].copy().astype(np.float64)

    # Converged starting state
    flux_norm_conv = output_dict["flux_norm"][last].copy().astype(np.float64)
    mask_conv = output_dict["mask"][last].copy().astype(np.int32)

    # Pre-allocate output arrays (overwritten each repeat, never reset)
    flux_total = np.zeros(n_grid, dtype=np.float64)
    chi_sq_err = np.array([0.0], dtype=np.float64)
    lcfs_r = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_z = np.zeros(n_lcfs_max, dtype=np.float64)
    lcfs_n = np.array([0], dtype=np.int32)
    coef = np.zeros(n_coef, dtype=np.float64)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)
    lcfs_err_code = np.array([0], dtype=np.int32)
    lapack_info = np.array([0], dtype=np.int32)
    meas_model = np.zeros(n_meas, dtype=np.float64)
    r_mag = np.array([0.0], dtype=np.float64)
    z_mag = np.array([0.0], dtype=np.float64)
    mag_flux = np.array([0.0], dtype=np.float64)
    r_cen = np.array([0.0], dtype=np.float64)
    z_cen = np.array([0.0], dtype=np.float64)
    xpt_r = np.zeros(n_xpt_max, dtype=np.float64)
    xpt_z = np.zeros(n_xpt_max, dtype=np.float64)
    xpt_flux = np.zeros(n_xpt_max, dtype=np.float64)
    xpt_n = np.array([0], dtype=np.int32)
    xpt_div = np.array([0], dtype=np.int32)
    timing_buf = np.zeros(T_NTIMERS, dtype=np.float64)

    results = np.zeros((N_REPEATS, T_NTIMERS), dtype=np.float64)

    def ptr(arr, ctype):
        return arr.ctypes.data_as(ctypes.POINTER(ctype))

    for i in range(N_REPEATS):
        # Re-copy per-call inputs that rtgsfit may modify in place
        flux_norm = flux_norm_conv.copy()
        mask = mask_conv.copy()
        mp = meas_pcs_fixed.copy()
        cc = coil_curr_fixed.copy()

        lib.rtgsfit_timing_reset()
        lib.rtgsfit(
            ptr(mp, ctypes.c_double),
            ptr(cc, ctypes.c_double),
            ptr(flux_norm, ctypes.c_double),
            ptr(mask, ctypes.c_int32),
            ptr(flux_total, ctypes.c_double),
            ptr(chi_sq_err, ctypes.c_double),
            ptr(lcfs_r, ctypes.c_double),
            ptr(lcfs_z, ctypes.c_double),
            ptr(lcfs_n, ctypes.c_int32),
            ptr(coef, ctypes.c_double),
            ptr(flux_boundary, ctypes.c_double),
            ptr(plasma_current, ctypes.c_double),
            ptr(lcfs_err_code, ctypes.c_int32),
            ptr(lapack_info, ctypes.c_int32),
            ptr(meas_model, ctypes.c_double),
            ctypes.c_int32(n_meas),
            ptr(r_mag, ctypes.c_double),
            ptr(z_mag, ctypes.c_double),
            ptr(mag_flux, ctypes.c_double),
            ptr(r_cen, ctypes.c_double),
            ptr(z_cen, ctypes.c_double),
            ptr(xpt_r, ctypes.c_double),
            ptr(xpt_z, ctypes.c_double),
            ptr(xpt_flux, ctypes.c_double),
            ctypes.c_int32(n_xpt_max),
            ptr(xpt_n, ctypes.c_int32),
            ptr(xpt_div, ctypes.c_int32),
        )
        lib.rtgsfit_timing_get(
            ptr(timing_buf, ctypes.c_double),
            ctypes.c_int(T_NTIMERS),
        )
        results[i, :] = timing_buf

    return results


@pytest.mark.parametrize("pulse_num,time", test_cases)
def test_timing_rtgsfit(pulse_num, time):
    """
    Timing benchmark: compile, warm-up, then N_REPEATS calls to rtgsfit().
    Uploads per-section stats to GitHub Gist and prints a summary table.
    """
    cfg = config_loader.load_and_prepare_config(
        pulse_num=pulse_num,
        run_name=config_loader.next_test_run_name(52_000_000 + pulse_num),
    )
    cfg["time"] = time
    cfg["rt_timing"] = True

    # One-time setup (same as integration test)
    rtgsfit_compile_setup.initialise_rtgsfit_node(cfg)
    rtgsfit_compile_setup.compile_rtgsfit(cfg)
    # Warm-up convergence run — also saves output_dict used as starting state
    replay_rtgsfit.replay_rtgsfit(cfg)

    timing_us = _run_timing_benchmark(cfg)

    labels = [
        "meas_prep", "basis", "meas_matrix", "weighting", "ls_fit",
        "source", "model_meas", "chi2", "poisson", "coil_flux",
        "vessel_flux", "xpts_and_axis", "xpts_sort", "limiter",
        "lcfs", "inside", "normalise", "total",
    ]

    # Verify basic sanity
    assert timing_us.shape == (N_REPEATS, T_NTIMERS)
    assert np.all(timing_us >= 0.0)
    assert np.all(timing_us[:, -1] >= timing_us[:, :-1].sum(axis=1) * 0.9), \
        "total timer is unexpectedly smaller than sum of parts"

    # Print summary table
    print(f"\nRTGSFIT timing  pulse={pulse_num}  t={time}s  N={N_REPEATS}")
    print(f"{'section':<16}  {'mean':>8}  {'median':>8}  {'std':>8}  {'min':>8}  {'max':>8}  {'p95':>8}  {'p99':>8}  (us)")
    print("-" * 90)
    for i, label in enumerate(labels):
        col = timing_us[:, i]
        print(
            f"{label:<16}  {col.mean():8.1f}  {np.median(col):8.1f}  "
            f"{col.std():8.1f}  {col.min():8.1f}  {col.max():8.1f}  "
            f"{np.percentile(col, 95):8.1f}  {np.percentile(col, 99):8.1f}"
        )

    # Build CSV rows for this case
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    commit = _git_commit()
    host = socket.gethostname()
    csv_rows = ""
    for i, label in enumerate(labels):
        col = timing_us[:, i]
        csv_rows += (
            f"{ts},{commit},{host},{pulse_num},{time},{label},"
            f"{col.mean():.3f},{np.median(col):.3f},{col.std():.3f},"
            f"{col.min():.3f},{col.max():.3f},"
            f"{np.percentile(col, 95):.3f},{np.percentile(col, 99):.3f}\n"
        )

    _upload_to_gist(csv_rows)
