# Example 01: full pipeline with a mock PCS

[example_01_run_full_pipeline_pcs_mock.ipynb](example_01_run_full_pipeline_pcs_mock.ipynb)
uses **RTGSFIT 99000230 / PRESHOT10** and recorded inputs from **ST40 pulse 14926**:

1. Read the stored preshot; no GSFit preshot build or database write.
2. Run the Python-free C build from MDSplus and check the compiled preshot checksum.
3. Read calibrated TEPCS inputs and call the C solver once per sample, buffering every output.
4. After the loop, write the native NetCDF-4 capture using the shared C writer.
5. Read the capture, compute offline diagnostics, save local NPZ results and check their readback.
6. Plot the results; optionally write/read back a private MDSplus results run.

Python mocks the PCS caller; this is not the PCS/Simulink integration or a real-time benchmark.
Samples are 1 ms apart, not the 50-microsecond PCS cycle. Available PSS calibration and correction
settings are used; missing corrections are reported rather than claiming exact historical PCS inputs.
PRESHOT10 supports pulses from 12551 onward.

## Files and database access

| Item | Behaviour |
|---|---|
| `RTGSFIT` 99000230 / `PRESHOT10` | Read only; never rebuilt or overwritten |
| `c/cache/example_01_pcs_mock/<UTC>/` | New directory per run containing native `.nc` and processed `.npz` files |
| `WRITE_MDSPLUS = False` | Default: no database writes |
| Optional result destination | Private pulse 52014926 / `PCS_MOCK01`; choose your own allocated private pulse and unused run |
| `BEST` | Never changed |

With `WRITE_MDSPLUS = True`, the notebook checks that the destination run is free, writes only
processed results, and compares every node read back. It refuses the shared preshot pulse and `BEST`.
It does not publish a new preshot. The network C-generator test still reads its separate private
fixture, 52000230 / `TEST04`.

The C build uses the checkout's `c/` by default (or `RTGSFIT_C_DIR` if set).
Only one build/replay may use a C directory at a time. For an independent run, point
`RTGSFIT_C_DIR` to a private copy of the C sources before starting the notebook kernel.

## Requirements

- TE network access to `RTGSFIT`, `TEPCS` and `PSS` on `smaug`.
- The RT-GSFit Python package, `mdsthin`, Matplotlib and Jupyter.
  From the repository root: `uv pip install -e ".[dev]" jupyter`.
  **GSFit is not required**: the preshot is already stored.
- A C11 compiler, GNU make, `cmp` (diffutils), OpenBLAS with CBLAS/LAPACKE,
  NetCDF-C with NetCDF-4 support, and the MDSplus C development/runtime libraries.
- Only for optional MDSplus result writes: `standard_utility`, database write access and its matching
  `RTGSFIT.csv` (checked against `rtgsfit.data_dictionary`).

## Run

Open the notebook, select the RT-GSFit environment, restart the kernel and run all cells.
Use this directory as the working directory so `import helpers` works.
The first cell selects one OpenBLAS thread before NumPy or the solver is loaded.

Alternatively, from the repository root:

```sh
mkdir -p c/cache/example_01_pcs_mock
cp examples/example_01_run_full_pipeline_pcs_mock/preshot_to_rtgsfit.svg c/cache/example_01_pcs_mock/
.venv/bin/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 \
    --output-dir c/cache/example_01_pcs_mock \
    examples/example_01_run_full_pipeline_pcs_mock/example_01_run_full_pipeline_pcs_mock.ipynb
```

This writes a separate executed copy without replacing the sample outputs in the source notebook.
Network reads and compilation dominate runtime; no GSFit inverse solve is performed.
The notebook prints the exact output paths.

## Saved example outputs

As a one-off exception, this notebook retains an approved run for ST40 pulse 14926 with PRESHOT10,
including its seven plots and printed summaries, so it can be read directly on GitHub.
These are saved results, not a live execution; timings are illustrative, not PCS qualification.
Machine-specific checkout paths in the saved text are shown relative to the repository root.
The raw NetCDF capture and NPZ results remain local and are not committed.
This exception does not change the output policy for other notebooks or future runs.

[helpers.py](helpers.py) holds checks and plots; [preshot_to_rtgsfit.svg](preshot_to_rtgsfit.svg)
shows the data flow. See the [manual](../../documentation/manual.md#76-native-netcdf-captures-replay-and-pcs)
for the capture contract and command-line post-processing.
