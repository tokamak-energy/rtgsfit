# Example 02: build the ST40 PCS preshot (`st40_tepcs`)

[example_02_build_preshot_st40_tepcs.ipynb](example_02_build_preshot_st40_tepcs.ipynb) rebuilds the ST40
campaign preshot **RTGSFIT 99000230 / PRESHOT04**, which the PCS compiles, from GSFit with the packaged
`st40_tepcs` recipe:

1. Show the recipe: GSFit settings, grid, PCS sensors and coil circuits.
2. Build the preshot with GSFit (`build_preshot`).
3. Inspect it: the inputs the PCS supplies, the grid, and the coil and passive Greens functions.
4. Compare it, node by node, with the stored PRESHOT04.
5. Save it as a local NPZ and a C preshot file, compile `librtgsfit.so` from the file, and check the
   compiled checksum.
6. Optionally write it to a new private MDSplus run, read it back, and compile that run with C only,
   as the PCS does.
7. Describe how a new campaign preshot is published for the PCS.

The notebook never writes to 99000230. Publishing a campaign preshot changes what the PCS compiles, so
it is a deployment step agreed with the PCS team, not a side effect of this example.
[Example 01](../example_01_run_full_pipeline_pcs_mock/README.md) runs the C solver on a pulse with PRESHOT04.

## Files and database access

| Item | Behaviour |
|---|---|
| `RTGSFIT` 99000230 / `PRESHOT04` | Read only, for the comparison |
| `c/cache/example_02_build_preshot/<UTC>/` | New directory per run, with the preshot as `.npz` and as a C preshot file (`.bin`) |
| `WRITE_MDSPLUS = False` | Default: no database writes |
| Optional destination | Private pulse 52014922 / `EXAMPLE02`; choose your own allocated private pulse and an unused run |
| Existing runs and `BEST` | Never overwritten or changed |

With `WRITE_MDSPLUS = True`, the notebook checks that the run is free before building and again
before writing, refuses 99000230 and `BEST`, and compares every node read back.

## Reproducing PRESHOT04 exactly

The comparison accepts round-off in floating-point nodes (1e-12 of each node's largest value) and
reports whether the build is bit-for-bit identical. Only a bit-for-bit identical build has the same C
checksum, and so compiles to the same library, as PRESHOT04.

The Poisson LU factors come from a threaded LAPACK factorisation, whose rounding can depend on the
machine and the number of threads. PRESHOT04 was built with the default number of OpenBLAS threads, so
the notebook does not set `OPENBLAS_NUM_THREADS`. With `OPENBLAS_NUM_THREADS=1`, the LU factors differ
by about 3e-16 of their size and the checksum changes.

## Requirements

- TE network access to `RTGSFIT` and the ST40 machine description on `smaug`.
- The RT-GSFit Python package with its TE extras, which include
  [GSFit](https://github.com/tokamak-energy/gsfit) and `mdsthin`, plus `shapely`, Matplotlib and Jupyter.
  From the repository root: `uv pip install -e ".[dev,te]" shapely jupyter`.
- A C11 compiler, GNU make, `cmp` (diffutils), OpenBLAS with CBLAS/LAPACKE, NetCDF-C, and the MDSplus
  C development/runtime libraries.
- Only for the optional MDSplus write: `standard_utility`, database write access and its matching
  `RTGSFIT.csv` (checked against `rtgsfit.data_dictionary`).

The C build uses the checkout's `c/` by default (or `RTGSFIT_C_DIR` if set).
Only one build/replay may use a C directory at a time. For an independent run, point
`RTGSFIT_C_DIR` to a private copy of the C sources before starting the notebook kernel.

## Run

Open the notebook, select the RT-GSFit environment, restart the kernel and run all cells.
Use this directory as the working directory so `import helpers` works.

Alternatively, from the repository root:

```sh
mkdir -p c/cache/example_02_build_preshot
.venv/bin/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 \
    --output-dir c/cache/example_02_build_preshot \
    examples/example_02_build_preshot_st40_tepcs/example_02_build_preshot_st40_tepcs.ipynb
```

This writes an executed copy to `c/cache/example_02_build_preshot/`, beside each run's timestamped
output directory. The whole notebook takes about a minute, mostly in the GSFit build; compiling
takes seconds. The committed notebook has no outputs.

[helpers.py](helpers.py) holds the checks, the comparison, the GSFit geometry and the plots. The
[manual](../../documentation/manual.md#5-the-preshot) describes the preshot and its storage.
