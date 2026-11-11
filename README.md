# RT-GSFit: Real-Time Grad-Shafranov Fit

RT-GSFit is a real-time tokamak equilibrium reconstruction code. It reconstructs the plasma
equilibrium within 1 ms of receiving live data from the magnetic sensors and PF coil power supplies.
It was deployed for shape control on the ST40 tokamak in the November 2025 experimental campaign.

<p align="center">
  <img src="https://raw.githubusercontent.com/aleksyprok/rtgsfit_media/refs/heads/main/SVGs/vessel_geometry_14757_t_100.svg" alt="ST40 Vessel Geometry" style="width:80%; height:auto;">
  <br>
  <em>Figure 1: Snapshot of the magnetic field during an ST40 pulse from the November 2025 campaign. RT-GSFit was used in real time to control the gap between the MCT/MCB (merging compression top/bottom) coils and the last closed flux surface (plasma boundary).</em>
</p>

RT-GSFit has two parts:

* **A C library** (`c/`) called once per Plasma Control System (PCS) cycle. It fits the plasma and
  passive currents to the measurements, solves the free-boundary Grad-Shafranov problem, and finds
  the magnetic axis, the x-points and the plasma boundary. It does not allocate memory. Every array
  size is a compile-time constant.
* **A Python package** (`python/rtgsfit`) that does everything else.
  * Before a pulse, it computes the **preshot**: all the expensive quantities that do not depend on
    the measurements, such as Greens functions and LU factors. It stores them in MDSplus and
    compiles them into the C library.
  * After a pulse, it **replays** the C library on the recorded measurements, exactly as the PCS
    ran it, and stores the results.
  * Offline post-processing adds fields, geometry and fit diagnostics after replay, or from a
    complete captured sequence of inputs/outputs. Exported flux is Wb, including the analytic case.

It has the same layout as its post-shot counterpart [GSFit](https://github.com/tokamak-energy/gsfit):
* Hierarchical settings (JSON with a shared JSON Schema)
* `database_readers/` and `database_writers/`
* results stored with GSFit's IMAS-style node names.

The method follows [J.-M. Moret et al., "Tokamak equilibrium reconstruction code LIUQE and its real
time implementation", Fusion Eng. Des. 91 (2015)](https://doi.org/10.1016/j.fusengdes.2014.09.019).
Coordinates are [COCOS 13](https://doi.org/10.1016/j.cpc.2012.09.010).
For every bundled case, including the analytic test, `psi` is **total poloidal flux in Wb**:
`Delta* psi = -2 pi mu0 R j_phi`, with
`B_R = -d(psi)/dZ / (2 pi R)` and `B_Z = d(psi)/dR / (2 pi R)`.
The analytic reference, synthetic flux-loop signals, preshot and C solver all use this same
normalisation. `PRESHOT.INFO:PSI_TO_WB` is 1 for every newly built preshot; the conversion factor
is retained only for reading older per-radian analytic archives.
See the [flux conventions](documentation/manual.md#2-conventions).

## Documentation

* [Manual](documentation/manual.md): covers the algorithm, the glossary, the C API, the preshot, the
  Python package, settings and storage, verification, and the migration from v1.
* [Data dictionary](documentation/data_dictionary.md): every quantity stored in MDSplus and in local
  files. This file is generated.
* [Style guide](documentation/style_guide.md): naming and coding rules for C and Python.
* [Source and variable reference](documentation/code_reference.md): what each module does and what
  the working variables mean, including the `local` NPZ readers/writers.
* [Negative-current assessment](documentation/negative_current.md): proposed work and acceptance
  tests; negative-current topology support is not implemented yet.
* [Examples](examples/): Jupyter notebooks that run RT-GSFit step by step.
  [Example 01](examples/example_01_run_full_pipeline_st40/README.md) runs the full ST40 pipeline: it
  builds the preshot with GSFit, replays a pulse and stores both in MDSplus.

## What works where

Some data and packages are only available inside the Tokamak Energy (TE) network. In particular, the
MDSplus server and TE data need it. Database access uses the pure-Python `mdsthin` client, not the
native `MDSplus` Python package or its shared libraries. `standard_utility` is needed for tree
creation and writing, with its `mdsthin` transport explicitly selected.

| Feature | Off the TE network | Needs |
|---|---|---|
| Build the C library and run the C unit tests (`make -C c test`) | ✅ | gcc, OpenBLAS (CBLAS and LAPACK) |
| Analytic replay and its tests (`rtgsfit replay analytic`) | ✅ | numpy, scipy, contourpy |
| Local preshot and results files (`.npz`), C preshot from a local file | ✅ | |
| Native NetCDF capture and offline post-processing | ✅ | `netCDF4` plus NetCDF-C development/runtime library for the shared C writer; no compiler needed for readback |
| FreeGSNKE MAST-U-like replay (`rtgsfit replay freegsnke`) | ✅ | [GSFit](https://github.com/tokamak-energy/gsfit) and [FreeGSNKE](https://github.com/FusionComputingLab/freegsnke), with FreeGSNKE's example machine files |
| ST40 preshot (`rtgsfit preshot st40`) | ❌ | GSFit and ST40 data in MDSplus |
| ST40 replay (`rtgsfit replay st40`) | ❌ | `mdsthin` and access to the MAG/ST40 and RTGSFIT trees |
| PPAC post-shot replay of one pulse (`rtgsfit.RtgsfitReplay`, [manual section 7.5](documentation/manual.md#75-post-shot-replays-in-ppac)) | ❌ | GSFit, `mdsthin`, ST40 data in MDSplus; `standard_utility` to write |
| Build the PCS library from MDSplus (`make -C c SHOT=... RUN_NAME=...`) | ❌ | MDSplus |
| Write to the `RTGSFIT` MDSplus tree | ❌ | MDSplus, `standard_utility`, write access |
| ST40 tests and RT-GSFit vs GSFit (`pytest --te-network`) | ❌ | all of the above |

## Quick start

These commands and CI use Python 3.13 via [uv](https://docs.astral.sh/uv/), which downloads the
interpreter if needed. The package minimum remains Python 3.11; other virtual-environment tools
work too.

    uv venv --python 3.13 && . .venv/bin/activate
    uv pip install -e ".[dev]"       # add [freegsnke] or [te] for those readers
    make -C c test                   # C unit tests
    pytest                           # off-network tests; TE and FreeGSNKE tests are skipped
    pytest --figures -m figures      # draw what the analytic and FreeGSNKE tests compare: tests/figures/<UTC time>/
    rtgsfit replay analytic          # writes the results to c/cache/rtgsfit_analytic.npz

On the TE network:

    rtgsfit replay st40 --set replay/pulse=13343 --set replay/time_end=0.2
    pytest --te-network
    pytest --te-network --figures -m figures            # the figures of every suite, ST40 vs GSFit too
    pytest --te-network --timing tests/test_timing.py   # rtgsfit() timing: c/cache/timing/*.csv

The figures are only drawn when asked for, never in CI. Each run writes PNG and PDF files to a new
directory, so earlier figures are kept; see the manual's Verification section.

The ST40 default reads the current-schema preshot from `RTGSFIT` pulse 52014922, `RUN22`.
The older `RUN21` lacks `PRESHOT.FIT:DOF_UNITS` and cannot be used with the current schema.

Configurations are in `python/rtgsfit/settings/<name>/rtgsfit.json`; a path to a custom JSON
file also works. VS Code provides completion, descriptions and validation through the shared
schema. CLI overrides use JSON values, for example `--set build/timing=true` or
`--set 'writers/names=["local"]'`; string values such as paths can remain unquoted.

Read the results:

```python
import numpy as np
results = np.load("c/cache/rtgsfit_analytic.npz")
results["GLOBAL.IP"], results["PROFILES_2D.R_Z.PSI"]   # [N_TIME], [N_TIME, N_Z, N_R]
```

To emulate the future PCS logger, save a native capture before offline processing:

```sh
rtgsfit replay analytic --set capture/path=analytic.nc
rtgsfit post-process c/cache/analytic.nc
```

The second command reads the captured C inputs/outputs and embedded preshot, then writes
`c/cache/analytic_processed.npz` without rebuilding C, rerunning the solver or reading a database.
Replay and PCS use the same C writer and versioned capture contract; legacy v1 NetCDF files are not reinterpreted.
See [the native capture specification](documentation/manual.md#76-native-netcdf-captures-replay-and-pcs).

## Repository layout

    c/src/                C library: rtgsfit.h is the public interface; the active preshot is generated here
    c/tools/              write_preshot.c, which generates the preshot C files from MDSplus or a file
    c/tests/              C unit tests
    c/build/              compiled library, objects and generator (not committed)
    c/cache/              local writer outputs (not committed)
    python/rtgsfit/       Python package: settings/, database_readers/, database_writers/, replay, preshot
    tests/                pytest: every integration test is a replay
    tests/figures/        diagnostic figures from pytest --figures, one directory per run (not committed)
    examples/             Jupyter notebooks that run RT-GSFit step by step, committed without outputs
    documentation/        manual, data dictionary, style guide

Private analysis notebooks and one-off replays live in a separate repository,
`rtgsfit_workbench`, so that this repository contains only the code, its tests and its examples.

## What a fresh clone contains

The repository includes the C/Python sources, JSON configurations and schema, documentation,
and test fixtures. It does **not** include a machine's generated C constants or a compiled solver.

| File or directory | In Git? | How it is obtained |
|---|---|---|
| `c/src/preshot.h` | Yes | Hand-maintained declarations and compile-time checks |
| `c/tools/write_preshot.c` | Yes | The C generator that writes the two files below |
| `c/src/preshot_sizes.h` | No | Generated compile-time size enums |
| `c/src/preshot_data.c` | No | Generated numerical arrays; replaces the old large constants source |
| `c/build/` | No | `librtgsfit.so`, objects and the generator |
| Replay/preshot `.npz` files | No, except committed test references | Generated by local writers |
| `c/cache/` | No | Local writer outputs |

The two generated files in `c/src/` are gitignored, just as the old generated constants source was.
They do not exist immediately after cloning. Do not edit them.

For an off-network first run:

```sh
uv venv --python 3.13 && . .venv/bin/activate
uv pip install -e ".[dev]"
make -C c test
rtgsfit replay analytic
```

There is one active preshot: the last one built. The replay writes it into `c/src/` and compiles
`c/build/librtgsfit.so`, inside this repository whatever the working directory. Make recompiles only
what changed, so replaying the same preshot again does not recompile; another preshot replaces it.
Run one build at a time. The library holds a CRC-32 of its preshot, and Python refuses a library
built from a different one.

Relative local-writer filenames, including those in the bundled settings, go under `c/cache/`; use an
absolute filename for an explicit export elsewhere. `RTGSFIT_CACHE_DIR` remains an optional override.
C unit tests need no preshot.

The C build itself needs no Python. The generator `c/tools/write_preshot.c` reads the preshot
either from MDSplus or from a file that Python writes:

```sh
make -C c SHOT=99000230 RUN_NAME=RUN19   # campaign preshot, after regeneration in the v2 schema
rtgsfit write-c-preshot analytic          # a preshot that is only in Python (analytic, FreeGSNKE, ...)
make -C c                                 # compile the active preshot
```

The PCS example uses preshot `RTGSFIT` 99000230, `RUN19`. Existing v1 data under that run must be
regenerated in the v2 schema first. The bundled `st40` replay preset still uses the separate
52014922 `RUN22` reference. The notebook and generator test use private pulse 52000230 `TEST03`,
not the shared 99000230 pulse. Tutorial writes use private development destinations; adapt them
to your own allocated pulse range.

`write-c-preshot` writes `c/build/preshot.bin` and runs `make PRESHOT_FILE=c/build/preshot.bin preshot`.
Replace `analytic` with `st40` on the TE network to read the configured stored ST40 preshot.
To calculate a new ST40 preshot using GSFit instead, add `--set preshot/reader=gsfit`.
These commands generate local files; they do not write to the shared database.
The generated files are machine/preshot-specific, not universal constants.

## Status

RT-GSFit v2 is a rewrite of the code that ran on ST40. The reconstructed equilibria are preserved:
the tests compare v2 with v1, with GSFit, with FreeGSNKE and with an analytic solution. All names
and the PCS interface have changed. The [manual](documentation/manual.md#11-migrating-from-v1) lists
every change.

The code is public, but some setup is still specific to ST40. If you would like to use RT-GSFit on
another machine, please get in touch.

## License

See [LICENSE](LICENSE).
