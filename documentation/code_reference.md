# Source and variable reference

Start with the [manual](manual.md) for the algorithm and the generated
[data dictionary](data_dictionary.md) for every stored variable, its shape, units and description.
Names follow GSFit's object-then-property order: `current_centroid_r`, not `r_current_centroid`.

## Directories

| Directory | Purpose |
|---|---|
| `c/src/` | Allocation-free reconstruction and its checked public/private interfaces |
| `c/tests/` | Independent numerical and invalid-buffer tests; fixtures in `data/` |
| `c/tools/` | `write_preshot.c`: generates the preshot C files from MDSplus or a preshot file, without Python |
| `python/rtgsfit/` | Importable Python package; `python/` separates languages, `rtgsfit/` defines the import name |
| `python/rtgsfit/settings/` | JSON loader, shared schema and named configurations; not Python runtime state |
| `python/rtgsfit/database_readers/` | Machine/preshot builders and sources of recorded or synthetic input signals |
| `python/rtgsfit/database_writers/` | Serialisation of the same dictionary to disk, generated C, or the TE database |
| `tests/` | Python unit tests and end-to-end replays; network tests opt in |
| `documentation/` | Manual, conventions, generated dictionary and historical derivations |
| `.github/workflows/` | CI validation |
| `.vscode/` | Editor formatting configuration |

`local` means **a local NPZ file**, not a local database service. The reader loads a preshot or
results archive; the writer stores a pickle-free NPZ archive with exactly the dictionary keys
used by the database writer. It needs no TE connection.

## Python modules

All paths in this table are relative to `python/rtgsfit/`.

| File or package | Responsibility |
|---|---|
| `__init__.py`, `py.typed` | Package version and typed-package marker |
| `__main__.py` | CLI commands: replay, native-capture post-processing, preshot, dictionary and CSV generation |
| `constants.py` | Physical constants, fixed model sizes and model coefficient labels |
| `data_dictionary.py` | Single source of stored paths, units, shapes and generated CSV/manual tables |
| `preshot.py` | Assemble/validate preshot matrices, initial state and generated C contract |
| `library.py` | Build the active C library with make, reload it when rebuilt, check ctypes buffers, run one iteration, decode status; `isolated_build()` compiles in a private temporary copy of `c/` |
| `replay.py` | Orchestrate input reading, iteration, offline collection, then output writers |
| `ppac.py` | `RtgsfitReplay`: PPAC's post-shot replay of one ST40 pulse (GSFit preshot, MAG `BEST` replay over the PCS window, MDSplus write to `RUN01_POST`) |
| `capture.py` | ctypes bridge to the shared atomic C NetCDF writer, plus versioned solver-free readback |
| `results.py` | Collect complete captured slices and validate the result dictionary |
| `post_processing.py` | Offline fit diagnostics, poloidal fields, geometry and export-unit conversion |
| `boundary.py` | Contour tracing and boundary shape; never part of the real-time solver |
| `settings/__init__.py` | Strict JSON parsing, schema validation, effective settings serialisation and overrides |
| `settings/settings.schema.json` | Types, allowed values, descriptions and units shared by runtime validation and VS Code |
| `settings/analytic/rtgsfit.json` | Self-contained circular analytic replay |
| `settings/st40/rtgsfit.json` | TE-network replay and GSFit preshot builder |
| `settings/freegsnke/rtgsfit.json` | FreeGSNKE example machine settings |
| `database_readers/__init__.py` | Reader interfaces, registration and `Signals` contract |
| `database_readers/analytic/__init__.py` | Synthetic inputs and analytic preshot assembly |
| `database_readers/analytic/solution.py` | Independent analytic equilibrium reference |
| `database_readers/analytic/greens.py` | Large-aspect-ratio filament response functions |
| `database_readers/gsfit/__init__.py` | Use GSFit machine descriptions to build a preshot |
| `database_readers/freegsnke/__init__.py` | Adapt FreeGSNKE machine data via the GSFit builder |
| `database_readers/local/__init__.py` | Read a local NPZ archive |
| `database_readers/rtgsfit_mdsplus/__init__.py` | Read RTGSFIT preshots through mdsthin, including explicit legacy-schema mode |
| `database_readers/st40_mdsplus/__init__.py` | Read PCS-equivalent ST40 magnetic and PSU signals through mdsthin |
| `database_writers/__init__.py` | Writer interfaces and dispatch |
| `database_writers/local/__init__.py` | Write the same dictionary as a local NPZ archive |
| `database_writers/c_source/__init__.py` | Write the preshot file for the C generator and make it the active preshot; CRC-32 checksum |
| `database_writers/rtgsfit_mdsplus/__init__.py` | Check CSV consistency and write via standard_utility with mdsthin explicitly selected |

## C modules

Each implementation has a matching header. Header comments define buffer sizes, units and return
values; `argument_check.h`, `grid.h` and `preshot.h` are header-only contracts.

| Stem in `c/src/` | Responsibility |
|---|---|
| `rtgsfit` | Wire the stages of one Picard update; sole consumer of compiled preshot constants |
| `measurement_processing` | Replace sensors, subtract known coil response and weight the fit target |
| `current_basis` | Form the three current-density basis functions from the previous good state |
| `least_squares` | Allocation-free LAPACK SVD solve and classification of raw LAPACK info |
| `finite_difference` | Checked grid derivatives and interpolation |
| `magnetic_nulls` | Locate/classify critical points and test core visibility |
| `plasma_boundary` | Select boundary flux, flood-fill the core and normalise flux |
| `poisson_solver` | Banded-LU fixed/free-boundary solves, Hagenow boundary construction |
| `timing` | Optional per-stage timing; no output-file or database operations |
| `argument_check` | Shared array, exact-length and minimum-length checks |
| `grid` | Grid structure, dimensions, checked coordinate buffers and index helpers |
| `preshot` | Compile-time size relations and extern declarations of generated arrays |

## Working variables not stored as independent nodes

| Name | Meaning |
|---|---|
| `fit_matrix` | Weighted response, `[N_FIT_TARGET, N_FIT_DOF]`, column-major in C |
| `fit_target` | Weighted measured constraints minus coil response, followed by zero regularisation targets |
| `fit_dof` | Unknown plasma amplitudes followed by passive-current amplitudes; these are the degrees of freedom |
| `fit_prediction` | `fit_matrix @ fit_dof`, in weighted target-row units |
| `fit_matrix_lapack`, `fit_target_lapack` | Disposable copies: LAPACK overwrites its inputs |
| `singular_values`, `rank`, `least_squares_work` | SVD diagnostics and caller-owned LAPACK scratch storage |
| `basis` | Three current-density basis arrays; constructed from the previous `psi_norm` and `mask` |
| `constraint_coil` | Known coil response at each fitted magnetic constraint |
| `constraint_measured` | Sensor-replaced values actually used as constraints; length `N_CONSTRAINT` |
| `meas` | Raw PCS measurements, including extra replacement signals; length `N_MEAS` |
| `psi_plasma` | Plasma contribution before adding known coil and fitted passive flux |
| `limiter_psi` | Candidate boundary flux on visible limiter points |
| `work` | Static C workspace, reused on each call; not shared safely between simultaneous callers |
| `i_*`, `j_*` | Loop indices; `n_*` are checked run-time counts, `N_*` compile-time counts |

`chi_sq_magnetic` and `chi_sq_total` are sums, not reduced chi-square statistics. The total adds
regularisation penalties to the magnetic sum. Their offline family/per-row breakdown is stored
alongside the original C scalars.
