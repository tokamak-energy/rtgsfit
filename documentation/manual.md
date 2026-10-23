# RT-GSFit manual

RT-GSFit (Real-Time Grad-Shafranov Fit) reconstructs the plasma equilibrium of a tokamak, in well
under a millisecond, from magnetic measurements and PF coil currents. It ran in the ST40 Plasma
Control System (PCS) for shape control in the November 2025 campaign. This manual describes
RT-GSFit v2: the algorithm, the names used for every quantity, the C library, the Python package
that builds the preshot and replays pulses, and how data is stored.

Contents

1. [Overview](#1-overview)
2. [Conventions](#2-conventions)
3. [Glossary](#3-glossary)
4. [The algorithm](#4-the-algorithm)
5. [The preshot](#5-the-preshot)
6. [The C library](#6-the-c-library)
7. [The Python package](#7-the-python-package)
8. [Settings](#8-settings)
9. [Storage: MDSplus and local files](#9-storage-mdsplus-and-local-files)
10. [Verification](#10-verification)
11. [Migrating from v1](#11-migrating-from-v1)
12. [Future work](#12-future-work)
13. [References](#13-references)

The names of every stored quantity are listed in [data_dictionary.md](data_dictionary.md), which is
generated from `python/rtgsfit/data_dictionary.py`. The coding style is in
[style_guide.md](style_guide.md).

---

## 1. Overview

RT-GSFit is split in two, because a real-time code must do as little as possible while the plasma
is running:

| | Preshot (before the pulse) | Real time (every PCS cycle) |
|---|---|---|
| Language | Python (`python/rtgsfit`) | C (`c/src`) |
| Inputs | Machine description, sensor geometry, settings | `meas`, `coil_current`, last `psi_norm` and `mask` |
| Work | Greens functions, LU factors of the Grad-Shafranov operator, weighted response matrices, ... | One least-squares fit, one Poisson solve, null and boundary search |
| Output | The preshot: `PRESHOT.*` in MDSplus, compiled into C as constants | The equilibrium of this cycle |
| Allocation | Any | None: all memory is static |

The preshot is computed by a **database reader** (GSFit for ST40, FreeGSNKE + GSFit for the
MAST-U-like example, or the built-in analytic model), stored by a **database writer** (MDSplus,
local `.npz` or C source), and then compiled into `librtgsfit.so`, or into the PCS, as `const`
arrays and `enum` sizes. RT-GSFit has the same layout as [GSFit](https://github.com/tokamak-energy/gsfit):
XML settings, `database_readers/` and `database_writers/`, and results stored with GSFit's
IMAS-style node names.

A **replay** runs the compiled C library from Python on recorded (or synthetic) measurements,
exactly as the PCS would, and stores every input and output. The tests are replays (section 10).

RT-GSFit iterates *in time*: each call to `rtgsfit()` does one Picard iteration, starting from the
`psi_norm` and `mask` of the previous call. At the PCS cycle rate the equilibrium changes little
from one call to the next, so one iteration per cycle is enough to track it. A replay with
`replay/repeat > 1` calls `rtgsfit()` several times per time slice, to study convergence.

## 2. Conventions

**Coordinates.** Right-handed cylindrical (R, phi, Z), COCOS 13 [Sauter and Medvedev 2013]. Flux
`psi` is in weber for the machine cases (ST40, FreeGSNKE); the analytic case uses weber per radian
(see its settings). The poloidal angle increases counterclockwise from the outboard mid-plane.

**Units.** SI everywhere: m, A, Wb, T, A/m², s. Units are stored with every node.

**Grid.** The grid is rectangular with `N_R` columns at `GRID:R` and `N_Z` rows at `GRID:Z`,
uniformly spaced by `GRID:D_R` and `GRID:D_Z`. Grid arrays are stored flat, row-major `[N_Z, N_R]`:
grid point `(i_r, i_z)` has the flat index `i_z * N_R + i_r`. In Python,
`psi.reshape(N_Z, N_R)` gives the 2D array, and the results already store it as `[N_TIME, N_Z, N_R]`.

**Edge order.** The grid edge has `N_EDGE = 2 N_R + 2 N_Z - 4` points. They start at the
bottom-left corner `(R_min, Z_min)`, go up the left side, right along the top, down the right side
and left along the bottom: clockwise in the (R, Z) plane. v1 called the edge `LTRB`
(left-top-right-bottom); the order is unchanged.

**Matrices.** Written `[rows, columns]` and stored row-major in C and Python, except `a_matrix`,
which is column-major because LAPACK expects that.

**Sizes.** Every array size has a name `N_...` (section 6.3), and every C function takes an array
together with its length.

## 3. Glossary

The names follow the least-squares notation of LAPACK and Wikipedia, `a_matrix x_vec ≈ b_vec`
(issue #68). X-points are `xpt`, so that `x` always means the unknowns.

| Symbol (code name) | Size | Definition |
|---|---|---|
| `meas` | `N_MEAS` | PCS measurements passed to `rtgsfit()`, in the order of `PRESHOT.MEAS:NAMES`. They include signals used only to replace faulty sensors. (v1: `meas_pcs`) |
| `coil_signal` | len(`COIL:SIGNALS`) | PCS coil signals, in the order of `PRESHOT.COIL:SIGNALS` (PSU currents, Rogowski wires). Python only. |
| `coil_current` | `N_COIL` | PF coil circuit currents: `coil_current = COIL:MATRIX coil_signal`. (v1: `coil_curr`) |
| `constraint_measured` | `N_CONSTRAINT` | Constraint values from the measurements: `CONSTRAINT:REPLACEMENT meas`. Flux loops, then Bp probes, then Rogowski coils. |
| `constraint_reconstructed` | `N_CONSTRAINT` | The same constraints computed from the fitted equilibrium. |
| `b_vec` | `N_B_VEC` | Right-hand side of the least-squares problem: the weighted constraints minus the coil contribution, then `N_REG` zeros for the regularisation rows. (v1: `meas`) |
| `a_matrix` | `[N_B_VEC, N_X_VEC]` | The weighted response of every row of `b_vec` to every entry of `x_vec`. (v1: `g_coef_meas_w`) |
| `x_vec` | `N_X_VEC` | The unknowns: `N_PLASMA` plasma degrees of freedom, then `N_PASSIVE` passive degrees of freedom. (v1: `coef`) |
| `b_pred` | `N_B_VEC` | `a_matrix x_vec`, the prediction of `b_vec`. (v1: `meas_model`) |
| `basis` | `[N_PLASMA, N_GRID]` | Plasma current basis functions (section 4.2). |
| `j_phi` | `N_GRID` | Toroidal plasma current density. |
| `psi` | `N_GRID` | Total poloidal flux: coils + passives + plasma. (v1: `flux_total`) |
| `psi_plasma` | `N_GRID` | Flux of the plasma current alone. (v1: `flux_pls`) |
| `psi_norm` | `N_GRID` | Normalised flux: 0 on the magnetic axis, 1 on the boundary and outside the plasma. (v1: `flux_norm`) |
| `mask` | `N_GRID` | 1 inside the plasma core, 0 outside. |
| `axis_psi`, `axis_r`, `axis_z` | scalar | The magnetic axis: the o-point with the largest `psi`. (v1: `mag_axis_flux`, `r_mag_axis`, `z_mag_axis`) |
| `xpt_r`, `xpt_z`, `xpt_psi`, `xpt_n` | `N_XPT_MAX` | X-points, sorted by descending `psi`; only the first `xpt_n` are valid. (v1: `xpt_flux`) |
| `limiter_psi` | scalar | The largest `psi` on the limiter points visible from the axis. |
| `boundary_psi` | scalar | Flux on the plasma boundary. (v1: `flux_boundary`, `lcfs_flux`) |
| `xpt_diverted` | scalar | 1 if an x-point sets the boundary, 0 if the plasma is limited. |
| `ip` | scalar | Plasma current. (v1: `plasma_current`) |
| `chi_sq_magnetic`, `chi_sq_total` | scalar | `sum (b_vec - b_pred)²` over the constraint rows, and over every row. (v1: `chi_sq_err`) |
| `edge` | `N_EDGE` | The grid points on the edge of the grid. (v1: `LTRB`) |
| `G_...` | | Greens functions: `X:G_Y` is the response of X to unit Y, e.g. `GRID:G_COIL` is `psi` on the grid per ampere in each coil. |

## 4. The algorithm

The method follows LIUQE [Moret et al. 2015]. One call of `rtgsfit()` does the 12 steps below
(they are numbered in `c/src/rtgsfit.c` and timed separately, section 6.4).

### 4.1 Constraints and `b_vec`

Faulty sensors are replaced by combinations of other signals (for example a broken Rogowski coil by
`I = V / R` from neighbouring flux loop voltages). The replacement is the matrix
`CONSTRAINT:REPLACEMENT [N_CONSTRAINT, N_MEAS]`; an identity row keeps a sensor as it is:

    constraint_measured = CONSTRAINT:REPLACEMENT meas

The coils' contribution is known, so it is subtracted, and every row is weighted by
`FIT:WEIGHT = 1 / uncertainty`:

    b_vec[i] = FIT:WEIGHT[i] (constraint_measured[i] - (CONSTRAINT:G_COIL coil_current)[i])   for i < N_CONSTRAINT
    b_vec[i] = 0                                                                              for the N_REG regularisation rows

### 4.2 Current model

Inside the plasma core (`mask` = 1 and `GRID:MASK_LIMITER` = 1) the toroidal current density is

    j_phi = x_vec[0] basis_0 + x_vec[1] basis_1 + x_vec[2] basis_2

| `x_vec` name | Basis function | Meaning |
|---|---|---|
| `P_PRIME` | `(1 - psi_norm) R` | p' term of the Grad-Shafranov equation, linear in `psi_norm` |
| `FF_PRIME` | `(1 - psi_norm) / (mu0 R)` | FF' term |
| `DELTA_Z` | `∂psi_norm / ∂Z` | A rigid vertical shift of the current; stabilises the vertical position |

and zero outside. The `psi_norm` and `mask` are those of the previous call, so the problem is
linear in `x_vec`. This is the transpose of `T_yg` in eq. 61 of Moret et al.

The passive structures (vessel, passive plates, ...) have `N_PASSIVE` degrees of freedom: currents
in individual passive elements, or eigenmodes of the vessel such as `IVC_EIG_01`. They are the last
entries of `x_vec`. Regularisation rows (`N_REG`) penalise passive currents that the measurements
do not constrain.

### 4.3 The least-squares fit

    a_matrix = [ basis FIT:A_GRID | FIT:A_PASSIVE^T ]        [N_B_VEC, N_X_VEC]

`FIT:A_GRID [N_GRID, N_B_VEC]` is the weighted response of `b_vec` to the current in each grid cell,
including the cell area `d_r d_z`, so the plasma columns are a single matrix product. The passive
columns do not depend on the plasma and are fixed at preshot. Then

    minimise |a_matrix x_vec - b_vec|²

is solved with LAPACK `dgelss` (SVD). The workspace is static and generous enough for LAPACK's
blocked algorithm, so the result is identical to `LAPACKE_dgelss()` but with no allocation.

### 4.4 Plasma current and residual

    j_phi = sum_k x_vec[k] basis_k,     ip = d_r d_z sum j_phi
    current_centroid_r = sum(j_phi R) / sum(j_phi)   (only if ip > SETTINGS:IP_MIN)
    b_pred = a_matrix x_vec
    chi_sq_magnetic = sum over constraint rows of (b_vec - b_pred)²
    constraint_reconstructed = CONSTRAINT:G_COIL coil_current + b_pred / FIT:WEIGHT

### 4.5 Poisson solve: the plasma flux

The plasma flux satisfies

    Δ* psi_plasma = -mu0 R j_phi,    Δ* = R ∂/∂R (1/R ∂/∂R) + ∂²/∂Z²

with free-boundary conditions (psi produced by `j_phi` alone). Following von Hagenow and Lackner
(Moret et al. section 4.2):

1. Solve with `psi_plasma = 0` on the edge (a fixed-boundary solve).
2. Its normal derivative on the edge is a surface current. The edge Greens function
   `EDGE:G_EDGE` turns that into the true `psi_plasma` on the edge.
3. Solve again with that Dirichlet boundary condition.

The finite-difference matrix (5-point stencil, identity rows on the edge, scaled by `d_z²`) is
LU-factorised at preshot, `P A = L U`. Its factors are banded with bandwidth `N_R`, so they are
stored compactly as `POISSON:LOWER_BAND [N_GRID, N_R]`, `POISSON:UPPER_BAND [N_GRID, N_R + 2]`
(column 0 is `1 / diag(U)`) and `POISSON:PERMUTATION`. Each solve is two banded triangular
solves. `POISSON:RHS_SCALE = mu0 R d_z²` turns `j_phi` into the right-hand side.

### 4.6 Total flux

    psi = GRID:G_COIL coil_current + GRID:G_PASSIVE x_vec[N_PLASMA:] + psi_plasma

### 4.7 Magnetic axis and x-points

Every grid cell inside the limiter is examined from each of its four corners: `psi` is expanded to
second order about the corner, and the stationary point of the expansion is accepted if it lies in
the cell. A positive Hessian determinant is an o-point, a negative one an x-point;
`|det| > SETTINGS:TOLERANCE` is required. The o-point with the largest `psi` is the magnetic axis.

X-points hidden behind another x-point, as seen from the axis, are removed: x-point `j` is removed
unless it is on the core side of every other x-point `i`, where "core side" is the half-plane
through `i` perpendicular to the line from `i` to the axis. The rest are sorted by descending `psi`.

### 4.8 Plasma boundary

`limiter_psi` is the largest `psi`, interpolated bilinearly (`LIMITER:GRID_INDEX`,
`LIMITER:GRID_WEIGHT`), over the limiter points that are on the core side of every x-point. Then

    diverted_psi = XPT_FRACTION xpt_psi[0] + (1 - XPT_FRACTION) axis_psi
    diverted if diverted_psi > limiter_psi:  boundary_psi = diverted_psi,  else boundary_psi = limiter_psi

`SETTINGS:XPT_FRACTION` (e.g. 0.999) sets the boundary just inside the separatrix, so the flood
fill does not leak out through the x-point.

### 4.9 Plasma core and normalised flux

The core is flood-filled (4-connected) from the grid point nearest the axis, over points with
`psi > boundary_psi`, inside the limiter, and on the core side of every x-point. Then

    psi_norm = (psi - axis_psi) / (boundary_psi - axis_psi)   inside the core, 1 elsewhere

`psi_norm` and `mask` are the input of the next call. If any step fails, they are left unchanged,
so the next call restarts from the last good equilibrium.

### 4.10 Status

`rtgsfit()` returns a bitwise OR of these flags; 0 is success. Values 1 to 32 are those of v1's
`lcfs_err_code`.

| Value | Name | Meaning |
|---|---|---|
| 0 | `RTGSFIT_OK` | Success |
| 1 | `RTGSFIT_ERR_AXIS_EQUALS_BOUNDARY` | `\|boundary_psi - axis_psi\| < SETTINGS:TOLERANCE` |
| 2 | `RTGSFIT_ERR_NO_AXIS` | No magnetic axis found (for example, no plasma current) |
| 4 | `RTGSFIT_ERR_TOO_MANY_XPTS` | `N_XPT_MAX` x-points found before the search ended |
| 8 | `RTGSFIT_ERR_TOO_MANY_OPTS` | Reserved (v1 compatibility) |
| 16 | `RTGSFIT_ERR_AXIS_OUTSIDE_CORE` | The grid point nearest the axis is not in the core |
| 32 | `RTGSFIT_ERR_BOUNDARY_ABOVE_AXIS` | `boundary_psi > axis_psi` |
| 64 | `RTGSFIT_ERR_ARGUMENT_SIZE` | An array length does not match the preshot, or a null pointer |
| 128 | `RTGSFIT_ERR_LAPACK` | `dgelss` returned non-zero `lapack_info` |

They are stored at `GLOBAL:STATUS`; `rtgsfit.library.describe_status()` turns a value into names.

### 4.11 Post-processing (Python only)

The replay adds what the PCS does not need: the boundary outline (the `psi = boundary_psi` contour
around the axis, with `contourpy`), and from it the geometric axis, minor radius, elongation and
triangularities (`BOUNDARY.*`). These are not computed in C.

## 5. The preshot

The preshot is every quantity that does not depend on the measurements. It is defined node by node
in `data_dictionary.py` and stored at `PRESHOT.<GROUP>:<NODE>`; in C the same quantity is the
constant `<GROUP>_<NODE>`.

| Group | Contents |
|---|---|
| `SIZES` | Array sizes; `enum` constants in C |
| `GRID` | Grid, limiter mask, coil and passive Greens functions on the grid |
| `EDGE` | Edge points, and the edge Greens function for the free-boundary condition |
| `LIMITER` | Limiter points and their interpolation stencils |
| `POISSON` | LU factors of the discrete Grad-Shafranov operator |
| `CONSTRAINT` | Sensor replacement and coil Greens functions of the constraints |
| `FIT` | Weights and weighted response matrices of the least-squares problem |
| `SETTINGS` | `XPT_FRACTION`, `TOLERANCE`, `IP_MIN` |
| `COIL`, `MEAS` | Names and the coil matrix; used to build the inputs of `rtgsfit()` |
| `INITIAL` | `psi_norm` and `mask` for the first call |
| `INFO` | The settings XML and the reader that made the preshot |

Machine-specific readers provide only the physics (`preshot.PreshotInputs`: Greens functions,
weights, geometry); `preshot.assemble_preshot()` derives everything else (LU factors, weighted
matrices, interpolation stencils, ...) in one place, and `Preshot` checks every node against the
data dictionary and every size relation (the same relations as the `_Static_assert`s in
`c/src/preshot.h`).

**Making a preshot for ST40** (TE network):

    python -m rtgsfit preshot st40 --set writers/names=rtgsfit_mdsplus \
        --set writers/rtgsfit_mdsplus/enabled=true --set writers/rtgsfit_mdsplus/run_name=RUN22

This runs GSFit with the `st40_setup_for_rtgsfit` settings (`preshot/builder = gsfit`) and writes
`PRESHOT` to the `RTGSFIT` tree. The PCS build then reads it:

    make -C c SHOT=52014922 RUN_NAME=RUN22

which calls `python -m rtgsfit write-c-preshot --pulse 52014922 --run-name RUN22`, writing
`c/preshot/preshot_sizes.h` and `c/preshot/preshot_data.c`, and compiles `c/build/librtgsfit.so`.

ST40 preshots in `RTGSFIT` pulse 52014922: `RUN19` was written by an early v1, `RUN20` by v1, and
`RUN21` by v2. v1 and v2 trees are not interchangeable (the node names changed, section 11).

## 6. The C library

### 6.1 Layout

| File | Contents |
|---|---|
| `rtgsfit.h` | The public interface: `rtgsfit()`, status flags, sizes, timers |
| `rtgsfit.c` | The 12 steps of section 4; the only file that includes the preshot |
| `preshot.h` | `extern` declarations of every preshot constant, and `_Static_assert`s on the sizes |
| `grid.h` | `struct rtgsfit_grid` and grid helpers |
| `argument_check.h` | `RTGSFIT_CHECK_SIZE` and `RTGSFIT_CHECK_MIN_SIZE` |
| `current_basis.c` | Section 4.2 |
| `least_squares.c` | Section 4.3: LAPACK `dgelss` without allocation |
| `poisson_solver.c` | Section 4.5 |
| `finite_difference.c` | ∂/∂Z on the grid, and the normal derivative on the edge |
| `magnetic_nulls.c` | Section 4.7 |
| `plasma_boundary.c` | Sections 4.8 and 4.9 |
| `timing.c` | Per-section CPU timers (`TIMING=1`) |

Every module except `rtgsfit.c` takes all its data as arguments, so it can be unit-tested on any
grid (`c/tests`).

### 6.2 Interface

```c
int32_t rtgsfit(const double meas[], int32_t n_meas,
                const double coil_current[], int32_t n_coil_current,
                double psi_norm[], int32_t n_psi_norm,     /* in and out */
                int32_t mask[], int32_t n_mask,            /* in and out */
                double psi[], int32_t n_psi,
                double j_phi[], int32_t n_j_phi,
                double x_vec[], int32_t n_x_vec,
                double constraint_measured[], int32_t n_constraint_measured,
                double constraint_reconstructed[], int32_t n_constraint_reconstructed,
                double xpt_r[], int32_t n_xpt_r,
                double xpt_z[], int32_t n_xpt_z,
                double xpt_psi[], int32_t n_xpt_psi,
                int32_t *xpt_n, int32_t *xpt_diverted,
                double *ip, double *axis_psi, double *axis_r, double *axis_z,
                double *boundary_psi, double *current_centroid_r, double *current_centroid_z,
                double *chi_sq_magnetic, double *chi_sq_total, int32_t *lapack_info);
```

Every array is followed by its length, which is checked against the preshot on every call; a
mismatch returns `RTGSFIT_ERR_ARGUMENT_SIZE` before anything is written. The full description of
each argument is in `c/src/rtgsfit.h`.

### 6.3 Sizes

The sizes are `enum` constants in the generated `preshot_sizes.h`, so they are true compile-time
constants and every array has a fixed size (issue #97; `const int` is not a constant expression in
C). `rtgsfit_get_sizes()` returns them in the order of `enum rtgsfit_size_index`, so a caller can
check that the library was built from the preshot it expects.

### 6.4 Timing

Built with `TIMING=1`, `rtgsfit_timing_get()` returns the thread CPU time of each step, summed since
`rtgsfit_timing_reset()`. A replay stores them per call at `TIMING:DURATION`. Typical medians on an
x86 workstation (one OpenBLAS thread):

| Case | Grid | `rtgsfit()` total | Poisson solve |
|---|---|---|---|
| ST40 | 33 × 65 | 0.76 ms | 0.39 ms |
| FreeGSNKE MAST-U-like | 65 × 129 | 5.5 ms | |

### 6.5 Real-time properties

* No memory is allocated: all intermediate arrays are static. As a result `rtgsfit()` is not
  re-entrant; call it from one thread at a time.
* On x86, `rtgsfit()` sets flush-to-zero and denormals-are-zero in MXCSR, because denormal
  arithmetic is very slow.
* The library links to OpenBLAS (CBLAS and LAPACK). Pin it to one thread
  (`OPENBLAS_NUM_THREADS=1`) on a real-time core. If NumPy in the same process uses MKL, set
  `MKL_THREADING_LAYER=SEQUENTIAL` to avoid conflicts between the two thread pools.

### 6.6 Building

    make -C c                                  # c/build/librtgsfit.so from c/preshot
    make -C c SHOT=52014922 RUN_NAME=RUN21     # write c/preshot from MDSplus first (TE network)
    make -C c test                             # C unit tests; no preshot needed
    make -C c DEBUG=1 TIMING=1 ARCH_FLAGS=...

The code is C11 and compiles without warnings with
`-Wall -Wextra -Wpedantic -Wshadow -Wvla -Wstrict-prototypes`. `-Wvla` guarantees that no
variable-length arrays are used.

**Integrating into the PCS.** Compile `c/src/*.c` and `c/preshot/preshot_data.c` with
`c/preshot` and `c/src` on the include path, link OpenBLAS, and call `rtgsfit()` once per cycle with
`psi_norm` and `mask` kept between calls and initialised from `PRESHOT.INITIAL`. The v2 signature
differs from v1 (section 11), so the Simulink wrapper must be updated.

## 7. The Python package

```
python/rtgsfit/
├── __main__.py           command line interface
├── settings/             XML settings: analytic/, st40/, freegsnke/
├── data_dictionary.py    every stored quantity; generates RTGSFIT.csv and data_dictionary.md
├── preshot.py            Preshot, assemble_preshot(), validation
├── library.py            builds and loads librtgsfit.so (ctypes); Rtgsfit, describe_status()
├── replay.py             read preshot -> read signals -> call rtgsfit() -> collect -> write
├── results.py            Results, and collect_results()
├── boundary.py           post-processed boundary outline and shape (section 4.11)
├── database_readers/     analytic, local, gsfit, freegsnke, rtgsfit_mdsplus, st40_mdsplus
└── database_writers/     local, c_source, rtgsfit_mdsplus
```

### 7.1 Command line

    rtgsfit list-settings
    rtgsfit replay analytic
    rtgsfit replay st40 --set replay/pulse=13343 --set replay/time_end=0.2
    rtgsfit preshot st40
    rtgsfit write-c-preshot --pulse 52014922 --run-name RUN21
    rtgsfit write-csv path/to/standard_utility/csv_files/RTGSFIT.csv
    rtgsfit write-data-dictionary documentation/data_dictionary.md

`rtgsfit` is `python -m rtgsfit`. `--set PATH=VALUE` overrides any setting; the type is that of the
existing setting, or given explicitly as `PATH:TYPE=VALUE` for a new one (`int`, `float`, `bool`,
`str`, `float_array`, `str_array`).

### 7.2 From Python

```python
from rtgsfit.settings import load_settings
from rtgsfit.replay import replay

settings = load_settings("st40").copy_with({"replay/pulse": 13343, "replay/time_end": 0.2})
results = replay(settings)                  # no writers are run
results["GLOBAL.IP"]                         # [N_TIME]
results["PROFILES_2D.R_Z.PSI"]               # [N_TIME, N_Z, N_R]
results.nested["CONSTRAINTS"]["BP_PROBE"]["ALL"]["RECONSTRUCT"]
```

`replay(settings, preshot=..., signals=...)` replaces the readers, for example to replay modified
measurements. `run_replay(settings)` also runs the writers. Each preshot is compiled once into
`$RTGSFIT_CACHE_DIR` (default `~/.cache/rtgsfit`), keyed by the preshot content, the C source and
the build options.

### 7.3 Database readers

A reader provides a preshot (`PreshotReader.read_preshot`), measurements (`SignalsReader.read_signals`),
or both.

| Reader | Preshot | Signals | Needs |
|---|---|---|---|
| `analytic` | yes | yes | nothing: a circular large-aspect-ratio plasma with an exact solution |
| `local` | yes | | a preshot `.npz` written by the `local` writer |
| `rtgsfit_mdsplus` | yes | | the `RTGSFIT` MDSplus tree (TE network) |
| `gsfit` | yes | yes | GSFit, and the data its settings read |
| `st40_mdsplus` | | yes | ST40 `MAG` and `PSU` trees (TE network) |
| `freegsnke` | yes | yes | FreeGSNKE and GSFit: synthetic MAST-U-like equilibrium |

`st40_mdsplus` reads exactly the signals the PCS sees, by their PCS names (`PSI_FLOOP_007`,
`I_SOL_PSU`, ...). Signals are linearly interpolated onto the replay times; the PSU currents use
their own time base.

### 7.4 Database writers

| Writer | Preshot | Results | Notes |
|---|---|---|---|
| `local` | `.npz` | `.npz` | keys are data dictionary keys, e.g. `GLOBAL.IP`; no pickling |
| `c_source` | `preshot_sizes.h`, `preshot_data.c` | | only rewrites files whose content changed |
| `rtgsfit_mdsplus` | `PRESHOT` | everything | needs `writers/rtgsfit_mdsplus/enabled = true` |

A results file or tree contains the whole preshot, the settings XML, the code version and the raw
`meas` and `coil_signal`, so every replay can be reproduced from its own output.

## 8. Settings

Settings are XML files in the style of the IMAS `code.parameters`: nested elements, each leaf with a
`type` (`int`, `float`, `bool`, `str`, `int_array`, `float_array`, `str_array`) and optional `unit`
and `description`. `load_settings(name)` loads `python/rtgsfit/settings/<name>/rtgsfit.xml`, or any
path. Settings are stored with every preshot (`PRESHOT.INFO:SETTINGS`) and every replay
(`INPUT:SETTINGS`).

| Section | Contents |
|---|---|
| `preshot` | `reader`, optional `builder`, and the boundary settings that go into `PRESHOT.SETTINGS` |
| `signals` | `reader` and its options |
| `replay` | `pulse`, `time_start`, `time_end`, `time_step`, `repeat`, `stop_on_status`, `stop_below_ip` |
| `build` | `timing`, `arch_flags` |
| `writers` | `names`, and one section per writer |
| `analytic`, `gsfit`, `freegsnke` | Reader-specific options: grid, sensors, coils, initial guess |

## 9. Storage: MDSplus and local files

Results use GSFit's IMAS-style layout (`GLOBAL:IP`, `GLOBAL.MAG_AXIS:R`, `PROFILES_2D.R_Z:PSI`,
`CONSTRAINTS.FLUX_LOOP.ALL:MEASURED`, ...), so scripts written for GSFit also read RT-GSFit. The
`RTGSFIT` tree is defined by `standard_utility/csv_files/RTGSFIT.csv`, which is generated by
`python -m rtgsfit write-csv` and must be kept identical to it. A local `.npz` has exactly the same
nodes (keys such as `PROFILES_2D.R_Z.PSI`). Every node is listed in
[data_dictionary.md](data_dictionary.md).

Node names are at most 12 characters (an MDSplus limit). The tests check that the CSV, the C header
and the data dictionary agree.

## 10. Verification

The tests are replays: the same code path as `rtgsfit replay`.

| Test | Where | What it checks |
|---|---|---|
| C unit tests | `c/tests`, `tests/test_c_core.py` | Least squares vs LAPACKE, Poisson vs a dense solve, finite differences, nulls, boundary, flood fill, argument checks |
| Analytic replay | `tests/test_analytic_replay.py` | `psi` within 1 %, `b_theta` and `j_phi` vs the exact solution; Ip, axis, circular boundary; monotonic convergence; zero current gives `NO_AXIS`; bit-for-bit agreement with v1 |
| Data dictionary, settings, I/O, CLI | `tests/test_*.py` | Names, CSV, generated docs, round trips, preshot validation, enum sizes |
| FreeGSNKE | `tests/test_freegsnke.py` | Ip 599.0 vs 600 kA; `boundary_psi` 0.18165 vs 0.18178 Wb; axis within one cell; diverted |
| ST40 | `tests/test_st40.py` (TE network) | RUN21 equals a fresh GSFit build; a pulse replay; and vs GSFit for pulses 13343, 13345, 13346 at 30 and 100 ms |

ST40 vs GSFit (13343 at 30 ms): `psi` within 2e-3 Wb, flux loops 5.7e-4 Wb, Bp probes 4.3e-4 T,
Rogowski coils 3.3 A, IVC eigenmode currents 1.9 A. The vs-GSFit test uses GSFit's coil currents,
because the PCS PSU signals differ from them by up to 1 %, which alone moves the Bp probes by more
than 1e-3 T at 100 ms.

The analytic solution and the large-aspect-ratio mutual inductances used by the analytic reader are
derived in [analytic_solution.pdf](analytic_solution/analytic_solution.pdf) and
[mutual_inductance.pdf](mutual_inductance/mutual_inductance.pdf) (LaTeX sources alongside).

The analytic regression against v1 needs one detail: v1 wrote the analytic sensor positions to text
files with 6 decimals and read them back, so its measurements differ by about 1e-7. The fit is
ill-conditioned along `P_PRIME` and `DELTA_Z`, which amplifies this to 0.03 % in `x_vec`. The test
rebuilds v1's rounded measurements, and then agrees with v1 to 1e-4 on every output.

## 11. Migrating from v1

### 11.1 Repository

| v1 | v2 |
|---|---|
| `src/` | `c/src/` |
| `src/constants.c` (generated) | `c/preshot/preshot_sizes.h` + `preshot_data.c` (generated) |
| `utility/mdsplus_const_to_file.py` | `python -m rtgsfit write-c-preshot` |
| `tests/rtgsfit_verify_analytic` | `tests/test_analytic_replay.py` + `database_readers/analytic` |
| `tests/rtgsfit_vs_gsfit` | `tests/test_st40.py` + `database_readers/gsfit` |
| `replay_rtgsfit_st40` (GitLab) | `rtgsfit replay st40`; private analysis in `rtgsfit_workbench` |
| PRESHOT node written by `replay_rtgsfit_st40` | `rtgsfit preshot st40` |

### 11.2 Preshot constants

| v1 C constant / MDSplus node | v2 node (C name `GROUP_NODE`) |
|---|---|
| `N_FLUX_LOOPS` / `N_F_LOOPS` | `SIZES:N_FLUX_LOOP` |
| `N_BP_PROBES` | `SIZES:N_BP_PROBE` |
| `N_ROGOWSKI_COILS` / `N_ROG_COILS` | `SIZES:N_ROGOWSKI` |
| `N_LTRB` | `SIZES:N_EDGE` |
| `N_COEF` | `SIZES:N_X_VEC` |
| `N_PLS` | `SIZES:N_PLASMA` |
| `N_VESS` | `SIZES:N_PASSIVE` |
| `N_SENS_PCS` | `SIZES:N_MEAS` |
| v1 `N_MEAS` (constraints + regularisation) | `SIZES:N_B_VEC` |
| `N_LIMIT` | `SIZES:N_LIMITER` |
| `N_LCFS_MAX` | removed (the outline is post-processed) |
| `FRAC` | `SETTINGS:XPT_FRACTION` |
| `THRESH` | `SETTINGS:TOLERANCE` |
| `DR`, `DZ`, `R_VEC`, `Z_VEC` | `GRID:D_R`, `GRID:D_Z`, `GRID:R`, `GRID:Z` |
| `R_GRID`, `Z_GRID` | `GRID:R_FLAT`, `GRID:Z_FLAT` |
| `INV_R_MU0` | `GRID:INV_MU0_R` |
| `MASK_LIM` | `GRID:MASK_LIMITER` |
| `G_GRID_COIL`, `G_GRID_VESSEL` | `GRID:G_COIL`, `GRID:G_PASSIVE` |
| `INV_R_LTRB_MU0` / `INV_R_L_MU0` | `EDGE:INV_MU0_R` |
| `G_LTRB` | `EDGE:G_EDGE` |
| `LIMIT_R`, `LIMIT_Z` | `LIMITER:R`, `LIMITER:Z` |
| `LIMIT_IDX`, `LIMIT_WEIGHT` / `LIMIT_W` | `LIMITER:GRID_INDEX`, `LIMITER:GRID_WEIGHT` |
| `R_MU0_DZ2` | `POISSON:RHS_SCALE` |
| `LOWER_BAND`, `UPPER_BAND`, `PERM_IDX` | `POISSON:LOWER_BAND`, `POISSON:UPPER_BAND`, `POISSON:PERMUTATION` |
| `SENS_NAMES` | `CONSTRAINT:NAMES` |
| `SENSOR_REPLACEMENT_MATRIX` / `SENS_REP_MAT` | `CONSTRAINT:REPLACEMENT` |
| `G_MEAS_COIL` / `GREENS.MEAS_COIL` | `CONSTRAINT:G_COIL` |
| `WEIGHT` | `FIT:WEIGHT` |
| `G_GRID_MEAS_WEIGHT` / `GREENS.GRID_MEAS_W` | `FIT:A_GRID` |
| `G_COEF_MEAS_WEIGHT` / `GREENS.COEF_MEAS_W` | `FIT:A_PASSIVE` (passive rows only) |
| `COEF_NAMES` | `FIT:X_VEC_NAMES` |
| `MEAS_NAMES` | `MEAS:NAMES` |
| `COIL_NAMES`, `COIL_SIGNALS`, `COIL_MATRIX` | `COIL:NAMES`, `COIL:SIGNALS`, `COIL:MATRIX` |
| `INITIAL_COND.FLUX_NORM`, `INITIAL_COND.MASK` | `INITIAL:PSI_NORM`, `INITIAL:MASK` |

### 11.3 `rtgsfit()` arguments

| v1 | v2 |
|---|---|
| `meas_pcs` | `meas` |
| `coil_curr` | `coil_current` |
| `flux_norm`, `flux_total` | `psi_norm`, `psi` |
| `coef` | `x_vec` |
| `meas_model`, `n_meas_model` | `constraint_reconstructed` (`b_pred` is internal) |
| | `j_phi`, `constraint_measured`, `chi_sq_total` (new) |
| `chi_sq_err` | `chi_sq_magnetic` |
| `flux_boundary`, `plasma_current` | `boundary_psi`, `ip` |
| `r_mag_axis`, `z_mag_axis`, `mag_axis_flux` | `axis_r`, `axis_z`, `axis_psi` |
| `r_cur_centroid`, `z_cur_centroid` | `current_centroid_r`, `current_centroid_z` |
| `xpt_flux`, `xpt_arrays_size` | `xpt_psi`, and one length per array |
| `lcfs_r`, `lcfs_z`, `lcfs_n` | removed: v1 already returned zeros |
| `lcfs_err_code` (output) | the return value |
| `lapack_dgelss_info` | `lapack_info` |

### 11.4 Behaviour changes

* Every array argument has a length, checked on every call; sizes are `enum` constants.
* The status is returned rather than written to `lcfs_err_code`, and includes
  `ARGUMENT_SIZE` and `LAPACK`. On failure, `psi_norm` and `mask` are unchanged (as in v1);
  `boundary_psi` is still reported if it was found.
* The PSU2COIL correction of the replay is dropped: replays use the PSU signals the PCS sees.
* The boundary outline and shape are computed in Python after the replay.
* The analytic test uses exact sensor positions (v1 rounded them to 6 decimals).

## 12. Future work

* **Faster Poisson solver** (PR #105): a sine transform in Z and a tridiagonal solve per mode in R,
  storing three small stencil arrays instead of the banded LU factors. It would roughly halve the
  ST40 call time, make the preshot `POISSON` group 100 times smaller, and remove the dense LU
  factorisation that dominates preshot time on large grids (15 s for 65 × 129). It changes the
  `PRESHOT.POISSON` nodes, so it needs a new preshot run.
* More plasma degrees of freedom (higher-order p' and FF'), and other constraints (MSE, pressure).
* Gaps and the boundary outline in C, if the PCS needs them.

## 13. References

* J.-M. Moret et al., "Tokamak equilibrium reconstruction code LIUQE and its real time
  implementation", Fusion Eng. Des. 91 (2015) 1. doi:10.1016/j.fusengdes.2014.09.019
* O. Sauter and S. Yu. Medvedev, "Tokamak coordinate conventions: COCOS", Comput. Phys. Commun. 184
  (2013) 293. doi:10.1016/j.cpc.2012.09.010
* GSFit: <https://github.com/tokamak-energy/gsfit>
* FreeGSNKE: <https://github.com/FusionComputingLab/freegsnke>
* IMAS Data Dictionary: <https://imas-data-dictionary.readthedocs.io/>
