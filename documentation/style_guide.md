# RT-GSFit style guide

The rules below are followed everywhere in `c/` and `python/`; new code must follow them too.
Formatting is automatic: run the formatters before committing, and CI checks it.
Install `.[dev]` for the pinned clang-format 21.1.8 used by CI.
The VS Code workspace enables format-on-save for C and Python; install the configured clang-format
and Ruff extensions. C formatting uses the executable installed into `.venv`.

    .venv/bin/clang-format -i c/src/*.[ch] c/tests/*.[ch]
    ruff format python tests && ruff check python tests

## 1. Names (C and Python)

Names follow the [manual](manual.md) and the [data dictionary](data_dictionary.md). If a quantity is
stored, its variable has the name used there. The rules are those of GSFit
(`gsfit/documentation/gsfit_coding_standards.md`).

**Hierarchical.** The name starts with the object and ends with the property, so related variables
sort together and the shared prefix says they belong together:

| Correct | Wrong |
|---|---|
| `axis_r`, `axis_z`, `axis_psi` | `r_axis`, `r_mag`, `psi_a` |
| `xpt_r`, `xpt_z`, `xpt_psi`, `xpt_n` | `r_xpt`, `xpt_flux`, `n_xpt` (for the count found) |
| `boundary_psi`, `limiter_psi` | `psi_b`, `lcfs_flux` |
| `current_centroid_r` | `r_cur_centroid` |

**Sizes and loops.** A size is `n_<thing>` (C compile-time sizes: `N_<THING>`). A loop over it uses
`i_<thing>`, then `j_<thing>`, `k_<thing>` for nested loops over the same thing:

```c
for (int32_t i_z = 0; i_z < N_Z; i_z++)
{
    for (int32_t i_r = 0; i_r < N_R; i_r++)
    {
        const int32_t i_grid = i_z * N_R + i_r;
        /* ... */
    }
}
```

Never loop to a literal or to a size of something else. A variable whose name is a quantity has the
length of that quantity: `meas` has `N_MEAS` entries, `fit_target` has `N_FIT_TARGET`. This is why v1's `meas`
(which had `N_FIT_TARGET` entries) became `fit_target`.

**Least squares.** `fit_matrix fit_dof ≈ fit_target`, with prediction `fit_prediction`.
Only `fit_dof` contains degrees of freedom; the target contains constraints, not a
magnetic field or unknown coefficients. These descriptive names supersede issue #68's `a/x/b`
notation. X-points are always `xpt`.

**Flux.** `psi`, never `flux`: `psi`, `psi_norm`, `psi_plasma`, `axis_psi`.

**Greens functions.** `g_<to>_<from>` locally, `<GROUP>:G_<FROM>` when stored. For example
`GRID:G_COIL` is the response of `psi` on the grid to unit current in each coil.

**Units** are SI and are not part of the name; they are in the data dictionary and in comments.

**Spelling.** British English in prose (`normalise`, `neighbour`), as in the existing code.

## 2. C

**Standard.** C11, compiled with
`-Wall -Wextra -Wpedantic -Wshadow -Wvla -Wstrict-prototypes -Warray-bounds` and no warnings.

**Formatting** is `.clang-format`: Allman braces (opening braces on their own lines), 4-space
indent, 100 columns, braces on every `if`/`for`/`while`, one parameter per line when a signature
does not fit on one line. This follows the ST40 PCS style guide. Hierarchical names follow GSFit.
Python retains its existing Ruff 120-column limit; PCS's C whitespace rules do not apply to Python.

**Arrays and sizes.**

* Every array argument is followed by its length, and the function checks it before using the array:

  ```c
  int32_t magnetic_nulls_find(const struct rtgsfit_grid *grid,
                              const double psi[],
                              const int32_t n_psi,
                              ...)
  {
      RTGSFIT_CHECK_SIZE(n_psi, grid->n_grid);
  ```

  `RTGSFIT_CHECK_SIZE` and `RTGSFIT_CHECK_MIN_SIZE` (`argument_check.h`) return
  `RTGSFIT_ERR_ARGUMENT_SIZE`. Array parameters are written `double x[]`, not `double *x`;
  `double *x` is for a pointer to a single output value.
  This applies to private helpers too. Use `RTGSFIT_CHECK_ARRAY` to reject nonempty null buffers
  and negative lengths. Buffer lengths inside operator structures must also be checked.
* Compile-time sizes are `enum` constants (issue #97), never `const int` or `#define`:
  `const int` is not a constant expression in C, so arrays sized by it would be variable-length.
* No variable-length arrays (`-Wvla`) and no `malloc` in the real-time path. Workspaces are
  `static` and sized by `enum` constants.
* Check size relations at compile time with `_Static_assert` where possible (`preshot.h`).

**Functions.**

* Public functions are named `<module>_<verb>[_<object>]`: `magnetic_nulls_find`,
  `plasma_boundary_flood_fill`, `poisson_solve_free_boundary`. Functions private to a file are
  `static`, and need no prefix.
* A function that can fail returns `int32_t` status (`enum rtgsfit_status`); results go through
  pointer arguments. Outputs come after inputs.
* Inputs are `const`. Local variables that do not change are `const`. Declare variables where they
  are first used, and loop indices in the `for`.
  A value parameter such as `const int32_t n_meas` cannot be reassigned inside the function.
  This top-level `const` does not change the caller's ABI; `const double meas[]` additionally
  prevents modification of the pointed-to values through that parameter.
* Integers that are stored or cross the API are `int32_t` (MDSplus and Simulink use 32 bits);
  `double` for every real number.

**Modules.** One header per `.c` file, with an include guard `RTGSFIT_<MODULE>_H`. Only `rtgsfit.c`
includes the preshot (`preshot.h`); every other module takes its data as arguments, so it can be
unit-tested on any grid. A new module gets a test in `c/tests`.

**Preshot constants** are named `<GROUP>_<NODE>`, the C form of `PRESHOT.<GROUP>:<NODE>`, e.g.
`FIT_A_GRID` for `PRESHOT.FIT:A_GRID`. They are generated; never edit `c/preshot/*`.

**Comments.** Every public function has a comment in its header that says what it does, its inputs,
outputs and return value, with units. In the `.c` file, comment why, not what. Use `/* */`.
Reference equations by paper and number (e.g. Moret et al. 2015, eq. 61).

## 3. Python

**Version and tools.** Python ≥ 3.11. `ruff format` and `ruff check` with the configuration in
`pyproject.toml` (120 columns; imports one per line, sorted).

**Types.** Every function is annotated. Arrays are `numpy.typing.NDArray[np.float64]` (aliases such as
`FloatArray`); use `from __future__ import annotations`.

**Structure.**

* Follow the GSFit layout: `settings/<name>/rtgsfit.json`, `database_readers/<name>/__init__.py`,
  `database_writers/<name>/__init__.py`. A reader or writer is found by the name of its directory.
* Readers and writers that need Tokamak Energy or optional packages (`mdsthin`,
  `gsfit`, `freegsnke`) import them inside the function that uses them, so that
  `import rtgsfit` works anywhere.
* A stored quantity is defined once, in `data_dictionary.py`. Do not write node names as strings
  elsewhere unless they are read from the data dictionary or are keys of `Results`.
* Settings are ordinary JSON, validated against `settings/settings.schema.json`; put units,
  descriptions and allowed values in that shared schema, not in repeated per-setting wrappers.
  Do not hard-code a machine's numbers in code. Archive effective settings through `Settings.json_text`.
* Errors are exceptions with a message that says what to do; do not `print` and continue.

**Docstrings.** A module docstring that says what the module is for. A one-line docstring for each
public function, longer when the behaviour is not obvious. Use backticks for code names.

## 4. Tests

* A test is a replay, or a unit test of one C module or Python function. Tests that need the TE
  network are marked `te_network`; tests that need GSFit and FreeGSNKE are marked `freegsnke`. Both
  are skipped by default.
* Compare to an independent reference (an analytic solution, GSFit, FreeGSNKE, LAPACK), and state
  where each tolerance comes from.
* Every bug fix gets a test that would have failed before.

## 5. Documentation

* Change the manual when you change behaviour, and the data dictionary when you change a stored
  quantity. Then regenerate the generated files, which the tests check:

      python -m rtgsfit write-data-dictionary documentation/data_dictionary.md
      python -m rtgsfit write-csv <standard_utility>/standard_utility/csv_files/RTGSFIT.csv

* The CSV is the definition of the MDSplus tree: commit it to `standard_utility` together with the
  change here.
