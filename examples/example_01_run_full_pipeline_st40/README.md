# Example 01: the full RT-GSFit pipeline for ST40

[example_01_run_full_pipeline_st40.ipynb](example_01_run_full_pipeline_st40.ipynb) runs every step of
RT-GSFit v2 for ST40 pulse 14926, one function call at a time, and plots what each step produces:

1. build the preshot with GSFit and store it in MDSplus;
2. read it back, write it as C source and compile `librtgsfit.so`;
3. read the measurements of the pulse and replay it, calling `rtgsfit()` once per millisecond as the
   PCS does;
4. plot the results, store them in MDSplus and read them back.

It shows how the parts of RT-GSFit fit together, with links to the relevant sections of the
[manual](../../documentation/manual.md).

## What it writes

The notebook writes two runs to the `RTGSFIT` tree on the MDSplus server `smaug`:

| `RTGSFIT` pulse | Run | Contents | Description (`\RTGSFIT::TOP.TEST03:HELP`) |
|---|---|---|---|
| 52000230 | `TEST03` | The preshot | Testing example_01 for rtgsfit_v2: preshot built by GSFit for ST40 pulse 14926 |
| 52014926 | `TEST03` | The replay of ST40 pulse 14926 | Testing example_01 for rtgsfit_v2: replay of ST40 pulse 14926 with the preshot from RTGSFIT 52000230 TEST03 |

Both destinations are in private development pulse space, not the shared PCS preshot at
99000230 `RUN19`. Use your own allocated private pulse range when adapting this example.
The network C-generator test reads 52000230 `TEST03`; it does not write to MDSplus.

The notebook stops before writing anything if either run already exists. To run it again, change
`RUN_NAME` in section 1, or set `OVERWRITE = True` to write over the existing runs. It never changes
`BEST`: `link_best` is false.

Section 5 makes this preshot the active one: it writes `c/src/preshot_sizes.h` and
`c/src/preshot_data.c` (gitignored) and compiles `c/build/librtgsfit.so`, inside this repository. It
prints the paths and inspects those exact files. Only one preshot is active at a time, so running
another replay or build afterwards replaces them.

## Files

| File | Contents |
|---|---|
| `example_01_run_full_pipeline_st40.ipynb` | The pipeline, with explanations. It is committed without outputs. |
| `preshot_to_rtgsfit.svg` | Data-flow diagram showing how the preshot is compiled and used during replay. |
| `helpers.py` | Checks of the environment and the MDSplus runs, the replay summary and the plots |

## Requirements

* The Tokamak Energy network, and write access to the `RTGSFIT` tree on `smaug`.
* The Python package with the TE extras (`gsfit` with `st40_database`, `mdsthin` and `standard_utility`,
  from TE's GitLab), and Jupyter. From the repository root:

      uv pip install -e ".[te,dev]" jupyter

  uv does not read `pip.conf`: add TE's GitLab package index to `~/.config/uv/uv.toml` as
  `extra-index-url = ["..."]`, or set `UV_EXTRA_INDEX_URL`.
* `standard_utility/csv_files/RTGSFIT.csv` must be identical to the output of `rtgsfit write-csv`.
  The notebook checks this before writing. If the released `standard_utility` is older than the data
  dictionary, install it from a checkout that has the new CSV:
  `uv pip install --no-deps -e path/to/standard_utility`.
* gcc, make and OpenBLAS (CBLAS and LAPACK), to compile `librtgsfit.so`.

## Running it

Open the notebook in Jupyter or VS Code and run all the cells. The working directory must be this
directory, which is the default, so that `import helpers` works. Or, from the command line:

    cd examples/example_01_run_full_pipeline_st40
    mkdir -p ../../c/cache/example_01
    cp preshot_to_rtgsfit.svg ../../c/cache/example_01/
    jupyter nbconvert --to notebook --execute --output-dir ../../c/cache/example_01 example_01_run_full_pipeline_st40.ipynb

This writes the executed copy to `c/cache/example_01/` inside the repository and leaves the notebook
here unchanged. The diagram is copied alongside it. A run takes about a minute and a half.
Building the preshot with GSFit takes about 40 s
and each MDSplus write about 20 s; the replay itself takes under a second.

Clear the outputs before committing a change to the notebook, because they contain ST40 data:

    jupyter nbconvert --clear-output --inplace example_01_run_full_pipeline_st40.ipynb

Section 11 of the notebook gives the equivalent `rtgsfit` command lines, and the commands for the
PCS build.
