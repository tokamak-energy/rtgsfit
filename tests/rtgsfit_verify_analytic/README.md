## Description
rtgsfit_verify_analytic sets up and runs analytic verification tests for RTGSFIT. Input data based on a known analytic solution to the Grad-Shafranov equation for large aspect-ratio Tokamaks. Runs RTGSFIT and checks that the numerical solution converges to the analytic as the number of iterations increases.

## Installation

```
uv venv --python 3.13
source .venv/bin/activate
uv pip install -e .
```

## Running Modules

To run a module such as `generate_constants_c` as a script, use Python’s `-m` flag:

```bash
python -m rtgsfit_verify_analytic.generate_constants_c
```

This works from **any directory** as long as the package has been installed in editable mode (e.g. via `uv pip install -e .`).

## Generating plots

First you need to make sure you have created the `output_dict.npy`. The easiest way is to simply run `pytest`, this has the added benefit of checking the code has been installed ok.

You can make plots of the current, magnetic field and $\psi$ values using
```
python -m rtgsfit_verify_analytic.plot.current
```
and replacing `current` with `b_field` and `psi` respectively.

## Maths

The derivation of the analytic solution as well as the formulas we used for the mutual induction and self inductance are the in `latex` directory.