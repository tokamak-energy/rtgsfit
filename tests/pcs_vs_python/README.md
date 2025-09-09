# pcs_vs_python

`pcs_vs_python` is a sub-repository within the `tests/` directory of the `rtgsfit` project.  
It provides tests to verify that results we get from running `rtgsfit` using python agree with the
results we get by running through the ST40 Plasma Control System on Simulink.

## Installation

```bash
uv venv --python 3.13
source .venv/bin/activate
uv pip install -e .
uv pip install /home/alex.prokopyszyn/my_mdsplus/python/MDSplus/.
uv pip install "numpy<2"
```

Note that the same virtual environment used in `rtgsfit_vs_gsfit` should work here.

## Repository Structure

- **`data/`**  
  Initially contains only `default_config.yaml`. This directory will also data produced by the scripts in `src/pcs_vs_python/replay/`.

- **`plots/`**  
  Initially empty. This directory will store plots produced by the scripts in `src/pcs_vs_python/plot/`.

- **`src/pcs_vs_python/`**  
  Contains the core Python modules:

  - `config_loader.py`  
    Defines `load_and_prepare_config()`, which loads `default_config.yaml` from the `data/` directory and returns a configuration dictionary (`cfg`) with details such as pulse number, run name, resolution, time snapshots etc.

  - `preshot.py`  
    Contains routines to setup the "PRESHOT" node and compile the code.

- **`src/pcs_vs_python/examples/`**
  Contains example python scripts which run RT-GSFit and possibly GSFit over a range of pulses. They run the full pipeline including the "PRESHOT" node setup, compling the code, running the code and plotting the results.

- **`src/pcs_vs_python/plot/`**  
  Contains modules for plotting data produced by `replay/rtgsfit.py` and `replay/gsfit.py`.
  Also contains the `rtgsfit_pred_meas.py` module as this is needed to calculate the predicted measurments
  from the `coef` array outputted by RT-GSFit.

- **`src/pcs_vs_python/`**
  Contains the moddules to replay RT-GSFit and GSFit and save the results onto MDS+.
  Also contains the `prep_meas_coil_curr.py` to provide the measuremtns and PF coil currents at every time step.
   
- **`pyproject.toml`**  
  Contains project metadata including Python package dependencies.

## Tests

To run the test you first need to run the `preshot.py` code to generate the MDS+ node for pulse number 99,000,230 and also compile the code.
Then run the PCS code through Simulink which will then upload data to MDS+.
The run
```bash
python tests
```
and this will then run the RT-GSFit C code using python routines and then check these agree with
results obtained through Simulink.