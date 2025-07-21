### `rtgsfit-vs-gsfit`

This is a sub-repository located in the `tests/` directory of the `rtgsfit` project.  
It provides tests to verify that the results from `rtgsfit` are consistent with those from `gsfit` for ST40 pulses.

## Installation

```
uv venv --python 3.13
source .venv/bin/activate
uv pip install -e .
uv pip install /home/alex.prokopyszyn/my_mdsplus/python/MDSplus/.
uv pip install "numpy<2"
```