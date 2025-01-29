import matplotlib.pyplot as plt
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
from netCDF4 import Dataset
import numpy as np

# Variables to write to MDSplus
pulseNo = 12_050
pulseNo_write = pulseNo + 30_000_000
run_name = "RUN01"
run_description = "test writing"

# Load data from netCDF file
data_file_name = f"/home/filip.janky/ops/rtgsfit/tests/rtgsfit_results_{pulseNo}.nc"
with Dataset(data_file_name, "r") as nc:
    print("Variables in netCDF file:")
    for var in nc.variables:
        print(var)
    psi = np.array(nc.variables["psi_total"][:])
    psi_b = np.array(nc.variables["psi_b"][:])
    mask = np.array(nc.variables["mask"][:])
    time = np.array(nc.variables["time"][:])
    r = np.array(nc.variables["r"][:])
    z = np.array(nc.variables["z"][:])


# Create a nested dictionary to store data
results = NestedDict()

## Store results
results["TWO_D"]["PSI"] = psi
results["TWO_D"]["MASK"] = mask
results["TWO_D"]["RGRID"] = r
results["TWO_D"]["ZGRID"] = z

results["TIME"] = time
results["GLOBAL"]["IP"] = 0.0 * time # BUXTON: need to fix this
results["GLOBAL"]["PSI_B"] = psi_b

# Write to MDSplus
util.create_script_nodes(
    script_name="RTGSFIT",
    pulseNo_write=pulseNo_write,
    pulseNo_cal=None,
    run_name=run_name,
    run_info=run_description,
    workflows=None,
    link_best=False,
)
util.write_script_data(
    script_name="RTGSFIT",
    pulseNo_write=pulseNo_write,
    data_to_write=results.to_dictionary(),
    pulseNo_cal=None,
    run_name=run_name,
    run_description=run_description,
)
