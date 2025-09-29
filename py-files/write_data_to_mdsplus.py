import matplotlib.pyplot as plt
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
from netCDF4 import Dataset
import numpy as np
import os

# username = "peter.buxton"
# username = "filip.janky"
username = os.getlogin()

if username == "filip.janky":
    # Variables to write to MDSplus
    pulseNo = 13349
    pulseNo_write = pulseNo + 30_000_000
    run_name = "RUN01"
    run_description = "test writing"
    #data_file_name = f"/home/filip.janky/ops/rtgsfit/results/rtgsfit_results_{pulseNo}.nc"
    data_file_name = f"/home/filip.janky/ops/pcs/model/ST40PCS/results/rtgsfit_results_{pulseNo}.nc"
elif username == "peter.buxton":
    # Variables to write to MDSplus
    pulseNo = 12_050
    pulseNo_write = pulseNo + 11_000_000
    run_name = "RUN01"
    run_description = "test writing"
    data_file_name = f"/home/peter.buxton/0_Version_Controlled/rtgsfit_github/results/rtgsfit_results_{pulseNo}.nc"

# Load data from netCDF file
with Dataset(data_file_name, "r") as nc:
    print("Variables in netCDF file:")
    for var in nc.variables:
        print(var)
    psi = np.array(nc.variables["flux_total"][:])
    psi_b = np.array(nc.variables["flux_boundary"][:])
#    mask = np.array(nc.variables["mask"][:])
    time = np.array(nc.variables["time"][:])
    plasma_current = np.array(nc.variables["plasma_current"][:])
    r = np.array(nc.variables["r"][:])
    z = np.array(nc.variables["z"][:])
    chi_sq_err = np.array(nc.variables["chi_sq_err"][:])


# Create a nested dictionary to store data
results = NestedDict()

## Store results
results["TWO_D"]["PSI"] = psi
results["TWO_D"]["MASK"] = mask
results["TWO_D"]["RGRID"] = r
results["TWO_D"]["ZGRID"] = z
results["TIME"] = time
results["GLOBAL"]["IP"] = plasma_current
results["GLOBAL"]["PSI_B"] = psi_b
results["GLOBAL"]["CHIT"] = chi_sq_err

# Write to MDSplus
print(results)
print(results["TIME"].shape)
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
