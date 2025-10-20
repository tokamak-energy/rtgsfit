import matplotlib.pyplot as plt
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
from netCDF4 import Dataset
import numpy as np
import sys


def write_data_to_mdsplus(
    pulseNo: int,
    run_name: str = "RUN01",
    run_description: str = "Standard run with default settings",
    settings_path: str = "default",
    write_to_mds: bool = True,
    pulseNo_write: int | None = None,
)-> None:
    data_file_name = f"/home/pcs.user/ops_double_threaded/st40pcs_dtacq/results/rtgsfit_results_{pulseNo}.nc"
    data_file_name = f"/home/filip.janky/ops/pcs/model/ST40PCS/results/rtgsfit_results_{pulseNo}.nc"

    # Load data from netCDF file
    with Dataset(data_file_name, "r") as nc:
        print("Variables in netCDF file new:")
        for var in nc.variables:
	        print(var)
        time = np.array(nc.variables["time"][:])
        flux_norm = np.array(nc.variables["flux_norm"][:])
        mask = np.array(nc.variables["mask"][:])
        chi_sq_err = np.array(nc.variables["chi_sq_err"][:])
        lcfs_n = np.array(nc.variables["lcfs_n"][:])
        flux_total = np.array(nc.variables["flux_total"][:])
        lcfs_r = np.array(nc.variables["lcfs_r"][:])
        lcfs_z = np.array(nc.variables["lcfs_z"][:])
        coef = np.array(nc.variables["coef"][:])
        flux_boundary = np.array(nc.variables["flux_boundary"][:])
        plasma_current = np.array(nc.variables["plasma_current"][:])
        lcfs_err_code = np.array(nc.variables["lcfs_err_code"][:])
        lapack_dgelss_info = np.array(nc.variables["lapack_dgelss_info"][:])
        meas_model = np.array(nc.variables["meas_model"][:])
        sensors_IN = np.array(nc.variables["sensors_IN"][:])
        I_PF_IN = np.array(nc.variables["I_PF_IN"][:])
        r = np.array(nc.variables["r"][:])
        z = np.array(nc.variables["z"][:])

        # Create a nested dictionary to store data
        results = NestedDict()

        ## Store results
        results["TIME"] = time
        results["TWO_D"]["PSI_N"] = flux_norm
        results["TWO_D"]["MASK"] = mask
        results["GLOBAL"]["CHIT"] = chi_sq_err
        results["P_BOUNDARY"]["NBND"] = lcfs_n
        results["TWO_D"]["PSI"] = flux_total
        results["P_BOUNDARY"]["RBND"] = lcfs_r
        results["P_BOUNDARY"]["ZBND"] = lcfs_z
        #results[""][""] = coef
        results["GLOBAL"]["PSI_B"] = flux_boundary
        results["GLOBAL"]["IP"] = plasma_current
        results["GLOBAL"]["LCFS_ERR"] = lcfs_err_code
        results["GLOBAL"]["DGELSS_INFO"] = lapack_dgelss_info
        #results[""][""] = meas_model
        #results[""][""] = sensors_IN
        results["TWO_D"]["RGRID"] = r
        results["TWO_D"]["ZGRID"] = z

        util.create_script_nodes(
            script_name="RTGSFIT",
            pulseNo_write=pulseNo + 30000000,
            pulseNo_cal=None,
            run_name=run_name,
            run_info=run_description,
            workflows=None,
            link_best=False,
        )
        util.write_script_data(
            script_name="RTGSFIT",
            pulseNo_write=pulseNo + 30000000,
            data_to_write=results.to_dictionary(),
            pulseNo_cal=None,
            run_name=run_name,
            run_description=run_description,
        )


args = sys.argv
pulseNo = int(args[1])
print(args[0])
print(args[1])
print(len(args))
if len(args) > 2:
    run_name = args[2]
    write_data_to_mdsplus(pulseNo, run_name)
else:
    write_data_to_mdsplus(pulseNo)

