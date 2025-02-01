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
    pulseNo = 12_050
    pulseNo_write = pulseNo + 30_000_000
    run_name = "RUN01"
    run_description = "test writing"
    data_file_name = (
        f"/home/filip.janky/ops/rtgsfit/results/rtgsfit_results_{pulseNo}.nc"
    )
elif username == "peter.buxton":
    # Variables to write to MDSplus
    pulseNo = 12_050
    pulseNo_write = pulseNo + 11_000_000
    run_name = "RUN05"
    run_description = "Current version of RTGSFIT"
    data_file_name = f"/home/peter.buxton/0_Version_Controlled/rtgsfit_github/results/rtgsfit_results_{pulseNo}.nc"

# Load data from netCDF file
with Dataset(data_file_name, "r") as nc:
    print("Variables in netCDF file:")
    for var in nc.variables:
        print(var)
    psi = np.array(nc.variables["flux_total"][:])
    psi_b = np.array(nc.variables["flux_boundary"][:])
    mask = np.array(nc.variables["mask"][:])
    time = np.array(nc.variables["time"][:])
    plasma_current = np.array(nc.variables["plasma_current"][:])
    r = np.array(nc.variables["r"][:])
    z = np.array(nc.variables["z"][:])
    flux_loops_measured = np.array(nc.variables["flux_loops_measured"][:])
    bp_probes_measured = np.array(nc.variables["bp_probes_measured"][:])
    rogowski_coils_measured = np.array(nc.variables["rogowski_coils_measured"][:])


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

# NEED TO FIX THIS FLOW. Perhaps read these in from the *.mat binary file?
flux_loops_names = [
    "PSI_FLOOP_001",
    "PSI_FLOOP_002",
    "PSI_FLOOP_003",
    "PSI_FLOOP_004",
    "PSI_FLOOP_005",
    "PSI_FLOOP_006",
    "PSI_FLOOP_007",
    "PSI_FLOOP_008",
    "PSI_FLOOP_009",
    "PSI_FLOOP_010",
    "PSI_FLOOP_011",
    "PSI_FLOOP_012",
    "PSI_FLOOP_013",
    "PSI_FLOOP_014",
    "PSI_FLOOP_015",
    "PSI_FLOOP_016",
    "PSI_FLOOP_017",
    "PSI_FLOOP_018",
    "PSI_FLOOP_019",
    "PSI_FLOOP_020",
    "PSI_FLOOP_021",
    "PSI_FLOOP_023",
    "PSI_FLOOP_024",
    "PSI_FLOOP_025",
    "PSI_FLOOP_026",
    "PSI_FLOOP_027",
    "PSI_FLOOP_029",
    "PSI_FLOOP_201",
    "PSI_FLOOP_212",
]
bp_probes_names = [
    "B_BPPROBE_101",
    "B_BPPROBE_102",
    "B_BPPROBE_104",
    "B_BPPROBE_105",
    "B_BPPROBE_109",
    "B_BPPROBE_112",
    "B_BPPROBE_113",
    "B_BPPROBE_114",
    "B_BPPROBE_115",
    "B_BPPROBE_116",
    "B_BPPROBE_118",
    "B_BPPROBE_119",
    "B_BPPROBE_120",
    "B_BPPROBE_121",
    "B_BPPROBE_122",
    "B_BPPROBE_124",
    "B_BPPROBE_125",
    "B_BPPROBE_126",
    "B_BPPROBE_127",
    "B_BPPROBE_128",
    "B_BPPROBE_130",
    "B_BPPROBE_131",
    "B_BPPROBE_132",
    "B_BPPROBE_134",
]
rogowski_coils_names = [
    "I_ROG_INIVC000",
    "I_ROG_BVLT",
    "I_ROG_BVLB",
    "I_ROG_GASBFLT",
    "I_ROG_GASBFLB",
    "I_ROG_HFSPSRT",
    "I_ROG_HFSPSRB",
    "I_ROG_DIVPSRT",
    "I_ROG_DIVPSRB",
]

results["CONSTRAINTS"]["FLOOP"]["MVALUE"] = flux_loops_measured
results["CONSTRAINTS"]["FLOOP"]["CVALUE"] = 0.0 * flux_loops_measured
results["CONSTRAINTS"]["FLOOP"]["NAME"] = np.array(flux_loops_names)

results["CONSTRAINTS"]["BPPROBE"]["MVALUE"] = bp_probes_measured
results["CONSTRAINTS"]["BPPROBE"]["CVALUE"] = 0.0 * bp_probes_measured
results["CONSTRAINTS"]["BPPROBE"]["NAME"] = np.array(bp_probes_names)

results["CONSTRAINTS"]["ROG"]["MVALUE"] = rogowski_coils_measured
results["CONSTRAINTS"]["ROG"]["CVALUE"] = 0.0 * rogowski_coils_measured
results["CONSTRAINTS"]["ROG"]["NAME"] = np.array(rogowski_coils_names)

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

print(f"Written to pulseNo={pulseNo_write}, run_name={run_name}")
