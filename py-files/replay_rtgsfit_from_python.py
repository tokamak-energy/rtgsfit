import ctypes
import numpy as np
from st40_database import GetData
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
import os

username = os.getlogin()
if username == "peter.buxton":
    pulseNo = 12_050
    pulseNo_write = pulseNo + 11_000_000
    run_name = "RUN09"
    run_description = "floop / (2.0 * pi)"
elif username == "alex.prokopyszyn":
    pulseNo = 12_050
    pulseNo_write = pulseNo + 52_000_000
    run_name = "RUN05"
    run_description = "test run"

# Load the shared library
rtgsfit_lib = ctypes.CDLL(os.path.dirname(os.path.abspath(__file__)) + '/../lib/librtgsfit.so')

# Define the argument types for the rtgsfit function
rtgsfit_lib.rtgsfit.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # meas
    ctypes.POINTER(ctypes.c_double),  # coil_curr
    ctypes.POINTER(ctypes.c_double),  # flux_norm
    ctypes.POINTER(ctypes.c_int),     # mask
    ctypes.POINTER(ctypes.c_double),  # flux_total
    ctypes.POINTER(ctypes.c_double),  # error
    ctypes.POINTER(ctypes.c_double),  # lcfs_r
    ctypes.POINTER(ctypes.c_double),  # lcfs_z
    ctypes.POINTER(ctypes.c_int),     # lcfs_n
    ctypes.POINTER(ctypes.c_double),  # coef
    ctypes.POINTER(ctypes.c_double),  # flux_boundary
    ctypes.POINTER(ctypes.c_double)   # plasma_current
]
# Define the return type for the rtgsfit function
rtgsfit_lib.rtgsfit.restype = ctypes.c_int

mag = GetData(pulseNo, "MAG#BEST")

times_to_reconstruct = np.arange(16e-3, 150e-3, 1e-3)
# times_to_reconstruct = np.array([80.0e-3])
n_time = len(times_to_reconstruct)

psi_list = [
    "FLOOP.L001.PSI",
    "FLOOP.L002.PSI",
    "FLOOP.L003.PSI",
    "FLOOP.L004.PSI",
    "FLOOP.L005.PSI",
    "FLOOP.L006.PSI",
    "FLOOP.L007.PSI",
    "FLOOP.L008.PSI",
    "FLOOP.L009.PSI",
    "FLOOP.L010.PSI",
    "FLOOP.L011.PSI",
    "FLOOP.L012.PSI",
    "FLOOP.L013.PSI",
    "FLOOP.L014.PSI",
    "FLOOP.L015.PSI",
    "FLOOP.L016.PSI",
    "FLOOP.L017.PSI",
    "FLOOP.L018.PSI",
    "FLOOP.L019.PSI",
    "FLOOP.L020.PSI",
    "FLOOP.L021.PSI",
    "FLOOP.L023.PSI",
    "FLOOP.L024.PSI",
    "FLOOP.L025.PSI",
    "FLOOP.L026.PSI",
    "FLOOP.L027.PSI",
    "FLOOP.L029.PSI",
    # "FLOOP.L101.PSI",
    # "FLOOP.L106.PSI",
    "FLOOP.L201.PSI",
    "FLOOP.L212.PSI",
]
n_psi = len(psi_list)
psi_data = np.zeros((n_time, n_psi))
time_experimental = mag.get("TIME")
for i_psi in range(0, n_psi):
    psi_experimental = mag.get(psi_list[i_psi])
    psi_recon_time = np.interp(times_to_reconstruct, time_experimental, psi_experimental)
    psi_data[:, i_psi] = psi_recon_time


bp_probes_names = [
    "BPPROBE.P101.B",
    "BPPROBE.P102.B",
    "BPPROBE.P104.B",
    "BPPROBE.P105.B",
    "BPPROBE.P109.B",
    "BPPROBE.P112.B",
    "BPPROBE.P113.B",
    "BPPROBE.P114.B",
    "BPPROBE.P115.B",
    "BPPROBE.P116.B",
    "BPPROBE.P118.B",
    "BPPROBE.P119.B",
    "BPPROBE.P120.B",
    "BPPROBE.P121.B",
    "BPPROBE.P122.B",
    "BPPROBE.P124.B",
    "BPPROBE.P125.B",
    "BPPROBE.P126.B",
    "BPPROBE.P127.B",
    "BPPROBE.P128.B",
    "BPPROBE.P130.B",
    "BPPROBE.P131.B",
    "BPPROBE.P132.B",
    "BPPROBE.P134.B",
]
n_bp_probe = len(bp_probes_names)
bp_probe_data = np.zeros((n_time, n_bp_probe))
time_experimental = mag.get("TIME")
for i_bp_probe in range(0, n_bp_probe):
    bp_probe_experimental = mag.get(bp_probes_names[i_bp_probe])
    bp_probe_recon_time = np.interp(times_to_reconstruct, time_experimental, bp_probe_experimental)
    bp_probe_data[:, i_bp_probe] = bp_probe_recon_time


rogowski_coils_names = [
    "ROG.INIVC000.I",
    "ROG.BVLT.I",
    "ROG.BVLB.I",
    "ROG.GASBFLT.I",
    "ROG.GASBFLB.I",
    "ROG.HFSPSRT.I",
    "ROG.HFSPSRB.I",
    "ROG.DIVPSRT.I",
    "ROG.DIVPSRB.I",
]
n_rogowski_coil = len(rogowski_coils_names)
rogowski_coil_data = np.zeros((n_time, n_rogowski_coil))
time_experimental = mag.get("TIME")
for i_rogowski_coil in range(0, n_rogowski_coil):
    rogowski_coil_experimental = mag.get(rogowski_coils_names[i_rogowski_coil])
    rogowski_coil_recon_time = np.interp(times_to_reconstruct, time_experimental, rogowski_coil_experimental)
    rogowski_coil_data[:, i_rogowski_coil] = rogowski_coil_recon_time

rogowski_coil_data[:, 0] = rogowski_coil_data[:, 0] - 8.0e3


## PF coil currents
psu2coil = GetData(pulseNo, "PSU2COIL#RUN01")

# Order:
# SOL, MC, DIV, BVL, BVUT, BVUB, PSH
i_pf_rtgsfit_order = np.vstack((
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.SOL.I")),   # SOL
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.MC.I")),    # MC
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.DIV.I")),   # DIV
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.BVL.I")),   # BVL
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.BVUT.I")),  # BVUT
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.BVUB.I")),  # BVUB
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.PSH.I"))    # PSH
)).T


# import matplotlib.pyplot as plt
# plt.plot(i_pf_rtgsfit_order)
# plt.savefig("i_pf_rtgsfit_order.png")


# Load initialisation data
flux_norm = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/flux_norm.txt', dtype=np.float64)
mask = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/mask.txt', dtype=np.float64)
psi_total = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/psi_total.txt', dtype=np.float64)


psi_output = np.zeros((n_time, 65, 33))
psi_b_output = np.zeros((n_time))
plasma_current_output = np.zeros((n_time))
n_psi_n = 100
psi_n = np.linspace(0.0, 1.0, n_psi_n)
p_prime_coeff = np.zeros((n_time))
ff_prime_coeff = np.zeros((n_time))
psin_dz_coeff = np.zeros((n_time))
gs_error = np.zeros((n_time))
p_prime_output = np.zeros((n_time, n_psi_n))
ff_prime_output = np.zeros((n_time, n_psi_n))

for i_time in range(0, n_time):
    # for i_iter in range(0, 30):
    meas = np.concatenate((psi_data[i_time, :], bp_probe_data[i_time, :], rogowski_coil_data[i_time, :], np.array([0.0, 0.0])))
    coil_curr = i_pf_rtgsfit_order[i_time, :] * 0.0 + 1.0
    coil_curr[0] = i_pf_rtgsfit_order[i_time, 0]
    coil_curr[1] = i_pf_rtgsfit_order[i_time, 1]
    coil_curr[2] = i_pf_rtgsfit_order[i_time, 2]
    coil_curr[3] = i_pf_rtgsfit_order[i_time, 3]
    coil_curr[4] = i_pf_rtgsfit_order[i_time, 4]
    coil_curr[5] = i_pf_rtgsfit_order[i_time, 5]
    coil_curr[6] = i_pf_rtgsfit_order[i_time, 6]

    # flux_norm = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    # mask = np.array([1, 0, 1], dtype=np.int32)
    # flux_total = np.array([1.0], dtype=np.float64)
    error = np.array([0.1], dtype=np.float64)
    lcfs_r = np.zeros(500, dtype=np.float64)
    lcfs_z = np.zeros(500, dtype=np.float64)
    lcfs_n = np.zeros(500, dtype=np.int32)
    coef = np.zeros(15, dtype=np.float64)
    flux_boundary = np.array([0.0], dtype=np.float64)
    plasma_current = np.array([0.0], dtype=np.float64)

    # Call the rtgsfit function
    result = rtgsfit_lib.rtgsfit(
        meas.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coil_curr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        flux_norm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        psi_total.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        error.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        lcfs_r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        lcfs_z.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        lcfs_n.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        flux_boundary.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        plasma_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )

    print(f"Result={result}")

    psi_output[i_time, :, :] = np.reshape(psi_total, (65, 33))
    psi_b_output[i_time] = flux_boundary[0]
    plasma_current_output[i_time] = plasma_current[0]

    p_prime_coeff[i_time] = coef[0]
    ff_prime_coeff[i_time] = coef[1]
    psin_dz_coeff[i_time] = coef[2]

    gs_error[i_time] = error[0]

    p_prime_boundary_condition = -coef[0]
    p_prime_output[i_time, :] = coef[0] * (1.0 - psi_n) / (2.0 * np.pi)# + p_prime_boundary_condition

    ff_prime_boundary_condition = -coef[1]
    ff_prime_output[i_time, :] = coef[1] * (1.0 - psi_n) / (2.0 * np.pi)# + ff_prime_boundary_condition

    # coef[0] = p_prime
    # coef[1] = ff_prime
    # coef[2-13] = ivc & coil cases dof
    # coef[14]


# Create a nested dictionary to store data
results = NestedDict()

## Store results
results["TWO_D"]["PSI"] = psi_output
# results["TWO_D"]["MASK"] = mask_output
results["TWO_D"]["RGRID"] = np.linspace(110.000e-3, 1.070e0, 33)
results["TWO_D"]["ZGRID"] = np.linspace(-960.000e-3, 960.000e-3, 65)

results["PROFILES"]["RHO"]["PSI_N"] = psi_n
results["PROFILES"]["RHO"]["P_PRIME"] = p_prime_output
results["PROFILES"]["RHO"]["FF_PRIME"] = ff_prime_output
results["PROFILES"]["RHO"]["P_PRIME_COEFF"] = p_prime_coeff
results["PROFILES"]["RHO"]["FF_PRIME_COEFF"] = ff_prime_coeff
results["PROFILES"]["RHO"]["PSIN_DZ_COEFF"] = psin_dz_coeff


results["TIME"] = times_to_reconstruct
results["GLOBAL"]["IP"] = plasma_current_output
results["GLOBAL"]["PSI_B"] = psi_b_output
results["GLOBAL"]["GS_ERROR"] = gs_error


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
