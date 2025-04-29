import ctypes
import numpy as np
from st40_database import GetData
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
import matplotlib.pyplot as plt
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
    run_name = "GJ_EIG12345"
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

# time_to_reconstruct = 0.126
time_to_reconstruct = 16e-3

psi_list = ["PSI_FLOOP_001 ","PSI_FLOOP_002 ","PSI_FLOOP_004 ","PSI_FLOOP_005 ","PSI_FLOOP_006 ","PSI_FLOOP_007 ","PSI_FLOOP_008 ","PSI_FLOOP_009 ","PSI_FLOOP_010 ","PSI_FLOOP_011 ","PSI_FLOOP_012 ","PSI_FLOOP_013 ","PSI_FLOOP_014 ","PSI_FLOOP_015 ","PSI_FLOOP_016 ","PSI_FLOOP_017 ","PSI_FLOOP_018 ","PSI_FLOOP_019 ","PSI_FLOOP_020 ","PSI_FLOOP_021 ","PSI_FLOOP_023 ","PSI_FLOOP_024 ","PSI_FLOOP_025 ","PSI_FLOOP_026 ","PSI_FLOOP_027 ","PSI_FLOOP_029 ","PSI_FLOOP_101 ","PSI_FLOOP_106"]
def convert_psi_list_format(inpsi_list):
    outpsi_list = []
    for item in inpsi_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"FLOOP.L{number_part}.PSI"
        outpsi_list.append(new_format)
    return outpsi_list
psi_list = convert_psi_list_format(psi_list)
n_psi = len(psi_list)
psi_data = np.zeros(n_psi)
time_experimental = mag.get("TIME")
for i_psi in range(0, n_psi):
    psi_experimental = mag.get(psi_list[i_psi])
    psi_recon_time = np.interp(time_to_reconstruct, time_experimental, psi_experimental)
    psi_data[i_psi] = psi_recon_time
print("Psi data")
print(psi_data)


bp_probes_names = ["B_BPPROBE_101 ","B_BPPROBE_102 ","B_BPPROBE_103 ","B_BPPROBE_104 ","B_BPPROBE_105 ","B_BPPROBE_109 ","B_BPPROBE_112 ","B_BPPROBE_113 ","B_BPPROBE_114 ","B_BPPROBE_115 ","B_BPPROBE_118 ","B_BPPROBE_119 ","B_BPPROBE_120 ","B_BPPROBE_121 ","B_BPPROBE_122 ","B_BPPROBE_124 ","B_BPPROBE_125 ","B_BPPROBE_126 ","B_BPPROBE_127 ","B_BPPROBE_128 ","B_BPPROBE_131 ","B_BPPROBE_134 "]
def convert_bp_list_format(inbp_list):
    outbp_list = []
    for item in inbp_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"BPPROBE.P{number_part}.B"
        outbp_list.append(new_format)
    return outbp_list
bp_probes_names = convert_bp_list_format(bp_probes_names)
n_bp_probe = len(bp_probes_names)
bp_probe_data = np.zeros(n_bp_probe)
time_experimental = mag.get("TIME")
for i_bp_probe in range(0, n_bp_probe):
    bp_probe_experimental = mag.get(bp_probes_names[i_bp_probe])
    bp_probe_recon_time = np.interp(time_to_reconstruct, time_experimental, bp_probe_experimental)
    bp_probe_data[i_bp_probe] = bp_probe_recon_time
print("BP Probe data")
print(bp_probe_data)


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
rogowski_coil_data = np.zeros(n_rogowski_coil)
time_experimental = mag.get("TIME")
for i_rogowski_coil in range(0, n_rogowski_coil):
    rogowski_coil_experimental = mag.get(rogowski_coils_names[i_rogowski_coil])
    rogowski_coil_recon_time = np.interp(time_to_reconstruct, time_experimental, rogowski_coil_experimental)
    rogowski_coil_data[i_rogowski_coil] = rogowski_coil_recon_time
rogowski_coil_data[0] = rogowski_coil_data[0] - 8.0e3
print("Rogowski coil data")
print(rogowski_coil_data)

## PF coil currents
psu2coil = GetData(pulseNo, "PSU2COIL#RUN01")

# Order:
# SOL, MC, DIV, BVL, BVUT, BVUB, PSH
i_pf_rtgsfit_order = np.array([
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.SOL.I")),   # SOL
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.MC.I")),    # MC
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.DIV.I")),   # DIV
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.BVL.I")),   # BVL
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.BVUT.I")),  # BVUT
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.BVUB.I")),  # BVUB
    np.interp(time_to_reconstruct, time_experimental, psu2coil.get("PF.PSH.I"))    # PSH
])



# Load initialisation data
flux_norm = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/flux_norm.txt', dtype=np.float64)
mask = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/mask.txt', dtype=np.float64)
psi_total = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/psi_total.txt', dtype=np.float64)

n_iter = 100
psi_output = np.zeros((n_iter, 65, 33))
psi_b_output = np.zeros((n_iter))
plasma_current_output = np.zeros((n_iter))
n_psi_n = 100
psi_n = np.linspace(0.0, 1.0, n_psi_n)
p_prime_coeff = np.zeros((n_iter))
ff_prime_coeff = np.zeros((n_iter))
psin_dz_coeff = np.zeros((n_iter))
gs_error = np.zeros((n_iter))
p_prime_output = np.zeros((n_iter, n_psi_n))
ff_prime_output = np.zeros((n_iter, n_psi_n))

for i_iter in range(0, n_iter):

    meas = np.concatenate((psi_data, bp_probe_data, rogowski_coil_data, np.array([0.0, 0.0])))
    coil_curr = i_pf_rtgsfit_order
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

    psi_output[i_iter, :, :] = np.reshape(psi_total, (65, 33))
    psi_b_output[i_iter] = flux_boundary[0]
    plasma_current_output[i_iter] = plasma_current[0]

    p_prime_coeff[i_iter] = coef[0]
    ff_prime_coeff[i_iter] = coef[1]
    psin_dz_coeff[i_iter] = coef[2]

    gs_error[i_iter] = error[0]

    p_prime_boundary_condition = -coef[0]
    p_prime_output[i_iter, :] = coef[0] * (1.0 - psi_n) / (2.0 * np.pi)# + p_prime_boundary_condition

    ff_prime_boundary_condition = -coef[1]
    ff_prime_output[i_iter, :] = coef[1] * (1.0 - psi_n) / (2.0 * np.pi)# + ff_prime_boundary_condition



import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
# Make a video of psi_output at each iteration from a series of png files
for i_iter in range(0, n_iter):
    psi_i = psi_output[i_iter, :, :]
    fig, ax = plt.subplots()
    cs = ax.contour(psi_i,
                    levels=20)
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('psi')
    ax.set_title(f"Iteration {i_iter}")
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_aspect('equal')
    os.makedirs(f"plots/psi", exist_ok=True)
    fig.savefig(f"plots/psi/psi_{i_iter:03d}.png")
    plt.close(fig)

    
