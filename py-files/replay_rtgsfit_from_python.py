import ctypes
import numpy as np
from st40_database import GetData
import standard_utility as util  # type: ignore
from diagnostics_analysis_base import NestedDict
import os

username = os.getlogin()
if username == "peter.buxton":
    pulseNo = 13_343
    pulseNo_write = pulseNo + 11_000_000
    run_name = "TEST03"
    run_description = "floop / (2.0 * pi)"
elif username == "alex.prokopyszyn":
    pulseNo = 12_050
    pulseNo_write = pulseNo + 52_000_000
    run_name = "GJ_EIG12345"
    run_description = "test run"
# pulseNo = 13349 # sed_tag:pulse_num_replay_range_of_pulses
# pulseNo_write = 52013349 # sed_tag:pulse_num_write_replay_range_of_pulses
# run_name = 'SCAN3' # sed_tag:run_name_replay_range_of_pulses
# run_description = 'Part of a scan of pulses to test the RT-GSFit code. Starting at a later time, namely, 60e-3 s. Also changed MC current from equal to MCT to average of MCB and MCT. Now using pulse number 99000230 instead of 12001000 in the RT-GSFit Matlab setup code also using elmag_run=RUN18 instead of RUN16. Where the aim is to fix the issue of assuming the MCB and MCT currents are equal when they are actually of opposite sign in the latter part of the pulse.' # sed_tag:run_description_replay_range_of_pulses

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

times_to_reconstruct = np.arange(60e-3, 150e-3, 1e-3)
# times_to_reconstruct = np.array([80.0e-3])
n_time = len(times_to_reconstruct)


# psi_list = ["PSI_FLOOP_001 ","PSI_FLOOP_002 ","PSI_FLOOP_004 ","PSI_FLOOP_005 ","PSI_FLOOP_006 ","PSI_FLOOP_007 ","PSI_FLOOP_008 ","PSI_FLOOP_009 ","PSI_FLOOP_010 ","PSI_FLOOP_011 ","PSI_FLOOP_012 ","PSI_FLOOP_013 ","PSI_FLOOP_014 ","PSI_FLOOP_015 ","PSI_FLOOP_016 ","PSI_FLOOP_017 ","PSI_FLOOP_018 ","PSI_FLOOP_019 ","PSI_FLOOP_020 ","PSI_FLOOP_021 ","PSI_FLOOP_023 ","PSI_FLOOP_024 ","PSI_FLOOP_025 ","PSI_FLOOP_026 ","PSI_FLOOP_027 ","PSI_FLOOP_029 ","PSI_FLOOP_101 ","PSI_FLOOP_106"]
psi_list = ["PSI_FLOOP_001 ","PSI_FLOOP_002 ","PSI_FLOOP_004 ","PSI_FLOOP_005 ","PSI_FLOOP_006 ","PSI_FLOOP_007 ","PSI_FLOOP_008 ","PSI_FLOOP_009 ","PSI_FLOOP_010 ","PSI_FLOOP_011 ","PSI_FLOOP_012 ","PSI_FLOOP_013 ","PSI_FLOOP_014 ","PSI_FLOOP_015 ","PSI_FLOOP_016 ","PSI_FLOOP_017 ","PSI_FLOOP_018 ","PSI_FLOOP_019 ","PSI_FLOOP_020 ","PSI_FLOOP_021 ","PSI_FLOOP_023 ","PSI_FLOOP_024 ","PSI_FLOOP_025 ","PSI_FLOOP_026 ","PSI_FLOOP_027 ","PSI_FLOOP_029 "]
def convert_psi_list_format(inpsi_list):
    outpsi_list = []
    for item in inpsi_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"FLOOP.L{number_part}.PSI"
        outpsi_list.append(new_format)
    return outpsi_list
psi_list = convert_psi_list_format(psi_list)
n_psi = len(psi_list)
psi_data = np.zeros((n_time, n_psi))
time_experimental = mag.get("TIME")
for i_psi in range(0, n_psi):
    psi_experimental = mag.get(psi_list[i_psi])
    psi_recon_time = np.interp(times_to_reconstruct, time_experimental, psi_experimental)
    psi_data[:, i_psi] = psi_recon_time


# bp_probes_names = ["B_BPPROBE_101 ","B_BPPROBE_102 ","B_BPPROBE_103 ","B_BPPROBE_104 ","B_BPPROBE_105 ","B_BPPROBE_109 ","B_BPPROBE_112 ","B_BPPROBE_113 ","B_BPPROBE_114 ","B_BPPROBE_115 ","B_BPPROBE_118 ","B_BPPROBE_119 ","B_BPPROBE_120 ","B_BPPROBE_121 ","B_BPPROBE_122 ","B_BPPROBE_124 ","B_BPPROBE_125 ","B_BPPROBE_126 ","B_BPPROBE_127 ","B_BPPROBE_128 ","B_BPPROBE_131 ","B_BPPROBE_134 "]
bp_probes_names = ["B_BPPROBE_101 ","B_BPPROBE_102 ","B_BPPROBE_103 ","B_BPPROBE_104 ","B_BPPROBE_105 ","B_BPPROBE_109 ","B_BPPROBE_112 ","B_BPPROBE_113 ","B_BPPROBE_114 ","B_BPPROBE_115 ","B_BPPROBE_118 ","B_BPPROBE_119 ","B_BPPROBE_121 ","B_BPPROBE_122 ","B_BPPROBE_124 ","B_BPPROBE_125 ","B_BPPROBE_126 ","B_BPPROBE_127 ","B_BPPROBE_128 ","B_BPPROBE_131 ","B_BPPROBE_134 "]
def convert_bp_list_format(inbp_list):
    outbp_list = []
    for item in inbp_list:
        number_part = item.split('_')[-1].strip()
        new_format = f"BPPROBE.P{number_part}.B"
        outbp_list.append(new_format)
    return outbp_list
bp_probes_names = convert_bp_list_format(bp_probes_names)
n_bp_probe = len(bp_probes_names)
bp_probe_data = np.zeros((n_time, n_bp_probe))
time_experimental = mag.get("TIME")
for i_bp_probe in range(0, n_bp_probe):
    bp_probe_experimental = mag.get(bp_probes_names[i_bp_probe])
    bp_probe_recon_time = np.interp(times_to_reconstruct, time_experimental, bp_probe_experimental)
    bp_probe_data[:, i_bp_probe] = bp_probe_recon_time

# ,"I_ROG_INIVC000","I_ROG_BVLT ","I_ROG_BVLB ","V_FLOOP_024 ","V_FLOOP_022 ","V_FLOOP_020 ","V_FLOOP_007 ","V_FLOOP_009 ","V_FLOOP_011 "
rogowski_coils_names = [
    "ROG.INIVC000.I",
    "ROG.BVLT.I",
    "ROG.BVLB.I",
    # "ROG.GASBFLT.I",
    # "ROG.GASBFLB.I",
    # "ROG.HFSPSRT.I",
    # "ROG.HFSPSRB.I",
    # "ROG.DIVPSRT.I",
    # "ROG.DIVPSRB.I",
]
n_rogowski_coil = len(rogowski_coils_names)
rogowski_coil_data = np.zeros((n_time, n_rogowski_coil))
time_experimental = mag.get("TIME")
for i_rogowski_coil in range(0, n_rogowski_coil):
    rogowski_coil_experimental = mag.get(rogowski_coils_names[i_rogowski_coil])
    rogowski_coil_recon_time = np.interp(times_to_reconstruct, time_experimental, rogowski_coil_experimental)
    rogowski_coil_data[:, i_rogowski_coil] = rogowski_coil_recon_time
# Make the last two columns zero as 
# rogowski_coil_data[:, -2:] = np.zeros((n_time, 2))
# rogowski_coil_data[:, 0] = rogowski_coil_data[:, 0] - 8.0e3

v_floop_names = [
    "FLOOP.L024.V",
    "FLOOP.L022.V",
    "FLOOP.L020.V",
    "FLOOP.L007.V",
    "FLOOP.L009.V",
    "FLOOP.L011.V"
]
n_v_floop = len(v_floop_names)
v_floop_data = np.zeros((n_time, n_v_floop))
time_experimental = mag.get("TIME")
for i_v_floop in range(0, n_v_floop):
    v_floop_experimental = mag.get(v_floop_names[i_v_floop])
    v_floop_recon_time = np.interp(times_to_reconstruct, time_experimental, v_floop_experimental)
    v_floop_data[:, i_v_floop] = v_floop_recon_time


## PF coil currents
psu2coil = GetData(pulseNo, "PSU2COIL#RUN05")

# Order:
# SOL, MC, DIV, BVL, BVUT, BVUB, PSH
i_pf_rtgsfit_order = np.vstack((
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.SOL.I")),   # SOL
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.MCT.I")),    # MCT
    np.interp(times_to_reconstruct, time_experimental, psu2coil.get("PF.MCT.I")),    # MCB
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
import mdsthin
conn = mdsthin.Connection("smaug")
conn.openTree("RTGSFIT", pulseNo_write)
flux_norm = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT.INITIAL_COND:FLUX_NORM").data()
mask = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT.INITIAL_COND:MASK").data()
psi_total = conn.get(f"\\RTGSFIT::TOP.{run_name}.PRESHOT.INITIAL_COND:PSI_TOTAL").data()
conn.disconnect()

# flux_norm = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/flux_norm.txt', dtype=np.float64)
# mask = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/mask.txt', dtype=np.float64)
# psi_total = np.loadtxt('/home/alex.prokopyszyn/GitHub/rtgsfit/data/psi_total.txt', dtype=np.float64)


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
    meas = np.concatenate((psi_data[i_time, :], bp_probe_data[i_time, :], rogowski_coil_data[i_time, :], v_floop_data[i_time, :]))
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

    print(f"meas_python={meas}")
    print(f"meas_python.shape={meas.shape}")

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

print("Writing to MDSplus")
print("Creating script nodes")
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
print("Writing data to MDSplus")
util.write_script_data(
    script_name="RTGSFIT",
    pulseNo_write=pulseNo_write,
    data_to_write=results.to_dictionary(),
    pulseNo_cal=None,
    run_name=run_name,
    run_description=run_description,
)
print("Finished writing to MDSplus")