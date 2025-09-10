import mdsthin
import numpy as np

class OutputDataRTGSFITClass:
    """
    Class to hold output data for RTGSFIT.
    We will save all the data as a dictionary in the data_dict attribute.
    """

    def __init__(self, cfg):
    
        self.cfg = cfg

        self.mu_0 = cfg["mu_0"]

        with mdsthin.Connection('smaug') as conn:
            conn.openTree("RTGSFIT", cfg["pulse_num_preshot"])
            r_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:R_VEC").data()
            self.z_vec = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:Z_VEC").data()
            n_lcfs_max = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_LCFS_MAX").data()
            n_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_COIL").data()
            n_f_loops = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_F_LOOPS").data()
            n_bp_probes = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_BP_PROBES").data()
            n_rog_coils = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_ROG_COILS").data()
            coil_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:COIL_NAMES").data()
            meas_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:MEAS_NAMES").data()
            self.weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:WEIGHT").data()
            sens_rep_mat = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:SENS_REP_MAT").data()
            self.coef_names = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:COEF_NAMES").data()
            g_grid_meas_weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT.GREENS:GRID_MEAS_W").data()
            g_coef_meas_weight = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT.GREENS:COEF_MEAS_W").data()
            g_meas_coil = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT.GREENS:MEAS_COIL").data()
            self.n_grid = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_GRID").data()
            n_meas = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_MEAS").data()
            n_coef = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_COEF").data()
            self.r_flat = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:R_GRID").data()
            self.n_pls = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_PLS").data()
            self.n_r = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_R").data()
            self.n_z = conn.get(f"\\RTGSFIT::TOP.{cfg["run_name_preshot"]}.PRESHOT:N_Z").data()

        n_t = cfg["n_t"]
        self.flux_loop_indices = []
        self.bp_probe_indices = []
        self.rog_coil_indices = []
        for i, meas_name in enumerate(meas_names):
            if meas_name.startswith("L") and meas_name[1:].isdigit():
                self.flux_loop_indices.append(i)
            elif meas_name.startswith("P") and meas_name[1:].isdigit():
                self.bp_probe_indices.append(i)
            elif meas_name.startswith("IVC"):
                continue
            elif meas_name.startswith("OVC"):
                continue
            else:
                self.rog_coil_indices.append(i)
        self.ivc_indices = []
        for i, coef_name in enumerate(self.coef_names):
            if coef_name.startswith("eig_"):
                self.ivc_indices.append(i)
            elif coef_name == "OVC":
                self.ovc_index = i

        self.g_grid_meas_weight = g_grid_meas_weight.reshape((self.n_grid, n_meas))
        self.g_coef_meas_weight = g_coef_meas_weight.reshape((n_coef, n_meas))
        self.g_meas_coil = g_meas_coil.reshape((n_meas, n_coil)).T
        n_sens = int(np.sqrt(len(sens_rep_mat)))
        self.sens_rep_mat = sens_rep_mat.reshape((n_sens, n_sens))

        self.data_dict = {}
        self.data_dict["GLOBAL"] = {}
        self.data_dict["GLOBAL"]["IP"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["GLOBAL"]["PSI_A"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["GLOBAL"]["PSI_B"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["GLOBAL"]["CHIT"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["GLOBAL"]["DELTA_Z"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["CONSTRAINTS"] = {}
        self.data_dict["CONSTRAINTS"]["COIL"] = {}
        self.data_dict["CONSTRAINTS"]["COIL"]["MVALUE"] = np.zeros((n_t, n_coil), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["COIL"]["NAME"] = coil_names
        self.data_dict["CONSTRAINTS"]["FLOOP"] = {}
        self.data_dict["CONSTRAINTS"]["FLOOP"]["CVALUE"] = np.zeros((n_t, n_f_loops), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["FLOOP"]["MVALUE"] = np.zeros((n_t, n_f_loops), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["FLOOP"]["NAME"] = meas_names[self.flux_loop_indices]
        self.data_dict["CONSTRAINTS"]["FLOOP"]["WEIGHT"] = self.weight[self.flux_loop_indices]
        self.data_dict["CONSTRAINTS"]["BPPROBE"] = {}
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["CVALUE"] = np.zeros((n_t, n_bp_probes), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["MVALUE"] = np.zeros((n_t, n_bp_probes), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["NAME"] = meas_names[self.bp_probe_indices]
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["WEIGHT"] = self.weight[self.bp_probe_indices]
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"] = {}
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["CVALUE"] = np.zeros((n_t, n_rog_coils), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["MVALUE"] = np.zeros((n_t, n_rog_coils), dtype=np.float64)
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["NAME"] = meas_names[self.rog_coil_indices]
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["WEIGHT"] = self.weight[self.rog_coil_indices]
        self.data_dict["PASSIVES"] = {}
        self.data_dict["PASSIVES"]["OVC"] = {}
        self.data_dict["PASSIVES"]["OVC"]["DOF"] = {}
        self.data_dict["PASSIVES"]["OVC"]["DOF"]["CONSTANT_J"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["PASSIVES"]["IVC"] = {}
        self.data_dict["PASSIVES"]["IVC"]["DOF"] = {}
        for i in range(len(self.ivc_indices)):
            self.data_dict["PASSIVES"]["IVC"]["DOF"][f"EIG_{i:02d}"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["PROFILES"] = {}
        self.data_dict["PROFILES"]["SOURCE_FUN"] = {}
        self.data_dict["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"] = {}
        self.data_dict["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["P_PRIME_DOF"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["FF_PRIM_DOF"] = np.zeros(n_t, dtype=np.float64)
        self.data_dict["P_BOUNDARY"] = {}
        self.data_dict["P_BOUNDARY"]["NBND"] = np.zeros(n_t, dtype=np.int32)
        self.data_dict["P_BOUNDARY"]["RBND"] = np.zeros((n_t, n_lcfs_max), dtype=np.float64)
        self.data_dict["P_BOUNDARY"]["ZBND"] = np.zeros((n_t, n_lcfs_max), dtype=np.float64)
        self.data_dict["TWO_D"] = {}
        self.data_dict["TWO_D"]["PSI"] = np.zeros((n_t, self.n_z, self.n_r), dtype=np.float64)
        self.data_dict["TWO_D"]["MASK"] = np.zeros((n_t, self.n_z, self.n_r), dtype=np.int32)
        self.data_dict["TWO_D"]["RGRID"] = r_vec
        self.data_dict["TWO_D"]["ZGRID"] = self.z_vec
        self.data_dict["TIME"] = np.linspace(cfg["t_min"], cfg["t_max"], n_t, dtype=np.float64)

    def update_mvalues(self,
                       iteration : int,
                       meas_pcs : np.ndarray,
                       coil_curr : np.ndarray):
        """
        Update the MVALUE fields in the data_dict based on the current measurements.
        """

        meas_no_reg = meas_pcs @ self.sens_rep_mat.T
        self.data_dict["CONSTRAINTS"]["FLOOP"]["MVALUE"][iteration] = meas_no_reg[self.flux_loop_indices]
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["MVALUE"][iteration] = meas_no_reg[self.bp_probe_indices]
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["MVALUE"][iteration] = meas_no_reg[self.rog_coil_indices]
        self.data_dict["CONSTRAINTS"]["COIL"]["MVALUE"][iteration] = coil_curr

    def update_cvalues(self,
                       iteration : int,
                       flux_norm : np.ndarray,
                       mask : np.ndarray,
                       flux_total : np.ndarray,
                       error : np.ndarray,
                       lcfs_r : np.ndarray,
                       lcfs_z : np.ndarray,
                       lcfs_n : np.ndarray,
                       coef : np.ndarray,
                       flux_boundary : np.ndarray,
                       plasma_current : np.ndarray):

        """
        Update the CVALUE fields in the data_dict based on the predicted values.
        """

        self.data_dict["GLOBAL"]["IP"][iteration] = plasma_current[0]

        # Calculate PSI_A using
        # flux_norm = (psi_a - flux_total) / (psi_a - psi_b)
        # Rearranging gives:
        # psi_a = (flux_norm * psi_b - flux_total) / (flux_norm - 1)
        min_index = np.argmin(flux_norm)
        psi_a = (flux_norm[min_index] * flux_boundary[0] - flux_total[min_index]) / (flux_norm[min_index] - 1)
        self.data_dict["GLOBAL"]["PSI_A"][iteration] = psi_a
    
        self.data_dict["GLOBAL"]["PSI_B"][iteration] = flux_boundary[0]

        self.data_dict["GLOBAL"]["CHIT"][iteration] = error[0]

        pls2_idx = np.where(self.coef_names == "pls2")[0][0]
        self.data_dict["GLOBAL"]["DELTA_Z"][iteration] = coef[pls2_idx]
       
        # Calculate g_pls_grid
        g_pls_grid = np.zeros((self.n_pls, self.n_grid), dtype=np.float64)
        g_pls_grid[0, :] = (1 - flux_norm) * self.r_flat
        g_pls_grid[1, :] = (1 - flux_norm)  / (self.r_flat * self.mu_0)
        dflux_norm_dz = np.gradient(flux_norm.reshape((self.n_z, self.n_r)), self.z_vec, axis=0)
        g_pls_grid[2, :] = dflux_norm_dz.flatten()

        # Calculate first n_pls rows of g_coef_meas_weight
        self.g_coef_meas_weight[:self.n_pls, :] = g_pls_grid @ self.g_grid_meas_weight

        # Now calculate the predicted measurements from the plasma current
        pred_meas = coef @ self.g_coef_meas_weight
        pred_meas = pred_meas / self.weight

        # Now need to calculate the contribution from the PF coils
        coil_curr = self.data_dict["CONSTRAINTS"]["COIL"]["MVALUE"][iteration]
        pred_meas += coil_curr @ self.g_meas_coil

        self.data_dict["CONSTRAINTS"]["FLOOP"]["CVALUE"][iteration] = pred_meas[self.flux_loop_indices]
        self.data_dict["CONSTRAINTS"]["BPPROBE"]["CVALUE"][iteration] = pred_meas[self.bp_probe_indices]
        self.data_dict["CONSTRAINTS"]["ROGOWSKI"]["CVALUE"][iteration] = pred_meas[self.rog_coil_indices]

        for i, ivc_idx in enumerate(self.ivc_indices):
            self.data_dict["PASSIVES"]["IVC"]["DOF"][f"EIG_{i:02d}"][iteration] = coef[ivc_idx]
        self.data_dict["PASSIVES"]["OVC"]["DOF"]["CONSTANT_J"][iteration] = coef[self.ovc_index]

        pls0_idx = np.where(self.coef_names == "pls0")[0][0]
        self.data_dict["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["P_PRIME_DOF"][iteration] = coef[pls0_idx]
        pls1_idx = np.where(self.coef_names == "pls1")[0][0]
        self.data_dict["PROFILES"]["SOURCE_FUN"]["LIUQE_POLY"]["FF_PRIM_DOF"][iteration] = coef[pls1_idx]

        self.data_dict["P_BOUNDARY"]["NBND"][iteration] = lcfs_n[0]
        self.data_dict["P_BOUNDARY"]["RBND"][iteration] = lcfs_r
        self.data_dict["P_BOUNDARY"]["ZBND"][iteration] = lcfs_z
        self.data_dict["TWO_D"]["PSI"][iteration] = flux_total.reshape((self.n_z, self.n_r))
        self.data_dict["TWO_D"]["MASK"][iteration] = mask.reshape((self.n_z, self.n_r))
