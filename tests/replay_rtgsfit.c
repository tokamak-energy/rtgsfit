#include <netcdf.h>
#include <stdlib.h>
#include "mds_tools.h"
#include "replay_rtgsfit.h"
#include "utils.h"

int count_lines(char *, uint32_t *);

int main(int argc, char *argv[]) {
  /* local vars */
  int i, ret;
  char *server = "smaug";
  char *tree = "tepcs";
  char str[100];
  char sig_name[200];
  double lcfs_r[N_LCFS_MAX];
  double lcfs_z[N_LCFS_MAX];
  int lcfs_n = 0;
  double coef[N_COEF];
  double flux_boundary = 0;
  double plasma_current = 0;

  /* Convert input argument to unsigned long long with base 10; */
  // check_pulseNo_validity(argv[1], filename, &pulseNo);
  // check_input_float_validity(argv[2], filename, &t_start);
  // check_input_float_validity(argv[3], filename, &rtgsfit_step);
  // check_input_float_validity(argv[4], filename, &t_end);
  pulseNo = atoi(argv[1]);
  t_start = atof(argv[2]);
  rtgsfit_step = atof(argv[3]);
  t_end = atof(argv[4]);
  printf("t_start %f t_step %f\n", t_start, rtgsfit_step);
  if ((t_start < -1.1) || (t_start > 0.5)) {
    snprintf(error_msg, sizeof(error_msg), "There is no plasma in selected time"
      " %f\n", t_start);
    errorExit(error_msg);
  }
  if ((rtgsfit_step < 0.00005) || (rtgsfit_step > 0.1)) {
    snprintf(error_msg, sizeof(error_msg), "rtgsfit_step is not set correctly "
      "%f\n", rtgsfit_step);
    errorExit(error_msg);
  }
  if (t_end < t_start) {
    snprintf(error_msg, sizeof(error_msg), "t_end %f is smaller than t_start "
      "%f\n", t_end, t_start);
    errorExit(error_msg);
  }

  uint32_t sensors_idx_num, pf_coils_idx_num;
  if (count_lines(fn_sensors_idx, &sensors_idx_num) != 0) {
    errorExit("Cannot count lines");
  }

  if (N_MEAS != sensors_idx_num) {
    snprintf(error_msg, sizeof(error_msg), "Incorrect amount of sensors index "
    "%d %u\n", N_MEAS, sensors_idx_num);
    errorExit(error_msg);
  }
  int sensors_idx[sensors_idx_num];
  if (get_int_from_file(fn_sensors_idx, sensors_idx_num, sensors_idx) != 0) {
    errorExit("Cannot get array for coil sensors index.\n");
  }
  // for (i = 0; i < sensors_idx_num; i++){
  //   printf("sensors_idx %d\n", sensors_idx[i]);
  // }
  if (count_lines(fn_pf_coil_idx, &pf_coils_idx_num) != 0) {
    errorExit("Cannot count lines");
  }
  if (N_COIL != pf_coils_idx_num) {
    snprintf(error_msg, sizeof(error_msg), "Incorrect amount of pf_coils "
      "%d %u\n", N_COIL, pf_coils_idx_num);
    errorExit(error_msg);
  }
  int32_t coil_idx[N_COIL];
  if (get_int_from_file(fn_pf_coil_idx, N_COIL, coil_idx) != 0) {
    errorExit("Cannot get array for coil current index.\n");
  }
  // for (i = 0; i < N_COIL; i++){
  //   printf("coil_idx %d\n", coil_idx[i]);
  // }
  double meas[N_MEAS];
  if (get_double_from_file(fn_meas, N_MEAS, meas) != 0) {
    errorExit("Cannot get array for meas data.\n");
  }
  double coil_curr[N_COIL];
  if (get_double_from_file(fn_coil_curr, N_COIL, coil_curr) != 0) {
    errorExit("Cannot get array for coil current.\n");
  }
  double flux_norm[N_GRID];
  if (get_double_from_file(fn_flux_norm, N_GRID, flux_norm) != 0) {
    errorExit("Cannot get array for normalised flux.\n");
  }
  int mask[N_GRID];
  if (get_int_from_file(fn_mask, N_GRID, mask) != 0) {
    errorExit("Cannot get array for mask.\n");
  }
  double psi_total[N_GRID];
  if (get_double_from_file(fn_psi_total, N_GRID, psi_total) != 0) {
    errorExit("Cannot get array for total flux.\n");
  }

  double error = 0;

  /* use get_signal_length to get size of sig_name */
  snprintf(sig_name, signal_name_len, "%s.CH%03d", pcs1, 1);
  len_pcs1_sign = get_signal_length(server, pulseNo, tree, sig_name);
  if (len_pcs1_sign < 1) {
    errorExit("Error retrieving length of sig_name pcs1 channel 1.");
  }
  /* use get_signal_length to get size of sig_name */
  snprintf(sig_name, signal_name_len, "%s.CH%03d", pcs2, 1);
  len_pcs2_sign = get_signal_length(server, pulseNo, tree, sig_name);
  if (len_pcs2_sign < 1) {
    errorExit("Error retrieving length of sig_name pcs2 channel 1.");
  }
  /* We assume that pcs1 and pcs2 collects the same amount of data */
  if (len_pcs1_sign != len_pcs2_sign) {
    errorExit("Different data lengths. This is not implemented.");
  }
  N_iter = len_pcs1_sign;
  /* use calloc() to allocate memory for all the pcs1 channels + time_axis */
  pcs1_data_len = len_pcs1_sign * (num_pcs1_channels + 1);
  pcs1_data = (float *)calloc(pcs1_data_len, sizeof(float));
  if (!pcs1_data) {
    errorExit("Failed to allocate memory block for pcs1 data.");
  }
  ret = get_dtacqbox(server, pulseNo, tree, pcs1, num_pcs1_channels,
    len_pcs1_sign, pcs1_data, &size_pcs1_data);
  if (ret == EXIT_FAILURE) {
    errorExit("pcs1 data not returned");
  }

  /* get pcs2 data */
  /* use calloc() to allocate memory for all the pcs2 channels + time_axis */
  pcs2_data_len = len_pcs2_sign * (num_pcs2_channels + 1);
  pcs2_data = (float *)calloc(pcs2_data_len, sizeof(float));
  if (!pcs2_data) {
    errorExit("Failed to allocate memory block for pcs2 data.");
  }
  ret = get_dtacqbox(server, pulseNo, tree, pcs2, num_pcs2_channels,
    len_pcs2_sign, pcs2_data, &size_pcs2_data);
  if (ret == EXIT_FAILURE) {
    errorExit("pcs2 data not returned");
  }

  /* get decimate parameter */
  ret = get_value(server, pulseNo, tree, decm_sig, &decimate);
  if(ret == EXIT_FAILURE) {
    errorExit("Failed to get decimate value.");
  }
  dec = (int)decimate;
  dec_N_iter = N_iter/dec;
  /* decimate pcs1 data and its time axis */
  pcs1_dec_data = (float *)calloc(pcs1_data_len/dec, sizeof(float));
  pcs1_dec_time = (float *)calloc(dec_N_iter, sizeof(float));
  ret = downsample(pcs1_data, pcs1_data_len, dec, pcs1_dec_data, &pcs1_dec_len);
  if (ret == EXIT_FAILURE) {
    errorExit("replay_st40pcs: data downsample doesn't work");
  }

  memcpy(pcs1_dec_time, pcs1_dec_data, sizeof(float) * dec_N_iter);
  printf("here\n");
  /* decimate pcs2 data and its time axis */
  pcs2_dec_data = (float *)calloc(pcs2_data_len/dec, sizeof(float));
  pcs2_dec_time = (float *)calloc(dec_N_iter, sizeof(float));
  ret = downsample(pcs2_data, pcs2_data_len, dec, pcs2_dec_data, &pcs2_dec_len);
  if (ret == EXIT_FAILURE) {
    errorExit("replay_st40pcs: data downsample doesn't work");
  }
  memcpy(pcs2_dec_time, pcs2_dec_data, sizeof(float) * dec_N_iter);
  /* Combine only measurement data from pcs1 and pcs2 into pcs_dec_data excluding
    time axis */
  float *pcs_dec_data = (float *)calloc((num_pcs1_channels + num_pcs1_channels) *
    dec_N_iter, sizeof(float));
  memcpy(pcs_dec_data, &pcs1_dec_data[dec_N_iter], sizeof(float) * dec_N_iter *
    num_pcs1_channels);
  memcpy(&pcs_dec_data[dec_N_iter * num_pcs1_channels], &pcs2_dec_data[dec_N_iter],
    sizeof(float) * dec_N_iter * num_pcs2_channels);
  /* ST40PCS part */
  /* timing vectors */
  struct timespec t0, t1, t2, t3;

  long *t21 = calloc((size_t)N_iter, sizeof(long));
  if (t21 == NULL) {
    errorExit("test_st40pcs: Error! memory t21 not allocated.");
  }
  long *t32 = calloc((size_t)N_iter, sizeof(long));
  if (t32 == NULL) {
    errorExit("test_st40pcs_mds: Error! memory t32 not allocated.");
  }

  float data_step = (pcs1_dec_time[dec_N_iter - 1] - pcs1_dec_time[0]) / dec_N_iter;
  /* t_start is referenced to plasma breakdown */
  long int t_start_idx = (long int) ((t_start - pcs1_dec_time[0]) / data_step);
  long int t_end_idx = (long int) ((t_end - pcs1_dec_time[0]) / data_step);
  uint32_t rtgsfit_step_idx = (uint32_t) (rtgsfit_step / data_step);
  printf("Starting check t_start %ld t_end %ld step %u\n", t_start_idx, t_end_idx, rtgsfit_step_idx);
  uint32_t n_iter = (uint32_t) ((t_end_idx - t_start_idx) / (rtgsfit_step_idx));
  printf("iters %u\n", n_iter);
  long *t10 = (long *)calloc((size_t)n_iter, sizeof(long));
  if (t10 == NULL) {
    errorExit("test_st40pcs: Error! memory t10 not allocated.");
  }

  char *fn_coef_out = "../data/results_coef_";
  char *fn_lcfs_r = "../data/results_lcfs_r_";
  char *fn_lcfs_z = "../data/results_lcfs_z_";

  // Dynamically allocate memory for "output" arrays
  double *time_out = (double *)malloc(n_iter * sizeof(double));
  double *flux_boundary_out = (double *)malloc(n_iter * sizeof(double));
  double *plasma_current_out = (double *)malloc(n_iter * sizeof(double));
  double *flux_total_out = (double *) malloc(n_iter * sizeof(double) * N_GRID);
  double *flux_loops_measured_out = (double *) malloc(n_iter * sizeof(double) * 29);
  double *bp_probes_measured_out = (double *) malloc(n_iter * sizeof(double) * 24);
  double *rogowski_coils_measured_out = (double *) malloc(n_iter * sizeof(double) * 9);
  int *mask_out = (int *) malloc(n_iter * sizeof(int) * N_GRID);

  // Loop over time
  /****************************************************************************/
  int iter = 0;
  for (long int idx = t_start_idx; idx < (t_end_idx - rtgsfit_step_idx);
        idx += rtgsfit_step_idx, ++iter) {
    printf("t_actual %f, t_end %lf, t_start_idx %ld t_end_idx %ld step %u idx %ld\n",
    pcs1_dec_time[idx], pcs1_dec_time[t_end_idx], t_start_idx, t_end_idx, rtgsfit_step_idx, idx);
    for (i = 0; i < sensors_idx_num; ++i) {
      if (sensors_idx[i] >= 0) {
        meas[i] = (double)(pcs_dec_data[dec_N_iter * sensors_idx[i] + idx]) * gains[sensors_idx[i]];
      } else {
        meas[i] = 0;
      }
    }
    for (i = 0; i < N_COIL; ++i) {
      coil_curr[i] = (double)(pcs_dec_data[dec_N_iter * coil_idx[i] + idx]) *
      gains[coil_idx[i]];
    }
    clock_gettime(CLOCK_MONOTONIC, &t0);
    ret = rtgsfit(meas, coil_curr, flux_norm, mask, psi_total, &error, lcfs_r, lcfs_z,
        &lcfs_n, coef, &flux_boundary, &plasma_current);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (ret) {
      printf("lapack returns %d", ret);
    }
    t10[iter] = (long)
        ((t1.tv_sec * 1e9 + t1.tv_nsec - t0.tv_sec * 1e9 - t0.tv_nsec) / 1e3);

    // Store results for writing to NetCDF later
    memcpy(&flux_total_out[iter * N_GRID], psi_total, N_GRID * sizeof(double));
    memcpy(&mask_out[iter * N_GRID], mask, N_GRID * sizeof(int));
    memcpy(&flux_loops_measured_out[iter * 29], meas, 29 * sizeof(double));
    memcpy(&bp_probes_measured_out[iter * 24], &meas[29], 24 * sizeof(double));
    memcpy(&rogowski_coils_measured_out[iter * 9], &meas[29+24], 9 * sizeof(double));
    time_out[iter] = (double)pcs1_dec_time[idx];
    flux_boundary_out[iter] = flux_boundary;
    plasma_current_out[iter] = plasma_current;

    printf("idx steps %ld iter %d n_iter %u\n", (idx - t_start_idx)/rtgsfit_step_idx, iter, n_iter);
  }
  /****************************************************************************/
  // End of the loop
  stash_st40runner("510", filename, sizeof(long), (size_t)iter, (void *)t10);
    printf("about to write to netcdf 01\n");

  int ncid; // NetCDF file ID
  char netcdf_file_name[100];
  snprintf(netcdf_file_name, sizeof(netcdf_file_name),
    "../results/rtgsfit_results_%d.nc", pulseNo);

  // Create the NetCDF file
  if (nc_create(netcdf_file_name, NC_CLOBBER, &ncid) != NC_NOERR) {
      perror("Failed to create NetCDF file");
      exit(1);
  }

  // Define base dimensions
  int dimid_time, dimid_z, dimid_r, dimid_flux_loops, dimid_bp_probes, dimid_rogowski_coils;
  nc_def_dim(ncid, "n_time", n_iter, &dimid_time);
  nc_def_dim(ncid, "n_z", N_Z, &dimid_z);
  nc_def_dim(ncid, "n_r", N_R, &dimid_r);
  nc_def_dim(ncid, "n_flux_loops", 29, &dimid_flux_loops);
  nc_def_dim(ncid, "n_bp_probes", 24, &dimid_bp_probes);
  nc_def_dim(ncid, "n_rogowski_coils", 9, &dimid_rogowski_coils);

  // Define compound dimensions
  int varid_psi, varid_mask, varid_time, varid_psi_b, varid_plasma_current, varid_r, varid_z, varid_flux_loops, varid_bp_probes, varid_rogowski_coils;
  int dimids_time_z_r[3] = {dimid_time, dimid_z, dimid_r};
  int dimids_time[1] = {dimid_time};
  int dimids_r[1] = {dimid_r};
  int dimids_z[1] = {dimid_z};
  int dimids_flux_loops[2] = {dimid_time, dimid_flux_loops};
  int dimids_bp_probes[2] = {dimid_time, dimid_bp_probes};
  int dimids_rogowski_coils[2] = {dimid_time, dimid_rogowski_coils};

  // Define variables
  nc_def_var(ncid, "flux_total", NC_DOUBLE, 3, dimids_time_z_r, &varid_psi);
  nc_def_var(ncid, "mask", NC_INT, 3, dimids_time_z_r, &varid_mask);
  nc_def_var(ncid, "time", NC_DOUBLE, 1, dimids_time, &varid_time);
  nc_def_var(ncid, "flux_boundary", NC_DOUBLE, 1, dimids_time, &varid_psi_b);
  nc_def_var(ncid, "plasma_current", NC_DOUBLE, 1, dimids_time, &varid_plasma_current);
  nc_def_var(ncid, "r", NC_DOUBLE, 1, dimids_r, &varid_r);
  nc_def_var(ncid, "z", NC_DOUBLE, 1, dimids_z, &varid_z);
  nc_def_var(ncid, "flux_loops_measured", NC_DOUBLE, 2, dimids_flux_loops, &varid_flux_loops);
  nc_def_var(ncid, "bp_probes_measured", NC_DOUBLE, 2, dimids_bp_probes, &varid_bp_probes);
  nc_def_var(ncid, "rogowski_coils_measured", NC_DOUBLE, 2, dimids_rogowski_coils, &varid_rogowski_coils);

  // End definitions
  nc_enddef(ncid);

  // Write data into NetCDF file
  nc_put_var_double(ncid, varid_psi, flux_total_out);
  nc_put_var_int(ncid, varid_mask, mask_out);
  nc_put_var_double(ncid, varid_time, time_out);
  nc_put_var_double(ncid, varid_psi_b, flux_boundary_out);
  nc_put_var_double(ncid, varid_plasma_current, plasma_current_out);
  nc_put_var_double(ncid, varid_r, R_VEC);
  nc_put_var_double(ncid, varid_z, Z_VEC);
  nc_put_var_double(ncid, varid_flux_loops, flux_loops_measured_out);
  nc_put_var_double(ncid, varid_bp_probes, bp_probes_measured_out);
  nc_put_var_double(ncid, varid_rogowski_coils, rogowski_coils_measured_out);

  // Close the NetCDF file
  nc_close(ncid);

  printf("finished writing to netcdf\n");

  /* free the dynamically allocated memory when done */
  free((void *)pcs1_data);
  free((void *)pcs1_dec_data);
  free((void *)pcs1_time);
  free((void *)pcs1_dec_time);
  free((void *)pcs2_data);
  free((void *)pcs2_dec_data);
  free((void *)pcs2_dec_time);
  free((void *)pcs2_time);
  free((void *)pcs_dec_data);
  free((void *)t10);
  // free((void *)t21);
  // free((void *)t32);

  free((void *)flux_total_out);
  free((void *)flux_boundary_out);
  free((void *)mask_out);
  free((void *)time_out);

  /* done */
  return EXIT_SUCCESS;
}

int count_lines(char *file_name, uint32_t *lines_count) {
  *lines_count = 0;
  char character;
  p_fid = fopen(file_name, "r");
  if (p_fid == NULL) {
    snprintf(error_msg, sizeof(error_msg), "File does not exist %s.", file_name);
    errorExit(error_msg);
  }
  while(!feof(p_fid)) {
    character = fgetc(p_fid);
    if (character == '\n') {
      (*lines_count)++;
    }
  }
  fclose(p_fid);
  return 0;
}

int get_int_from_file(char *file_name, uint32_t array_size, int32_t *array) {
  int i;
  p_fid = fopen(file_name, "r");
  if (p_fid == NULL) {
    snprintf(error_msg, sizeof(error_msg), "File does not exist %s.", file_name);
    errorExit(error_msg);
  }
  for (i = 0; i < array_size; ++i) {
    fscanf(p_fid, "%d", &array[i]);
  }
  fclose(p_fid);
  return 0;
}

int get_double_from_file(char *file_name, uint32_t array_size, double *array) {
  int i;
  p_fid = fopen(file_name, "r");
  if (p_fid == NULL) {
    snprintf(error_msg, sizeof(error_msg), "File does not exist %s.", file_name);
    errorExit(error_msg);
  }
  for (i = 0; i < array_size; ++i) {
    fscanf(p_fid, "%lf", &array[i]);
  }
  fclose(p_fid);
  return 0;
}
