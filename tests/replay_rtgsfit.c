#include "mds_tools.h"
#include "replay_rtgsfit.h"
#include <stdlib.h>

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

  /* Convert input argument to unsigned long long with base 10; */
  // check_pulseNo_validity(argv[1], filename, &pulseNo);
  // check_input_float_validity(argv[2], filename, &t_start);
  // check_input_float_validity(argv[3], filename, &rtgsfit_step);
  // check_input_float_validity(argv[4], filename, &t_end);
  pulseNo = atoi(argv[1]);
  t_start= atof(argv[2]);
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
  uint32_t istep = 0;

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
  long int idx = 0;
  uint32_t n_iter = (uint32_t) ((t_end_idx - t_start_idx) / (rtgsfit_step_idx));
  printf("iters %u\n", n_iter);
  long *t10 = (long *)calloc((size_t)n_iter, sizeof(long));
  if (t10 == NULL) {
    errorExit("test_st40pcs: Error! memory t10 not allocated.");
  }

  char *fn_psi_total_out = "../data/psi_total_out_1p2";
  for (idx = t_start_idx; idx < t_end_idx;  idx += rtgsfit_step_idx) {
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
      coil_curr[i] = (double)(pcs_dec_data[dec_N_iter * coil_idx[i] + idx]) * gains[coil_idx[i]];
    }
    clock_gettime(CLOCK_MONOTONIC, &t0);
    rtgsfit(meas, coil_curr, flux_norm, mask, psi_total, &error, lcfs_r, lcfs_z,
        &lcfs_n, coef);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    snprintf(str, sizeof(str), "%s2_%01.4f", fn_psi_total_out, pcs1_dec_time[idx]);
    printf("%s_%01.3f time %lf\n", fn_psi_total_out, pcs1_dec_time[idx], istep * data_step);
    p_fid = fopen(str, "w");
    for (i = 0; i < N_GRID; ++i) {
      fprintf(p_fid, "%lf\n", psi_total[i]);
    }
    fclose(p_fid);
    // // t10[istep] = (long)
    // //   ((t1.tv_sec * 1e9 + t1.tv_nsec - t0.tv_sec * 1e9 - t0.tv_nsec) / 1e3);
    // //   printf("here time %ld\n", t10[istep]);
    // istep++;
  }
  // }
  printf("data are stored in %s_YY\n", fn_psi_total_out);
  // stash_st40runner_float("pcs1_ch001", filename,
  //   sizeof(float), (size_t)dec_N_iter, &pcs_dec_data[dec_N_iter * sensors_idx[29]]);

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
  free((void *)t21);
  free((void *)t32);

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
