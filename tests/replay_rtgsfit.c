#include "mds_tools.h"
#include "replay_rtgsfit.h"

int main(int argc, char *argv[]) {
  /* local vars */

  char *endptr, *input_str;
  int sigdesc;     /* signal descriptor */
  int timedesc;    /* descriptor for timebase */
  int i, iter, len, ret;
  char msg[msg_size];
  server = "smaug";
  tree_name = "tepcs";
  char *decm_sig = "acq2106_032.decimate";
  float decimate;
  int dec;
  char signal_name[signal_name_len];
  char *fn_sensors_idx = "../data/sensor_index.txt";
  char *fn_pf_coil_idx = "../data/pfcoil_index.txt";
  char *fn_meas = "../data/meas.txt";
  char *fn_coil_curr = "../data/coil_curr.txt";
  char *fn_flux_norm = "../data/flux_norm.txt";
  char *fn_mask = "../data/mask.txt";
  char *fn_psi_total = "../data/psi_total.txt";
  FILE *p_fid;
  uint32_t sensors_num, pf_coils_num;

  /* Convert input argument to unsigned long long with base 10; */
  check_pulseNo_validity(argv[1], filename, &pulseNo);
  check_input_float_validity(argv[2], filename, &t_start);

  uint32_t sensors_idx_num, pf_coils_idx_num;
  if (count_lines(fn_sensors_idx, &sensors_idx_num) != 0) {
    errorExit("Cannot count lines");
  }

  if (N_MEAS != sensors_idx_num) {
    snprintf(error_msg, sizeof(error_msg), "Incorrect amount of sensors index "
    "%d %d\n", N_MEAS, sensors_idx_num);
    errorExit(error_msg);
  }
  int sensors_idx[sensors_idx_num];
  if (get_int_from_file(fn_sensors_idx, sensors_idx_num, sensors_idx) != 0) {
    errorExit("Cannot get array for coil sensors index.\n");
  }
  if (count_lines(fn_pf_coil_idx, &pf_coils_idx_num) != 0) {
    errorExit("Cannot count lines");
  }
  if (N_COIL != pf_coils_idx_num) {
    snprintf(error_msg, sizeof(error_msg), "Incorrect amount of pf_coils "
      "%d %d\n", N_COIL, pf_coils_idx_num);
    errorExit(error_msg);
  }
  int32_t coil_idx[N_COIL];
  if (get_int_from_file(fn_pf_coil_idx, N_COIL, coil_idx) != 0) {
    errorExit("Cannot get array for coil current index.\n");
  }
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

  /* use get_signal_length to get size of signal */
  snprintf(signal_name, signal_name_len, "%s.CH%03u", pcs1, 1);
  len_pcs1_sign = get_signal_length(server, tree_name, pulseNo, signal_name);
  if (len_pcs1_sign < 1) {
    errorExit("Error retrieving length of signal pcs1 channel 1.");
  }
  /* use get_signal_length to get size of signal */
  snprintf(signal_name, signal_name_len, "%s.CH%03u", pcs2, 1);
  len_pcs2_sign = get_signal_length(server, tree_name, pulseNo, signal_name);
  if (len_pcs2_sign < 1) {
    errorExit("Error retrieving length of signal pcs2 channel 1.");
  }
  /* We assume that pcs1 and pcs2 collects the same amount of data */
  if (len_pcs1_sign != len_pcs2_sign) {
    errorExit("Different data lengths. This is not implemented.");
  }
  N_iter = len_pcs1_sign;
  /* use calloc() to allocate memory for all the pcs1 channels + time_axis */
  pcs1_data_len = len_pcs1_sign * num_pcs1_channels;
  pcs1_data = (float *)calloc(pcs1_data_len, sizeof(float));
  if (!pcs1_data) {
    errorExit("Failed to allocate memory block for pcs1 data.");
  }
  ret = get_dtacqbox(pulseNo, server, tree_name, pcs1, num_pcs1_channels,
    len_pcs1_sign, pcs1_data, &size_pcs1_data);
  if (ret == EXIT_FAILURE) {
    errorExit("pcs1 data not returned");
  }
  /* Get acq2106_032 time */
  float *pcs1_time = (float *)calloc(len_pcs1_sign, sizeof(float));
  size_t pcs1_time_length;
  snprintf(signal_name, signal_name_len, "%s", pcs1);
  ret = get_time(pulseNo, server, tree_name, signal_name, len_pcs1_sign,
    pcs1_time, &pcs1_time_length);

  /* get pcs2 data */
  /* use calloc() to allocate memory for all the pcs2 channels + time_axis */
  pcs2_data_len = len_pcs2_sign * num_pcs2_channels;
  pcs2_data = (float *)calloc(pcs2_data_len, sizeof(float));
  if (!pcs2_data) {
    errorExit("Failed to allocate memory block for pcs2 data.");
  }
  ret = get_dtacqbox(pulseNo, server, tree_name, pcs2, num_pcs2_channels,
    len_pcs2_sign, pcs2_data, &size_pcs2_data);
  if (ret == EXIT_FAILURE) {
    errorExit("pcs2 data not returned");
  }
  float *pcs2_time = (float *)calloc(len_pcs2_sign, sizeof(float));
  size_t pcs2_time_length;
  snprintf(signal_name, signal_name_len, "%s", pcs2);
  ret = get_time(pulseNo, server, tree_name, signal_name, len_pcs2_sign,
    pcs2_time, &pcs2_time_length);

  /* get decimate parameter */
  ret = get_value(
    pulseNo, server, tree_name, decm_sig, &decimate);
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

  /* get TEPCS time from PCS1 */
  /* use get_signal_length to get size of signal */
  snprintf(signal_name, signal_name_len, "%s:PCS_TIME", pcs1);
  int pcs1_dec_time_length;
  pcs1_dec_time_length = get_signal_length(server, tree_name, pulseNo, signal_name);
  if (pcs1_dec_time_length < 1) {
    errorExit("Error retrieving length of signal pcs1 channel 1.");
  }
  size_t ret_pcs1_dec_time_length;
  ret = get_signal(pulseNo, server, tree_name, signal_name,
  pcs1_dec_time_length, pcs1_dec_time, &ret_pcs1_dec_time_length);

  /* decimate pcs2 data and its time axis */
  pcs2_dec_data = (float *)calloc(pcs2_data_len/dec, sizeof(float));
  pcs2_dec_time = (float *)calloc(dec_N_iter, sizeof(float));
  ret = downsample(pcs2_data, pcs2_data_len, dec, pcs2_dec_data, &pcs2_dec_len);
  if (ret == EXIT_FAILURE) {
    errorExit("replay_st40pcs: data downsample doesn't work");
  }

  /* get TEPCS time from PCS1 */
  /* use get_signal_length to get size of signal */
  snprintf(signal_name, signal_name_len, "%s:PCS_TIME", pcs2);
  int pcs2_dec_time_length;
  pcs2_dec_time_length = get_signal_length(server, tree_name, pulseNo, signal_name);
  if (pcs2_dec_time_length < 1) {
    errorExit("Error retrieving length of signal pcs1 channel 1.");
  }
  size_t ret_pcs2_dec_time_length;
  ret = get_signal(pulseNo, server, tree_name, signal_name,
  pcs2_dec_time_length, pcs2_dec_time, &ret_pcs2_dec_time_length);

  /* ST40PCS part */
  /* timing vectors */
  struct timespec t0, t1, t2, t3;
  int chn, istep;
  /* All input channels 384 */
  /* Output channels without CALC 2*32 AO + 16 DO */
  short *outputs;
  outputs = calloc(num_of_out_channels * dec_N_iter, sizeof(short));
  if (outputs == NULL) {
    errorExit("test_st40pcs_mds: Cannot allocate memory.");
  }
  char channel_num[3];

  float *pcs_dec_data = (float *)calloc(pcs1_dec_len + pcs2_dec_len,
    sizeof(float));
  memcpy(pcs_dec_data, pcs1_dec_data, pcs1_dec_len * sizeof(float));
  memcpy(&pcs_dec_data[pcs1_dec_len], pcs2_dec_data,
     (pcs2_dec_len) * sizeof(float));
  char str[100];
  float t_step = (pcs1_dec_time[dec_N_iter - 1] - pcs1_dec_time[0]) / dec_N_iter;
  /* t_start is referenced to plasma breakdown */
  long int t_start_idx = (long int) ((t_start - pcs1_dec_time[0]) / t_step);
  long *t10 = calloc((size_t)N_iter, sizeof(long));
  if (t10 == NULL) {
    errorExit("test_st40pcs: Error! memory t10 not allocated.");
  }
  long *t21 = calloc((size_t)N_iter, sizeof(long));
  if (t21 == NULL) {
    errorExit("test_st40pcs: Error! memory t21 not allocated.");
  }
  long *t32 = calloc((size_t)N_iter, sizeof(long));
  if (t32 == NULL) {
    errorExit("test_st40pcs_mds: Error! memory t32 not allocated.");
  }
  for (i = 0; i < sensors_idx_num; ++i) {
    if (sensors_idx[i] >= 0) {
      meas[i] = (double)(pcs_dec_data[t_start_idx * (num_pcs1_channels + num_pcs2_channels) +
        sensors_idx[i]]);
    } else {
      meas[i] = 0;
    }
  }
  for (i = 0; i < N_COIL; ++i) {
    coil_curr[i] = (double)(pcs_dec_data[t_start_idx * (num_pcs1_channels + num_pcs2_channels) +
     coil_idx[i]]);
  }
  uint32_t N_steps = 1;
  stash_st40runner_float("pcs1_ch001", filename,
    sizeof(float), (size_t)dec_N_iter, &pcs_dec_data[dec_N_iter * sensors_idx[29]]);
  for (istep = 0; istep < N_steps;  ++istep) {
    clock_gettime(CLOCK_MONOTONIC, &t0);
    rtgsfit(meas, coil_curr, flux_norm, mask, psi_total, &error);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    t10[istep] = (long)
      ((t1.tv_sec * 1e9 + t1.tv_nsec - t0.tv_sec * 1e9 - t0.tv_nsec) / 1e3);
      printf("time %ld\n", t10[istep]);
  }
  char *fn_psi_total_out = "../data/psi_total_out_1p2.txt";
  snprintf(str, sizeof(str), "%s_%02d", fn_psi_total_out, 1);
  p_fid = fopen(str, "w");
  for (i = 0; i < N_GRID; ++i) {
    fprintf(p_fid, "%lf\n", psi_total[i]);
  }
  fclose(p_fid);
  // }
  printf("data are stored in %s\n", fn_psi_total_out);

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
  free((void *)outputs);
  free((void *)t10);
  free((void *)t21);
  free((void *)t32);

  /* done */
  return EXIT_SUCCESS;
}
