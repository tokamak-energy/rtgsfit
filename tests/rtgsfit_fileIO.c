#include <stdio.h>
#include "rtgsfit.h"
#include "constants.h"
#include "utils.h"

int main() {
  int i, iter;
  uint32_t lines_count;
  uint32_t sensors_num, pf_coils_num;
  FILE *p_fid;
  char *fn_sensors_idx = "../data/sensor_index.txt";
  char *fn_pf_coil_idx = "../data/pf_coil_index.txt";
  char *fn_meas = "../data/meas.txt";
  char *fn_coil_curr = "../data/coil_curr.txt";
  char *fn_flux_norm = "../data/flux_norm.txt";
  char *fn_mask = "../data/mask.txt";
  char *fn_psi_total = "../data/psi_total.txt";
  char *fn_error = "../data/error.txt";
  char character;

  if (count_lines(fn_sensors_idx, &sensors_num) != 0) {
    errorExit("Cannot count lines");
  }
  // if (N_MEAS != sensors_num) {
  //   snprintf(error_msg, sizeof(error_msg), "Incorrect amount of sensors index %d %d\n", N_MEAS, sensors_num);
  //   errorExit(error_msg);
  // }
  int sensors_idx[N_MEAS];
  if (get_int_from_file(fn_sensors_idx, N_MEAS, sensors_idx) != 0) {
    errorExit("Cannot get array for coil sensors index.\n");
  }
  if (count_lines(fn_pf_coil_idx, &pf_coils_num) != 0) {
    errorExit("Cannot count lines");
  }
  // if (N_COIL != pf_coils_num) {
  //   snprintf(error_msg, sizeof(error_msg), "Incorrect amount of pf_coils %d %d\n", N_COIL, pf_coils_num);
  //   errorExit(error_msg);
  // }
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

  char str[100];
  for (iter = 0; iter < 20; ++iter) {
    char *fn_psi_total_out = "../data/psi_total_out.txt";
    snprintf(str, sizeof(str), "%s_%02d", fn_psi_total_out, iter);
    rtgsfit(meas, coil_curr, flux_norm, mask, psi_total, &error);
    p_fid = fopen(str, "w");
    for (i = 0; i < N_GRID; ++i) {
      fprintf(p_fid, "%lf\n", psi_total[i]);
    }
    fclose(p_fid);
  }
  printf("I am working\n");
  return 0;
}
