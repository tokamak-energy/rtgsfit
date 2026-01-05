#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/constants.h"
#include "../src/find_plasma.h"

#define TOL 1e-6 // tolerance for comparing doubles

// CSV reader using N_R and N_Z from constants
int read_flux_csv(const char *filename, double *flux) {
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    perror("Error opening file");
    return -1;
  }

  char line[1024];
  int row = 0;

  while (fgets(line, sizeof(line), fp) && row < N_Z) {
    char *token = strtok(line, ",");
    int col = 0;
    while (token != NULL && col < N_R) {
      flux[row * N_R + col] = atof(token);
      token = strtok(NULL, ",");
      col++;
    }
    if (col != N_R) {
      fprintf(stderr, "Unexpected number of columns at row %d\n", row);
      fclose(fp);
      return -1;
    }
    row++;
  }

  if (row != N_Z) {
    fprintf(stderr, "Unexpected number of rows: %d\n", row);
    fclose(fp);
    return -1;
  }

  fclose(fp);
  return 0;
}

// Helper to compare doubles within a tolerance
int double_equal(double a, double b, double tol) { return fabs(a - b) < tol; }

int main() {
  double flux[N_R * N_Z];

  if (read_flux_csv("../test_data/flux_65x33.csv", flux) != 0) {
    return 1;
  }

  double opt_r[10], opt_z[10], opt_flux[10];
  int32_t opt_n = 0;

  double xpt_r[10], xpt_z[10], xpt_flux[10];
  int32_t xpt_n = 0;

  find_nulls(flux, opt_r, opt_z, opt_flux, &opt_n, xpt_r, xpt_z, xpt_flux,
             &xpt_n);

  // === Check counts ===
  if (opt_n != 3) {
    printf("FAILED: Expected 3 O-points, found %d\n", opt_n);
    return 1;
  }
  if (xpt_n != 3) {
    printf("FAILED: Expected 3 X-points, found %d\n", xpt_n);
    return 1;
  }

  // === Expected O-points ===
  double expected_opt_r[3] = {0.14585184486585667, 0.5504987458865702,
                              0.14588535485454304};
  double expected_opt_z[3] = {-0.3649616878952398, -0.00033762600061724556,
                              0.3660190742308108};

  // === Expected X-points ===
  double expected_xpt_r[3] = {0.3153374696769646, 0.14116510830539897,
                              0.31518440224014677};
  double expected_xpt_z[3] = {-0.5379372821229703, -0.0001838273485009369,
                              0.5378137974139481};

  // === Matching results ===
  int matched_O[3] = {0};
  int matched_X[3] = {0};

  // === Match O-points ===
  for (int i = 0; i < opt_n; i++) {
    int found = 0;
    for (int j = 0; j < 3; j++) {
      if (!matched_O[j] && double_equal(opt_r[i], expected_opt_r[j], TOL) &&
          double_equal(opt_z[i], expected_opt_z[j], TOL)) {
        matched_O[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("FAILED: O-point not expected: (%g, %g)\n", opt_r[i], opt_z[i]);
      return 1;
    }
  }

  // === Match X-points ===
  for (int i = 0; i < xpt_n; i++) {
    int found = 0;
    for (int j = 0; j < 3; j++) {
      if (!matched_X[j] && double_equal(xpt_r[i], expected_xpt_r[j], TOL) &&
          double_equal(xpt_z[i], expected_xpt_z[j], TOL)) {
        matched_X[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("FAILED: X-point not expected: (%g, %g)\n", xpt_r[i], xpt_z[i]);
      return 1;
    }
  }

  // === Success ===
  printf("Test PASSED\n");
  printf("Found O-points:\n");
  for (int i = 0; i < opt_n; i++) {
    printf("  (%g, %g)\n", opt_r[i], opt_z[i]);
  }

  printf("Found X-points:\n");
  for (int i = 0; i < xpt_n; i++) {
    printf("  (%g, %g)\n", xpt_r[i], xpt_z[i]);
  }

  return 0;
}