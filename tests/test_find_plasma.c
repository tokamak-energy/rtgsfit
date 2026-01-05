#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/constants.h"
#include "../src/find_plasma.h"

#define TOL 1e-6

// Maximum number of points this test is willing to compare.
// This is a test-only bound chosen to avoid variable length arrays (VLA) and
// heap allocation while keeping stack usage obvious and portable. If this limit
// is exceeded, the test should fail loudly rather than invoke undefined
// behavior.
#define MAX_TEST_POINTS 64

// -------- CSV reader using N_R and N_Z from constants --------
static int read_flux_csv(const char *filename, double *flux) {
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

// -------- Helpers --------
static int double_equal(double a, double b, double tol) {
  return fabs(a - b) < tol;
}

static void fail_count(const char *label, int got, int expected) {
  fprintf(stderr, "FAILED: Expected %d %s, found %d\n", expected, label, got);
}

static void print_points(const char *label, const double *r, const double *z,
                         int32_t n) {
  printf("%s (%d):\n", label, (int)n);
  for (int32_t i = 0; i < n; i++) {
    printf("  (%0.17g, %0.17g)\n", r[i], z[i]);
  }
}

static int match_points_unordered(const char *label, const double *got_r,
                                  const double *got_z, int32_t got_n,
                                  const double *exp_r, const double *exp_z,
                                  int32_t exp_n, double tol) {
  if (got_n != exp_n) {
    fail_count(label, (int)got_n, (int)exp_n);
    return 0;
  }

  // One-to-one bookkeeping for unordered point comparison.
  // matched[j] == 1 means expected point j has already been matched.
  // A fixed-size array is used intentionally (no VLAs / no malloc);
  // the runtime check below ensures we never write out of bounds.
  int matched[MAX_TEST_POINTS] = {0};

  if (exp_n > MAX_TEST_POINTS) {
    fprintf(stderr,
            "FAILED: %s: expected %d points, exceeds MAX_TEST_POINTS=%d\n",
            label, (int)exp_n, MAX_TEST_POINTS);
    return 0;
  }

  for (int32_t i = 0; i < got_n; i++) {
    int found = 0;
    for (int32_t j = 0; j < exp_n; j++) {
      if (!matched[j] && double_equal(got_r[i], exp_r[j], tol) &&
          double_equal(got_z[i], exp_z[j], tol)) {
        matched[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found) {
      fprintf(stderr, "FAILED: %s point not expected: (%0.17g, %0.17g)\n",
              label, got_r[i], got_z[i]);
      return 0;
    }
  }

  return 1;
}

// -------- find_nulls() test --------
static int test_find_nulls(void) {
  double flux[N_R * N_Z];
  if (read_flux_csv("../test_data/flux_65x33.csv", flux) != 0)
    return 0;

  double opt_r[MAX_TEST_POINTS];
  double opt_z[MAX_TEST_POINTS];
  double opt_flux[MAX_TEST_POINTS];
  int32_t opt_n = 0;

  double xpt_r[MAX_TEST_POINTS];
  double xpt_z[MAX_TEST_POINTS];
  double xpt_flux[MAX_TEST_POINTS];
  int32_t xpt_n = 0;

  (void)find_nulls(flux, opt_r, opt_z, opt_flux, &opt_n, xpt_r, xpt_z, xpt_flux,
                   &xpt_n);

  // expected
  const int32_t expected_opt_n = 3;
  const double expected_opt_r[3] = {0.14585184486585667, 0.5504987458865702,
                                    0.14588535485454304};
  const double expected_opt_z[3] = {
      -0.3649616878952398, -0.00033762600061724556, 0.3660190742308108};

  const int32_t expected_xpt_n = 3;
  const double expected_xpt_r[3] = {0.3153374696769646, 0.14116510830539897,
                                    0.31518440224014677};
  const double expected_xpt_z[3] = {-0.5379372821229703, -0.0001838273485009369,
                                    0.5378137974139481};

  // checks
  if (!match_points_unordered("O-points", opt_r, opt_z, opt_n, expected_opt_r,
                              expected_opt_z, expected_opt_n, TOL)) {
    print_points("Found O-points", opt_r, opt_z, opt_n);
    return 0;
  }

  if (!match_points_unordered("X-points", xpt_r, xpt_z, xpt_n, expected_xpt_r,
                              expected_xpt_z, expected_xpt_n, TOL)) {
    print_points("Found X-points", xpt_r, xpt_z, xpt_n);
    return 0;
  }

  return 1;
}

// -------- main --------
int main(void) {
  if (!test_find_nulls())
    return 1;
  printf("Test PASSED: find_nulls\n");
  return 0;
}
