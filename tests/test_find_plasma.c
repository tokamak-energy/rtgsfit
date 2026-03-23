#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/constants.h"
#include "../src/find_plasma.h"

#define TOL 1e-6

// Maximum number of points needed for for the
// opt and xpt arrays in the tests below.
#define MAX_NUM_TEST_POINTS 6

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

static int read_int_csv(const char *filename, int32_t *data) {
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
      data[row * N_R + col] = (int32_t)strtol(token, NULL, 10);
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
  int matched[MAX_NUM_TEST_POINTS] = {0};

  if (exp_n > MAX_NUM_TEST_POINTS) {
    fprintf(stderr,
            "FAILED: %s: expected %d points, exceeds MAX_NUM_TEST_POINTS=%d\n",
            label, (int)exp_n, MAX_NUM_TEST_POINTS);
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

static int assert_core_side_verbose(const char *case_label, double r, double z,
                                    double r_mag_axis, double z_mag_axis,
                                    double *xpt_r, double *xpt_z, int32_t xpt_n,
                                    int expected) {
  int got =
      is_core_side_of_xpoint(r, z, r_mag_axis, z_mag_axis, xpt_r, xpt_z, xpt_n);

  if (got != expected) {
    fprintf(stderr, "FAILED: %s\n", case_label);
    fprintf(stderr, "  point = (%0.17g, %0.17g)\n", r, z);
    fprintf(stderr, "  mag_axis = (%0.17g, %0.17g)\n", r_mag_axis, z_mag_axis);
    print_points("X-points", xpt_r, xpt_z, xpt_n);
    fprintf(stderr, "  expected is_core_side_of_xpoint = %d, got %d\n",
            expected, got);
    return 0;
  }

  return 1;
}

// -------- find_nulls() test --------
// Verifies that the C null-point finder reproduces the results of the
// reference Python implementation in the find_null_in_gradient_march
// repository (Tokamak Energy GitLab).
// We use data from pulse number 14,926 at time 0.1s.
static int test_find_nulls(void) {
  double flux[N_GRID];
  if (read_flux_csv("../test_data/flux_65x33.csv", flux) != 0) {
    fprintf(stderr,
            "FAILED: test_find_nulls: could not read flux data from '%s'\n",
            "../test_data/flux_65x33.csv");
    return 0;
  }

  double opt_r[MAX_NUM_TEST_POINTS];
  double opt_z[MAX_NUM_TEST_POINTS];
  double opt_flux[MAX_NUM_TEST_POINTS];
  int32_t opt_n = 0;

  double xpt_r[MAX_NUM_TEST_POINTS];
  double xpt_z[MAX_NUM_TEST_POINTS];
  double xpt_flux[MAX_NUM_TEST_POINTS];
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

// -------- filter_xpts() test --------
static int test_filter_xpts(void) {
  double xpt_r[MAX_NUM_TEST_POINTS];
  double xpt_z[MAX_NUM_TEST_POINTS];
  double xpt_flux[MAX_NUM_TEST_POINTS];
  int32_t xpt_n = 6;

  /* initialize only the active entries */
  xpt_r[0] = -0.5;
  xpt_z[0] = +0.0;
  xpt_r[1] = +0.5;
  xpt_z[1] = +0.0;
  xpt_r[2] = +1.0;
  xpt_z[2] = +0.0;
  xpt_r[3] = +0.0;
  xpt_z[3] = -0.5;
  xpt_r[4] = +0.0;
  xpt_z[4] = +0.5;
  xpt_r[5] = +0.0;
  xpt_z[5] = +1.0;

  xpt_flux[0] = -1.0;
  xpt_flux[1] = +0.0;
  xpt_flux[2] = +1.0;
  xpt_flux[3] = +2.0;
  xpt_flux[4] = +3.0;
  xpt_flux[5] = +4.0;

  double r_mag_axis = 0.0;
  double z_mag_axis = 0.0;

  filter_xpts(xpt_r, xpt_z, xpt_flux, &xpt_n, r_mag_axis, z_mag_axis);

  const int32_t expected_xpt_n = 4;
  const double expected_xpt_r[4] = {-0.5, +0.5, +0.0, +0.0};
  const double expected_xpt_z[4] = {+0.0, +0.0, -0.5, +0.5};
  const double expected_xpt_flux[4] = {-1.0, +0.0, +2.0, +3.0};

  if (!match_points_unordered("Filtered X-points", xpt_r, xpt_z, xpt_n,
                              expected_xpt_r, expected_xpt_z, expected_xpt_n,
                              TOL)) {
    print_points("Filtered X-points", xpt_r, xpt_z, xpt_n);
    return 0;
  }

  // Check that the flux values were correctly filtered in an unordered manner.
  for (int32_t i = 0; i < xpt_n; i++) {
    int found = 0;
    for (int32_t j = 0; j < expected_xpt_n; j++) {
      if (double_equal(xpt_flux[i], expected_xpt_flux[j], TOL)) {
        found = 1;
        break;
      }
    }
    if (!found) {
      fprintf(stderr,
              "FAILED: Filtered X-point flux mismatch at index %d: got "
              "%0.17g, not found in expected values\n",
              (int)i, xpt_flux[i]);
      return 0;
    }
  }

  // Check case where xpt_n is zero
  xpt_n = 0;
  filter_xpts(xpt_r, xpt_z, xpt_flux, &xpt_n, r_mag_axis, z_mag_axis);
  if (xpt_n != 0)
    return 0;

  return 1;
}

// -------- is_core_side_of_xpoint() test --------
static int test_is_core_side_of_xpoint(void) {

  double xpt_r[MAX_NUM_TEST_POINTS];
  double xpt_z[MAX_NUM_TEST_POINTS];
  int32_t xpt_n = 4;

  xpt_r[0] = -1.0;
  xpt_z[0] = +0.0;
  xpt_r[1] = +1.0;
  xpt_z[1] = +0.0;
  xpt_r[2] = +0.0;
  xpt_z[2] = +1.0;
  xpt_r[3] = +0.0;
  xpt_z[3] = -1.0;

  double r_mag_axis = 0.0;
  double z_mag_axis = 0.0;

  /* Core-side points */
  if (!assert_core_side_verbose("core-side point #1", 0.5, 0.5, r_mag_axis,
                                z_mag_axis, xpt_r, xpt_z, xpt_n, 1))
    return 0;

  if (!assert_core_side_verbose("core-side point #2", 0.99, 0.0, r_mag_axis,
                                z_mag_axis, xpt_r, xpt_z, xpt_n, 1))
    return 0;

  /* Non-core-side points */
  if (!assert_core_side_verbose("non-core-side point #1", -1.5, 0.0, r_mag_axis,
                                z_mag_axis, xpt_r, xpt_z, xpt_n, 0))
    return 0;

  if (!assert_core_side_verbose("non-core-side point #2", 1.5, 1.5, r_mag_axis,
                                z_mag_axis, xpt_r, xpt_z, xpt_n, 0))
    return 0;

  return 1;
}

// -------- flood_fill_plasma_core() test --------
static int test_flood_fill_plasma_core(void) {
  double flux[N_GRID];
  int32_t mask[N_GRID], expected_mask[N_GRID];
  int32_t xpt_n = 3, lcfs_err_code = 0;
  double xpt_r[3] = {0.3153374696769646, 0.14116510830539897,
                     0.31518440224014677};
  double xpt_z[3] = {-0.5379372821229703, -0.0001838273485009369,
                     0.5378137974139481};
  double r_mag_axis = 0.5504987458865702;
  double z_mag_axis = -0.00033762600061724556;
  double flux_boundary = -0.029051824398111502;

  if (read_flux_csv("../test_data/flux_65x33.csv", flux) != 0) {
    fprintf(stderr, "FAILED: test_flood_fill_plasma_core: flux CSV read failed "
                    "(see previous error output)\n");
    return 0;
  }

  if (read_int_csv("../test_data/mask_65x33.csv", expected_mask) != 0) {
    fprintf(stderr, "FAILED: test_flood_fill_plasma_core: mask CSV read failed "
                    "(see previous error output)\n");
    return 0;
  }

  lcfs_err_code = flood_fill_plasma_core(mask, flux, flux_boundary, r_mag_axis,
                                         z_mag_axis, xpt_r, xpt_z, xpt_n);
  if (lcfs_err_code != ERR_NONE) {
    fprintf(stderr,
            "FAILED: test_flood_fill_plasma_core: flood_fill_plasma_core "
            "returned error code %d\n",
            lcfs_err_code);
    return 0;
  }

  // expected mask (manually verified)
  for (int32_t i = 0; i < N_GRID; i++) {
    if ((int32_t)mask[i] != (int32_t)expected_mask[i]) {
      fprintf(stderr,
              "FAILED: test_flood_fill_plasma_core: mask mismatch at index %d: "
              "expected %d, got %d\n",
              (int)i, (int32_t)expected_mask[i], (int32_t)mask[i]);
      return 0;
    }
  }

  return 1;
}

// -------- sort_xpts() test --------
static int test_sort_xpts(void) {
  double xpt_r[MAX_NUM_TEST_POINTS];
  double xpt_z[MAX_NUM_TEST_POINTS];
  double xpt_flux[MAX_NUM_TEST_POINTS];
  int32_t xpt_n = 6;

  xpt_r[0] = -0.5;
  xpt_z[0] = +0.0;
  xpt_flux[0] = -1.0;
  xpt_r[1] = +0.5;
  xpt_z[1] = +0.0;
  xpt_flux[1] = +0.0;
  xpt_r[2] = +1.0;
  xpt_z[2] = +0.0;
  xpt_flux[2] = +1.0;
  xpt_r[3] = +0.0;
  xpt_z[3] = -0.5;
  xpt_flux[3] = +2.0;
  xpt_r[4] = +0.0;
  xpt_z[4] = +0.5;
  xpt_flux[4] = +3.0;
  xpt_r[5] = +0.0;
  xpt_z[5] = +1.0;
  xpt_flux[5] = +4.0;
  sort_xpts(xpt_r, xpt_z, xpt_flux, xpt_n);

  const double expected_xpt_r[6] = {+0.0, +0.0, +0.0, +1.0, +0.5, -0.5};
  const double expected_xpt_z[6] = {+1.0, +0.5, -0.5, +0.0, +0.0, +0.0};
  const double expected_xpt_flux[6] = {+4.0, +3.0, +2.0, +1.0, +0.0, -1.0};

  for (int32_t i = 0; i < xpt_n; i++) {
    if (!double_equal(xpt_r[i], expected_xpt_r[i], TOL) ||
        !double_equal(xpt_z[i], expected_xpt_z[i], TOL) ||
        !double_equal(xpt_flux[i], expected_xpt_flux[i], TOL)) {
      fprintf(stderr, "FAILED: test_sort_xpts: mismatch at index %d\n", (int)i);
      fprintf(stderr, "  got:    (%0.17g, %0.17g), flux=%0.17g\n", xpt_r[i],
              xpt_z[i], xpt_flux[i]);
      fprintf(stderr, "  expect: (%0.17g, %0.17g), flux=%0.17g\n",
              expected_xpt_r[i], expected_xpt_z[i], expected_xpt_flux[i]);
      return 0;
    }
  }

  return 1;
}

// -------- main --------
int main(void) {

  if (!test_find_nulls())
    return 1;
  printf("Test PASSED: find_nulls\n");
  if (!test_filter_xpts())
    return 1;
  printf("Test PASSED: filter_xpts\n");
  if (!test_is_core_side_of_xpoint())
    return 1;
  printf("Test PASSED: is_core_side_of_xpoint\n");
  if (!test_flood_fill_plasma_core())
    return 1;
  printf("Test PASSED: flood_fill_plasma_core\n");
  if (!test_sort_xpts())
    return 1;
  printf("Test PASSED: sort_xpts\n");

  return 0;
}
