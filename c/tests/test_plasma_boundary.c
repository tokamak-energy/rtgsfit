/*
 * Tests for magnetic_nulls.c and plasma_boundary.c.
 *
 * The psi data is ST40 pulse 14926 at t = 0.1 s on a 33 x 65 grid. The expected axis and
 * x-points come from the reference Python implementation (find_null_in_gradient_march), and
 * the expected core mask was verified by hand.
 */
#include <float.h>

#include "magnetic_nulls.h"
#include "plasma_boundary.h"
#include "test_utilities.h"

enum { N_R = 33, N_Z = 65, N_GRID = N_R * N_Z, N_XPT_MAX = 6 };

static double grid_r[N_R];
static double grid_z[N_Z];
static double psi[N_GRID];
static int32_t mask_limiter[N_GRID];
static int32_t expected_mask[N_GRID];

static const double expected_axis_r = 0.5504987458865702;
static const double expected_axis_z = -0.00033762600061724556;
static const double expected_xpt_r[3] = {0.3153374696769646, 0.14116510830539897,
                                         0.31518440224014677};
static const double expected_xpt_z[3] = {-0.5379372821229703, -0.0001838273485009369,
                                         0.5378137974139481};
static const double boundary_psi = -0.029051824398111502;

static struct rtgsfit_grid test_grid(void)
{
    const struct rtgsfit_grid grid = {N_R, N_Z, 0.030000000000000002, 0.03, grid_r, grid_z};
    return grid;
}

static int32_t load_data(void)
{
    return read_csv("data/grid_r_33.csv", grid_r, N_R) ||
           read_csv("data/grid_z_65.csv", grid_z, N_Z) ||
           read_csv("data/psi_65x33.csv", psi, N_GRID) ||
           read_csv_int("data/mask_limiter_65x33.csv", mask_limiter, N_GRID) ||
           read_csv_int("data/mask_65x33.csv", expected_mask, N_GRID);
}

/* 1 if (r, z) matches one of the expected points that has not been matched yet. */
static int32_t match_point(double r,
                           double z,
                           const double r_expected[],
                           const double z_expected[],
                           int32_t matched[],
                           int32_t n)
{
    for (int32_t i = 0; i < n; i++) {
        if (!matched[i] && fabs(r - r_expected[i]) < 1e-6 && fabs(z - z_expected[i]) < 1e-6) {
            matched[i] = 1;
            return 1;
        }
    }
    return 0;
}

static void test_find_nulls(void)
{
    const struct rtgsfit_grid grid = test_grid();
    double xpt_r[N_XPT_MAX], xpt_z[N_XPT_MAX], xpt_psi[N_XPT_MAX];
    int32_t xpt_n = -1, axis_found = 0;
    double axis_r = 0.0, axis_z = 0.0, axis_psi = 0.0;

    const int32_t status = magnetic_nulls_find(
        &grid, psi, N_GRID, mask_limiter, N_GRID, 1e-10, &axis_found, &axis_r, &axis_z, &axis_psi,
        xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX, xpt_psi, N_XPT_MAX, &xpt_n);
    CHECK(status == RTGSFIT_OK);
    CHECK(axis_found == 1);
    CHECK_CLOSE(axis_r, expected_axis_r, 1e-6);
    CHECK_CLOSE(axis_z, expected_axis_z, 1e-6);
    CHECK(xpt_n == 3);
    int32_t matched[3] = {0, 0, 0};
    for (int32_t i = 0; i < xpt_n && i < 3; i++) {
        CHECK(match_point(xpt_r[i], xpt_z[i], expected_xpt_r, expected_xpt_z, matched, 3));
    }
}

static void test_find_nulls_too_many_xpts(void)
{
    const struct rtgsfit_grid grid = test_grid();
    double xpt_r[2], xpt_z[2], xpt_psi[2];
    int32_t xpt_n = 0, axis_found = 0;
    double axis_r, axis_z, axis_psi;

    const int32_t status =
        magnetic_nulls_find(&grid, psi, N_GRID, mask_limiter, N_GRID, 1e-10, &axis_found, &axis_r,
                            &axis_z, &axis_psi, xpt_r, 2, xpt_z, 2, xpt_psi, 2, &xpt_n);
    CHECK(status == RTGSFIT_ERR_TOO_MANY_XPTS);
    CHECK(xpt_n == 2);
}

/* A saddle exactly on a grid vertex is shared by four cells, but must be found once. */
static void test_find_nulls_on_vertex(void)
{
    enum { N = 7 };
    double r[N], z[N], saddle_psi[N * N];
    int32_t saddle_mask[N * N];
    for (int32_t i = 0; i < N; i++) {
        r[i] = (double)i;
        z[i] = (double)i;
    }
    for (int32_t i_z = 0; i_z < N; i_z++) {
        for (int32_t i_r = 0; i_r < N; i_r++) {
            saddle_psi[i_z * N + i_r] =
                (r[i_r] - 3.0) * (r[i_r] - 3.0) - (z[i_z] - 3.0) * (z[i_z] - 3.0);
            saddle_mask[i_z * N + i_r] = 1;
        }
    }
    const struct rtgsfit_grid grid = {N, N, 1.0, 1.0, r, z};
    double xpt_r[2], xpt_z[2], xpt_psi[2];
    int32_t xpt_n = 0, axis_found = 0;
    double axis_r, axis_z, axis_psi;

    const int32_t status =
        magnetic_nulls_find(&grid, saddle_psi, N * N, saddle_mask, N * N, 1e-10, &axis_found,
                            &axis_r, &axis_z, &axis_psi, xpt_r, 2, xpt_z, 2, xpt_psi, 2, &xpt_n);
    CHECK(status == RTGSFIT_OK);
    CHECK(axis_found == 0);
    CHECK(xpt_n == 1);
    CHECK_CLOSE(xpt_r[0], 3.0, 1e-12);
    CHECK_CLOSE(xpt_z[0], 3.0, 1e-12);
}

static void test_find_nulls_argument_size(void)
{
    const struct rtgsfit_grid grid = test_grid();
    double xpt_r[N_XPT_MAX], xpt_z[N_XPT_MAX], xpt_psi[N_XPT_MAX];
    int32_t xpt_n = 0, axis_found = 0;
    double axis_r, axis_z, axis_psi;

    CHECK(magnetic_nulls_find(&grid, psi, N_GRID - 1, mask_limiter, N_GRID, 1e-10, &axis_found,
                              &axis_r, &axis_z, &axis_psi, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX,
                              xpt_psi, N_XPT_MAX, &xpt_n) == RTGSFIT_ERR_ARGUMENT_SIZE);
    CHECK(magnetic_nulls_find(&grid, psi, N_GRID, mask_limiter, N_GRID, 1e-10, &axis_found, &axis_r,
                              &axis_z, &axis_psi, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX - 1, xpt_psi,
                              N_XPT_MAX, &xpt_n) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_is_core_side(void)
{
    const double xpt_r[4] = {-1.0, 1.0, 0.0, 0.0};
    const double xpt_z[4] = {0.0, 0.0, 1.0, -1.0};

    CHECK(magnetic_nulls_is_core_side(0.5, 0.5, 0.0, 0.0, xpt_r, xpt_z, 4) == 1);
    CHECK(magnetic_nulls_is_core_side(0.99, 0.0, 0.0, 0.0, xpt_r, xpt_z, 4) == 1);
    CHECK(magnetic_nulls_is_core_side(-1.5, 0.0, 0.0, 0.0, xpt_r, xpt_z, 4) == 0);
    CHECK(magnetic_nulls_is_core_side(1.5, 1.5, 0.0, 0.0, xpt_r, xpt_z, 4) == 0);
}

static void test_filter_xpts(void)
{
    double xpt_r[N_XPT_MAX] = {-0.5, 0.5, 1.0, 0.0, 0.0, 0.0};
    double xpt_z[N_XPT_MAX] = {0.0, 0.0, 0.0, -0.5, 0.5, 1.0};
    double xpt_psi[N_XPT_MAX] = {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    int32_t keep[N_XPT_MAX];
    int32_t xpt_n = 6;

    CHECK(magnetic_nulls_filter_xpts(0.0, 0.0, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX, xpt_psi,
                                     N_XPT_MAX, &xpt_n, keep, N_XPT_MAX) == RTGSFIT_OK);
    /* (1, 0) is hidden behind (0.5, 0), and (0, 1) behind (0, 0.5); the order is preserved. */
    const double expected_r[4] = {-0.5, 0.5, 0.0, 0.0};
    const double expected_z[4] = {0.0, 0.0, -0.5, 0.5};
    const double expected_psi[4] = {-1.0, 0.0, 2.0, 3.0};
    CHECK(xpt_n == 4);
    for (int32_t i = 0; i < 4; i++) {
        CHECK_CLOSE(xpt_r[i], expected_r[i], 0.0);
        CHECK_CLOSE(xpt_z[i], expected_z[i], 0.0);
        CHECK_CLOSE(xpt_psi[i], expected_psi[i], 0.0);
    }

    xpt_n = 0;
    CHECK(magnetic_nulls_filter_xpts(0.0, 0.0, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX, xpt_psi,
                                     N_XPT_MAX, &xpt_n, keep, N_XPT_MAX) == RTGSFIT_OK);
    CHECK(xpt_n == 0);

    xpt_n = N_XPT_MAX + 1;
    CHECK(magnetic_nulls_filter_xpts(0.0, 0.0, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX, xpt_psi,
                                     N_XPT_MAX, &xpt_n, keep,
                                     N_XPT_MAX) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_sort_xpts(void)
{
    double xpt_r[4] = {1.0, 2.0, 3.0, 4.0};
    double xpt_z[4] = {-1.0, -2.0, -3.0, -4.0};
    double xpt_psi[4] = {0.1, 0.4, -0.2, 0.3};

    CHECK(magnetic_nulls_sort_xpts(xpt_r, 4, xpt_z, 4, xpt_psi, 4, 4) == RTGSFIT_OK);
    const double expected_psi[4] = {0.4, 0.3, 0.1, -0.2};
    const double expected_r[4] = {2.0, 4.0, 1.0, 3.0};
    for (int32_t i = 0; i < 4; i++) {
        CHECK_CLOSE(xpt_psi[i], expected_psi[i], 0.0);
        CHECK_CLOSE(xpt_r[i], expected_r[i], 0.0);
        CHECK_CLOSE(xpt_z[i], -expected_r[i], 0.0);
    }
}

static void test_flood_fill(void)
{
    const struct rtgsfit_grid grid = test_grid();
    static int32_t mask[N_GRID];
    static int32_t queue[N_GRID];

    const int32_t status = plasma_boundary_flood_fill(
        &grid, psi, N_GRID, mask_limiter, N_GRID, boundary_psi, expected_axis_r, expected_axis_z,
        expected_xpt_r, 3, expected_xpt_z, 3, 3, mask, N_GRID, queue, N_GRID);
    CHECK(status == RTGSFIT_OK);
    int32_t n_different = 0;
    for (int32_t i = 0; i < N_GRID; i++) {
        n_different += mask[i] != expected_mask[i];
    }
    CHECK(n_different == 0);
}

static void test_flood_fill_axis_outside_core(void)
{
    const struct rtgsfit_grid grid = test_grid();
    static int32_t mask[N_GRID];
    static int32_t queue[N_GRID];
    for (int32_t i = 0; i < N_GRID; i++) {
        mask[i] = 7;
    }
    /* With the boundary above every psi value, no point is in the core. */
    const int32_t status = plasma_boundary_flood_fill(
        &grid, psi, N_GRID, mask_limiter, N_GRID, DBL_MAX, expected_axis_r, expected_axis_z,
        expected_xpt_r, 3, expected_xpt_z, 3, 3, mask, N_GRID, queue, N_GRID);
    CHECK(status == RTGSFIT_ERR_AXIS_OUTSIDE_CORE);
    CHECK(mask[0] == 7); /* unchanged on failure */
}

static void test_select_boundary(void)
{
    const double xpt_psi[2] = {0.5, 0.2};
    double boundary_psi_out = 0.0;
    int32_t diverted = -1;

    /* Diverted: 0.9 * 0.5 + 0.1 * 1.0 = 0.55 > limiter 0.3. */
    CHECK(plasma_boundary_select(1.0, 0.3, xpt_psi, 2, 2, 0.9, 1e-10, &boundary_psi_out,
                                 &diverted) == RTGSFIT_OK);
    CHECK_CLOSE(boundary_psi_out, 0.55, 1e-15);
    CHECK(diverted == 1);

    /* Limited. */
    CHECK(plasma_boundary_select(1.0, 0.7, xpt_psi, 2, 2, 0.9, 1e-10, &boundary_psi_out,
                                 &diverted) == RTGSFIT_OK);
    CHECK_CLOSE(boundary_psi_out, 0.7, 0.0);
    CHECK(diverted == 0);

    CHECK(plasma_boundary_select(1.0, 1.0, xpt_psi, 2, 0, 0.9, 1e-10, &boundary_psi_out,
                                 &diverted) == RTGSFIT_ERR_AXIS_EQUALS_BOUNDARY);
    CHECK(plasma_boundary_select(1.0, 2.0, xpt_psi, 2, 0, 0.9, 1e-10, &boundary_psi_out,
                                 &diverted) == RTGSFIT_ERR_BOUNDARY_ABOVE_AXIS);
    CHECK(plasma_boundary_select(1.0, 0.0, xpt_psi, 2, 3, 0.9, 1e-10, &boundary_psi_out,
                                 &diverted) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_limiter_psi(void)
{
    /* Two limiter points on a 2-point "grid"; the second is hidden behind the x-point. */
    const double psi_values[2] = {1.0, 3.0};
    const double limiter_r[2] = {0.5, 2.0};
    const double limiter_z[2] = {0.0, 0.0};
    const int32_t index[8] = {0, 0, 1, 1, 1, 1, 1, 1};
    const double weight[8] = {0.25, 0.25, 0.25, 0.25, 1.0, 0.0, 0.0, 0.0};
    const double xpt_r[1] = {1.5};
    const double xpt_z[1] = {0.0};
    double limiter_psi = 0.0;

    CHECK(plasma_boundary_limiter_psi(psi_values, 2, limiter_r, 2, limiter_z, 2, index, 8, weight,
                                      8, 0.0, 0.0, xpt_r, 1, xpt_z, 1, 0,
                                      &limiter_psi) == RTGSFIT_OK);
    CHECK_CLOSE(limiter_psi, 3.0, 0.0);
    CHECK(plasma_boundary_limiter_psi(psi_values, 2, limiter_r, 2, limiter_z, 2, index, 8, weight,
                                      8, 0.0, 0.0, xpt_r, 1, xpt_z, 1, 1,
                                      &limiter_psi) == RTGSFIT_OK);
    CHECK_CLOSE(limiter_psi, 2.0, 0.0);
    CHECK(plasma_boundary_limiter_psi(psi_values, 2, limiter_r, 2, limiter_z, 2, index, 7, weight,
                                      8, 0.0, 0.0, xpt_r, 1, xpt_z, 1, 1,
                                      &limiter_psi) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_normalise(void)
{
    const double psi_values[3] = {2.0, 1.5, 1.0};
    const int32_t mask[3] = {1, 1, 0};
    const int32_t limiter[3] = {1, 0, 1};
    double psi_norm[3];

    CHECK(plasma_boundary_normalise(psi_values, 3, mask, 3, limiter, 3, 2.0, 1.0, psi_norm, 3) ==
          RTGSFIT_OK);
    CHECK_CLOSE(psi_norm[0], 0.0, 0.0);
    CHECK_CLOSE(psi_norm[1], 1.0, 0.0);
    CHECK_CLOSE(psi_norm[2], 1.0, 0.0);
}

int main(void)
{
    if (load_data() != 0) {
        return 1;
    }
    RUN_TEST(test_find_nulls);
    RUN_TEST(test_find_nulls_too_many_xpts);
    RUN_TEST(test_find_nulls_on_vertex);
    RUN_TEST(test_find_nulls_argument_size);
    RUN_TEST(test_is_core_side);
    RUN_TEST(test_filter_xpts);
    RUN_TEST(test_sort_xpts);
    RUN_TEST(test_flood_fill);
    RUN_TEST(test_flood_fill_axis_outside_core);
    RUN_TEST(test_select_boundary);
    RUN_TEST(test_limiter_psi);
    RUN_TEST(test_normalise);
    return test_failures == 0 ? 0 : 1;
}
