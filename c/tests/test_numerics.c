/*
 * Tests for finite_difference.c, current_basis.c, least_squares.c and poisson_solver.c.
 */
#include <lapacke.h>
#include <string.h>

#include "current_basis.h"
#include "finite_difference.h"
#include "least_squares.h"
#include "poisson_solver.h"
#include "test_utilities.h"

enum
{
    N_R = 6,
    N_Z = 7,
    N_GRID = N_R * N_Z,
    N_EDGE = 2 * N_R + 2 * N_Z - 4
};

static double grid_r[N_R];
static double grid_z[N_Z];

static struct rtgsfit_grid make_grid(void)
{
    const double d_r = 0.1;
    const double d_z = 0.15;
    for (int32_t i = 0; i < N_R; i++)
    {
        grid_r[i] = 0.5 + d_r * i;
    }
    for (int32_t i = 0; i < N_Z; i++)
    {
        grid_z[i] = -0.45 + d_z * i;
    }
    const struct rtgsfit_grid grid = {N_R, N_Z, d_r, d_z, grid_r, grid_z};
    return grid;
}

static void test_d_dz_is_exact_for_quadratics(void)
{
    const struct rtgsfit_grid grid = make_grid();
    double f[N_GRID], df_dz[N_GRID];
    for (int32_t i_z = 0; i_z < N_Z; i_z++)
    {
        for (int32_t i_r = 0; i_r < N_R; i_r++)
        {
            const double z = grid_z[i_z];
            f[grid_index(&grid, i_r, i_z)] = 3.0 * z * z - 2.0 * z + grid_r[i_r];
        }
    }
    CHECK(finite_difference_d_dz(&grid, f, N_GRID, df_dz, N_GRID) == RTGSFIT_OK);
    for (int32_t i_z = 0; i_z < N_Z; i_z++)
    {
        for (int32_t i_r = 0; i_r < N_R; i_r++)
        {
            CHECK_CLOSE(df_dz[grid_index(&grid, i_r, i_z)], 6.0 * grid_z[i_z] - 2.0, 1e-12);
        }
    }
    CHECK(finite_difference_d_dz(&grid, f, N_GRID, df_dz, N_GRID + 1) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_edge_d_dn_order_and_values(void)
{
    const struct rtgsfit_grid grid = make_grid();
    double f[N_GRID], df_dn[N_EDGE];
    for (int32_t i = 0; i < N_GRID; i++)
    {
        f[i] = 1.0 + i * i;
    }
    CHECK(finite_difference_edge_d_dn(&grid, f, N_GRID, df_dn, N_EDGE) == RTGSFIT_OK);
    const double side = grid.d_z / grid.d_r;
    const double top = grid.d_r / grid.d_z;
#define F(i_r, i_z) f[grid_index(&grid, (i_r), (i_z))]
    /* First point of each side, and the last point of the edge. */
    CHECK_CLOSE(df_dn[0], side * (2.0 * F(1, 0) - 0.5 * F(2, 0)), 1e-9);
    CHECK_CLOSE(df_dn[N_Z - 1], side * (2.0 * F(1, N_Z - 1) - 0.5 * F(2, N_Z - 1)), 1e-9);
    CHECK_CLOSE(df_dn[N_Z], top * (2.0 * F(1, N_Z - 2) - 0.5 * F(1, N_Z - 3)), 1e-9);
    CHECK_CLOSE(df_dn[N_Z + N_R - 2],
                side * (2.0 * F(N_R - 2, N_Z - 1) - 0.5 * F(N_R - 3, N_Z - 1)), 1e-9);
    CHECK_CLOSE(df_dn[2 * N_Z + N_R - 3], side * (2.0 * F(N_R - 2, 0) - 0.5 * F(N_R - 3, 0)), 1e-9);
    CHECK_CLOSE(df_dn[2 * N_Z + N_R - 2], top * (2.0 * F(N_R - 2, 1) - 0.5 * F(N_R - 2, 2)), 1e-9);
    CHECK_CLOSE(df_dn[N_EDGE - 1], top * (2.0 * F(1, 1) - 0.5 * F(1, 2)), 1e-9);
#undef F
    CHECK(finite_difference_edge_d_dn(&grid, f, N_GRID, df_dn, N_EDGE - 1) ==
          RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_current_basis(void)
{
    const struct rtgsfit_grid grid = make_grid();
    double psi_norm[N_GRID], r_flat[N_GRID], inv_r[N_GRID], basis[3 * N_GRID];
    int32_t mask[N_GRID], mask_limiter[N_GRID];
    for (int32_t i_z = 0; i_z < N_Z; i_z++)
    {
        for (int32_t i_r = 0; i_r < N_R; i_r++)
        {
            const int32_t i = grid_index(&grid, i_r, i_z);
            psi_norm[i] = 0.1 * i_z;
            r_flat[i] = grid_r[i_r];
            inv_r[i] = 1.0 / grid_r[i_r];
            mask[i] = i_r > 1;
            mask_limiter[i] = i_z > 0 && i_z < N_Z - 1;
        }
    }
    CHECK(current_basis_compute(&grid, psi_norm, N_GRID, mask, N_GRID, mask_limiter, N_GRID, r_flat,
                                N_GRID, inv_r, N_GRID, basis, 3 * N_GRID) == RTGSFIT_OK);
    for (int32_t i = 0; i < N_GRID; i++)
    {
        const int32_t inside = mask[i] && mask_limiter[i];
        CHECK_CLOSE(basis[i], inside ? (1.0 - psi_norm[i]) * r_flat[i] : 0.0, 1e-15);
        CHECK_CLOSE(basis[N_GRID + i], inside ? (1.0 - psi_norm[i]) * inv_r[i] : 0.0, 1e-15);
        CHECK_CLOSE(basis[2 * N_GRID + i], inside ? 0.1 / grid.d_z : 0.0, 1e-12);
    }
    CHECK(current_basis_compute(&grid, psi_norm, N_GRID, mask, N_GRID, mask_limiter, N_GRID, r_flat,
                                N_GRID, inv_r, N_GRID, basis,
                                2 * N_GRID) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

static void test_least_squares_matches_lapacke_dgelss(void)
{
    enum
    {
        N_B = 44,
        N_X = 7,
        N_WORK = 4096
    };
    double a[N_B * N_X], b[N_B], a_reference[N_B * N_X], b_reference[N_B];
    double s[N_X], s_reference[N_X], work[N_WORK];
    uint32_t state = 12345u;
    for (int32_t i = 0; i < N_B * N_X; i++)
    {
        state = state * 1664525u + 1013904223u;
        a[i] = (double)(state >> 8) / 16777216.0 - 0.5;
    }
    for (int32_t i = 0; i < N_B; i++)
    {
        state = state * 1664525u + 1013904223u;
        b[i] = (double)(state >> 8) / 16777216.0;
    }
    memcpy(a_reference, a, sizeof(a));
    memcpy(b_reference, b, sizeof(b));

    int32_t rank = 0, info = -1;
    CHECK(least_squares_solve(N_B, N_X, a, N_B * N_X, b, N_B, s, N_X, work, N_WORK, -1.0, &rank,
                              &info) == RTGSFIT_OK);
    lapack_int rank_reference = 0;
    const lapack_int info_reference =
        LAPACKE_dgelss(LAPACK_COL_MAJOR, N_B, N_X, 1, a_reference, N_B, b_reference, N_B,
                       s_reference, -1.0, &rank_reference);
    CHECK(info == 0 && info_reference == 0);
    CHECK(rank == N_X && rank_reference == N_X);
    for (int32_t i = 0; i < N_X; i++)
    {
        CHECK_CLOSE(b[i], b_reference[i], 0.0); /* bit-for-bit */
    }
    CHECK(least_squares_solve(N_B, N_X, a, N_B * N_X, b, N_B, s, N_X, work, 10, -1.0, &rank,
                              &info) == RTGSFIT_ERR_ARGUMENT_SIZE);
    CHECK(least_squares_solve(N_X - 1, N_X, a, (N_X - 1) * N_X, b, N_X - 1, s, N_X, work, N_WORK,
                              -1.0, &rank, &info) == RTGSFIT_ERR_ARGUMENT_SIZE);
    a[0] = NAN;
    CHECK(least_squares_solve(N_B, N_X, a, N_B * N_X, b, N_B, s, N_X, work, N_WORK, -1.0, &rank,
                              &info) == RTGSFIT_ERR_FIT_NONFINITE_INPUT);
    CHECK(info == 0 && rank == 0);
    a[0] = 0.0;
    b[0] = INFINITY;
    CHECK(least_squares_solve(N_B, N_X, a, N_B * N_X, b, N_B, s, N_X, work, N_WORK, -1.0, &rank,
                              &info) == RTGSFIT_ERR_FIT_NONFINITE_INPUT);
}

/*
 * Build the 5-point operator with identity edge rows, factorise it without pivoting (it is
 * diagonally dominant), store the factors in the banded layout, and check A psi = rhs.
 */
static void test_poisson_fixed_boundary(void)
{
    const struct rtgsfit_grid grid = make_grid();
    static double a[N_GRID][N_GRID], lu[N_GRID][N_GRID];
    static double lower_band[N_GRID * N_R], upper_band[N_GRID * (N_R + 2)];
    static int32_t permutation[N_GRID];
    const double ratio = (grid.d_z / grid.d_r) * (grid.d_z / grid.d_r);

    memset(a, 0, sizeof(a));
    for (int32_t i_z = 0; i_z < N_Z; i_z++)
    {
        for (int32_t i_r = 0; i_r < N_R; i_r++)
        {
            const int32_t i = grid_index(&grid, i_r, i_z);
            permutation[i] = i;
            if (i_r == 0 || i_r == N_R - 1 || i_z == 0 || i_z == N_Z - 1)
            {
                a[i][i] = 1.0;
                continue;
            }
            const double r = grid_r[i_r];
            a[i][i - 1] = ratio * r / (r - 0.5 * grid.d_r);
            a[i][i + 1] = ratio * r / (r + 0.5 * grid.d_r);
            a[i][i - N_R] = 1.0;
            a[i][i + N_R] = 1.0;
            a[i][i] = -(a[i][i - 1] + a[i][i + 1] + 2.0);
        }
    }
    memcpy(lu, a, sizeof(a));
    for (int32_t k = 0; k < N_GRID; k++)
    {
        for (int32_t i = k + 1; i < N_GRID; i++)
        {
            lu[i][k] /= lu[k][k];
            for (int32_t j = k + 1; j < N_GRID; j++)
            {
                lu[i][j] -= lu[i][k] * lu[k][j];
            }
        }
    }
    memset(lower_band, 0, sizeof(lower_band));
    memset(upper_band, 0, sizeof(upper_band));
    for (int32_t i = 0; i < N_GRID; i++)
    {
        for (int32_t j = 0; j < N_R; j++)
        {
            const int32_t col = i - N_R + j;
            lower_band[i * N_R + j] = (col >= 0) ? lu[i][col] : 0.0;
        }
        upper_band[i * (N_R + 2)] = 1.0 / lu[i][i];
        for (int32_t j = 1; j < N_R + 2; j++)
        {
            upper_band[i * (N_R + 2) + j] = (i + j < N_GRID) ? lu[i][i + j] : 0.0;
        }
    }

    double edge_inv_mu0_r[N_EDGE], edge_g_edge[N_EDGE * N_EDGE];
    memset(edge_inv_mu0_r, 0, sizeof(edge_inv_mu0_r));
    memset(edge_g_edge, 0, sizeof(edge_g_edge));
    const struct poisson_operator op = {
        grid,   lower_band,     N_GRID * N_R, upper_band,  N_GRID * (N_R + 2), permutation,
        N_GRID, edge_inv_mu0_r, N_EDGE,       edge_g_edge, N_EDGE * N_EDGE};

    double rhs[N_GRID], psi[N_GRID], work[2 * N_GRID + 2 * N_EDGE];
    for (int32_t i = 0; i < N_GRID; i++)
    {
        rhs[i] = sin(0.7 * i) + 0.1 * i;
    }
    CHECK(poisson_solve_fixed_boundary(&op, rhs, N_GRID, psi, N_GRID, work, N_GRID) == RTGSFIT_OK);
    for (int32_t i = 0; i < N_GRID; i++)
    {
        double a_psi = 0.0;
        for (int32_t j = 0; j < N_GRID; j++)
        {
            a_psi += a[i][j] * psi[j];
        }
        CHECK_CLOSE(a_psi, rhs[i], 1e-11);
    }

    /* With a zero edge Green's function the free-boundary solution has psi = 0 on the edge. */
    CHECK(poisson_solve_free_boundary(&op, rhs, N_GRID, psi, N_GRID, work,
                                      poisson_work_size(&grid)) == RTGSFIT_OK);
    CHECK_CLOSE(psi[0], 0.0, 0.0);
    CHECK_CLOSE(psi[N_GRID - 1], 0.0, 0.0);
    CHECK(poisson_solve_free_boundary(&op, rhs, N_GRID, psi, N_GRID, work,
                                      poisson_work_size(&grid) - 1) == RTGSFIT_ERR_ARGUMENT_SIZE);
    struct poisson_operator invalid_op = op;
    invalid_op.n_lower_band--;
    CHECK(poisson_solve_fixed_boundary(&invalid_op, rhs, N_GRID, psi, N_GRID, work, N_GRID) ==
          RTGSFIT_ERR_ARGUMENT_SIZE);
    invalid_op = op;
    invalid_op.permutation = NULL;
    CHECK(poisson_solve_fixed_boundary(&invalid_op, rhs, N_GRID, psi, N_GRID, work, N_GRID) ==
          RTGSFIT_ERR_ARGUMENT_SIZE);
}

int main(void)
{
    RUN_TEST(test_d_dz_is_exact_for_quadratics);
    RUN_TEST(test_edge_d_dn_order_and_values);
    RUN_TEST(test_current_basis);
    RUN_TEST(test_least_squares_matches_lapacke_dgelss);
    RUN_TEST(test_poisson_fixed_boundary);
    return test_failures == 0 ? 0 : 1;
}
