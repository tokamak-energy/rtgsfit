/*
 * File: test_poisson_fast.c
 * -------------------------
 * Unit test for the fast Poisson solver (src/poisson_fast.c) and its
 * integration in poisson_solver().
 *
 * Built against test_data/constants.c (ST40-shaped grid, LU with pivoting):
 *   cd src && cp ../test_data/constants.c . && make SHOT=0 RUN_NAME=no_mds test_poisson_fast && ./test_poisson_fast
 * checks that the fast solver is selected and reproduces solve_tria() and the
 * banded Hagenow flow to rounding, and that it satisfies the LIUQE eq. (45)
 * stencil independently of solve_tria().
 *
 * Built with -DEXPECT_FALLBACK against test_data/constants_poisson_fallback.c
 * (a small grid whose operator is deliberately not separable):
 *   make SHOT=0 RUN_NAME=no_mds test_poisson_fallback && ./test_poisson_fallback
 * checks that initialisation rejects the operator and that poisson_solver()
 * then gives exactly the banded LU result.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "poisson_fast.h"
#include "poisson_solver.h"
#include "solve_tria.h"

#define TOL 1.0e-10

static int n_fail = 0;

#define CHECK(cond, ...)                                             \
    do                                                               \
    {                                                                \
        if (cond)                                                    \
        {                                                            \
            printf("  ok   : " __VA_ARGS__);                         \
        }                                                            \
        else                                                         \
        {                                                            \
            printf("  FAIL : " __VA_ARGS__);                         \
            n_fail++;                                                \
        }                                                            \
        printf("\n");                                                \
    } while (0)

static int is_boundary(int i, int j)
{
    return (i == 0 || i == N_Z - 1 || j == 0 || j == N_R - 1);
}

/* kind 0: smooth source, zero boundary; 1: everything nonzero; 2: all zero. */
static void fill_rhs(int kind, double* b)
{
    for (int i = 0; i < N_Z; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            const int p = i * N_R + j;
            if (kind == 2)
            {
                b[p] = 0.0;
            }
            else if (kind == 1)
            {
                b[p] = cos(0.37 * (double)p) + 0.1 * sin(1.3 * (double)i) * cos(0.7 * (double)j);
            }
            else
            {
                b[p] = is_boundary(i, j) ? 0.0
                     : sin(3.1 * (double)i / (double)N_Z) * cos(2.3 * (double)j / (double)N_R) + 0.5;
            }
        }
    }
}

#ifndef EXPECT_FALLBACK
static double rel_dev(const double* x, const double* y, int n)
{
    double max_x = 0.0, max_d = 0.0;
    for (int p = 0; p < n; p++)
    {
        if (fabs(x[p]) > max_x) max_x = fabs(x[p]);
        if (fabs(x[p] - y[p]) > max_d) max_d = fabs(x[p] - y[p]);
    }
    return (max_x > 0.0) ? max_d / max_x : max_d;
}

static int all_finite(const double* x, int n)
{
    for (int p = 0; p < n; p++)
    {
        if (!isfinite(x[p])) return 0;
    }
    return 1;
}

/* Residual of the LIUQE eq. (45) stencil (full cylindrical form, as generated
 * by GSFit for ST40) for the solution x of A x = b, relative to max |b|. */
static double stencil_residual(const double* x, const double* b)
{
    const double k = (DZ / DR) * (DZ / DR);
    double max_r = 0.0, max_b = 0.0;
    for (int i = 1; i < N_Z - 1; i++)
    {
        for (int j = 1; j < N_R - 1; j++)
        {
            const int p = i * N_R + j;
            const double r = R_VEC[j];
            const double a_j = k * r / (r + DR / 2.0);
            const double b_j = k * r / (r - DR / 2.0);
            const double c_j = 2.0 + a_j + b_j;
            const double res = x[p + N_R] + x[p - N_R] + a_j * x[p + 1] + b_j * x[p - 1] - c_j * x[p] - b[p];
            if (fabs(res) > max_r) max_r = fabs(res);
        }
    }
    for (int p = 0; p < N_GRID; p++)
    {
        if (fabs(b[p]) > max_b) max_b = fabs(b[p]);
    }
    return (max_b > 0.0) ? max_r / max_b : max_r;
}
#endif

int main(void)
{
    double* b1 = malloc((size_t)N_GRID * sizeof(double));
    double* b2 = malloc((size_t)N_GRID * sizeof(double));
    double* x1 = malloc((size_t)N_GRID * sizeof(double));
    double* x2 = malloc((size_t)N_GRID * sizeof(double));

    const int method = poisson_solver_init();
    const int status = poisson_fast_status();
    printf("grid %d x %d, method %d, fast status %d, max deviation at init %.3e\n",
           N_R, N_Z, method, status, poisson_fast_max_deviation());

#ifdef EXPECT_FALLBACK
    CHECK(method == POISSON_SOLVER_METHOD_LU, "banded LU selected for a non-separable operator");
    CHECK(status == POISSON_FAST_ERR_Z_COUPLING, "init reports the z-coupling structure error (status %d)", status);
    for (int kind = 0; kind < 3; kind++)
    {
        fill_rhs(kind, b1);
        memcpy(b2, b1, (size_t)N_GRID * sizeof(double));
        poisson_solver(b1, x1);
        poisson_solver_lu(b2, x2);
        CHECK(memcmp(x1, x2, (size_t)N_GRID * sizeof(double)) == 0 &&
              memcmp(b1, b2, (size_t)N_GRID * sizeof(double)) == 0,
              "poisson_solver() is bitwise identical to poisson_solver_lu() (rhs kind %d)", kind);
    }
#else
    CHECK(method == POISSON_SOLVER_METHOD_FAST, "fast solver selected");
    CHECK(status == POISSON_FAST_OK, "init status ok");
    CHECK(poisson_fast_max_deviation() >= 0.0 && poisson_fast_max_deviation() < TOL,
          "init cross-check deviation %.3e < %.0e", poisson_fast_max_deviation(), TOL);

    for (int kind = 0; kind < 3; kind++)
    {
        fill_rhs(kind, b1);
        solve_tria(b1, x1);
        poisson_fast_solve(b1, x2);
        const double dev = rel_dev(x1, x2, N_GRID);
        CHECK(dev < TOL && all_finite(x2, N_GRID), "single solve matches solve_tria (rhs kind %d): %.3e", kind, dev);
        if (kind == 2)
        {
            int all_zero = 1;
            for (int p = 0; p < N_GRID; p++)
            {
                if (x2[p] != 0.0) all_zero = 0;
            }
            CHECK(all_zero, "zero right-hand side gives an exactly zero solution");
        }
        else
        {
            const double res = stencil_residual(x2, b1);
            CHECK(res < TOL, "fast solution satisfies the LIUQE stencil independently (rhs kind %d): %.3e", kind, res);
        }

        fill_rhs(kind, b1);
        memcpy(b2, b1, (size_t)N_GRID * sizeof(double));
        poisson_solver_lu(b1, x1);
        poisson_fast_poisson_solver(b2, x2);
        const double dev_x = rel_dev(x1, x2, N_GRID);
        const double dev_b = rel_dev(b1, b2, N_GRID);
        CHECK(dev_x < TOL && dev_b < TOL && all_finite(x2, N_GRID),
              "Hagenow flow matches banded flow (rhs kind %d): flux %.3e, boundary %.3e", kind, dev_x, dev_b);
    }

    /* Re-initialisation and lazy initialisation. */
    CHECK(poisson_solver_init() == POISSON_SOLVER_METHOD_FAST, "re-initialisation selects the fast solver again");
    poisson_fast_free();
    CHECK(poisson_fast_status() == POISSON_FAST_NOT_INITIALISED, "poisson_fast_free() resets the status");
    CHECK(poisson_solver_init() == POISSON_SOLVER_METHOD_FAST, "initialisation after free works");
    fill_rhs(1, b1);
    memcpy(b2, b1, (size_t)N_GRID * sizeof(double));
    poisson_solver(b1, x1);
    poisson_solver_lu(b2, x2);
    CHECK(rel_dev(x2, x1, N_GRID) < TOL, "poisson_solver() dispatches to the fast solver: %.3e", rel_dev(x2, x1, N_GRID));
#endif

    free(b1); free(b2); free(x1); free(x2);
    if (n_fail == 0)
    {
        printf("PASSED\n");
        return 0;
    }
    printf("FAILED: %d check(s)\n", n_fail);
    return 1;
}
