/*
 * RT-GSFit: one iteration of the real-time equilibrium reconstruction.
 *
 * This is the only file that includes the preshot. Every algorithm module takes its data as
 * explicit (array, length) arguments, so it can be tested with any grid.
 */
#include "rtgsfit.h"

#include <cblas.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "current_basis.h"
#include "grid.h"
#include "least_squares.h"
#include "magnetic_nulls.h"
#include "measurement_processing.h"
#include "plasma_boundary.h"
#include "poisson_solver.h"
#include "preshot.h"
#include "timing.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define RTGSFIT_HAS_FLUSH_TO_ZERO 1
#endif

/* Workspace for LAPACK dgelss; generous enough for its optimal, blocked, path. */
enum
{
    LEAST_SQUARES_WORK_SIZE =
        128 * (N_FIT_TARGET + N_FIT_DOF) + N_FIT_DOF * N_FIT_DOF + N_FIT_TARGET * N_FIT_DOF,
    POISSON_WORK_SIZE = 2 * N_GRID + 2 * N_EDGE,
};

/*
 * All intermediate arrays. They are static rather than on the stack, because for large grids
 * they would not fit on a real-time thread's stack, and the real-time path must not allocate.
 * As a consequence `rtgsfit()` is not re-entrant: call it from one thread at a time.
 */
static struct
{
    double constraint_coil[N_CONSTRAINT];        /* constraint values from the coils alone */
    double basis[N_PLASMA * N_GRID];             /* plasma current basis functions */
    double fit_matrix[N_FIT_DOF * N_FIT_TARGET]; /* column-major [N_FIT_TARGET x N_FIT_DOF] */
    double fit_matrix_lapack[N_FIT_DOF * N_FIT_TARGET]; /* copy destroyed by LAPACK */
    double fit_target[N_FIT_TARGET];
    double fit_target_lapack[N_FIT_TARGET]; /* copy overwritten by LAPACK */
    double fit_prediction[N_FIT_TARGET];
    double singular_values[N_FIT_DOF];
    double least_squares_work[LEAST_SQUARES_WORK_SIZE];
    double poisson_rhs[N_GRID];
    double psi_plasma[N_GRID];
    double psi_passive[N_GRID];
    double poisson_work[POISSON_WORK_SIZE];
    int32_t xpt_keep[N_XPT_MAX];
    int32_t flood_fill_queue[N_GRID];
} work;

static void enable_flush_to_zero(void)
{
#ifdef RTGSFIT_HAS_FLUSH_TO_ZERO
    /* Denormals are very slow on x86. FTZ is MXCSR bit 15 and DAZ is bit 6. */
    const unsigned int denormals_are_zero = 1u << 6;
    _mm_setcsr(_mm_getcsr() | _MM_FLUSH_ZERO_ON | denormals_are_zero);
#endif
}

static struct rtgsfit_grid preshot_grid(void)
{
    const struct rtgsfit_grid grid = {N_R, N_Z, GRID_D_R, GRID_D_Z, GRID_R, GRID_Z};
    return grid;
}

static struct poisson_operator preshot_poisson_operator(void)
{
    const struct poisson_operator op = {preshot_grid(),
                                        POISSON_LOWER_BAND,
                                        N_GRID * N_R,
                                        POISSON_UPPER_BAND,
                                        N_GRID * (N_R + 2),
                                        POISSON_PERMUTATION,
                                        N_GRID,
                                        EDGE_INV_MU0_R,
                                        N_EDGE,
                                        EDGE_G_EDGE,
                                        N_EDGE * N_EDGE};
    return op;
}

int32_t rtgsfit_get_sizes(int32_t sizes[], const int32_t n_sizes)
{
    RTGSFIT_CHECK_SIZE(n_sizes, RTGSFIT_SIZE_COUNT);
    RTGSFIT_CHECK_ARRAY(sizes, n_sizes);
    sizes[RTGSFIT_SIZE_N_R] = N_R;
    sizes[RTGSFIT_SIZE_N_Z] = N_Z;
    sizes[RTGSFIT_SIZE_N_GRID] = N_GRID;
    sizes[RTGSFIT_SIZE_N_EDGE] = N_EDGE;
    sizes[RTGSFIT_SIZE_N_LIMITER] = N_LIMITER;
    sizes[RTGSFIT_SIZE_N_COIL] = N_COIL;
    sizes[RTGSFIT_SIZE_N_MEAS] = N_MEAS;
    sizes[RTGSFIT_SIZE_N_FLUX_LOOP] = N_FLUX_LOOP;
    sizes[RTGSFIT_SIZE_N_BP_PROBE] = N_BP_PROBE;
    sizes[RTGSFIT_SIZE_N_ROGOWSKI] = N_ROGOWSKI;
    sizes[RTGSFIT_SIZE_N_CONSTRAINT] = N_CONSTRAINT;
    sizes[RTGSFIT_SIZE_N_REG] = N_REG;
    sizes[RTGSFIT_SIZE_N_FIT_TARGET] = N_FIT_TARGET;
    sizes[RTGSFIT_SIZE_N_PLASMA] = N_PLASMA;
    sizes[RTGSFIT_SIZE_N_PASSIVE] = N_PASSIVE;
    sizes[RTGSFIT_SIZE_N_FIT_DOF] = N_FIT_DOF;
    sizes[RTGSFIT_SIZE_N_XPT_MAX] = N_XPT_MAX;
    return RTGSFIT_OK;
}

/*
 * Build the column-major fit_matrix [N_FIT_TARGET x N_FIT_DOF]. Column k is the weighted response of
 * every row of fit_target to unit fit_dof[k]:
 *   plasma columns   basis FIT_A_GRID (FIT_A_GRID already contains FIT_WEIGHT and d_r d_z)
 *   passive columns  FIT_A_PASSIVE, fixed at preshot
 */
static void build_fit_matrix(void)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_PLASMA, N_FIT_TARGET, N_GRID, 1.0,
                work.basis, N_GRID, FIT_A_GRID, N_FIT_TARGET, 0.0, work.fit_matrix, N_FIT_TARGET);
    if (N_PASSIVE > 0)
    {
        memcpy(&work.fit_matrix[N_PLASMA * N_FIT_TARGET], FIT_A_PASSIVE,
               sizeof(double) * (size_t)N_PASSIVE * N_FIT_TARGET);
    }
}

/* j_phi = sum_k fit_dof[k] basis_k, together with its integral and centroid. */
static int32_t compute_plasma_current(const double fit_dof[],
                                      const int32_t n_fit_dof,
                                      double j_phi[],
                                      const int32_t n_j_phi,
                                      double *ip,
                                      double *current_centroid_r,
                                      double *current_centroid_z)
{
    RTGSFIT_CHECK_SIZE(n_fit_dof, N_FIT_DOF);
    RTGSFIT_CHECK_SIZE(n_j_phi, N_GRID);
    RTGSFIT_CHECK_ARRAY(fit_dof, n_fit_dof);
    RTGSFIT_CHECK_ARRAY(j_phi, n_j_phi);
    cblas_dgemv(CblasRowMajor, CblasTrans, N_PLASMA, N_GRID, 1.0, work.basis, N_GRID, fit_dof, 1,
                0.0, j_phi, 1);

    double j_phi_sum = 0.0;
    double j_phi_r_sum = 0.0;
    double j_phi_z_sum = 0.0;
    for (int32_t i = 0; i < N_GRID; i++)
    {
        j_phi_sum += j_phi[i];
        j_phi_r_sum += j_phi[i] * GRID_R_FLAT[i];
        j_phi_z_sum += j_phi[i] * GRID_Z_FLAT[i];
    }
    *ip = j_phi_sum * GRID_D_R * GRID_D_Z;
    /* The centroid is meaningless when there is (almost) no plasma current. */
    if (*ip > SETTINGS_IP_MIN)
    {
        *current_centroid_r = j_phi_r_sum / j_phi_sum;
        *current_centroid_z = j_phi_z_sum / j_phi_sum;
    }
    else
    {
        *current_centroid_r = 0.0;
        *current_centroid_z = 0.0;
    }
    return RTGSFIT_OK;
}

/* fit_prediction = fit_matrix fit_dof, the chi squared values, and the reconstructed constraints. */
static int32_t compute_fit_residual(const double fit_dof[],
                                    const int32_t n_fit_dof,
                                    double constraint_reconstructed[],
                                    const int32_t n_constraint_reconstructed,
                                    double *chi_sq_magnetic,
                                    double *chi_sq_total)
{
    RTGSFIT_CHECK_SIZE(n_fit_dof, N_FIT_DOF);
    RTGSFIT_CHECK_SIZE(n_constraint_reconstructed, N_CONSTRAINT);
    RTGSFIT_CHECK_ARRAY(fit_dof, n_fit_dof);
    RTGSFIT_CHECK_ARRAY(constraint_reconstructed, n_constraint_reconstructed);
    cblas_dgemv(CblasRowMajor, CblasTrans, N_FIT_DOF, N_FIT_TARGET, 1.0, work.fit_matrix,
                N_FIT_TARGET, fit_dof, 1, 0.0, work.fit_prediction, 1);

    double chi_sq = 0.0;
    for (int32_t i = 0; i < N_CONSTRAINT; i++)
    {
        const double residual = work.fit_target[i] - work.fit_prediction[i];
        chi_sq += residual * residual;
    }
    *chi_sq_magnetic = chi_sq;
    for (int32_t i = N_CONSTRAINT; i < N_FIT_TARGET; i++)
    {
        const double residual = work.fit_target[i] - work.fit_prediction[i];
        chi_sq += residual * residual;
    }
    *chi_sq_total = chi_sq;

    for (int32_t i = 0; i < N_CONSTRAINT; i++)
    {
        constraint_reconstructed[i] =
            work.constraint_coil[i] + work.fit_prediction[i] / FIT_WEIGHT[i];
    }
    return RTGSFIT_OK;
}

/* psi = psi_coil + (psi_plasma + psi_passive) on the grid. */
static int32_t compute_total_psi(const double coil_current[],
                                 const int32_t n_coil_current,
                                 const double x_passive[],
                                 const int32_t n_x_passive,
                                 double psi[],
                                 const int32_t n_psi)
{
    RTGSFIT_CHECK_SIZE(n_coil_current, N_COIL);
    RTGSFIT_CHECK_SIZE(n_x_passive, N_PASSIVE);
    RTGSFIT_CHECK_SIZE(n_psi, N_GRID);
    RTGSFIT_CHECK_ARRAY(coil_current, n_coil_current);
    RTGSFIT_CHECK_ARRAY(x_passive, n_x_passive);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_COIL, 1.0, GRID_G_COIL, N_COIL, coil_current,
                1, 0.0, psi, 1);
    if (N_PASSIVE > 0)
    {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_PASSIVE, 1.0, GRID_G_PASSIVE, N_PASSIVE,
                    x_passive, 1, 0.0, work.psi_passive, 1);
        for (int32_t i = 0; i < N_GRID; i++)
        {
            psi[i] += work.psi_plasma[i] + work.psi_passive[i];
        }
    }
    else
    {
        for (int32_t i = 0; i < N_GRID; i++)
        {
            psi[i] += work.psi_plasma[i];
        }
    }
    return RTGSFIT_OK;
}

int32_t rtgsfit(const double meas[],
                const int32_t n_meas,
                const double coil_current[],
                const int32_t n_coil_current,
                double psi_norm[],
                const int32_t n_psi_norm,
                int32_t mask[],
                const int32_t n_mask,
                double psi[],
                const int32_t n_psi,
                double j_phi[],
                const int32_t n_j_phi,
                double fit_dof[],
                const int32_t n_fit_dof,
                double constraint_measured[],
                const int32_t n_constraint_measured,
                double constraint_reconstructed[],
                const int32_t n_constraint_reconstructed,
                double xpt_r[],
                const int32_t n_xpt_r,
                double xpt_z[],
                const int32_t n_xpt_z,
                double xpt_psi[],
                const int32_t n_xpt_psi,
                int32_t *xpt_n,
                int32_t *xpt_diverted,
                double *ip,
                double *axis_psi,
                double *axis_r,
                double *axis_z,
                double *boundary_psi,
                double *current_centroid_r,
                double *current_centroid_z,
                double *chi_sq_magnetic,
                double *chi_sq_total,
                int32_t *lapack_info)
{
    enable_flush_to_zero();
    const uint64_t start_total_ns = timing_now_ns();
    uint64_t start_ns = start_total_ns;

    RTGSFIT_CHECK_SIZE(n_meas, N_MEAS);
    RTGSFIT_CHECK_SIZE(n_coil_current, N_COIL);
    RTGSFIT_CHECK_SIZE(n_psi_norm, N_GRID);
    RTGSFIT_CHECK_SIZE(n_mask, N_GRID);
    RTGSFIT_CHECK_SIZE(n_psi, N_GRID);
    RTGSFIT_CHECK_SIZE(n_j_phi, N_GRID);
    RTGSFIT_CHECK_SIZE(n_fit_dof, N_FIT_DOF);
    RTGSFIT_CHECK_SIZE(n_constraint_measured, N_CONSTRAINT);
    RTGSFIT_CHECK_SIZE(n_constraint_reconstructed, N_CONSTRAINT);
    RTGSFIT_CHECK_SIZE(n_xpt_r, N_XPT_MAX);
    RTGSFIT_CHECK_SIZE(n_xpt_z, N_XPT_MAX);
    RTGSFIT_CHECK_SIZE(n_xpt_psi, N_XPT_MAX);
    RTGSFIT_CHECK_ARRAY(meas, n_meas);
    RTGSFIT_CHECK_ARRAY(coil_current, n_coil_current);
    RTGSFIT_CHECK_ARRAY(psi_norm, n_psi_norm);
    RTGSFIT_CHECK_ARRAY(mask, n_mask);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(j_phi, n_j_phi);
    RTGSFIT_CHECK_ARRAY(fit_dof, n_fit_dof);
    RTGSFIT_CHECK_ARRAY(constraint_measured, n_constraint_measured);
    RTGSFIT_CHECK_ARRAY(constraint_reconstructed, n_constraint_reconstructed);
    RTGSFIT_CHECK_ARRAY(xpt_r, n_xpt_r);
    RTGSFIT_CHECK_ARRAY(xpt_z, n_xpt_z);
    RTGSFIT_CHECK_ARRAY(xpt_psi, n_xpt_psi);
    if (!xpt_n || !xpt_diverted || !ip || !axis_psi || !axis_r || !axis_z || !boundary_psi ||
        !current_centroid_r || !current_centroid_z || !chi_sq_magnetic || !chi_sq_total ||
        !lapack_info)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }

    const struct rtgsfit_grid grid = preshot_grid();
    const struct poisson_operator poisson = preshot_poisson_operator();
    int32_t status = RTGSFIT_OK;

    *xpt_n = 0;
    *xpt_diverted = 0;
    *axis_psi = 0.0;
    *axis_r = 0.0;
    *axis_z = 0.0;
    *boundary_psi = 0.0;
    *lapack_info = 0;

    /* 1. Constraints and fit_target. */
    status = measurement_processing_prepare(
        meas, n_meas, coil_current, n_coil_current, CONSTRAINT_REPLACEMENT, N_CONSTRAINT * N_MEAS,
        CONSTRAINT_G_COIL, N_CONSTRAINT * N_COIL, FIT_WEIGHT, N_FIT_TARGET, constraint_measured,
        n_constraint_measured, work.constraint_coil, N_CONSTRAINT, work.fit_target, N_FIT_TARGET);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_CONSTRAINTS, start_ns);

    /* 2. Current basis functions, from the previous psi_norm and mask. */
    status = current_basis_compute(&grid, psi_norm, N_GRID, mask, N_GRID, GRID_MASK_LIMITER, N_GRID,
                                   GRID_R_FLAT, N_GRID, GRID_INV_R, N_GRID, work.basis,
                                   N_PLASMA * N_GRID);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_BASIS, start_ns);

    /* 3. fit_matrix (fit_target is already weighted). */
    build_fit_matrix();
    start_ns = timing_accumulate(RTGSFIT_TIMER_FIT_MATRIX, start_ns);

    /* 4. Least squares: fit_matrix fit_dof = fit_target. */
    memcpy(work.fit_matrix_lapack, work.fit_matrix, sizeof(work.fit_matrix));
    memcpy(work.fit_target_lapack, work.fit_target, sizeof(work.fit_target));
    int32_t rank = 0;
    status = least_squares_solve(N_FIT_TARGET, N_FIT_DOF, work.fit_matrix_lapack,
                                 N_FIT_DOF * N_FIT_TARGET, work.fit_target_lapack, N_FIT_TARGET,
                                 work.singular_values, N_FIT_DOF, work.least_squares_work,
                                 LEAST_SQUARES_WORK_SIZE, -1.0, &rank, lapack_info);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    memcpy(fit_dof, work.fit_target_lapack, sizeof(double) * N_FIT_DOF);
    start_ns = timing_accumulate(RTGSFIT_TIMER_LEAST_SQUARES, start_ns);

    /* 5. Plasma current. */
    status = compute_plasma_current(fit_dof, n_fit_dof, j_phi, n_j_phi, ip, current_centroid_r,
                                    current_centroid_z);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_CURRENT, start_ns);

    /* 6. Residual. */
    status = compute_fit_residual(fit_dof, n_fit_dof, constraint_reconstructed,
                                  n_constraint_reconstructed, chi_sq_magnetic, chi_sq_total);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_RESIDUAL, start_ns);

    /* 7. Plasma flux; RHS_SCALE includes 2 pi for total flux and the d_z^2 stencil scaling. */
    for (int32_t i = 0; i < N_GRID; i++)
    {
        work.poisson_rhs[i] = -POISSON_RHS_SCALE[i] * j_phi[i];
    }
    status = poisson_solve_free_boundary(&poisson, work.poisson_rhs, N_GRID, work.psi_plasma,
                                         N_GRID, work.poisson_work, POISSON_WORK_SIZE);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_POISSON, start_ns);

    /* 8. Total psi. */
    status =
        compute_total_psi(coil_current, n_coil_current, &fit_dof[N_PLASMA], N_PASSIVE, psi, n_psi);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_EXTERNAL_PSI, start_ns);

    /* 9. Magnetic axis and x-points. */
    int32_t axis_found = 0;
    status = magnetic_nulls_find(&grid, psi, N_GRID, GRID_MASK_LIMITER, N_GRID, SETTINGS_TOLERANCE,
                                 &axis_found, axis_r, axis_z, axis_psi, xpt_r, N_XPT_MAX, xpt_z,
                                 N_XPT_MAX, xpt_psi, N_XPT_MAX, xpt_n);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    if (!axis_found)
    {
        status = RTGSFIT_ERR_NO_AXIS;
        goto finish;
    }
    status = magnetic_nulls_filter_xpts(*axis_r, *axis_z, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX,
                                        xpt_psi, N_XPT_MAX, xpt_n, work.xpt_keep, N_XPT_MAX);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    status =
        magnetic_nulls_sort_xpts(xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX, xpt_psi, N_XPT_MAX, *xpt_n);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_NULLS, start_ns);

    /* 10. Boundary: limited or diverted. */
    double limiter_psi = 0.0;
    status = plasma_boundary_limiter_psi(psi, N_GRID, LIMITER_R, N_LIMITER, LIMITER_Z, N_LIMITER,
                                         LIMITER_GRID_INDEX, N_LIMITER * PLASMA_BOUNDARY_N_INTERP,
                                         LIMITER_GRID_WEIGHT, N_LIMITER * PLASMA_BOUNDARY_N_INTERP,
                                         *axis_r, *axis_z, xpt_r, N_XPT_MAX, xpt_z, N_XPT_MAX,
                                         *xpt_n, &limiter_psi);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    status = plasma_boundary_select(*axis_psi, limiter_psi, xpt_psi, N_XPT_MAX, *xpt_n,
                                    SETTINGS_XPT_FRACTION, SETTINGS_TOLERANCE, boundary_psi,
                                    xpt_diverted);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_BOUNDARY, start_ns);

    /* 11. Plasma core mask. */
    status = plasma_boundary_flood_fill(
        &grid, psi, N_GRID, GRID_MASK_LIMITER, N_GRID, *boundary_psi, *axis_r, *axis_z, xpt_r,
        N_XPT_MAX, xpt_z, N_XPT_MAX, *xpt_n, mask, N_GRID, work.flood_fill_queue, N_GRID);
    if (status != RTGSFIT_OK)
    {
        goto finish;
    }
    start_ns = timing_accumulate(RTGSFIT_TIMER_FLOOD_FILL, start_ns);

    /* 12. Normalised flux, the starting point of the next iteration. */
    status = plasma_boundary_normalise(psi, N_GRID, mask, N_GRID, GRID_MASK_LIMITER, N_GRID,
                                       *axis_psi, *boundary_psi, psi_norm, N_GRID);
    timing_accumulate(RTGSFIT_TIMER_NORMALISE, start_ns);

finish:
    timing_accumulate(RTGSFIT_TIMER_TOTAL, start_total_ns);
    return status;
}
