/*
 * Linear least squares:  minimise |fit_matrix fit_dof - fit_target|_2  with LAPACK dgelss (SVD).
 *
 * Descriptive fit names avoid confusing the target with a magnetic field:
 *   fit_matrix [n_fit_target x n_fit_dof]  column-major, so column k is contiguous
 *   fit_target [n_fit_target]  weighted, coil-subtracted constraints; regularisation zeros
 *   fit_dof    [n_fit_dof]     unknown DOF values in physical units, not weighted
 *
 * Weights scale the target rows and matching matrix rows, not the unknowns. LAPACK overwrites
 * its target buffer with the solution; that buffer no longer holds a target after the solve.
 */
#ifndef RTGSFIT_LEAST_SQUARES_H
#define RTGSFIT_LEAST_SQUARES_H

#include <stdint.h>

/* Smallest workspace dgelss accepts; a larger one lets LAPACK use blocked algorithms. */
static inline int64_t least_squares_min_work_size(const int32_t n_fit_target,
                                                  const int32_t n_fit_dof)
{
    const int32_t n_min = n_fit_target < n_fit_dof ? n_fit_target : n_fit_dof;
    const int32_t n_max = n_fit_target > n_fit_dof ? n_fit_target : n_fit_dof;
    const int64_t two_n_min = 2 * (int64_t)n_min;
    return 3 * (int64_t)n_min + (two_n_min > n_max ? two_n_min : n_max);
}

/*
 * Solve the least-squares problem in place, without allocating memory.
 *
 * fit_matrix[n_fit_matrix] (n_fit_matrix == n_fit_target * n_fit_dof) is destroyed.
 * fit_target[n_fit_target] is overwritten: its first n_fit_dof entries become fit_dof.
 * singular_values[n_singular_values] receives min(n_fit_target, n_fit_dof) singular values.
 * work[n_work] is scratch space with n_work >= least_squares_min_work_size(). If LAPACK's optimal
 * size is not larger than n_work then the optimal size is used, and the result is identical to
 * LAPACKE_dgelss().
 * rcond: singular values below rcond * s_max are treated as zero; rcond < 0 means machine
 * precision.
 *
 * Returns RTGSFIT_OK, RTGSFIT_ERR_ARGUMENT_SIZE, RTGSFIT_ERR_FIT_NONFINITE_INPUT, or
 * RTGSFIT_ERR_LAPACK OR-ed with the corresponding fit-failure category.
 * `info` preserves LAPACK's exact info; 0 on success or pre-LAPACK nonfinite-input rejection.
 */
int32_t least_squares_solve(const int32_t n_fit_target,
                            const int32_t n_fit_dof,
                            double fit_matrix[],
                            const int32_t n_fit_matrix,
                            double fit_target[],
                            const int32_t n_fit_target_array,
                            double singular_values[],
                            const int32_t n_singular_values,
                            double work[],
                            const int32_t n_work,
                            double rcond,
                            int32_t *rank,
                            int32_t *info);

#endif /* RTGSFIT_LEAST_SQUARES_H */
