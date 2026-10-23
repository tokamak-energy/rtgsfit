/*
 * Linear least squares:  minimise |a_matrix x_vec - b_vec|_2  with LAPACK dgelss (SVD).
 *
 * The naming follows the LAPACK documentation and issue #68:
 *   a_matrix [n_b_vec x n_x_vec]  column-major, so column k is contiguous
 *   b_vec    [n_b_vec]
 *   x_vec    [n_x_vec]
 */
#ifndef RTGSFIT_LEAST_SQUARES_H
#define RTGSFIT_LEAST_SQUARES_H

#include <stdint.h>

/* Smallest workspace dgelss accepts; a larger one lets LAPACK use blocked algorithms. */
static inline int32_t least_squares_min_work_size(int32_t n_b_vec, int32_t n_x_vec)
{
    const int32_t n_min = n_b_vec < n_x_vec ? n_b_vec : n_x_vec;
    const int32_t n_max = n_b_vec > n_x_vec ? n_b_vec : n_x_vec;
    const int32_t two_n_min = 2 * n_min;
    return 3 * n_min + (two_n_min > n_max ? two_n_min : n_max);
}

/*
 * Solve the least-squares problem in place, without allocating memory.
 *
 * a_matrix[n_a_matrix] (n_a_matrix == n_b_vec * n_x_vec) is destroyed.
 * b_vec[n_b_vec] is overwritten: its first n_x_vec entries become x_vec.
 * singular_values[n_singular_values] receives min(n_b_vec, n_x_vec) singular values.
 * work[n_work] is scratch space with n_work >= least_squares_min_work_size(). If LAPACK's optimal
 * size is not larger than n_work then the optimal size is used, and the result is identical to
 * LAPACKE_dgelss().
 * rcond: singular values below rcond * s_max are treated as zero; rcond < 0 means machine
 * precision.
 *
 * Returns RTGSFIT_OK or RTGSFIT_ERR_ARGUMENT_SIZE. `info` is LAPACK's info, and is 0 on success.
 */
int32_t least_squares_solve(int32_t n_b_vec,
                            int32_t n_x_vec,
                            double a_matrix[],
                            int32_t n_a_matrix,
                            double b_vec[],
                            int32_t n_b_vec_array,
                            double singular_values[],
                            int32_t n_singular_values,
                            double work[],
                            int32_t n_work,
                            double rcond,
                            int32_t *rank,
                            int32_t *info);

#endif /* RTGSFIT_LEAST_SQUARES_H */
