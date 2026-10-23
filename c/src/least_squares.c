#include "least_squares.h"

#include <lapacke.h>

#include "argument_check.h"

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
                            int32_t *info)
{
    if (n_b_vec < 1 || n_x_vec < 1 || n_b_vec < n_x_vec) {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    RTGSFIT_CHECK_SIZE(n_a_matrix, n_b_vec * n_x_vec);
    RTGSFIT_CHECK_SIZE(n_b_vec_array, n_b_vec);
    RTGSFIT_CHECK_SIZE(n_singular_values, n_x_vec);
    RTGSFIT_CHECK_MIN_SIZE(n_work, least_squares_min_work_size(n_b_vec, n_x_vec));

    lapack_int lapack_rank = 0;
    lapack_int lapack_info =
        LAPACKE_dgelss_work(LAPACK_COL_MAJOR, n_b_vec, n_x_vec, 1, a_matrix, n_b_vec, b_vec,
                            n_b_vec, singular_values, rcond, &lapack_rank, work, -1);
    if (lapack_info == 0) {
        const double optimal = work[0];
        const lapack_int n_use =
            (optimal > 0.0 && optimal <= (double)n_work) ? (lapack_int)optimal : (lapack_int)n_work;
        lapack_info =
            LAPACKE_dgelss_work(LAPACK_COL_MAJOR, n_b_vec, n_x_vec, 1, a_matrix, n_b_vec, b_vec,
                                n_b_vec, singular_values, rcond, &lapack_rank, work, n_use);
    }
    *rank = (int32_t)lapack_rank;
    *info = (int32_t)lapack_info;
    return RTGSFIT_OK;
}
