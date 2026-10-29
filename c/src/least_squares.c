#include "least_squares.h"

#include <lapacke.h>
#include <math.h>

#include "argument_check.h"

int32_t rtgsfit_status_from_lapack_info(const int32_t info)
{
    if (info == 0)
    {
        return RTGSFIT_OK;
    }
    if (info == LAPACK_WORK_MEMORY_ERROR || info == LAPACK_TRANSPOSE_MEMORY_ERROR)
    {
        return RTGSFIT_ERR_LAPACK | RTGSFIT_ERR_FIT_RESOURCE_FAILURE;
    }
    return RTGSFIT_ERR_LAPACK |
           (info < 0 ? RTGSFIT_ERR_FIT_INVALID_ARGUMENT : RTGSFIT_ERR_FIT_NO_CONVERGENCE);
}

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
                            int32_t *info)
{
    if (n_fit_target < 1 || n_fit_dof < 1 || n_fit_target < n_fit_dof)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    RTGSFIT_CHECK_SIZE(n_fit_matrix, (int64_t)n_fit_target * n_fit_dof);
    RTGSFIT_CHECK_SIZE(n_fit_target_array, n_fit_target);
    RTGSFIT_CHECK_SIZE(n_singular_values, n_fit_dof);
    RTGSFIT_CHECK_MIN_SIZE(n_work, least_squares_min_work_size(n_fit_target, n_fit_dof));
    RTGSFIT_CHECK_ARRAY(fit_matrix, n_fit_matrix);
    RTGSFIT_CHECK_ARRAY(fit_target, n_fit_target_array);
    RTGSFIT_CHECK_ARRAY(singular_values, n_singular_values);
    RTGSFIT_CHECK_ARRAY(work, n_work);
    if (!rank || !info)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    *rank = 0;
    *info = 0;
    for (int32_t i_value = 0; i_value < n_fit_matrix; i_value++)
    {
        if (!isfinite(fit_matrix[i_value]))
        {
            return RTGSFIT_ERR_FIT_NONFINITE_INPUT;
        }
    }
    for (int32_t i_value = 0; i_value < n_fit_target_array; i_value++)
    {
        if (!isfinite(fit_target[i_value]))
        {
            return RTGSFIT_ERR_FIT_NONFINITE_INPUT;
        }
    }

    lapack_int lapack_rank = 0;
    lapack_int lapack_info = LAPACKE_dgelss_work(LAPACK_COL_MAJOR, n_fit_target, n_fit_dof, 1,
                                                 fit_matrix, n_fit_target, fit_target, n_fit_target,
                                                 singular_values, rcond, &lapack_rank, work, -1);
    if (lapack_info == 0)
    {
        const double optimal = work[0];
        const lapack_int n_use =
            (optimal > 0.0 && optimal <= (double)n_work) ? (lapack_int)optimal : (lapack_int)n_work;
        lapack_info = LAPACKE_dgelss_work(LAPACK_COL_MAJOR, n_fit_target, n_fit_dof, 1, fit_matrix,
                                          n_fit_target, fit_target, n_fit_target, singular_values,
                                          rcond, &lapack_rank, work, n_use);
    }
    *rank = (int32_t)lapack_rank;
    *info = (int32_t)lapack_info;
    return rtgsfit_status_from_lapack_info(*info);
}
