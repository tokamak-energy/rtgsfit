#include "measurement_processing.h"

#include <cblas.h>

#include "argument_check.h"

int32_t measurement_processing_prepare(const double meas[],
                                       const int32_t n_meas,
                                       const double coil_current[],
                                       const int32_t n_coil_current,
                                       const double replacement[],
                                       const int32_t n_replacement,
                                       const double g_coil[],
                                       const int32_t n_g_coil,
                                       const double weight[],
                                       const int32_t n_weight,
                                       double constraint_measured[],
                                       const int32_t n_constraint_measured,
                                       double constraint_coil[],
                                       const int32_t n_constraint_coil,
                                       double fit_target[],
                                       const int32_t n_fit_target)
{
    const int32_t n_constraint = n_constraint_measured;
    RTGSFIT_CHECK_MIN_SIZE(n_meas, 1);
    RTGSFIT_CHECK_MIN_SIZE(n_coil_current, 1);
    RTGSFIT_CHECK_MIN_SIZE(n_constraint, 1);
    RTGSFIT_CHECK_SIZE(n_replacement, (int64_t)n_constraint * n_meas);
    RTGSFIT_CHECK_SIZE(n_g_coil, (int64_t)n_constraint * n_coil_current);
    RTGSFIT_CHECK_SIZE(n_constraint_coil, n_constraint);
    RTGSFIT_CHECK_MIN_SIZE(n_fit_target, n_constraint);
    RTGSFIT_CHECK_SIZE(n_weight, n_fit_target);
    RTGSFIT_CHECK_ARRAY(meas, n_meas);
    RTGSFIT_CHECK_ARRAY(coil_current, n_coil_current);
    RTGSFIT_CHECK_ARRAY(replacement, n_replacement);
    RTGSFIT_CHECK_ARRAY(g_coil, n_g_coil);
    RTGSFIT_CHECK_ARRAY(weight, n_weight);
    RTGSFIT_CHECK_ARRAY(constraint_measured, n_constraint_measured);
    RTGSFIT_CHECK_ARRAY(constraint_coil, n_constraint_coil);
    RTGSFIT_CHECK_ARRAY(fit_target, n_fit_target);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_constraint, n_meas, 1.0, replacement, n_meas, meas,
                1, 0.0, constraint_measured, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_constraint, n_coil_current, 1.0, g_coil,
                n_coil_current, coil_current, 1, 0.0, constraint_coil, 1);
    for (int32_t i = 0; i < n_constraint; i++)
    {
        fit_target[i] = (constraint_measured[i] - constraint_coil[i]) * weight[i];
    }
    for (int32_t i = n_constraint; i < n_fit_target; i++)
    {
        fit_target[i] = 0.0;
    }
    return RTGSFIT_OK;
}
