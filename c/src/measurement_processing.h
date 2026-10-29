#ifndef RTGSFIT_MEASUREMENT_PROCESSING_H
#define RTGSFIT_MEASUREMENT_PROCESSING_H

#include <stdint.h>

#include "rtgsfit.h"

/*
 * Apply GSFit-style sensor replacement, subtract the known coil response, and weight the
 * constraints. Matrices are row-major. fit_target has n_constraint rows followed by zero-valued
 * regularisation rows; weight has n_fit_target entries. All lengths count elements, not bytes.
 *
 * This module has no preshot dependency, so preparation can later move to the PCS caller.
 * Returns RTGSFIT_ERR_ARGUMENT_SIZE before writing outputs if any buffer is invalid.
 */
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
                                       const int32_t n_fit_target);

#endif /* RTGSFIT_MEASUREMENT_PROCESSING_H */
