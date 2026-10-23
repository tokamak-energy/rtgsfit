/*
 * Magnetic nulls: points where grad(psi) = 0. The o-point with the largest psi is the magnetic
 * axis; the saddle points are x-points.
 */
#ifndef RTGSFIT_MAGNETIC_NULLS_H
#define RTGSFIT_MAGNETIC_NULLS_H

#include <stdint.h>

#include "grid.h"

/*
 * Locate the magnetic axis and the x-points of psi.
 *
 * Each grid cell is examined from each of its four corners that lie inside the limiter. psi is
 * expanded to second order about the corner, and the stationary point of that expansion is
 * accepted if it lies inside the cell. If several corners give a null, the one closest to its
 * corner is used. A positive Hessian determinant means an o-point; a negative one, an x-point.
 * The Hessian determinant must exceed `tolerance` in magnitude.
 *
 * psi[n_psi], mask_limiter[n_mask_limiter] are grid arrays.
 * xpt_r, xpt_z, xpt_psi all have length n_xpt_max; the first *xpt_n entries are written.
 * *axis_found is 1 if an o-point was found, and then axis_r, axis_z, axis_psi are set.
 *
 * Returns RTGSFIT_ERR_TOO_MANY_XPTS, and stops the search, as soon as n_xpt_max x-points have
 * been found.
 */
int32_t magnetic_nulls_find(const struct rtgsfit_grid *grid,
                            const double psi[],
                            int32_t n_psi,
                            const int32_t mask_limiter[],
                            int32_t n_mask_limiter,
                            double tolerance,
                            int32_t *axis_found,
                            double *axis_r,
                            double *axis_z,
                            double *axis_psi,
                            double xpt_r[],
                            int32_t n_xpt_r,
                            double xpt_z[],
                            int32_t n_xpt_z,
                            double xpt_psi[],
                            int32_t n_xpt_psi,
                            int32_t *xpt_n);

/*
 * 1 if the point (r, z) is on the same side as the magnetic axis of every x-point, where "side"
 * is the half-plane through the x-point perpendicular to the line joining it to the axis.
 * xpt_r and xpt_z have length n_xpt; only the first xpt_n entries are used.
 */
int32_t magnetic_nulls_is_core_side(double r,
                                    double z,
                                    double axis_r,
                                    double axis_z,
                                    const double xpt_r[],
                                    const double xpt_z[],
                                    int32_t xpt_n);

/*
 * Remove x-points that are hidden behind another x-point, as seen from the magnetic axis:
 * x-point j is removed if it is not on the core side of some other x-point i. The remaining
 * x-points keep their order and *xpt_n is reduced.
 *
 * keep[n_keep] is integer workspace with n_keep >= *xpt_n.
 */
int32_t magnetic_nulls_filter_xpts(double axis_r,
                                   double axis_z,
                                   double xpt_r[],
                                   int32_t n_xpt_r,
                                   double xpt_z[],
                                   int32_t n_xpt_z,
                                   double xpt_psi[],
                                   int32_t n_xpt_psi,
                                   int32_t *xpt_n,
                                   int32_t keep[],
                                   int32_t n_keep);

/* Sort the first xpt_n x-points by descending psi. */
int32_t magnetic_nulls_sort_xpts(double xpt_r[],
                                 int32_t n_xpt_r,
                                 double xpt_z[],
                                 int32_t n_xpt_z,
                                 double xpt_psi[],
                                 int32_t n_xpt_psi,
                                 int32_t xpt_n);

#endif /* RTGSFIT_MAGNETIC_NULLS_H */
