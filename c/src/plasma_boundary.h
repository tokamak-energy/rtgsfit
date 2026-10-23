/*
 * Plasma boundary: the limiter psi, the choice between a limited and a diverted boundary, the
 * plasma core mask, and the normalised flux.
 */
#ifndef RTGSFIT_PLASMA_BOUNDARY_H
#define RTGSFIT_PLASMA_BOUNDARY_H

#include <stdint.h>

#include "grid.h"

/* Number of grid points used to interpolate psi onto a limiter point (bilinear). */
enum { PLASMA_BOUNDARY_N_INTERP = 4 };

/*
 * The largest psi over the limiter points that are on the core side of every x-point.
 *
 * limiter_grid_index and limiter_grid_weight are row-major [n_limiter][PLASMA_BOUNDARY_N_INTERP].
 * Row i gives the grid indices of the four grid points around limiter point i, and their
 * bilinear interpolation weights.
 * If no limiter point qualifies, *limiter_psi = -DBL_MAX.
 */
int32_t plasma_boundary_limiter_psi(const double psi[],
                                    int32_t n_psi,
                                    const double limiter_r[],
                                    int32_t n_limiter_r,
                                    const double limiter_z[],
                                    int32_t n_limiter_z,
                                    const int32_t limiter_grid_index[],
                                    int32_t n_index,
                                    const double limiter_grid_weight[],
                                    int32_t n_weight,
                                    double axis_r,
                                    double axis_z,
                                    const double xpt_r[],
                                    int32_t n_xpt_r,
                                    const double xpt_z[],
                                    int32_t n_xpt_z,
                                    int32_t xpt_n,
                                    double *limiter_psi);

/*
 * Choose psi on the plasma boundary. The boundary is limited (boundary_psi = limiter_psi)
 * unless the highest x-point flux, moved a fraction (1 - xpt_fraction) towards the axis, is
 * larger, in which case the plasma is diverted.
 *
 * xpt_psi[n_xpt_psi] must already be sorted by descending psi (only the first xpt_n are used).
 *
 * Returns RTGSFIT_ERR_AXIS_EQUALS_BOUNDARY if |boundary_psi - axis_psi| < tolerance, and
 * RTGSFIT_ERR_BOUNDARY_ABOVE_AXIS if boundary_psi > axis_psi.
 */
int32_t plasma_boundary_select(double axis_psi,
                               double limiter_psi,
                               const double xpt_psi[],
                               int32_t n_xpt_psi,
                               int32_t xpt_n,
                               double xpt_fraction,
                               double tolerance,
                               double *boundary_psi,
                               int32_t *xpt_diverted);

/*
 * Flood fill the plasma core, starting from the grid point nearest the magnetic axis. A grid
 * point is in the core if psi > boundary_psi, it is inside the limiter, it is on the core side
 * of every x-point, and it is 4-connected to the starting point.
 *
 * mask[n_mask] is only changed on success. queue[n_queue] is workspace with n_queue == n_grid.
 * Returns RTGSFIT_ERR_AXIS_OUTSIDE_CORE if the starting point is not in the core.
 */
int32_t plasma_boundary_flood_fill(const struct rtgsfit_grid *grid,
                                   const double psi[],
                                   int32_t n_psi,
                                   const int32_t mask_limiter[],
                                   int32_t n_mask_limiter,
                                   double boundary_psi,
                                   double axis_r,
                                   double axis_z,
                                   const double xpt_r[],
                                   int32_t n_xpt_r,
                                   const double xpt_z[],
                                   int32_t n_xpt_z,
                                   int32_t xpt_n,
                                   int32_t mask[],
                                   int32_t n_mask,
                                   int32_t queue[],
                                   int32_t n_queue);

/*
 * psi_norm = (psi - axis_psi) / (boundary_psi - axis_psi) inside the plasma core
 * (mask && mask_limiter), and 1 elsewhere. All arrays are grid arrays of length n.
 */
int32_t plasma_boundary_normalise(const double psi[],
                                  int32_t n_psi,
                                  const int32_t mask[],
                                  int32_t n_mask,
                                  const int32_t mask_limiter[],
                                  int32_t n_mask_limiter,
                                  double axis_psi,
                                  double boundary_psi,
                                  double psi_norm[],
                                  int32_t n_psi_norm);

#endif /* RTGSFIT_PLASMA_BOUNDARY_H */
