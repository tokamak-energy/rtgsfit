/*
 * Finite-difference derivatives on the (R, Z) grid.
 */
#ifndef RTGSFIT_FINITE_DIFFERENCE_H
#define RTGSFIT_FINITE_DIFFERENCE_H

#include <stdint.h>

#include "grid.h"

/*
 * df_dz = d(f)/dZ, using second-order central differences in the interior and second-order
 * one-sided differences on the bottom and top rows.
 *
 * f[n_f] and df_dz[n_df_dz] are grid arrays, so n_f == n_df_dz == n_r * n_z.
 */
int32_t finite_difference_d_dz(const struct rtgsfit_grid *grid,
                               const double f[],
                               const int32_t n_f,
                               double df_dz[],
                               const int32_t n_df_dz);

/*
 * The scaled inward normal derivative of f on the edge of the grid. It is used by the
 * von Hagenow free-boundary method (Moret et al. 2015, section 4.2). f must be zero on the edge.
 *
 * The edge is traversed clockwise from the bottom-left corner (the "edge order", see the manual):
 *   left   (R = R_min, Z increasing, both corners)
 *   top    (Z = Z_max, R increasing, no corners)
 *   right  (R = R_max, Z decreasing, both corners)
 *   bottom (Z = Z_min, R decreasing, no corners)
 *
 * Each value is the one-sided second-order difference (2 f_1 - f_2 / 2), multiplied by
 * d_z / d_r on the left and right sides and by d_r / d_z on the top and bottom sides.
 *
 * f[n_f] is a grid array; df_dn[n_df_dn] has n_df_dn == 2 n_r + 2 n_z - 4.
 */
int32_t finite_difference_edge_d_dn(const struct rtgsfit_grid *grid,
                                    const double f[],
                                    const int32_t n_f,
                                    double df_dn[],
                                    const int32_t n_df_dn);

#endif /* RTGSFIT_FINITE_DIFFERENCE_H */
