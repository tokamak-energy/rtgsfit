/*
 * File: find_plasma.c
 * -------------------
 * Replacement for find_x_points.c from RT-GSFit (work in progress).
 *
 * Purpose:
 *   This file implements functions to locate null points in the magnetic field
 *   which correspond to points where the gradient of the poloidal flux is zero.
 *   The original find_x_points.c calculates the coordinates of the last closed
 *   flux surface (LCFS). In this replacement, we only compute the value of ψ at
 *   the boundary, avoiding unnecessary calculations.
 *
 * Notes:
 *   - We refer to the region inside the LCFS as the "plasma region". While
 *     plasma also exists outside the LCFS, the pressure and temperature are
 *     significantly higher within it. This is therefore the region of primary
 *     interest for equilibrium reconstruction.
 *   - Replacement is ongoing; current implementation may be incomplete.
 *
 * Author: Alex Prokopyszyn
 * Created: 2025-Dec-08
 */
#include "find_plasma.h"

#include <math.h>
#include <stdio.h>

#include "constants.h"

#define N_R_MIN_1 (N_R - 1)
#define N_R_PLS_1 (N_R + 1)
#define N_Z_MIN_1 (N_Z - 1)
#define N_Z_PLS_1 (N_Z + 1)

// find_nulls:
//   Finds null points of the magnetic field which correspond to points where ∇ψ
//   = 0. Also characterises them as either x-points or o-points based on the
//   Hessian matrix.
//
//   The numerical procedure closely follows the scheme described in Section 4.3
//   of Moret et al. (2015) (the LIUQE equilibrium reconstruction paper).
//   However, whereas the LIUQE implementation uses a 6-point interpolation
//   stencil, we employ a 9-point interpolation scheme.
//
// Parameters:
//   flux      - ψ 2D input array explicitly stored as a contiguous block
//               (i.e., 1D array)
//   opt_r     - output array for the o-point R coordinates
//   opt_z     - output array for the o-point Z coordinates
//   opt_flux  - output array for flux at o-points
//   opt_n     - pointer which will hold the number of o-points found
//   xpt_r     - output array for the x-point R coordinates
//   xpt_z     - output array for the x-point Z coordinates
//   xpt_flux  - output array for flux at the x-points
//   xpt_n     - pointer which will hold the number of x-points found
int find_nulls(double *flux, double *opt_r, double *opt_z, double *opt_flux,
               int32_t *opt_n, double *xpt_r, double *xpt_z, double *xpt_flux,
               int32_t *xpt_n) {
  *opt_n = 0;
  *xpt_n = 0;
  for (int32_t i_row = 1; i_row < N_Z_MIN_1; i_row++) {
    int32_t row = i_row * N_R;
    for (int32_t i_col = 1; i_col < N_R_MIN_1; i_col++) {
      int32_t idx = row + i_col;
      if (!MASK_LIM[idx]) {
        continue;
      }
      int32_t idx_rp = idx + 1;
      int32_t idx_rm = idx - 1;
      int32_t idx_zp = idx + N_R;
      int32_t idx_zm = idx - N_R;
      int32_t idx_rp_zp = idx_rp + N_R;
      int32_t idx_rp_zm = idx_rp - N_R;
      int32_t idx_rm_zp = idx_rm + N_R;
      int32_t idx_rm_zm = idx_rm - N_R;
      double a = 0.5 * (flux[idx_rp] - flux[idx_rm]);
      double b = 0.5 * (flux[idx_zp] - flux[idx_zm]);
      double c = flux[idx_rp] - 2.0 * flux[idx] + flux[idx_rm];
      double d = flux[idx_zp] - 2.0 * flux[idx] + flux[idx_zm];
      double e = 0.25 * (flux[idx_rp_zp] - flux[idx_rm_zp] - flux[idx_rp_zm] +
                         flux[idx_rm_zm]);
      double denom = (c * d - e * e);
      if (fabs(denom) < 1e-14)
        continue;
      double inv_denom = 1.0 / denom;
      double dr_norm = (b * e - a * d) * inv_denom; // (dr / ΔR)
      double dz_norm = (a * e - b * c) * inv_denom; // (dz / ΔZ)
      if (fabs(dr_norm) <= 0.5 && fabs(dz_norm) <= 0.5) {
        double null_r = R_VEC[i_col] + dr_norm * DR;
        double null_z = Z_VEC[i_row] + dz_norm * DZ;
        double hess_det = c * d - e * e;
        double flux_at_null = flux[idx] + a * dr_norm + b * dz_norm +
                              0.5 * c * dr_norm * dr_norm +
                              0.5 * d * dz_norm * dz_norm +
                              e * dr_norm * dz_norm;
        if (fabs(hess_det) < 1e-14)
          continue;
        if (hess_det > 0.0) {
          // o-point
          opt_r[*opt_n] = null_r;
          opt_z[*opt_n] = null_z;
          opt_flux[*opt_n] = flux_at_null;
          (*opt_n)++;
          if (*opt_n >= N_XPT_MAX) {
            return 1024;
          }
        } else if (hess_det < 0.0) {
          // x-point
          xpt_r[*xpt_n] = null_r;
          xpt_z[*xpt_n] = null_z;
          xpt_flux[*xpt_n] = flux_at_null;
          (*xpt_n)++;
          if (*xpt_n >= N_XPT_MAX) {
            return 512;
          }
        }
      }
    }
  }
  return 0;
}

// filter_xpts:
// loop over every x-point and check if they are behind another x-point relative
// to the o-point. If they are then remove them. See Fig. 2 of Moret et al.
// (2015) for reference.
//
// The ith x-point is located at (xpt_r[i], xpt_z[i]) and the
// jth x-point is located at (xpt_r[j], xpt_z[j]).
//
// The vector v points from the ith x-point to the jth x-point is
// v = (xpt_r[j] - xpt_r[i], xpt_z[j] - xpt_z[i]).
// The vector w that points from the ith x-point to the magnetic axis is
// w = (r_mag_axis - xpt_r[i], z_mag_axis - xpt_z[i]).
//
// If the dot product of v and w is negative, then the jth x-point
// is not considered for the flux calculation.
void filter_xpts(double *xpt_r, double *xpt_z, int32_t *xpt_n,
                 double r_mag_axis, double z_mag_axis) {

  // keep[i] = 1 if x-point i is kept, 0 if removed
  int keep[N_XPT_MAX];
  for (int i = 0; i < *xpt_n; ++i)
    keep[i] = 1;

  for (int i = 0; i < *xpt_n; ++i) {

    // w = mag_axis - xpt_i
    const double wi_r = r_mag_axis - xpt_r[i];
    const double wi_z = z_mag_axis - xpt_z[i];

    for (int j = 0; j < *xpt_n; ++j) {
      if (j == i)
        continue;

      // v = xpt_j - xpt_i
      const double vij_r = xpt_r[j] - xpt_r[i];
      const double vij_z = xpt_z[j] - xpt_z[i];

      const double dot = vij_r * wi_r + vij_z * wi_z;

      // If dot < 0, j is behind i relative to the axis => remove j
      if (dot < 0) {
        keep[j] = 0;
      }
    }
  }

  // Flush the kept x-points to the front of the arrays.
  int k = 0;
  for (int i = 0; i < *xpt_n; ++i) {
    if (keep[i]) {
      xpt_r[k] = xpt_r[i];
      xpt_z[k] = xpt_z[i];
      ++k;
    }
  }

  *xpt_n = k;
}

// find_mask:
//   Determines which grid points are inside the last closed flux surface (LCFS)
//   based on the value of ψ at the LCFS boundary.
//   Uses a flood-fill algorithm.
int32_t find_mask(double *flux_total, double z_mag_axis, double *lcfs_r,
                  double *lcfs_z, int32_t lcfs_n, int32_t *mask) {
  int32_t i_grid, i_lcfs, i_lcfs_next;
  double cross_product;
  int32_t inside;

  for (i_grid = 0; i_grid < N_GRID; i_grid++) {
    if (!MASK_LIM[i_grid]) {
      continue;
    }

    inside = 1;

    for (i_lcfs = 0; i_lcfs < lcfs_n; i_lcfs++) {
      i_lcfs_next = (i_lcfs + 1) % lcfs_n;

      cross_product = (lcfs_r[i_lcfs_next] - lcfs_r[i_lcfs]) *
                          (Z_GRID[i_grid] - lcfs_z[i_lcfs]) -
                      (lcfs_z[i_lcfs_next] - lcfs_z[i_lcfs]) *
                          (R_GRID[i_grid] - lcfs_r[i_lcfs]);

      if (cross_product > 0.0) {
        inside = 0;
        break;
      }
    }

    if (inside) {
      mask[i_grid] = 1;
    } else {
      mask[i_grid] = 0;
    }
  }

  return 0;
}