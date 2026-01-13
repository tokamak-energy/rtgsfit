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
 *   - We refer to the region inside the LCFS as the "plasma region" or the
 * "plasma core". While plasma also exists outside the LCFS, the pressure and
 * temperature are significantly higher within it. This is therefore the region
 * of primary interest for equilibrium reconstruction.
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
    int32_t i_row_n_r = i_row * N_R;
    for (int32_t i_col = 1; i_col < N_R_MIN_1; i_col++) {
      int32_t idx = i_row_n_r + i_col;
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
      double denom = c * d - e * e;
      if (fabs(denom) < THRESH)
        continue;
      double inv_denom = 1.0 / denom;
      double dr_norm = (b * e - a * d) * inv_denom; // (dr / ΔR)
      double dz_norm = (a * e - b * c) * inv_denom; // (dz / ΔZ)
      if (fabs(dr_norm) <= 0.5 && fabs(dz_norm) <= 0.5) {
        double null_r = R_VEC[i_col] + dr_norm * DR;
        double null_z = Z_VEC[i_row] + dz_norm * DZ;
        double hess_det = denom;
        double flux_at_null = flux[idx] + a * dr_norm + b * dz_norm +
                              0.5 * c * dr_norm * dr_norm +
                              0.5 * d * dz_norm * dz_norm +
                              e * dr_norm * dz_norm;
        // Redundant check (already handled by fabs(denom) above)
        // if (fabs(hess_det) < THRESH)
        //     continue;
        if (hess_det > 0.0) {
          // o-point
          opt_r[*opt_n] = null_r;
          opt_z[*opt_n] = null_z;
          opt_flux[*opt_n] = flux_at_null;
          (*opt_n)++;
          if (*opt_n >= N_XPT_MAX) {
            return ERR_NUM_OPTS;
          }
        } else if (hess_det < 0.0) {
          // x-point
          xpt_r[*xpt_n] = null_r;
          xpt_z[*xpt_n] = null_z;
          xpt_flux[*xpt_n] = flux_at_null;
          (*xpt_n)++;
          if (*xpt_n >= N_XPT_MAX) {
            return ERR_NUM_XPTS;
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
//
// Parameters:
//   xpt_r       - array of x-point R coordinates (with size N_XPT_MAX)
//   xpt_z       - array of x-point Z coordinates (with size N_XPT_MAX)
//   xpt_n       - pointer which holds the number of x-points found
//   r_mag_axis  - R coordinate of the magnetic axis
//   z_mag_axis  - Z coordinate of the magnetic axis
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

// is_core_side_of_xpoint:
// Similar to filter_xpts, but for a single grid point.
// Returns 1 if the grid point is on the magnetic axis side of all x-points
// and 0 otherwise.
//
// Parameters:
//   r_grid      - R coordinate of the grid point
//   z_grid      - Z coordinate of the grid point
//   r_mag_axis  - R coordinate of the magnetic axis
//   z_mag_axis  - Z coordinate of the magnetic axis
//   xpt_r       - array of x-point R coordinates (with size N_XPT_MAX)
//   xpt_z       - array of x-point Z coordinates (with size N_XPT_MAX)
//   xpt_n       - number of x-points
int is_core_side_of_xpoint(double r_grid, double z_grid, double r_mag_axis,
                           double z_mag_axis, double *xpt_r, double *xpt_z,
                           int32_t xpt_n) {

  for (int32_t i_xpt = 0; i_xpt < xpt_n; i_xpt++) {
    // v = grid_point - xpt_i
    const double vi_r = r_grid - xpt_r[i_xpt];
    const double vi_z = z_grid - xpt_z[i_xpt];

    // w = mag_axis - xpt_i
    const double wi_r = r_mag_axis - xpt_r[i_xpt];
    const double wi_z = z_mag_axis - xpt_z[i_xpt];

    const double dot = vi_r * wi_r + vi_z * wi_z;

    // If dot < 0, grid point is behind xpt_i relative to the axis => filter it
    if (dot < 0) {
      return 0; // dont include this grid point as part of the plasma core
    }
  }

  return 1; // include this grid point as part of the plasma core
}

// Converts a 2D grid index (column i_r, row i_z)
// into a row-major 1D array index.
static inline int32_t grid_idx(int32_t i_r, int32_t i_z) {
  // R varies fastest
  return i_r + N_R * i_z;
}

// Finds the index of the element in a 1D array `vec` (length `n`) 
// that is closest to the value `x`.
static int32_t nearest_index_1d(const double *vec, int32_t n, double x) {
  int32_t best_i = 0;
  double best_d = fabs(vec[0] - x);
  for (int32_t i = 1; i < n; i++) {
    const double d = fabs(vec[i] - x);
    if (d < best_d) {
      best_d = d;
      best_i = i;
    }
  }
  return best_i;
}

// is_core_candidate:
// Checks if the grid point at (i_r, i_z) is a candidate for inclusion
// in the plasma core region.
// A candidate must satisfy:
// 1. ψ > ψ_boundary
// 2. Be on the magnetic axis side of all x-points
// 3. MASK_LIM[idx] is true
// Parameters:
//   i_r         - R index of the grid point
//   i_z         - Z index of the grid point
//   R_VEC       - array of R coordinates (with length N_R)
//   Z_VEC       - array of Z coordinates (with length N_Z)
//   flux_total  - array of total flux values (with length N_GRID)
//   flux_boundary- flux value at the LCFS boundary
//   r_mag_axis  - R coordinate of the magnetic axis
//   z_mag_axis  - Z coordinate of the magnetic axis
//   xpt_r       - array of x-point R coordinates (with size N_XPT_MAX)
//   xpt_z       - array of x-point Z coordinates (with size N_XPT_MAX)
//   xpt_n       - number of x-points
static inline int is_core_candidate(int32_t i_r, int32_t i_z,
                                    const double *flux_total,
                                    double flux_boundary, double r_mag_axis,
                                    double z_mag_axis, double *xpt_r,
                                    double *xpt_z, int32_t xpt_n) {
  const int32_t idx = grid_idx(i_r, i_z);

  const double psi = flux_total[idx];

  // psi > flux_boundary in plasma core
  if (psi <= flux_boundary)
    return 0;

  if (!MASK_LIM[idx])
    return 0;

  const double r = R_VEC[i_r];
  const double z = Z_VEC[i_z];
  return is_core_side_of_xpoint(r, z, r_mag_axis, z_mag_axis, xpt_r, xpt_z,
                                xpt_n);
}

// flood_fill_plasma_core:
// Flood-fill algorithm to identify the plasma core region inside the LCFS.
// The mask array is updated in-place to mark grid cells inside the plasma core.
// We also make sure the private flux region is excluded by using the dot
// product trick above to make sure the plasma core is not behind any x-points
// relative to the magnetic axis.
// returns 0 on success, nonzero on error
int flood_fill_plasma_core(int32_t *mask, double *flux_total,
                           double flux_boundary, double r_mag_axis,
                           double z_mag_axis, double *xpt_r, double *xpt_z,
                           int32_t xpt_n) {

  int32_t queue[N_GRID];
  int32_t head = 0, tail = 0;

  // Seed near magnetic axis
  const int32_t i_r0 = nearest_index_1d(R_VEC, N_R, r_mag_axis);
  const int32_t i_z0 = nearest_index_1d(Z_VEC, N_Z, z_mag_axis);

  // Strict check: seed must satisfy core conditions
  if (!is_core_candidate(i_r0, i_z0, flux_total, flux_boundary, r_mag_axis,
                         z_mag_axis, xpt_r, xpt_z, xpt_n)) {
    return ERR_AXIS_OUT_CORE;
  }

  // Reset mask
  for (int32_t i = 0; i < N_GRID; i++) {
    mask[i] = 0;
  }

  // Push seed
  const int32_t seed_idx = grid_idx(i_r0, i_z0);
  mask[seed_idx] = 1;
  queue[tail++] = seed_idx;

  // 4-neighbour BFS (breadth-first search) flood fill
  while (head < tail) {

    const int32_t idx = queue[head++];
    const int32_t i_z = idx / N_R;
    const int32_t i_r = idx - i_z * N_R;

    const int32_t nbr_r[4] = {i_r - 1, i_r + 1, i_r, i_r};
    const int32_t nbr_z[4] = {i_z, i_z, i_z - 1, i_z + 1};

    for (int k = 0; k < 4; k++) {
      const int32_t i_r_nbr = nbr_r[k];
      const int32_t i_z_nbr = nbr_z[k];

      if (i_r_nbr < 0 || i_r_nbr >= N_R)
        continue;
      if (i_z_nbr < 0 || i_z_nbr >= N_Z)
        continue;

      const int32_t idx_nbr = grid_idx(i_r_nbr, i_z_nbr);
      if (mask[idx_nbr])
        continue;

      if (is_core_candidate(i_r_nbr, i_z_nbr, flux_total, flux_boundary,
                            r_mag_axis, z_mag_axis, xpt_r, xpt_z, xpt_n)) {
        mask[idx_nbr] = 1;
        queue[tail++] = idx_nbr;
      }
    }
  }

  return 0; // success
}
