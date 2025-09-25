#include "find_x_point.h"
#include "constants.h"
#include "gradient.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#define N_PATCH_CORNERS  4 // number of corners in a cell (rectangle)
#define N_PATCH_PAIRS    6 // maximum number of unique unordered connections between up to 4 edge crossings
#define N_R_MIN_1 (N_R - 1)
#define N_R_PLS_1 (N_R + 1)
#define N_Z_MIN_1 (N_Z - 1)
#define N_Z_PLS_1 (N_Z + 1)
#define N_R_TIMES_2 (N_R * 2 - 1)
#define N_Z_TIMES_2 (N_Z * 2 - 1)
#define MIN(aa,bb) ((aa)<=(bb)?(aa):(bb))
#define MAX(aa,bb) ((aa)>=(bb)?(aa):(bb))
#define ERR_MID_BDRY 0b000001
#define ERR_COL_START 0b000010
#define ERR_COL_END 0b000100
#define ERR_COL_START_END 0b001000
#define ERR_VERT_BDRY 0b010000
#define ERR_MASK_INDEX 0b100000

/// Find crossings along the edges of a patch.
/// 
/// Evaluates the array `grad_patch` on each edge of the patch and determines where it
/// changes sign (i.e., crosses zero). The function records:
///   - the total number of crossings (0–4),
///   - the coordinates of each crossing.
void find_zero_on_edge(
    double* grad_patch,
    double* cross_row,
    double* cross_col,
    int32_t* count) {
  const int32_t rowNum[N_PATCH_CORNERS] = {0, 0, 1, 1};
  const int32_t colNum[N_PATCH_CORNERS] = {0, 1, 0, 1};
  const int32_t nn_c[N_PATCH_CORNERS] = {1, 1, 0, 0}; // next node column
  const int32_t nn_r[N_PATCH_CORNERS] = {0, 1, 0, 1}; // next node row
  const int32_t nn_ind[N_PATCH_CORNERS] = {1, 3, 0, 2}; // next node index

  *count = 0;
  for (int ii = 0; ii < N_PATCH_CORNERS; ii++) {
    if (fabs(grad_patch[ii]) < THRESH) {
      cross_col[*count] = colNum[ii];
      cross_row[*count] = rowNum[ii];
      (*count)++;
    } else if (grad_patch[ii] * grad_patch[nn_ind[ii]] < 0) {
      int32_t diff_col = -colNum[ii] + nn_c[ii];
      int32_t diff_row = -rowNum[ii] + nn_r[ii];
      double off = fabs(grad_patch[ii]) / \
        (fabs(grad_patch[ii]) + fabs(grad_patch[nn_ind[ii]]));
      cross_col[*count] = colNum[ii] + diff_col * off;
      cross_row[*count] = rowNum[ii] + diff_row * off;
      (*count)++;
    }
  }
}

/// Generate all unique unordered index pairs between edge crossings in a patch.
///
/// unordered pairs of crossing indices. The results represent the line segments
/// that could connect crossing points of a gradient function across patch edges.
/// 
/// @param[in]  n_count  Number of edge crossings in the patch (0–4).
/// @param[out] n_comb   Number of unique unordered pairs generated.
/// @param[out] idx_a    Array of first indices of each pair (size up to N_PATCH_PAIRS).
/// @param[out] idx_b    Array of second indices of each pair (size up to N_PATCH_PAIRS).
void n_comb(
    int32_t n_count,
    int32_t *n_comb,
    int32_t *idx_a,
    int32_t *idx_b) {

  switch (n_count){
  case 0:
    *n_comb = 0;
    printf("Error in n_comb: this function should not be called when n_count=0.\n");
    break;
  case 1:
    *n_comb = 1;
    idx_a[0] = 0;
    idx_b[0] = 0;
    break;
  case 2:
    *n_comb = 1;
    idx_a[0] = 0;
    idx_b[0] = 1;
    break;
  case 3:
    *n_comb = 3;
    idx_a[0] = 0;
    idx_b[0] = 1;
    idx_a[1] = 0;
    idx_b[1] = 2;
    idx_a[2] = 1;
    idx_b[2] = 2;
    break;
  case 4:
    *n_comb = 6;
    idx_a[0] = 0;
    idx_b[0] = 1;
    idx_a[1] = 0;
    idx_b[1] = 2;
    idx_a[2] = 0;
    idx_b[2] = 3;
    idx_a[3] = 1;
    idx_b[3] = 2;
    idx_a[4] = 1;
    idx_b[4] = 3;
    idx_a[5] = 2;
    idx_b[5] = 3;
    break;
  default:
    *n_comb = 0;
    printf("Error in n_comb: unexpected n_count=%d. This should not be possible.\n", n_count);
    break;
  }
}

void find_xpoint_cross(
    double a_start_row,
    double a_end_row,
    double a_start_col,
    double a_end_col,
    double b_start_row,
    double b_end_row,
    double b_start_col,
    double b_end_col,
    double *col_off,
    double *row_off,
    bool *flag) {
  double a_diff_row = a_end_row - a_start_row;
  double a_diff_col = a_end_col - a_start_col;
  double b_diff_row = b_end_row - b_start_row;
  double b_diff_col = b_end_col - b_start_col;
  double det = -a_diff_col * b_diff_row + b_diff_col * a_diff_row;
  double ab_start_col = a_start_col - b_start_col;
  double ab_start_row = (a_start_row - b_start_row);
  double lambda_a = (b_diff_row * ab_start_col - b_diff_col * ab_start_row) / det;
  double lambda_b = (a_diff_row * ab_start_col - a_diff_col * ab_start_row) / det;

  *flag = false;
  if((lambda_a >= 0) && (lambda_a <= 1) && (lambda_b >= 0) && (lambda_b <= 1)) {
    *col_off = a_start_col + a_diff_col * lambda_a;
    *row_off = a_start_row + a_diff_row * lambda_a;
    *flag = true;
  }
}

void find_null_in_gradient_march(
    double* flux,
    double* opt_r,
    double* opt_z,
    double* opt_flux,
    int32_t* opt_n,
    double* xpt_r,
    double* xpt_z,
    double* xpt_flux,
    int32_t* xpt_n) {
  
  double hess_zz[N_GRID];
  double hess_rr[N_GRID];
  double hess_rz[N_GRID];
  double hess_det_at_null;
  double hess_rr_at_null;
  double hess_zz_at_null;
  double hess_rz_at_null;
  double grad_z[N_GRID];
  double grad_r[N_GRID];

  *opt_n = 0;
  *xpt_n = 0;
  // find hessian & gradients
  gradient_z(flux, grad_z);
  gradient_r(flux, grad_r);
  hessian_zz(flux, hess_zz);
  hessian_rr(flux, hess_rr);
  gradient_r(grad_z, hess_rz);
  for (int32_t i_row = 0; i_row < N_Z_MIN_1; i_row++) {
    for (int32_t i_col = 0; i_col < N_R_MIN_1; i_col++)  {
      double gr_patch[N_PATCH_CORNERS];
      double gz_patch[N_PATCH_CORNERS];
      double grs_row[N_PATCH_CORNERS];
      double gzs_row[N_PATCH_CORNERS];
      double grs_col[N_PATCH_CORNERS];
      double gzs_col[N_PATCH_CORNERS];
      int32_t idx = i_row * N_R + i_col;
      gr_patch[0] = grad_r[idx];
      gr_patch[1] = grad_r[idx + 1];
      gr_patch[2] = grad_r[idx + N_R];
      gr_patch[3] = grad_r[idx + N_R_PLS_1];
      gz_patch[0] = grad_z[idx];
      gz_patch[1] = grad_z[idx + 1];
      gz_patch[2] = grad_z[idx + N_R];
      gz_patch[3] = grad_z[idx + N_R_PLS_1];

      int32_t count_r;
      int32_t count_z;
      find_zero_on_edge(gr_patch, grs_row, grs_col, &count_r);
      find_zero_on_edge(gz_patch, gzs_row, gzs_col, &count_z);

      if ((count_r > 0) && (count_z > 0)) {
        int32_t idx_r_start[N_PATCH_PAIRS];
        int32_t idx_r_end[N_PATCH_PAIRS];
        int32_t idx_z_start[N_PATCH_PAIRS];
        int32_t idx_z_end[N_PATCH_PAIRS];
        int32_t n_comb_r;
        int32_t n_comb_z;
        n_comb(count_r, &n_comb_r, idx_r_start, idx_r_end);
        n_comb(count_z, &n_comb_z, idx_z_start, idx_z_end);

        for (int r_idx = 0; r_idx < n_comb_r; r_idx++) {
          double gr_row_start = grs_row[idx_r_start[r_idx]];
          double gr_row_end = grs_row[idx_r_end[r_idx]];
          double gr_col_start = grs_col[idx_r_start[r_idx]];
          double gr_col_end = grs_col[idx_r_end[r_idx]];
          for (int z_idx = 0; z_idx < n_comb_z; z_idx++) {
            double gz_row_start = gzs_row[idx_z_start[z_idx]];
            double gz_row_end = gzs_row[idx_z_end[z_idx]];
            double gz_col_start = gzs_col[idx_z_start[z_idx]];
            double gz_col_end = gzs_col[idx_z_end[z_idx]];

            if ((MIN(gz_row_start, gz_row_end) <= MAX(gr_row_start, gr_row_end)) &&
              (MIN(gr_row_start, gr_row_end) <= MAX(gz_row_start, gz_row_end)) &&
              (MIN(gz_col_start, gz_col_end) <= MAX(gr_col_start, gr_col_end)) &&
              (MIN(gr_col_start, gr_col_end) <= MAX(gz_col_start, gz_col_end))) {

              double col_off;
              double row_off;
              bool is_stnry;
              find_xpoint_cross(gz_row_start, gz_row_end, gz_col_start,
              gz_col_end, gr_row_start, gr_row_end, gr_col_start,
              gr_col_end, &col_off, &row_off, &is_stnry);

              if (is_stnry) {
                hess_rr_at_null = lin_intrp_2(hess_rr, idx, col_off, row_off);
                hess_zz_at_null = lin_intrp_2(hess_zz, idx, col_off, row_off);
                hess_rz_at_null = lin_intrp_2(hess_rz, idx, col_off, row_off);
                hess_det_at_null = hess_rr_at_null * hess_zz_at_null - \
                  hess_rz_at_null * hess_rz_at_null;

                if (hess_det_at_null > 0.0 && hess_rr_at_null < 0.0 && MASK_LIM[idx]) {
                  opt_r[*opt_n] = R_VEC[i_col] + DR * col_off;
                  opt_z[*opt_n] = Z_VEC[i_row] + DZ * row_off;
                  opt_flux[*opt_n] = lin_intrp_2(flux, idx, col_off, row_off);
                  (*opt_n)++;
                }
                else if (hess_det_at_null < 0.0 && MASK_LIM[idx]) {
                  xpt_r[*xpt_n] = R_VEC[i_col] + DR * col_off;
                  xpt_z[*xpt_n] = Z_VEC[i_row] + DZ * row_off;
                  xpt_flux[*xpt_n] = lin_intrp_2(flux, idx, col_off, row_off);
                  (*xpt_n)++;
                }
              }
            }
          }
        }
      }
    }
  }
}

double lin_intrp_2(
    double* flux,
    int32_t idx,
    double frac_dist_null_r,
    double frac_dist_null_z) {
  double flux_up, flux_down, flux_at_null;

  flux_up = (1.0 - frac_dist_null_r) * flux[idx + N_R] + \
    frac_dist_null_r * flux[idx + N_R_PLS_1];
  flux_down = (1.0 - frac_dist_null_r) * flux[idx] + frac_dist_null_r * flux[idx + 1];
  flux_at_null = (1.0 - frac_dist_null_z) * flux_down + frac_dist_null_z * flux_up;

  return flux_at_null;
}

double lin_intrp(
    double *flux,
    int32_t idx,
    double dist_to_null_r,
    double dist_to_null_z,
    double abs_dist_null_r,
    double abs_dist_null_z) {
  double flux_up, flux_down, flux_at_null;

  abs_dist_null_r = abs_dist_null_r / DR;
  abs_dist_null_z = abs_dist_null_z / DZ;
  if (dist_to_null_r > 0.0) {
    if (dist_to_null_z > 0.0) {
      flux_up = (1 - abs_dist_null_r) * flux[idx + N_R] + abs_dist_null_r * \
        flux[idx + N_R_PLS_1];
      flux_down = (1 - abs_dist_null_r) * flux[idx] + abs_dist_null_r * \
        flux[idx + 1];
      flux_at_null = (1 - abs_dist_null_z) * flux_down + abs_dist_null_z * flux_up;
    } else {
      flux_up = (1 - abs_dist_null_r) * flux[idx] + abs_dist_null_r * \
        flux[idx + 1];
      flux_down = (1 - abs_dist_null_r) * flux[idx - N_R] + abs_dist_null_r * \
        flux[idx - N_R_PLS_1];
      flux_at_null = (1 - abs_dist_null_z) * flux_up + abs_dist_null_z * flux_down;
    }
  } else {
    if (dist_to_null_z > 0.0) {
      flux_up = (1 - abs_dist_null_r) * flux[idx + N_R] + abs_dist_null_r * \
        flux[idx + N_R_MIN_1];
      flux_down = (1 - abs_dist_null_r) * flux[idx] + abs_dist_null_r * \
        flux[idx - 1];
      flux_at_null = (1 - abs_dist_null_z) * flux_down + abs_dist_null_z * flux_up;
    } else {
      flux_up = (1 - abs_dist_null_r)*flux[idx] + abs_dist_null_r * \
        flux[idx - 1];
      flux_down = (1 - abs_dist_null_r)*flux[idx - N_R] + abs_dist_null_r * \
        flux[idx - N_R_MIN_1];
      flux_at_null = (1 - abs_dist_null_z) * flux_up + abs_dist_null_z * flux_down;
    }
  }
  return flux_at_null;
}

void find_lcfs_rz(
    double *flux,
    double flux_lcfs,
    double *lcfs_r,
    double *lcfs_z,
    int32_t *lcfs_n) {
  double flux_offset[N_GRID];

  *lcfs_n = 0;

  // Offset the flux. The boundary is where `flux_offset = 0`
  for (int32_t idx = 0; idx < (N_GRID); idx++) {
    flux_offset[idx] = flux[idx] - flux_lcfs;
  }

  // Loop over z
  // We are excluding the last row because we will be comparing the flux at this gird point (r, z) to the flux at (r + d_r, z)
  for (int32_t i_row = 0; i_row < N_Z_MIN_1; i_row++) {
    // Loop over r
    // We are excluding the last column because we will be comparing the flux at this gird point (r, z) to the flux at (r, z + d_z)
    for (int32_t i_col = 0; i_col < N_R_MIN_1; i_col++) {
      int32_t idx = i_row * N_R + i_col;

      // If flux is small and within the vacuum vessel add grid point to lcfs
      // THRESH=1e-10 is really small; so it won't happen often; is there any point in doing this? or should we always do linear interpolation?
      if ((fabs(flux_offset[idx]) < THRESH)  && MASK_LIM[idx]) {
        lcfs_r[*lcfs_n] = R_VEC[i_col];
        lcfs_z[*lcfs_n] = Z_VEC[i_row];
        (*lcfs_n)++;
      } else {
        // Compare the flux at this grid point (r, z) to the flux at (r + d_r, z).
        // If there is a sign change then the lcfs must be between grid points
        // Also needs to be within the vacuum vessel
        if ((flux_offset[idx] * flux_offset[idx + 1] < 0)  && (MASK_LIM[idx] || MASK_LIM[idx + 1])) {
          // Calculate the ratio of the flux at this grid point to the flux at (R+dR, Z)
          // frac will be between 0.0 and 1.0; but can't be very close to 0.0 otherwise
          // it would have been picked up by the first test `fabs(flux[idx]) < THRESH`
          double off = fabs(flux_offset[idx]) / (fabs(flux_offset[idx]) + fabs(flux_offset[idx + 1]));
          if (((off <= 0.5) && MASK_LIM[idx]) || ((off >= 0.5) && MASK_LIM[idx + 1])) {
            lcfs_r[*lcfs_n] = R_VEC[i_col] + off * DR; // linear interpolation
            lcfs_z[*lcfs_n ] = Z_VEC[i_row];
            (*lcfs_n)++;
          }
        }
        // Compare the flux at this grid point (r, z) to the flux at (r, z + d_z).
        // If there is a sign change then the lcfs must be between grid points
        if ((flux_offset[idx] * flux_offset[idx + N_R] < 0) &&
             (MASK_LIM[idx] || MASK_LIM[idx + N_R])) {
          double off = fabs(flux_offset[idx]) / (fabs(flux_offset[idx]) + fabs(flux_offset[idx + N_R]));
          if (((off <= 0.5) && MASK_LIM[idx]) ||
              ((off >= 0.5) && MASK_LIM[idx + N_R])) {
            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n ] = Z_VEC[i_row] + off * DZ; // linear interpolation
            (*lcfs_n)++;
          }
        }
      }
    }
  }

  // Loop over z (excluding maximum z) with r at the last grid point
  // Checking the mesh points missed in the first loop; (max(r), max(z)) checked in next loop
  for (int32_t i_row = 0; i_row < N_Z - 1; i_row++) {
    int32_t idx = i_row * N_R + N_R - 1;
    if ((fabs(flux_offset[idx]) < THRESH) && MASK_LIM[idx]) {
      lcfs_r[*lcfs_n] = R_VEC[N_R_MIN_1];
      lcfs_z[*lcfs_n] = Z_VEC[i_row];
      (*lcfs_n)++;
    }
  }

  // Loop over r with z at the last grid point
  // Checking the mesh points missed in the first loop
  for (int32_t i_col = 0; i_col < N_R; i_col++) {
    int32_t idx = (N_Z - 1) * N_R + i_col;
    if ((fabs(flux_offset[idx]) < THRESH) && MASK_LIM[idx]) {
      lcfs_r[*lcfs_n] = R_VEC[i_col];
      lcfs_z[*lcfs_n] = Z_VEC[N_Z_MIN_1];
      (*lcfs_n)++;
    }
  }
}

int inside_lcfs(
    double r_opt,
    double z_opt,
    double *lcfs_r,
    double *lcfs_z,
    int32_t lcfs_n,
    int32_t *mask) {
  int32_t error = 0;
  double r_start = -DBL_MAX;
  double r_end = DBL_MAX;
  double r_tmp[N_R_TIMES_2];
  double z_tmp[N_Z_TIMES_2];
  int32_t n_mid_plane_bdry = 0;

  memset(mask, 0, N_GRID * sizeof(int));

  // Find the z grid point closest to the o-point (magnetic axis)
  double z_nearest = round((z_opt - Z_VEC[0]) / DZ) * DZ + Z_VEC[0];

  // Add r points to `r_tmp` which are at the same z as the o-point
  // Note: there should be at least the LFS and HFS points. If not, then there is an error
  // But there may be more points (around the solenoid)
  for (int32_t i_lcfs = 0; i_lcfs < lcfs_n; i_lcfs++) {
    if (fabs(lcfs_z[i_lcfs] - z_nearest) < THRESH) {
      r_tmp[n_mid_plane_bdry] = lcfs_r[i_lcfs];
      n_mid_plane_bdry++;
    }
  }

  // If we have not found the LFS and HFS boundary points, then there is an error
  // we cannot continue and should exit
  if (n_mid_plane_bdry < 2) {
    error |= ERR_MID_BDRY;
    return error;
  }

  // `r_start` = LFS; `r_end` = HFS
  for (int32_t i_count=0; i_count< n_mid_plane_bdry; i_count++) {
    if (r_tmp[i_count] < r_opt) {
      // For the LFS, we prefer points closer to the magnetic axis, i.e. not the solenoid
      if (r_tmp[i_count] > r_start) {
        r_start = r_tmp[i_count];
      }
    } else {
      if (r_tmp[i_count] < r_end) {
        r_end = r_tmp[i_count];
      }
    }
  }

  // `col_start` and `col_end` are the indexes of the LFS and HFS
  int32_t col_start = (int) ceil((r_start - R_VEC[0] ) / DR);
  int32_t col_end = (int) floor((r_end - R_VEC[0]) / DR);
  // Check that indexing is within bounds; will only happen under exotic condtions
  // such as the plasma extending outside the (R, Z) grid. But we cannot continue
  // under such conditions and should exit
  if ((col_start < 0) || (col_start > N_R)) {
    error |= ERR_COL_START;
  }
  if ((col_end < 0) || (col_end > N_R)) {
    error |= ERR_COL_END;
  }
  if (col_end < col_start) {
    error |= ERR_COL_START_END;
  }
  if (error) {
    return error;
  }

  // Assume plasma is widest at the magnetic axis -> if this is not true then the present algorithm
  // will miss part of the plasma boundary
  // Loop in r from LFS to HFS
  for (int32_t i_col = col_start; i_col <= col_end; i_col++) {
    // Loop over all possible boundary points and find the z's which could be the boundary
    // Note there should be at least 2 points (upper and lower),
    // but there may be more points (espescially around the private flux region)
    int32_t n_vert_bdry = 0;
    for (int32_t i_lcfs = 0; i_lcfs < lcfs_n; i_lcfs++) {
      if (fabs(lcfs_r[i_lcfs] - R_VEC[i_col]) < THRESH) {
        z_tmp[n_vert_bdry] = lcfs_z[i_lcfs];
        n_vert_bdry++;
      }
    }
    if (n_vert_bdry < 2) {
      error |= ERR_VERT_BDRY;
      return error;
    }

    // Find the two points which are vertically closest to the magnetic axis
    double min_diff = DBL_MAX;
    double second_min_diff = DBL_MAX;
    double z_closest = DBL_MAX; // initial value shouldn't matter
    double z_second_closest = DBL_MAX; // initial value shouldn't matter
    for (int32_t i_count = 0; i_count < n_vert_bdry; i_count++) {
      double diff = fabs(z_tmp[i_count] - z_opt);
      if (diff < min_diff) {
        second_min_diff = min_diff;
        z_second_closest = z_closest;
        z_closest = z_tmp[i_count];
        min_diff = diff;
      } else if (diff < second_min_diff) {
        second_min_diff = diff;
        z_second_closest = z_tmp[i_count];
      }
    }
    double z_start = MIN(z_closest, z_second_closest);
    double z_end = MAX(z_closest, z_second_closest);

    // March vertically and apply the mask
    int32_t row_start = (int) ceil((z_start - Z_VEC[0] ) / DZ);
    int32_t row_end = (int) floor((z_end - Z_VEC[0]) / DZ);
    if ((row_start >= 0) && (row_start <= row_end) && (row_end <= N_Z)) {
      for (int32_t i_row = row_start; i_row <= row_end; i_row++) {
        mask[i_row * N_R + i_col] = 1;
      }
    } else {
      // This error should not be reachable because it is handled by
      // `ERR_VERT_BDRY` error check (n_vert_bdry<2)
      error = ERR_MASK_INDEX;
      return error;
    }
  }
  return error;
}
