#include "magnetic_nulls.h"

#include <math.h>

/* Stationary point of the second-order Taylor expansion of psi about one grid point. */
struct null_candidate {
    double d_i_r; /* offset from the grid point, in units of d_r */
    double d_i_z; /* offset from the grid point, in units of d_z */
    double hessian_determinant;
    double psi;
};

/* Returns 1 and fills `candidate` unless the Hessian is (nearly) singular. */
static int32_t null_candidate_at(const struct rtgsfit_grid *grid,
                                 const double psi[],
                                 int32_t i_r,
                                 int32_t i_z,
                                 double tolerance,
                                 struct null_candidate *candidate)
{
    const int32_t n_r = grid->n_r;
    const int32_t i = grid_index(grid, i_r, i_z);

    /* Derivatives in grid units (d_r = d_z = 1). */
    const double psi_r = 0.5 * (psi[i + 1] - psi[i - 1]);
    const double psi_z = 0.5 * (psi[i + n_r] - psi[i - n_r]);
    const double psi_rr = psi[i + 1] - 2.0 * psi[i] + psi[i - 1];
    const double psi_zz = psi[i + n_r] - 2.0 * psi[i] + psi[i - n_r];
    const double psi_rz =
        0.25 * (psi[i + 1 + n_r] - psi[i - 1 + n_r] - psi[i + 1 - n_r] + psi[i - 1 - n_r]);

    const double determinant = psi_rr * psi_zz - psi_rz * psi_rz;
    if (fabs(determinant) < tolerance) {
        return 0;
    }
    const double inv_determinant = 1.0 / determinant;
    const double d_r = (psi_z * psi_rz - psi_r * psi_zz) * inv_determinant;
    const double d_z = (psi_r * psi_rz - psi_z * psi_rr) * inv_determinant;

    candidate->d_i_r = d_r;
    candidate->d_i_z = d_z;
    candidate->hessian_determinant = determinant;
    candidate->psi = psi[i] + psi_r * d_r + psi_z * d_z + 0.5 * psi_rr * d_r * d_r +
                     0.5 * psi_zz * d_z * d_z + psi_rz * d_r * d_z;
    return 1;
}

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
                            int32_t *xpt_n)
{
    const int32_t status = grid_check(grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    RTGSFIT_CHECK_SIZE(n_psi, grid_n_grid(grid));
    RTGSFIT_CHECK_SIZE(n_mask_limiter, grid_n_grid(grid));
    RTGSFIT_CHECK_MIN_SIZE(n_xpt_r, 1);
    RTGSFIT_CHECK_SIZE(n_xpt_z, n_xpt_r);
    RTGSFIT_CHECK_SIZE(n_xpt_psi, n_xpt_r);

    static const int32_t corner_d_i_r[4] = {0, 1, 0, 1};
    static const int32_t corner_d_i_z[4] = {0, 0, 1, 1};

    *axis_found = 0;
    *xpt_n = 0;

    /* Cell (i_r, i_z) has corners i_r..i_r+1 and i_z..i_z+1, all of which need four neighbours. */
    for (int32_t i_z = 1; i_z < grid->n_z - 2; i_z++) {
        for (int32_t i_r = 1; i_r < grid->n_r - 2; i_r++) {
            int32_t found = 0;
            double best_distance_sq = 0.0;
            double best_r = 0.0;
            double best_z = 0.0;
            struct null_candidate best = {0.0, 0.0, 0.0, 0.0};

            for (int32_t k = 0; k < 4; k++) {
                const int32_t corner_i_r = i_r + corner_d_i_r[k];
                const int32_t corner_i_z = i_z + corner_d_i_z[k];
                struct null_candidate candidate;

                if (!mask_limiter[grid_index(grid, corner_i_r, corner_i_z)]) {
                    continue;
                }
                if (!null_candidate_at(grid, psi, corner_i_r, corner_i_z, tolerance, &candidate)) {
                    continue;
                }
                /*
                 * Position of the null within the cell. The interval is half-open, [0, 1), so a
                 * null on a shared edge or vertex belongs to exactly one cell and is found once.
                 */
                const double u = (double)corner_d_i_r[k] + candidate.d_i_r;
                const double v = (double)corner_d_i_z[k] + candidate.d_i_z;
                if (u < 0.0 || u >= 1.0 || v < 0.0 || v >= 1.0) {
                    continue;
                }
                const double distance_sq =
                    candidate.d_i_r * candidate.d_i_r + candidate.d_i_z * candidate.d_i_z;
                if (!found || distance_sq < best_distance_sq) {
                    found = 1;
                    best_distance_sq = distance_sq;
                    best = candidate;
                    best_r = grid->r[corner_i_r] + candidate.d_i_r * grid->d_r;
                    best_z = grid->z[corner_i_z] + candidate.d_i_z * grid->d_z;
                }
            }
            if (!found) {
                continue;
            }

            if (best.hessian_determinant > 0.0) {
                if (!*axis_found || best.psi > *axis_psi) {
                    *axis_found = 1;
                    *axis_r = best_r;
                    *axis_z = best_z;
                    *axis_psi = best.psi;
                }
            } else if (best.hessian_determinant < 0.0) {
                xpt_r[*xpt_n] = best_r;
                xpt_z[*xpt_n] = best_z;
                xpt_psi[*xpt_n] = best.psi;
                (*xpt_n)++;
                if (*xpt_n >= n_xpt_r) {
                    return RTGSFIT_ERR_TOO_MANY_XPTS;
                }
            }
        }
    }
    return RTGSFIT_OK;
}

int32_t magnetic_nulls_is_core_side(double r,
                                    double z,
                                    double axis_r,
                                    double axis_z,
                                    const double xpt_r[],
                                    const double xpt_z[],
                                    int32_t xpt_n)
{
    for (int32_t i = 0; i < xpt_n; i++) {
        const double dot =
            (r - xpt_r[i]) * (axis_r - xpt_r[i]) + (z - xpt_z[i]) * (axis_z - xpt_z[i]);
        if (dot < 0.0) {
            return 0;
        }
    }
    return 1;
}

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
                                   int32_t n_keep)
{
    RTGSFIT_CHECK_SIZE(n_xpt_z, n_xpt_r);
    RTGSFIT_CHECK_SIZE(n_xpt_psi, n_xpt_r);
    if (*xpt_n < 0 || *xpt_n > n_xpt_r) {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    RTGSFIT_CHECK_MIN_SIZE(n_keep, *xpt_n);

    const int32_t n = *xpt_n;
    for (int32_t i = 0; i < n; i++) {
        keep[i] = 1;
    }
    /* Every x-point, including ones already removed, can hide the others. */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            if (j == i) {
                continue;
            }
            const double dot = (xpt_r[j] - xpt_r[i]) * (axis_r - xpt_r[i]) +
                               (xpt_z[j] - xpt_z[i]) * (axis_z - xpt_z[i]);
            if (dot < 0.0) {
                keep[j] = 0;
            }
        }
    }
    int32_t n_kept = 0;
    for (int32_t i = 0; i < n; i++) {
        if (keep[i]) {
            xpt_r[n_kept] = xpt_r[i];
            xpt_z[n_kept] = xpt_z[i];
            xpt_psi[n_kept] = xpt_psi[i];
            n_kept++;
        }
    }
    *xpt_n = n_kept;
    return RTGSFIT_OK;
}

static inline void swap(double *a, double *b)
{
    const double tmp = *a;
    *a = *b;
    *b = tmp;
}

int32_t magnetic_nulls_sort_xpts(double xpt_r[],
                                 int32_t n_xpt_r,
                                 double xpt_z[],
                                 int32_t n_xpt_z,
                                 double xpt_psi[],
                                 int32_t n_xpt_psi,
                                 int32_t xpt_n)
{
    RTGSFIT_CHECK_SIZE(n_xpt_z, n_xpt_r);
    RTGSFIT_CHECK_SIZE(n_xpt_psi, n_xpt_r);
    if (xpt_n < 0 || xpt_n > n_xpt_r) {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    /* Selection sort: n is small and the order of equal values is then reproducible. */
    for (int32_t i = 0; i < xpt_n - 1; i++) {
        int32_t i_max = i;
        for (int32_t j = i + 1; j < xpt_n; j++) {
            if (xpt_psi[j] > xpt_psi[i_max]) {
                i_max = j;
            }
        }
        if (i_max != i) {
            swap(&xpt_r[i], &xpt_r[i_max]);
            swap(&xpt_z[i], &xpt_z[i_max]);
            swap(&xpt_psi[i], &xpt_psi[i_max]);
        }
    }
    return RTGSFIT_OK;
}
