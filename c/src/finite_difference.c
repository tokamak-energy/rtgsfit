#include "finite_difference.h"

int32_t finite_difference_d_dz(
    const struct rtgsfit_grid *grid, const double f[], int32_t n_f, double df_dz[], int32_t n_df_dz)
{
    const int32_t status = grid_check(grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    RTGSFIT_CHECK_SIZE(n_f, grid_n_grid(grid));
    RTGSFIT_CHECK_SIZE(n_df_dz, grid_n_grid(grid));

    const int32_t n_r = grid->n_r;
    const int32_t n_z = grid->n_z;
    const double inv_2_d_z = 1.0 / (2.0 * grid->d_z);

    for (int32_t i_z = 0; i_z < n_z; i_z++) {
        for (int32_t i_r = 0; i_r < n_r; i_r++) {
            const int32_t i = grid_index(grid, i_r, i_z);
            if (i_z == 0) {
                df_dz[i] = (-3.0 * f[i] + 4.0 * f[i + n_r] - f[i + 2 * n_r]) * inv_2_d_z;
            } else if (i_z == n_z - 1) {
                df_dz[i] = (3.0 * f[i] - 4.0 * f[i - n_r] + f[i - 2 * n_r]) * inv_2_d_z;
            } else {
                df_dz[i] = (f[i + n_r] - f[i - n_r]) * inv_2_d_z;
            }
        }
    }
    return RTGSFIT_OK;
}

/* One-sided second-order difference towards the interior, without the 1/(grid spacing). */
static inline double inward_difference(double f_1, double f_2)
{
    return 2.0 * f_1 - 0.5 * f_2;
}

int32_t finite_difference_edge_d_dn(
    const struct rtgsfit_grid *grid, const double f[], int32_t n_f, double df_dn[], int32_t n_df_dn)
{
    const int32_t status = grid_check(grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    RTGSFIT_CHECK_SIZE(n_f, grid_n_grid(grid));
    RTGSFIT_CHECK_SIZE(n_df_dn, grid_n_edge(grid));

    const int32_t n_r = grid->n_r;
    const int32_t n_z = grid->n_z;
    const double d_z_over_d_r = grid->d_z / grid->d_r;
    const double d_r_over_d_z = grid->d_r / grid->d_z;
    int32_t i_edge = 0;

    /* Left side: i_r = 0, Z increasing, including both corners. */
    for (int32_t i_z = 0; i_z < n_z; i_z++) {
        df_dn[i_edge++] = d_z_over_d_r * inward_difference(f[grid_index(grid, 1, i_z)],
                                                           f[grid_index(grid, 2, i_z)]);
    }
    /* Top side: i_z = n_z - 1, R increasing, excluding the corners. */
    for (int32_t i_r = 1; i_r < n_r - 1; i_r++) {
        df_dn[i_edge++] = d_r_over_d_z * inward_difference(f[grid_index(grid, i_r, n_z - 2)],
                                                           f[grid_index(grid, i_r, n_z - 3)]);
    }
    /* Right side: i_r = n_r - 1, Z decreasing, including both corners. */
    for (int32_t i_z = n_z - 1; i_z >= 0; i_z--) {
        df_dn[i_edge++] = d_z_over_d_r * inward_difference(f[grid_index(grid, n_r - 2, i_z)],
                                                           f[grid_index(grid, n_r - 3, i_z)]);
    }
    /* Bottom side: i_z = 0, R decreasing, excluding the corners. */
    for (int32_t i_r = n_r - 2; i_r >= 1; i_r--) {
        df_dn[i_edge++] = d_r_over_d_z * inward_difference(f[grid_index(grid, i_r, 1)],
                                                           f[grid_index(grid, i_r, 2)]);
    }
    return RTGSFIT_OK;
}
