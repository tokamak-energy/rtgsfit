#include "plasma_boundary.h"

#include <float.h>
#include <math.h>

#include "magnetic_nulls.h"

int32_t plasma_boundary_limiter_psi(const double psi[],
                                    const int32_t n_psi,
                                    const double limiter_r[],
                                    const int32_t n_limiter_r,
                                    const double limiter_z[],
                                    const int32_t n_limiter_z,
                                    const int32_t limiter_grid_index[],
                                    const int32_t n_index,
                                    const double limiter_grid_weight[],
                                    const int32_t n_weight,
                                    double axis_r,
                                    double axis_z,
                                    const double xpt_r[],
                                    const int32_t n_xpt_r,
                                    const double xpt_z[],
                                    const int32_t n_xpt_z,
                                    int32_t xpt_n,
                                    double *limiter_psi)
{
    const int32_t n_limiter = n_limiter_r;
    RTGSFIT_CHECK_SIZE(n_limiter_z, n_limiter);
    RTGSFIT_CHECK_SIZE(n_index, (int64_t)n_limiter * PLASMA_BOUNDARY_N_INTERP);
    RTGSFIT_CHECK_SIZE(n_weight, (int64_t)n_limiter * PLASMA_BOUNDARY_N_INTERP);
    RTGSFIT_CHECK_SIZE(n_xpt_z, n_xpt_r);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(limiter_r, n_limiter_r);
    RTGSFIT_CHECK_ARRAY(limiter_z, n_limiter_z);
    RTGSFIT_CHECK_ARRAY(limiter_grid_index, n_index);
    RTGSFIT_CHECK_ARRAY(limiter_grid_weight, n_weight);
    RTGSFIT_CHECK_ARRAY(xpt_r, n_xpt_r);
    RTGSFIT_CHECK_ARRAY(xpt_z, n_xpt_z);
    if (!limiter_psi)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    if (xpt_n < 0 || xpt_n > n_xpt_r)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }

    double psi_max = -DBL_MAX;
    for (int32_t i = 0; i < n_limiter; i++)
    {
        int32_t is_core_side = 0;
        const int32_t status =
            magnetic_nulls_is_core_side(limiter_r[i], limiter_z[i], axis_r, axis_z, xpt_r, n_xpt_r,
                                        xpt_z, n_xpt_z, xpt_n, &is_core_side);
        if (status != RTGSFIT_OK)
        {
            return status;
        }
        if (!is_core_side)
        {
            continue;
        }
        double psi_i = 0.0;
        for (int32_t k = 0; k < PLASMA_BOUNDARY_N_INTERP; k++)
        {
            const int32_t i_grid = limiter_grid_index[i * PLASMA_BOUNDARY_N_INTERP + k];
            if (i_grid < 0 || i_grid >= n_psi)
            {
                return RTGSFIT_ERR_ARGUMENT_SIZE;
            }
            psi_i += limiter_grid_weight[i * PLASMA_BOUNDARY_N_INTERP + k] * psi[i_grid];
        }
        if (psi_i > psi_max)
        {
            psi_max = psi_i;
        }
    }
    *limiter_psi = psi_max;
    return RTGSFIT_OK;
}

int32_t plasma_boundary_select(double axis_psi,
                               double limiter_psi,
                               const double xpt_psi[],
                               const int32_t n_xpt_psi,
                               int32_t xpt_n,
                               double xpt_fraction,
                               double tolerance,
                               double *boundary_psi,
                               int32_t *xpt_diverted)
{
    RTGSFIT_CHECK_ARRAY(xpt_psi, n_xpt_psi);
    if (!boundary_psi || !xpt_diverted)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    if (xpt_n < 0 || xpt_n > n_xpt_psi)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }

    *boundary_psi = limiter_psi;
    *xpt_diverted = 0;
    if (xpt_n > 0)
    {
        /* Step slightly inside the separatrix so the flood fill does not leak through the x-point. */
        const double diverted_psi = xpt_fraction * xpt_psi[0] + (1.0 - xpt_fraction) * axis_psi;
        if (diverted_psi > limiter_psi)
        {
            *boundary_psi = diverted_psi;
            *xpt_diverted = 1;
        }
    }
    if (fabs(*boundary_psi - axis_psi) < tolerance)
    {
        return RTGSFIT_ERR_AXIS_EQUALS_BOUNDARY;
    }
    if (*boundary_psi > axis_psi)
    {
        return RTGSFIT_ERR_BOUNDARY_ABOVE_AXIS;
    }
    return RTGSFIT_OK;
}

static int32_t nearest_index(const double values[], int32_t n, double target, int32_t *index)
{
    RTGSFIT_CHECK_MIN_SIZE(n, 1);
    RTGSFIT_CHECK_ARRAY(values, n);
    int32_t i_best = 0;
    double distance_best = fabs(values[0] - target);
    for (int32_t i = 1; i < n; i++)
    {
        const double distance = fabs(values[i] - target);
        if (distance < distance_best)
        {
            distance_best = distance;
            i_best = i;
        }
    }
    *index = i_best;
    return RTGSFIT_OK;
}

/* Data shared by every core test during one flood fill. */
struct core_test
{
    const struct rtgsfit_grid *grid;
    const double *psi;
    const int32_t *mask_limiter;
    double boundary_psi;
    double axis_r;
    double axis_z;
    const double *xpt_r;
    int32_t n_xpt_r;
    const double *xpt_z;
    int32_t n_xpt_z;
    int32_t xpt_n;
};

static int32_t is_core(const struct core_test *test, int32_t i_r, int32_t i_z, int32_t *inside)
{
    *inside = 0;
    const int32_t i = grid_index(test->grid, i_r, i_z);
    if (test->psi[i] <= test->boundary_psi || !test->mask_limiter[i])
    {
        return RTGSFIT_OK;
    }
    return magnetic_nulls_is_core_side(test->grid->r[i_r], test->grid->z[i_z], test->axis_r,
                                       test->axis_z, test->xpt_r, test->n_xpt_r, test->xpt_z,
                                       test->n_xpt_z, test->xpt_n, inside);
}

int32_t plasma_boundary_flood_fill(const struct rtgsfit_grid *grid,
                                   const double psi[],
                                   const int32_t n_psi,
                                   const int32_t mask_limiter[],
                                   const int32_t n_mask_limiter,
                                   double boundary_psi,
                                   double axis_r,
                                   double axis_z,
                                   const double xpt_r[],
                                   const int32_t n_xpt_r,
                                   const double xpt_z[],
                                   const int32_t n_xpt_z,
                                   int32_t xpt_n,
                                   int32_t mask[],
                                   const int32_t n_mask,
                                   int32_t queue[],
                                   const int32_t n_queue)
{
    int32_t status = grid_check(grid);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    const int32_t n_grid = grid_n_grid(grid);
    RTGSFIT_CHECK_SIZE(n_psi, n_grid);
    RTGSFIT_CHECK_SIZE(n_mask_limiter, n_grid);
    RTGSFIT_CHECK_SIZE(n_mask, n_grid);
    RTGSFIT_CHECK_SIZE(n_queue, n_grid);
    RTGSFIT_CHECK_SIZE(n_xpt_z, n_xpt_r);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(mask_limiter, n_mask_limiter);
    RTGSFIT_CHECK_ARRAY(mask, n_mask);
    RTGSFIT_CHECK_ARRAY(queue, n_queue);
    RTGSFIT_CHECK_ARRAY(xpt_r, n_xpt_r);
    RTGSFIT_CHECK_ARRAY(xpt_z, n_xpt_z);
    if (xpt_n < 0 || xpt_n > n_xpt_r)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }

    const struct core_test test = {grid,  psi,     mask_limiter, boundary_psi, axis_r, axis_z,
                                   xpt_r, n_xpt_r, xpt_z,        n_xpt_z,      xpt_n};
    int32_t i_r_seed, i_z_seed, inside;
    status = nearest_index(grid->r, grid->n_r, axis_r, &i_r_seed);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    status = nearest_index(grid->z, grid->n_z, axis_z, &i_z_seed);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    status = is_core(&test, i_r_seed, i_z_seed, &inside);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    if (!inside)
    {
        return RTGSFIT_ERR_AXIS_OUTSIDE_CORE;
    }

    for (int32_t i = 0; i < n_grid; i++)
    {
        mask[i] = 0;
    }
    /* Breadth-first search; each grid point enters the queue at most once. */
    int32_t head = 0;
    int32_t tail = 0;
    const int32_t i_seed = grid_index(grid, i_r_seed, i_z_seed);
    mask[i_seed] = 1;
    queue[tail++] = i_seed;
    while (head < tail)
    {
        const int32_t i = queue[head++];
        const int32_t i_z = i / grid->n_r;
        const int32_t i_r = i - i_z * grid->n_r;
        const int32_t neighbour_i_r[4] = {i_r - 1, i_r + 1, i_r, i_r};
        const int32_t neighbour_i_z[4] = {i_z, i_z, i_z - 1, i_z + 1};

        for (int32_t k = 0; k < 4; k++)
        {
            const int32_t j_r = neighbour_i_r[k];
            const int32_t j_z = neighbour_i_z[k];
            if (j_r < 0 || j_r >= grid->n_r || j_z < 0 || j_z >= grid->n_z)
            {
                continue;
            }
            const int32_t j = grid_index(grid, j_r, j_z);
            if (!mask[j])
            {
                status = is_core(&test, j_r, j_z, &inside);
                if (status != RTGSFIT_OK)
                {
                    return status;
                }
                if (inside)
                {
                    mask[j] = 1;
                    queue[tail++] = j;
                }
            }
        }
    }
    return RTGSFIT_OK;
}

int32_t plasma_boundary_normalise(const double psi[],
                                  const int32_t n_psi,
                                  const int32_t mask[],
                                  const int32_t n_mask,
                                  const int32_t mask_limiter[],
                                  const int32_t n_mask_limiter,
                                  double axis_psi,
                                  double boundary_psi,
                                  double psi_norm[],
                                  const int32_t n_psi_norm)
{
    RTGSFIT_CHECK_SIZE(n_mask, n_psi);
    RTGSFIT_CHECK_SIZE(n_mask_limiter, n_psi);
    RTGSFIT_CHECK_SIZE(n_psi_norm, n_psi);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(mask, n_mask);
    RTGSFIT_CHECK_ARRAY(mask_limiter, n_mask_limiter);
    RTGSFIT_CHECK_ARRAY(psi_norm, n_psi_norm);

    const double inv_psi_difference = 1.0 / (boundary_psi - axis_psi);
    for (int32_t i = 0; i < n_psi; i++)
    {
        if (mask[i] && mask_limiter[i])
        {
            psi_norm[i] = (psi[i] - axis_psi) * inv_psi_difference;
        }
        else
        {
            psi_norm[i] = 1.0;
        }
    }
    return RTGSFIT_OK;
}
