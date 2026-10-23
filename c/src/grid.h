/*
 * The rectangular (R, Z) computational grid.
 *
 * Grid arrays are stored row-major with R varying fastest:
 *     index = i_z * n_r + i_r
 */
#ifndef RTGSFIT_GRID_H
#define RTGSFIT_GRID_H

#include <stdint.h>

#include "argument_check.h"
#include "rtgsfit.h"

/* Describes the grid. It does not own `r` and `z`; they normally point into the preshot. */
struct rtgsfit_grid {
    int32_t n_r;
    int32_t n_z;
    double d_r;
    double d_z;
    const double *r; /* [n_r] */
    const double *z; /* [n_z] */
};

static inline int32_t grid_n_grid(const struct rtgsfit_grid *grid)
{
    return grid->n_r * grid->n_z;
}

/* Number of points on the edge of the grid, with each corner counted once. */
static inline int32_t grid_n_edge(const struct rtgsfit_grid *grid)
{
    return 2 * grid->n_r + 2 * grid->n_z - 4;
}

static inline int32_t grid_index(const struct rtgsfit_grid *grid, int32_t i_r, int32_t i_z)
{
    return i_z * grid->n_r + i_r;
}

/* The grid needs at least 4 points in each direction for the finite-difference stencils. */
static inline int32_t grid_check(const struct rtgsfit_grid *grid)
{
    if (grid->n_r < 4 || grid->n_z < 4 || !(grid->d_r > 0.0) || !(grid->d_z > 0.0)) {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    return RTGSFIT_OK;
}

#endif /* RTGSFIT_GRID_H */
