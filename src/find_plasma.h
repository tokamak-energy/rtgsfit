// find_plasma.h
#ifndef FIND_PLASMA_H
#define FIND_PLASMA_H

#include <stdint.h>

// Finds null points: the highest-flux o-point (the magnetic axis) and all
// x-points. opt_r/opt_z/opt_flux are single values, and opt_n is set to 1 if an
// o-point was found and 0 otherwise.
// Returns 0 on success, nonzero on error.
int find_nulls(double *flux, double *opt_r, double *opt_z, double *opt_flux,
               int32_t *opt_n, double *xpt_r, double *xpt_z, double *xpt_flux,
               int32_t *xpt_n);

// Filters x-points: removes those "behind" other x-points relative to the axis.
void filter_xpts(double *xpt_r, double *xpt_z, double *xpt_flux, int32_t *xpt_n,
                 double r_mag_axis, double z_mag_axis);

// Sorts x-points in descending order of flux value.
void sort_xpts(double *xpt_r, double *xpt_z, double *xpt_flux, int32_t xpt_n);

// Returns 1 if (r_grid, z_grid) is on the magnetic-axis side of all x-points,
// else 0.
int is_core_side_of_xpoint(double r_grid, double z_grid, double r_mag_axis,
                           double z_mag_axis, double *xpt_r, double *xpt_z,
                           int32_t xpt_n);

// Flood-fill plasma core mask.
// mask[i] set to 1 inside core (connected to mag-axis seed), else 0.
// Returns 0 on success; otherwise returns an error flag/code (see constants.h).
int flood_fill_plasma_core(int32_t *mask, double *flux_total,
                           double flux_boundary, double r_mag_axis,
                           double z_mag_axis, double *xpt_r, double *xpt_z,
                           int32_t xpt_n);

#endif // FIND_PLASMA_H