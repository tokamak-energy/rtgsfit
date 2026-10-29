/*
 * The preshot: every quantity RT-GSFit needs that can be computed before the pulse.
 *
 * The sizes are defined as enum constants in the generated `preshot_sizes.h` and the data in
 * the generated `preshot_data.c`. Both are written by `rtgsfit.database_writers.c_source` from
 * the same preshot that is stored in MDSplus under PRESHOT. The C name `GROUP_NODE` corresponds
 * to the MDSplus node `PRESHOT.GROUP:NODE`. See the manual for the definition of each quantity.
 *
 * Matrices are stored row-major, `[rows][columns]`, and flattened.
 */
#ifndef RTGSFIT_PRESHOT_H
#define RTGSFIT_PRESHOT_H

#include <stdint.h>

#include "current_basis.h"
#include "plasma_boundary.h"
#include "preshot_sizes.h"

/* C forbids zero-length arrays, so optional dimensions are allocated at least one element. */
#define PRESHOT_AT_LEAST_ONE(n) ((n) > 0 ? (n) : 1)

/* GRID: the computational grid */
extern const double GRID_R[N_R];                  /* [m] */
extern const double GRID_Z[N_Z];                  /* [m] */
extern const double GRID_D_R;                     /* [m] */
extern const double GRID_D_Z;                     /* [m] */
extern const double GRID_R_FLAT[N_GRID];          /* [m] R at each grid point */
extern const double GRID_Z_FLAT[N_GRID];          /* [m] Z at each grid point */
extern const double GRID_INV_R[N_GRID];           /* [1/m] 1 / R */
extern const int32_t GRID_MASK_LIMITER[N_GRID];   /* 1 inside the limiter, 0 on the edge */
extern const double GRID_G_COIL[N_GRID * N_COIL]; /* psi per ampere of coil current */
extern const double GRID_G_PASSIVE[N_GRID * PRESHOT_AT_LEAST_ONE(N_PASSIVE)];

/* EDGE: the edge of the grid, in edge order */
extern const double EDGE_R[N_EDGE];
extern const double EDGE_Z[N_EDGE];
extern const double EDGE_INV_MU0_R[N_EDGE];
extern const double EDGE_G_EDGE[N_EDGE * N_EDGE];

/* LIMITER: the limiter points and their bilinear interpolation stencil on the grid */
extern const double LIMITER_R[N_LIMITER];
extern const double LIMITER_Z[N_LIMITER];
extern const int32_t LIMITER_GRID_INDEX[N_LIMITER * PLASMA_BOUNDARY_N_INTERP];
extern const double LIMITER_GRID_WEIGHT[N_LIMITER * PLASMA_BOUNDARY_N_INTERP];

/* POISSON: banded LU factors of the finite-difference Grad-Shafranov operator */
/* (2 pi / INFO:PSI_TO_WB) mu0 R d_z^2; rhs = -POISSON_RHS_SCALE * j_phi */
extern const double POISSON_RHS_SCALE[N_GRID];
extern const double POISSON_LOWER_BAND[N_GRID * N_R];
extern const double POISSON_UPPER_BAND[N_GRID * (N_R + 2)];
extern const int32_t POISSON_PERMUTATION[N_GRID];

/* CONSTRAINT: from the PCS measurements to the constraints */
extern const double CONSTRAINT_G_COIL[N_CONSTRAINT * N_COIL];
extern const double CONSTRAINT_REPLACEMENT[N_CONSTRAINT * N_MEAS];

/* FIT: the weighted least-squares problem, fit_matrix fit_dof = fit_target */
extern const double FIT_WEIGHT[N_FIT_TARGET];
extern const double FIT_A_GRID[N_GRID * N_FIT_TARGET];
extern const double FIT_A_PASSIVE[PRESHOT_AT_LEAST_ONE(N_PASSIVE) * N_FIT_TARGET];

/* SETTINGS: scalar algorithm settings */
extern const double SETTINGS_XPT_FRACTION;
extern const double SETTINGS_TOLERANCE;
extern const double SETTINGS_IP_MIN; /* [A] smallest ip with a current centroid */

_Static_assert(N_GRID == N_R * N_Z, "N_GRID must equal N_R * N_Z");
_Static_assert(N_EDGE == 2 * N_R + 2 * N_Z - 4, "N_EDGE must equal 2 N_R + 2 N_Z - 4");
_Static_assert(N_R >= 4 && N_Z >= 4, "the grid needs at least 4 points in each direction");
_Static_assert(N_CONSTRAINT == N_FLUX_LOOP + N_BP_PROBE + N_ROGOWSKI,
               "N_CONSTRAINT must equal N_FLUX_LOOP + N_BP_PROBE + N_ROGOWSKI");
_Static_assert(N_FIT_TARGET == N_CONSTRAINT + N_REG,
               "N_FIT_TARGET must equal N_CONSTRAINT + N_REG");
_Static_assert((int)N_PLASMA == (int)CURRENT_BASIS_N_DOF, "N_PLASMA must match the current basis");
_Static_assert(N_FIT_DOF == N_PLASMA + N_PASSIVE, "N_FIT_DOF must equal N_PLASMA + N_PASSIVE");
_Static_assert(N_FIT_TARGET >= N_FIT_DOF, "the least-squares problem must not be under-determined");
_Static_assert(N_LIMITER >= 1 && N_COIL >= 1 && N_MEAS >= 1 && N_XPT_MAX >= 1,
               "N_LIMITER, N_COIL, N_MEAS and N_XPT_MAX must be at least 1");

#endif /* RTGSFIT_PRESHOT_H */
