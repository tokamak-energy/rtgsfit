/*
 * RT-GSFit public interface.
 *
 * This is the only header a caller needs: the Plasma Control System, Simulink, or the Python
 * replay. All names follow `documentation/manual.md`. Every array argument is followed by its
 * length, and that length is checked against the compiled-in preshot sizes on every call.
 */
#ifndef RTGSFIT_H
#define RTGSFIT_H

#include <stdint.h>

/*
 * Status flags returned by `rtgsfit()`. They are powers of two, so several can be reported at
 * once with a bitwise OR. Values 1 to 32 are unchanged from RT-GSFit v1 (`lcfs_err_code`).
 */
enum rtgsfit_status
{
    RTGSFIT_OK = 0,
    RTGSFIT_ERR_AXIS_EQUALS_BOUNDARY = 1, /* |boundary_psi - axis_psi| < SETTINGS_TOLERANCE */
    RTGSFIT_ERR_NO_AXIS = 2,              /* no o-point (magnetic axis) was found */
    RTGSFIT_ERR_TOO_MANY_XPTS = 4,        /* N_XPT_MAX x-points were found before the scan ended */
    RTGSFIT_ERR_TOO_MANY_OPTS = 8,        /* reserved, kept for compatibility with v1 */
    RTGSFIT_ERR_AXIS_OUTSIDE_CORE = 16,   /* the grid point nearest the axis is not in the core */
    RTGSFIT_ERR_BOUNDARY_ABOVE_AXIS = 32, /* boundary_psi > axis_psi */
    RTGSFIT_ERR_ARGUMENT_SIZE = 64,       /* an array length does not match the preshot */
    RTGSFIT_ERR_LAPACK = 128,             /* LAPACK dgelss returned a non-zero `info` */
    RTGSFIT_ERR_FIT_INVALID_ARGUMENT = 256,
    RTGSFIT_ERR_FIT_NO_CONVERGENCE = 512,
    RTGSFIT_ERR_FIT_NONFINITE_INPUT = 1024,
    RTGSFIT_ERR_FIT_RESOURCE_FAILURE = 2048,
};

/** Classify LAPACK info as success, invalid argument or nonconvergence; retain raw info separately. */
int32_t rtgsfit_status_from_lapack_info(const int32_t info);

/*
 * Run one RT-GSFit iteration. It fits the plasma and passive currents to the measurements, solves
 * the free-boundary Poisson problem for psi, and then locates the magnetic axis, the x-points and
 * the plasma boundary.
 *
 * Inputs:
 *   meas                      PCS measurements, ordered as PRESHOT.MEAS:NAMES
 *   coil_current              PF coil currents [A], ordered as PRESHOT.COIL:NAMES
 *
 * Inputs and outputs (the result of one iteration is the starting guess of the next):
 *   psi_norm                  normalised poloidal flux on the grid: 0 on the axis, 1 on the boundary
 *   mask                      1 inside the plasma core, 0 elsewhere
 *
 * Outputs:
 *   psi                       total poloidal flux [Wb]; Wb/rad only with legacy analytic preshots
 *   j_phi                     toroidal plasma current density on the grid [A/m^2]
 *   fit_dof                   fitted DOF values, not weighted: plasma first, then passive
 *   constraint_measured       unweighted constraint values after sensor replacement
 *   constraint_reconstructed  unweighted constraints predicted by the fitted equilibrium
 *   xpt_r, xpt_z, xpt_psi     x-points sorted by descending psi; only the first `*xpt_n` are valid
 *   xpt_n                     number of valid x-points
 *   xpt_diverted              1 if an x-point sets the boundary, 0 if the plasma is limited
 *   ip                        plasma current [A]
 *   axis_psi, axis_r, axis_z  magnetic axis [psi units, m, m]
 *   boundary_psi              psi on the plasma boundary [psi units]
 *   current_centroid_r, current_centroid_z   current centroid [m]
 *   chi_sq_magnetic           weighted squared residual over the constraint rows of fit_target
 *   chi_sq_total              weighted squared residual over every row of fit_target
 *   lapack_info               `info` returned by LAPACK dgelss
 *
 * Returns a bitwise OR of `enum rtgsfit_status`. When the status is non-zero, `psi_norm` and
 * `mask` are left unchanged, so the next iteration restarts from the last good equilibrium.
 */
int32_t rtgsfit(const double meas[],
                const int32_t n_meas,
                const double coil_current[],
                const int32_t n_coil_current,
                double psi_norm[],
                const int32_t n_psi_norm,
                int32_t mask[],
                const int32_t n_mask,
                double psi[],
                const int32_t n_psi,
                double j_phi[],
                const int32_t n_j_phi,
                double fit_dof[],
                const int32_t n_fit_dof,
                double constraint_measured[],
                const int32_t n_constraint_measured,
                double constraint_reconstructed[],
                const int32_t n_constraint_reconstructed,
                double xpt_r[],
                const int32_t n_xpt_r,
                double xpt_z[],
                const int32_t n_xpt_z,
                double xpt_psi[],
                const int32_t n_xpt_psi,
                int32_t *xpt_n,
                int32_t *xpt_diverted,
                double *ip,
                double *axis_psi,
                double *axis_r,
                double *axis_z,
                double *boundary_psi,
                double *current_centroid_r,
                double *current_centroid_z,
                double *chi_sq_magnetic,
                double *chi_sq_total,
                int32_t *lapack_info);

/*
 * The compiled-in problem sizes. A caller can use them to allocate its buffers, and to check that
 * the library was built from the preshot it expects. `rtgsfit_get_sizes()` fills `sizes` in the
 * order of `enum rtgsfit_size_index`.
 */
enum rtgsfit_size_index
{
    RTGSFIT_SIZE_N_R = 0,
    RTGSFIT_SIZE_N_Z,
    RTGSFIT_SIZE_N_GRID,
    RTGSFIT_SIZE_N_EDGE,
    RTGSFIT_SIZE_N_LIMITER,
    RTGSFIT_SIZE_N_COIL,
    RTGSFIT_SIZE_N_MEAS,
    RTGSFIT_SIZE_N_FLUX_LOOP,
    RTGSFIT_SIZE_N_BP_PROBE,
    RTGSFIT_SIZE_N_ROGOWSKI,
    RTGSFIT_SIZE_N_CONSTRAINT,
    RTGSFIT_SIZE_N_REG,
    RTGSFIT_SIZE_N_FIT_TARGET,
    RTGSFIT_SIZE_N_PLASMA,
    RTGSFIT_SIZE_N_PASSIVE,
    RTGSFIT_SIZE_N_FIT_DOF,
    RTGSFIT_SIZE_N_XPT_MAX,
    RTGSFIT_SIZE_COUNT
};

/* Returns RTGSFIT_ERR_ARGUMENT_SIZE unless n_sizes == RTGSFIT_SIZE_COUNT. */
int32_t rtgsfit_get_sizes(int32_t sizes[], const int32_t n_sizes);

/*
 * CPU time spent in each section, in microseconds, summed since the last
 * `rtgsfit_timing_reset()`. The times are only recorded when the library was built with
 * `TIMING=1`; otherwise every value is zero.
 */
enum rtgsfit_timer
{
    RTGSFIT_TIMER_CONSTRAINTS = 0, /* sensor replacement and coil subtraction */
    RTGSFIT_TIMER_BASIS,           /* plasma current basis functions */
    RTGSFIT_TIMER_FIT_MATRIX,      /* plasma columns of the least-squares matrix */
    RTGSFIT_TIMER_LEAST_SQUARES,   /* LAPACK dgelss */
    RTGSFIT_TIMER_CURRENT,         /* j_phi, ip and current centroid */
    RTGSFIT_TIMER_RESIDUAL,        /* fit_prediction, chi squared and reconstructed constraints */
    RTGSFIT_TIMER_POISSON,         /* free-boundary Poisson solve for the plasma psi */
    RTGSFIT_TIMER_EXTERNAL_PSI,    /* coil and passive psi */
    RTGSFIT_TIMER_NULLS,           /* magnetic axis and x-points */
    RTGSFIT_TIMER_BOUNDARY,        /* limiter psi and boundary selection */
    RTGSFIT_TIMER_FLOOD_FILL,      /* plasma core mask */
    RTGSFIT_TIMER_NORMALISE,       /* psi_norm */
    RTGSFIT_TIMER_TOTAL,
    RTGSFIT_TIMER_COUNT
};

void rtgsfit_timing_reset(void);

/* Returns RTGSFIT_ERR_ARGUMENT_SIZE unless n_time_us == RTGSFIT_TIMER_COUNT. */
int32_t rtgsfit_timing_get(double time_us[], const int32_t n_time_us);

#endif /* RTGSFIT_H */
