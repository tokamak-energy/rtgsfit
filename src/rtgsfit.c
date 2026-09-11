#include "constants.h"
#include "gradient.h"
#include "rtgsfit.h"
#include "poisson_solver.h"
#include "find_x_point.h"
#include "find_plasma.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#define N_MEAS_NO_REG (N_BP_PROBES + N_FLUX_LOOPS + N_ROGOWSKI_COILS)

#ifdef ENABLE_RT_TIMING

static inline uint64_t thread_cpu_ns(void)
{
    struct timespec ts;
    int ret = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    assert(ret == 0);
    return (uint64_t)ts.tv_sec * 1000000000ull + ts.tv_nsec;
}

enum {
    T_MEAS_PREP = 0,
    T_BASIS,
    T_MEAS_MATRIX,
    T_WEIGHTING,
    T_LS_FIT,
    T_SOURCE,
    T_MODEL_MEAS,
    T_CHI2,
    T_POISSON,
    T_COIL_FLUX,
    T_VESSEL_FLUX,
    T_XPTS_AND_AXIS,
    T_XPTS_SORT,
    T_LIMITER,
    T_LCFS,
    T_INSIDE,
    T_NORMALISE,
    T_TOTAL,
    T_NTIMERS
};

static uint64_t timing_acc[T_NTIMERS];
static uint64_t timing_t0;

#define TSTART()        do { timing_t0 = thread_cpu_ns(); } while (0)
#define TACC(idx)       do { timing_acc[(idx)] += (thread_cpu_ns() - timing_t0); } while (0)

void rtgsfit_timing_reset(void)
{
    for (int i = 0; i < T_NTIMERS; i++) timing_acc[i] = 0;
}

void rtgsfit_timing_dump(void)
{
    static const char *names[T_NTIMERS] = {
        "meas_prep",
        "basis",
        "meas_matrix",
        "weighting",
        "ls_fit",
        "source",
        "model_meas",
        "chi2",
        "poisson",
        "coil_flux",
        "vessel_flux",
        "xpts_and_axis",
        "xpts_sort",
        "limiter",
        "lcfs",
        "inside",
        "normalise",
        "total"
    };

    for (int i = 0; i < T_NTIMERS; i++) {
        printf("%-14s : %10.3f us\n", names[i], (double)timing_acc[i] * 1e-3);
    }
}

/**
 * Copy timing counters into the provided buffer in microseconds.
 *
 * @param out_us Output buffer that receives timing values.
 * @param n Maximum number of entries to copy.
 */
void rtgsfit_timing_get(double *out_us, int n)
{
    int count = (n < T_NTIMERS) ? n : T_NTIMERS;
    for (int i = 0; i < count; i++) {
        out_us[i] = (double)timing_acc[i] * 1e-3;
    }
}

#else

/* When ENABLE_RT_TIMING is not defined, timing macros are intentional no-ops. */
#define TSTART()        do {} while (0)
#define TACC(idx)       do {} while (0)

#endif // ENABLE_RT_TIMING

int32_t max_idx(int32_t n_arr, double* arr)
{
    int i_arr;
    int i_max = 0;
    double arr_max = arr[0];

    for (i_arr=1; i_arr<n_arr; i_arr++)
    {
        if (arr_max < arr[i_arr])
        {
            arr_max = arr[i_arr];
            i_max = i_arr;
        }
    }
    return i_max;
}

void rm_coil_from_meas(
        const double* coil_curr,
        double* meas,
        double* meas_no_coil
        )
{
    int i_meas;
    // subtract PF (vessel) contributions from measurements
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_MEAS, N_COIL, 1.0, G_MEAS_COIL,
            N_COIL, coil_curr, 1, 0.0, meas_no_coil, 1);

    for (i_meas=0; i_meas<N_MEAS; i_meas++)
    {
        meas_no_coil[i_meas] = meas[i_meas] - meas_no_coil[i_meas];
    }
}

void make_basis(
        double* flux_norm,
        int* mask,
        double* basis
        )
{
    int i_grid;

    // could use gradient from previous iteration.  should apply mask
    gradient_z(flux_norm, &basis[2*N_GRID]);  // BUXTON: third column of basis which is for delta_z


    // could use 1 - flux_norm instead of flux_norm ?????
    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        if (mask[i_grid] && MASK_LIM[i_grid])
        {
            basis[i_grid] = (1 - flux_norm[i_grid]) * R_GRID[i_grid]; // BUXTON: probably p_prime
            // PROKOPYSZYN: Why not divide by R instead of (R * mu0)?
            basis[i_grid + N_GRID] = (1 -  flux_norm[i_grid]) * INV_R_MU0[i_grid];  // BUXTON: probably ff_prime
        }
        else
        {
            basis[i_grid] = 0.0;
            basis[i_grid + N_GRID] = 0.0;
            basis[i_grid + 2*N_GRID] = 0.0;
        }
    }
}

/**
 * @brief Calculates the flux on the limiter but excludes some of the
 * limit points based on the location of the x-points.
 *
 * The ith x-point is located at (xpt_r[i], xpt_z[i]) and the
 * jth limiter point is located at (LIMIT_R[j], LIMIT_Z[j]).
 * 
 * The vector v that points from the x-point to the limiter point is
 * v = (LIMIT_R[j] - xpt_r[i], LIMIT_Z[j] - xpt_z[i]).
 * The vector w that points from the x-point to the axis is
 * w = (r_mag_axis - xpt_r[i], z_mag_axis - xpt_z[i]).
 * 
 * If the dot product of v and w is negative, then the limiter point
 * is not considered for the flux calculation.
 *
 * @param flux_total Array of total flux values.
 * @param xpt_r Array of x-point R coordinates.
 * @param xpt_z Array of x-point Z coordinates.
 * @param xpt_n Number of x-points.
 * @param r_mag_axis R coordinate of the axis.
 * @param z_mag_axis Z coordinate of the axis.
 * @return The computed flux value on the limiter after x-point-filtering.
 */
double find_flux_on_limiter_xfiltered(double flux_total[],
                                      double xpt_r[],
                                      double xpt_z[],
                                      int32_t xpt_n,
                                      double r_mag_axis,
                                      double z_mag_axis)
{

    int skip;
    double flux_limit_max;

    flux_limit_max = -DBL_MAX;

    for (int32_t i_limit = 0; i_limit < N_LIMIT; i_limit++)
    {

        skip = 0;
        for (int32_t i_xpt = 0; i_xpt < xpt_n; i_xpt++)
        {
            double dot_product = (LIMIT_R[i_limit] - xpt_r[i_xpt]) * (r_mag_axis - xpt_r[i_xpt]) +
                          (LIMIT_Z[i_limit] - xpt_z[i_xpt]) * (z_mag_axis - xpt_z[i_xpt]);
            if (dot_product < 0.0)
            {
                skip = 1;
                break;
            }
        }

        if (skip) continue;

        double flux_limit = 0.0;
        for (int32_t i_intrp = 0; i_intrp < N_INTRP; i_intrp++)
        {
            int32_t idx = i_limit * N_INTRP + i_intrp;
            flux_limit += LIMIT_WEIGHT[idx] * flux_total[LIMIT_IDX[idx]];
        }
        if (flux_limit > flux_limit_max)
        {
            flux_limit_max = flux_limit;
        }
    }
    return flux_limit_max;
}

void normalise_flux(
        double* flux_total,
        double flux_lcfs,
        double flux_axis,
        int32_t* mask,
        double* flux_norm
        )
{
    double inv_flux_diff;
    inv_flux_diff = 1.0/(flux_lcfs - flux_axis);

    for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
    {
        if (mask[i_grid] && MASK_LIM[i_grid])
        {
            flux_norm[i_grid] = (flux_total[i_grid] - flux_axis) * inv_flux_diff;
        }
        else
        {
            flux_norm[i_grid] = 1.0;
        }
    }
}

void rtgsfit(
        const double* meas_pcs, // input
        const double* coil_curr, // input
        double* flux_norm, // input/output
        int32_t* mask, // input/output
        double* flux_total, // output
        double* chi_sq_err, // output
        double* lcfs_r, // output
        double* lcfs_z, // output
        int32_t* lcfs_n, // output
        double* coef, // output
        double* flux_boundary, // output
        double* plasma_current, // output
        int32_t *lcfs_err_code, // output
        int* lapack_dgelss_info, // output
        double *meas_model, // output
        const int32_t n_meas_model, // input
        double* r_mag_axis, // output
        double* z_mag_axis,  // output
        double* mag_axis_flux, // output
        double* r_cur_centroid, // output
        double* z_cur_centroid,  // output
        double* xpt_r, // output array
        double* xpt_z, // output array
        double* xpt_flux, // output array
        const int32_t xpt_arrays_size, // input
        int32_t* xpt_n, // output integer
        int32_t* xpt_diverted // output integer
        )
{
#ifdef ENABLE_RT_TIMING
    uint64_t t_total_0 = thread_cpu_ns();
#endif // ENABLE_RT_TIMING

    assert(n_meas_model == N_MEAS);
    assert(xpt_arrays_size == N_XPT_MAX);

    // N_MEAS includes the number of regularisations.
    // n_meas_no_reg is the number of measurements after the regularisations have been removed.
    // The meas array doesn't need the regularisations as we use meas_no_coil
    // when the LAPACKE_dgelss function is called.
    // meas = SENSOR_REPLACEMENT_MATRIX * meas_pcs
    // meas contains the post-processed measurements, but doesn't include the regularisation elements.
    // The regularisation elements are included in meas_no_coil.
    // meas_pcs contains the raw measurements from the PCS that need to be post-processed
    // using the SENSOR_REPLACEMENT_MATRIX.
    double meas[N_MEAS_NO_REG];

    TSTART();
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                N_MEAS_NO_REG, N_SENS_PCS,
                1.0, SENSOR_REPLACEMENT_MATRIX, N_SENS_PCS,
                meas_pcs, 1,
                0.0, meas, 1);
    TACC(T_MEAS_PREP);

    // will this be done during compilation?
    double g_coef_meas_w[N_COEF * N_MEAS];
    memcpy(g_coef_meas_w, G_COEF_MEAS_WEIGHT, sizeof(double) * N_MEAS * N_COEF);

    // subtract PF contributions from measurements
    // Note that this also sets the regularisation elements of meas_no_coil to zero.
    // meas_no_coil is the post-processed measurments with the PF coil contributions removed
    // and includes the regularisation elements, which can be thought of as fake Rogowski coil
    // measurements.
    double meas_no_coil[N_MEAS];

    TSTART();
    rm_coil_from_meas(coil_curr, meas, meas_no_coil);
    TACC(T_MEAS_PREP);

    // make basis
    // This makes the transpose of the T_{yg} matrix in eqn. (61) of the Moret et al. (2015)
    // LIUQE paper, without the Delta_R * Delta_Z factor.
    double g_pls_grid[N_PLS * N_GRID];

    TSTART();
    make_basis(flux_norm, mask, g_pls_grid);
    TACC(T_BASIS);

    // make meas-pls matrix
    // g_coef_meas_w = g_pls_grid * G_GRID_MEAS_WEIGHT
    TSTART();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N_PLS, N_MEAS, N_GRID,
                1.0, g_pls_grid, N_GRID,
                G_GRID_MEAS_WEIGHT, N_MEAS,
                0.0, g_coef_meas_w, N_MEAS);
    TACC(T_MEAS_MATRIX);

    // form meas vectors from measurements
    TSTART();
    for (int32_t i_meas = 0; i_meas < N_MEAS; i_meas++)
    {
        meas_no_coil[i_meas] *= WEIGHT[i_meas];
    }
    TACC(T_WEIGHTING);

    double meas_no_coil_cp[N_MEAS];
    double g_coef_meas_w_orig[N_COEF * N_MEAS];
    // copy measurment to coef due to overwritting in LAPACKE_dgelss
    memcpy(meas_no_coil_cp, meas_no_coil, sizeof(double) * N_MEAS);
    memcpy(g_coef_meas_w_orig, g_coef_meas_w, sizeof(double) * N_COEF * N_MEAS);

    // fit coeff or use dgelsd or  dgels or gelsy
    // BUXTON: g_coef_meas_w = "constraint_weights * fitting_matrix"
    // BUXTON: meas_no_coil_cp = "constraint_weights * s_measured - constraint_weights * constraint_values_from_coils"
    // BUXTON: GSFit.rs uses "dgelss" = same!!
    // "single_vals" = singular values, not used
    // "rcond" = -1 == machine precision
    lapack_int rank;
    double rcond = -1.0;
    double single_vals[N_COEF];

    TSTART();
    *lapack_dgelss_info = LAPACKE_dgelss(
        LAPACK_COL_MAJOR,
        N_MEAS,
        N_COEF,
        1,
        g_coef_meas_w,
        N_MEAS,
        meas_no_coil_cp,
        N_MEAS,
        single_vals,
        rcond,
        &rank
    );
    TACC(T_LS_FIT);

    // BUXTON: copy "meas_no_coil_cp" into "coef"
    memcpy(coef, meas_no_coil_cp, sizeof(double) * N_COEF);

    // apply coeff to find current
    // BUXTON: matrix-vector multiplication; result stored in "source"
    // BUXTON: "source = g_pls_grid * coef"
    // BUXTON: source = plasma current on (R, Z) grid
    double source[N_GRID];

    TSTART();
    cblas_dgemv(CblasRowMajor, CblasTrans,
                N_PLS, N_GRID,
                1.0, g_pls_grid, N_GRID,
                coef, 1,
                0.0, source, 1);

    // `source` is the current density in each grid cell;
    // plasma_current = sum(source) * d_area
    double source_sum = 0.0;
    *r_cur_centroid = 0.0;
    *z_cur_centroid = 0.0;
    for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++) {
        source_sum += source[i_grid];
        *r_cur_centroid += source[i_grid] * R_GRID[i_grid];
        *z_cur_centroid += source[i_grid] * Z_GRID[i_grid];
    }
    *plasma_current = source_sum * DR * DZ;
    // Divide r_cur_centroid, z_cur_centroid by source_sum to get centroid position
    // provided the source_sum is not too close to zero or negative.
    if (*plasma_current >= PLASMA_CURRENT_CUTOFF) {
        *r_cur_centroid /= source_sum;
        *z_cur_centroid /= source_sum;
    } else {
        *r_cur_centroid = 0.0;
        *z_cur_centroid = 0.0;
    }
    TACC(T_SOURCE);

    // modelled measurements
    // BUXTON: measurements
    // BUXTON: "meas_model = g_coef_meas_w_orig * coef"
    // double meas_model_arr[N_MEAS];
    TSTART();
    cblas_dgemv(CblasRowMajor, CblasTrans,
                N_COEF, N_MEAS,
                1.0, g_coef_meas_w_orig, N_MEAS,
                coef, 1,
                0.0, meas_model, 1);
    TACC(T_MODEL_MEAS);

    // find chi squared error between meas and model
    TSTART();
    *chi_sq_err = 0.0;
    for (int32_t i_meas = 0; i_meas < N_MEAS_NO_REG; i_meas++)
    {
        double diff = meas_no_coil[i_meas] - meas_model[i_meas];
        *chi_sq_err += diff * diff;
    }
    TACC(T_CHI2);

    // If the plasma current is below the cutoff, there is effectively no
    // plasma present. In that case `source` (the plasma current density on
    // the grid) is made up of values very close to zero, which drives the
    // Poisson solver into subnormal floating point arithmetic. Subnormal
    // arithmetic is significantly slower than normal arithmetic on most
    // platforms, so we exit early here and skip the Poisson solve (and all
    // subsequent flux/LCFS/x-point processing that depends on it).
    if (*plasma_current < PLASMA_CURRENT_CUTOFF) {
        *lcfs_err_code = ERR_LOW_PLASMA_CURRENT;
        return;
    }

    // convert current to RHS of eq
    for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
    {
        source[i_grid] *= -R_MU0_DZ2[i_grid];  // BUXTON: R_MU0_DZ2=mu0 * R * d_area
    }

    //  poisson solver -> psi_plasma
    // BUXTON: calculate psi_plasma
    double flux_pls[N_GRID];

    TSTART();
    poisson_solver(source, flux_pls);
    TACC(T_POISSON);

    // coil psi on grid: flux_total = G_GRID_COIL * coil_curr */
    TSTART();
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                N_GRID, N_COIL,
                1.0, G_GRID_COIL, N_COIL,
                coil_curr, 1,
                0.0, flux_total, 1);
    TACC(T_COIL_FLUX);

    // calculate vessel flux on grid
    if (N_VESS > 0)
    {
        double flux_vessel[N_GRID];

        TSTART();
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    N_GRID, N_VESS,
                    1.0, G_GRID_VESSEL, N_VESS,
                    &coef[N_PLS], 1,
                    0.0, flux_vessel, 1);

        for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
        {
            flux_total[i_grid] += flux_pls[i_grid] + flux_vessel[i_grid];
        }
        TACC(T_VESSEL_FLUX);
    }
    else
    {
        for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
        {
            flux_total[i_grid] += flux_pls[i_grid];
        }
    }

    // find x point & opt
    double opt_r;
    double opt_z;
    double opt_flux;

    *xpt_n = 0;
    int32_t opt_n = 0;

    TSTART();
    *lcfs_err_code = 0;
    *lcfs_err_code |= find_nulls(flux_total,
               &opt_r, &opt_z, &opt_flux, &opt_n,
               xpt_r, xpt_z, xpt_flux, xpt_n);
    if (*lcfs_err_code != 0) {
        return;
    }

    // Check if mag axis found
    if (opt_n == 0)
    {
        *lcfs_err_code = ERR_NO_AXIS;
        return;
    }

    *mag_axis_flux = opt_flux;
    *r_mag_axis = opt_r;
    *z_mag_axis = opt_z;

    // Filter x-points
    filter_xpts(xpt_r, xpt_z, xpt_flux, xpt_n, *r_mag_axis, *z_mag_axis);
    TACC(T_XPTS_AND_AXIS);

    // sort x-points in descending order of flux value
    TSTART();
    sort_xpts(xpt_r, xpt_z, xpt_flux, *xpt_n);
    TACC(T_XPTS_SORT);

    // limiter flux with x-point filtering
    TSTART();
    double lcfs_flux = find_flux_on_limiter_xfiltered(flux_total,
                                                      xpt_r, xpt_z, *xpt_n,
                                                      *r_mag_axis, *z_mag_axis);
    // By default the plasma is taken to be wall-limited, so the xpt_diverted flag is set to 0.
    *xpt_diverted = 0;
    TACC(T_LIMITER);

    // select xpt
    if (*xpt_n > 0)
    {
        // We assume the xpts have already been sorted in descending order of flux value, so the first xpt has the highest flux value.
        double xpt_flux_max = xpt_flux[0];
        xpt_flux_max = FRAC * xpt_flux_max + (1.0 - FRAC) * (*mag_axis_flux);
        if (xpt_flux_max > lcfs_flux)
        {
            lcfs_flux = xpt_flux_max;
            // The plasma is x-point limited, so set the xpt_diverted flag to 1.
            *xpt_diverted = 1;
        }
    }

    if (fabs(lcfs_flux - (*mag_axis_flux)) < THRESH) {
      // Don't call normalise_flux() if lcfs_flux is too close to mag_axis_flux
      // This avoids division by a very small number.
      *lcfs_err_code |= ERR_AX_EQ_BDRY;
      return;
    }

    if (lcfs_flux > (*mag_axis_flux)) {
      // lcfs_flux should never be greater than mag_axis_flux
      *lcfs_err_code |= ERR_BDRY_GT_AX;
      return;
    }

    // extract LCFS
    TSTART();
    // *lcfs_err_code |= find_lcfs_rz(flux_total, lcfs_flux, lcfs_r, lcfs_z, lcfs_n);
    // No longer calculating lcfs_r, lcfs_z as we don't use them.
    // Just set them to zero.
    for (int32_t i = 0; i < N_LCFS_MAX; i++) {
        lcfs_r[i] = 0.0;
        lcfs_z[i] = 0.0;
    }
    *lcfs_n = 0;
    TACC(T_LCFS);

    // inside LCFS mask
    TSTART();
    // *lcfs_err_code |= inside_lcfs(*r_mag_axis, *z_mag_axis,
    //                               lcfs_r, lcfs_z, *lcfs_n, mask);
// int flood_fill_plasma_core(int32_t *mask, double *flux_total,
//                            double flux_boundary, double r_mag_axis,
//                            double z_mag_axis, double *xpt_r, double *xpt_z,
//                            int32_t xpt_n);

    *lcfs_err_code |=
        flood_fill_plasma_core(mask, flux_total, lcfs_flux, *r_mag_axis,
                              *z_mag_axis, xpt_r, xpt_z, *xpt_n);
    if (*lcfs_err_code != 0) {
        return;
    }
    TACC(T_INSIDE);

    // normalise total psi
    TSTART();
    normalise_flux(flux_total, lcfs_flux, *mag_axis_flux, mask, flux_norm);
    TACC(T_NORMALISE);

    // Store psi_b for later
    *flux_boundary = lcfs_flux;

#ifdef ENABLE_RT_TIMING
    timing_acc[T_TOTAL] += (thread_cpu_ns() - t_total_0);
#endif

}
