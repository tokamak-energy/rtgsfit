#include "constants.h"
#include "gradient.h"
#include "rtgsfit.h"
#include "poisson_solver.h"
#include "find_x_point.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#define N_MEAS_NO_REG N_BP_PROBES + N_FLUX_LOOPS + N_ROGOWSKI_COILS

int32_t max_idx(
        int32_t n_arr,
        double* arr
        )
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
        double* coil_curr,
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

double find_flux_on_limiter(double* flux_total)
{

    int i_limit, i_intrp, idx;
    double flux_limit_max, flux_limit;

    flux_limit_max = -DBL_MAX;

    for (i_limit = 0; i_limit < N_LIMIT; i_limit++)
    {
        flux_limit = 0.0;
        for (i_intrp = 0; i_intrp < N_INTRP; i_intrp++)
        {
            idx = i_limit*N_INTRP + i_intrp;
            flux_limit += LIMIT_WEIGHT[idx] * flux_total[LIMIT_IDX[idx]];
        }
        if (flux_limit > flux_limit_max)
        {
            flux_limit_max = flux_limit;
        }
    }
    return flux_limit_max;
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
 * w = (axis_r - xpt_r[i], axis_z - xpt_z[i]).
 * 
 * If the dot product of v and w is negative, then the limiter point
 * is not considered for the flux calculation.
 *
 * @param flux_total Array of total flux values.
 * @param xpt_r Array of x-point R coordinates.
 * @param xpt_z Array of x-point Z coordinates.
 * @param xpt_n Number of x-points.
 * @param axis_r R coordinate of the axis.
 * @param axis_z Z coordinate of the axis.
 * @return The computed flux value on the limiter after x-point-filtering.
 */
double find_flux_on_limiter_xfiltered(double flux_total[],
                                      double xpt_r[],
                                      double xpt_z[],
                                      int xpt_n,
                                      double axis_r,
                                      double axis_z)
{

    int skip;
    double flux_limit_max;

    flux_limit_max = -DBL_MAX;

    for (int32_t i_limit = 0; i_limit < N_LIMIT; i_limit++)
    {

        // If dot product of (R_LIM[i_limit] - xpt_r, Z_LIM[i_limit] - xpt_z) and 
        // (axis_r - xpt_r, axis_z - xpt_z) is negative for any xpt, skip this limit point
        skip = 0;
        for (int32_t i_xpt = 0; i_xpt < xpt_n; i_xpt++)
        {
            double dot_product = (LIMIT_R[i_limit] - xpt_r[i_xpt]) * (axis_r - xpt_r[i_xpt]) +
                          (LIMIT_Z[i_limit] - xpt_z[i_xpt]) * (axis_z - xpt_z[i_xpt]);
            if (dot_product < 0.0)
            {
                skip = 1; // skip this limit point
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

    // psi norm has to be of total flux, as boundary is defined in terms of total flux !
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
        double* meas_pcs, // input
        double* coil_curr, // input
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
        int64_t *lapack_dgelss_info, // output
        double *meas_model, // output
        int32_t n_meas_model // input
        )
{
    assert(n_meas_model == N_MEAS);
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
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_MEAS_NO_REG, N_SENS_PCS, 1.0,
            SENSOR_REPLACEMENT_MATRIX, N_SENS_PCS, meas_pcs, 1, 0.0, meas, 1);

    // will this be done during compilation?
    double g_coef_meas_w[N_COEF * N_MEAS];
    memcpy(g_coef_meas_w, G_COEF_MEAS_WEIGHT, sizeof(double) * N_MEAS * N_COEF);

    // subtract PF contributions from measurements
    // Note that this also sets the regularisation elements of meas_no_coil to zero.
    // meas_no_coil is the post-processed measurments with the PF coil contributions removed
    // and includes the regularisation elements, which can be thought of as fake Rogowski coil
    // measurements.
    double meas_no_coil[N_MEAS];
    rm_coil_from_meas(coil_curr, meas, meas_no_coil);

    // make basis
    // This makes the transpose of the T_{yg} matrix in eqn. (61) of the Moret et al. (2015)
    // LIUQE paper, without the Delta_R * Delta_Z factor.
    double g_pls_grid[N_PLS * N_GRID];
    make_basis(flux_norm, mask, g_pls_grid);

    // make meas-pls matrix
    // g_coef_meas_w = g_pls_grid * G_GRID_MEAS_WEIGHT
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_PLS, N_MEAS, N_GRID,
            1.0, g_pls_grid, N_GRID, G_GRID_MEAS_WEIGHT, N_MEAS, 0.0,
            g_coef_meas_w, N_MEAS);

    // form meas vectors from measurements
    for (int32_t i_meas = 0; i_meas < N_MEAS; i_meas++)
    {
        meas_no_coil[i_meas] *= WEIGHT[i_meas];
    }

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
    *lapack_dgelss_info = (int64_t)LAPACKE_dgelss(
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

    // BUXTON: copy "meas_no_coil_cp" into "coef"
    memcpy(coef, meas_no_coil_cp, sizeof(double) * N_COEF);

    // apply coeff to find current
    // BUXTON: matrix-vector multiplication; result stored in "source"
    // BUXTON: "source = g_pls_grid * coef"
    // BUXTON: source = plasma current on (R, Z) grid
    double source[N_GRID];
    cblas_dgemv(CblasRowMajor, CblasTrans, N_PLS, N_GRID, 1.0, g_pls_grid,
            N_GRID, coef, 1, 0.0, source, 1);

    // `source` is the current density in each grid cell;
    // plasma_current = sum(source) * d_area
    double source_sum = 0.0;
    for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++) {
        source_sum += source[i_grid];
    }
    *plasma_current = source_sum * DR * DZ;


    // modelled measurements
    // BUXTON: measurements
    // BUXTON: "meas_model = g_coef_meas_w_orig * coef"
    // double meas_model_arr[N_MEAS];
    cblas_dgemv(CblasRowMajor, CblasTrans, N_COEF, N_MEAS, 1.0, g_coef_meas_w_orig,
            N_MEAS, coef, 1, 0.0, meas_model, 1);
    
    // find chi squared error between meas and model
    *chi_sq_err = 0.0;
    for (int32_t i_meas = 0; i_meas < N_MEAS; i_meas++)
    {
        *chi_sq_err = *chi_sq_err + (meas_no_coil[i_meas] -  meas_model[i_meas]) \
                    * (meas_no_coil[i_meas] -  meas_model[i_meas]);
    }

    // convert current to RHS of eq
    for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
    {
        source[i_grid] *= -R_MU0_DZ2[i_grid];  // BUXTON: R_MU0_DZ2=mu0 * R * d_area
    }

    //  poisson solver -> psi_plasma
    // BUXTON: calculate psi_plasma
    double flux_pls[N_GRID];
    poisson_solver(source, flux_pls);

    // calculate coil psi on grid
    // BUXTON: "flux_total = G_GRID_COIL * coil_curr"
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_COIL, 1.0, G_GRID_COIL,
            N_COIL, coil_curr, 1, 0.0, flux_total, 1);

    // calculate vessel flux on grid
    if (N_VESS > 0)
    {

        double flux_vessel[N_GRID];
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_VESS, 1.0, G_GRID_VESSEL,
                N_VESS, &coef[N_PLS], 1, 0.0, flux_vessel, 1);

        // calculate total flux = coil + plasma + vessel
        for (int32_t i_grid = 0; i_grid < N_GRID; i_grid++)
        {
            flux_total[i_grid] +=  flux_pls[i_grid] + flux_vessel[i_grid];
            // BUXTON: do I need to add delta_z*d(psi)/d(psi_n)
        }
    }
    else
    {
        for (int32_t i_grid=0; i_grid < N_GRID; i_grid++)
        {
            flux_total[i_grid] +=  flux_pls[i_grid];
        }
    }

    // flux value on limiter
    // lcfs_flux = find_flux_on_limiter(flux_total); // Using find_flux_on_limiter_xfiltered instead.

    // find x point & opt
    double xpt_r[N_XPT_MAX];
    double xpt_z[N_XPT_MAX];
    double xpt_flux[N_XPT_MAX];
    double opt_r[N_XPT_MAX];
    double opt_z[N_XPT_MAX];
    double opt_flux[N_XPT_MAX];

    int32_t xpt_n = 0;
    int32_t opt_n = 0;
    find_null_in_gradient_march(flux_total, opt_r, opt_z, opt_flux, &opt_n,
            xpt_r, xpt_z, xpt_flux, &xpt_n);

    // select opt
    int32_t i_opt = max_idx(opt_n, opt_flux);
    double axis_flux = opt_flux[i_opt];
    double axis_r = opt_r[i_opt];
    double axis_z = opt_z[i_opt];

    double lcfs_flux = find_flux_on_limiter_xfiltered(flux_total, xpt_r, xpt_z, xpt_n, axis_r, axis_z);

    // select xpt
    if (xpt_n > 0)
    {
        int32_t i_xpt = max_idx(xpt_n, xpt_flux);
        double xpt_flux_max = xpt_flux[i_xpt];
        xpt_flux_max = FRAC * xpt_flux_max + (1.0-FRAC)*axis_flux;
        if (xpt_flux_max > lcfs_flux)
        {
            lcfs_flux = xpt_flux_max;
        }
    }

    // extract LCFS
    find_lcfs_rz(flux_total, lcfs_flux, lcfs_r, lcfs_z, lcfs_n);

    // extract inside of LCFS
    // BUXTON: we think this might have an error??????
    *lcfs_err_code = inside_lcfs(axis_r, axis_z, lcfs_r, lcfs_z, *lcfs_n, mask);

    // normalise total psi
    normalise_flux(flux_total, lcfs_flux, axis_flux, mask, flux_norm);

    // store axis_r, axis_z and axis_flux in the meas_pcs array
    // meas_pcs[0] = axis_r;
    // meas_pcs[1] = axis_z;
    // meas_pcs[2] = axis_flux;

    // Store psi_b for later
    *flux_boundary = lcfs_flux;
}
