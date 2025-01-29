#include "constants.h"
#include <math.h>
#include <cblas.h>
#include "gradient.h"
#include "rtgsfit.h"
#include <lapacke.h>
#include "poisson_solver.h"
#include "find_x_point.h"
#include <stdio.h>
#include <float.h>
#include <string.h>


int max_idx(
        int n_arr,
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


void normalise_flux(
        double* flux_total,
        double flux_lcfs,
        double flux_axis,
        int* mask,
        double* flux_norm
        )
{
    double inv_flux_diff;
    int i_grid;
    inv_flux_diff = 1.0/(flux_lcfs - flux_axis);

    // psi norm has to be of total flux, as boundary is defined in terms of total flux !
    for (i_grid=0; i_grid<N_GRID; i_grid++)
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



int rtgsfit(
        double* meas,
        double* coil_curr,
        double* flux_norm,
        int* mask,
        double* flux_total,
        double* error,
        double* lcfs_r,
        double* lcfs_z,
        int* lcfs_n,
        double* coef,
        double* flux_boundary
        )
{
    double g_coef_meas_w[N_COEF*N_MEAS], g_coef_meas_w_orig[N_COEF*N_MEAS];
    double g_pls_grid[N_PLS*N_GRID];
    int info, i_grid, i_opt, i_xpt, i_meas;
    int rank;
    double rcond = -1.0;
    double xpt_flux_max;
    double single_vals[N_COEF], meas_no_coil_cp[N_MEAS];
    double source[N_GRID], meas_no_coil[N_MEAS], meas_model[N_MEAS];
    double flux_pls[N_GRID], flux_vessel[N_GRID];
    double lcfs_flux, axis_flux, axis_r, axis_z;
    double xpt_r[N_XPT_MAX], xpt_z[N_XPT_MAX], xpt_flux[N_XPT_MAX];
    double opt_r[N_XPT_MAX], opt_z[N_XPT_MAX], opt_flux[N_XPT_MAX];
    int xpt_n = 0;
    int opt_n = 0;

    double g_coef_meas_w_with_regularisation[N_COEF*(N_MEAS+5)];  // order: [(i_coef=0, i_meas=0), (i_coef=0, i_meas=1), (i_coef=0, i_meas=2), ...]
    double meas_no_coil_cp_with_regularisation[N_MEAS + 5];

    // will this be done during compilation?
    memcpy(g_coef_meas_w, G_COEF_MEAS_WEIGHT, sizeof(double)*N_MEAS*N_COEF);

    // subtract PF (vessel) contributions from measurements
    rm_coil_from_meas(coil_curr, meas, meas_no_coil);

    // make basis
    make_basis(flux_norm, mask, g_pls_grid);

    // make meas-pls matrix
    // g_coef_meas_w = g_pls_grid * G_GRID_MEAS_WEIGHT
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_PLS, N_MEAS, N_GRID,
            1.0, g_pls_grid, N_GRID, G_GRID_MEAS_WEIGHT, N_MEAS, 0.0,
            g_coef_meas_w, N_MEAS);

    // form meas vectors from measurements
    for (i_meas=0; i_meas<N_MEAS; i_meas++)
    {
        meas_no_coil[i_meas] *= WEIGHT[i_meas];
    }

    // copy measurment to coef due to overwritting in LAPACKE_dgelss
    memcpy(meas_no_coil_cp, meas_no_coil, sizeof(double)*N_MEAS);
    memcpy(g_coef_meas_w_orig, g_coef_meas_w, sizeof(double)*N_COEF*N_MEAS);

    // fit coeff or use dgelsd or  dgels or gelsy
    // BUXTON: g_coef_meas_w = "constraint_weights * fitting_matrix"
    // BUXTON: meas_no_coil_cp = "constraint_weights * s_measured - constraint_weights * constraint_values_from_coils"
    // BUXTON: GSFit.rs uses "dgelss" = same!!
    // "single_vals" = singular values, not used
    // "rcond" = -1 == machine precision

    // Add regularisations
    for (int i_coef = 0; i_coef < N_COEF; i_coef++) {
        for (int i_meas = 0; i_meas < N_MEAS; i_meas++) {
            g_coef_meas_w_with_regularisation[i_coef * (N_MEAS + 5) + i_meas] = g_coef_meas_w[i_coef*N_MEAS + i_meas];
        }
        for (int i_meas = N_MEAS; i_meas < N_MEAS + 5; i_meas++) {
            g_coef_meas_w_with_regularisation[i_coef * (N_MEAS + 5) + i_meas] = 0.0;
        }
    }





    // for (int i_meas = 0; i_meas < N_MEAS; i_meas++) {
    //     for (int i_coef = 0; i_coef < N_COEF; i_coef++) {
    //         g_coef_meas_w_with_regularisation[i_meas * (N_COEF) + i_coef] = g_coef_meas_w[i_meas * N_COEF + i_coef];
    //     }
    // }
    // for (int i_meas = N_MEAS; i_meas < N_MEAS + 5; i_meas++) {
    //     for (int i_coef = 0; i_coef < N_COEF; i_coef++) {
    //         g_coef_meas_w_with_regularisation[i_meas * (N_COEF) + i_coef] = 0.0;
    //     }
    // }


    g_coef_meas_w_with_regularisation[69*4-5] = 15.0 * 0.001 * 20.0;
    g_coef_meas_w_with_regularisation[69*4-5 + 69 + 1] = 15.0 * 0.001 * 20.0;
    g_coef_meas_w_with_regularisation[69*4-5 + 69*2 + 2] = 15.0 * 0.001 * 20.0;
    g_coef_meas_w_with_regularisation[69*4-5 + 69*3 + 3] = 15.0 * 0.001 * 20.0;
    g_coef_meas_w_with_regularisation[69*4-5 + 69*4 + 4] = 15.0 * 0.001 * 20.0;
    // g_coef_meas_w_with_regularisation[N_COEF * (N_MEAS + 2) - 1] = 15.0 * 0.001 * 20.0;
    // g_coef_meas_w_with_regularisation[N_COEF * (N_MEAS + 3) - 1] = 15.0 * 0.001 * 20.0;
    // g_coef_meas_w_with_regularisation[N_COEF * (N_MEAS + 4) - 1] = 15.0 * 0.001 * 20.0;
    // g_coef_meas_w_with_regularisation[N_COEF * (N_MEAS + 5) - 1] = 15.0 * 0.001 * 20.0;

    for (i_meas=0; i_meas<N_MEAS; i_meas++)
    {
        meas_no_coil_cp_with_regularisation[i_meas] = meas_no_coil_cp[i_meas];
    }
    for (i_meas=N_MEAS; i_meas<N_MEAS + 5; i_meas++)
    {
        meas_no_coil_cp_with_regularisation[i_meas] = 0.0;
    }

    // info = LAPACKE_dgelss(
    //   LAPACK_COL_MAJOR,
    //   N_MEAS,
    //   N_COEF,
    //   1,
    //   g_coef_meas_w,
    //   N_MEAS,
    //   meas_no_coil_cp,
    //   N_MEAS,
    //   single_vals,
    //   rcond,
    //   &rank
    // );
    info = LAPACKE_dgelss(
      LAPACK_COL_MAJOR,
      N_MEAS+5,
      N_COEF,
      1,
      g_coef_meas_w_with_regularisation,
      N_MEAS+5,
      meas_no_coil_cp_with_regularisation,
      N_MEAS+5,
      single_vals,
      rcond,
      &rank
    );

    // BUXTON: copy "meas_no_coil_cp" into "coef"
    memcpy(coef, meas_no_coil_cp_with_regularisation, sizeof(double) * N_COEF);

    // apply coeff to find current
    // BUXTON: matrix-vector multiplication; result stored in "source"
    // BUXTON: "source = g_pls_grid * coef"
    // BUXTON: source = plasma current on (R, Z) grid
    cblas_dgemv(CblasRowMajor, CblasTrans, N_PLS, N_GRID, 1.0, g_pls_grid,
            N_GRID, coef, 1, 0.0, source, 1);

    // Sum the variable "source" and print to screen; should be plasma current???
    double source_sum = 0.0;
    for (i_grid = 0; i_grid < N_GRID; i_grid++) {
        source_sum += source[i_grid];
    }
    printf("Sum of source: %lf\n", source_sum * 9e-4); // BUXTON: this is close to plasma current


    // modelled measurements
    // BUXTON: measurements
    // BUXTON: "meas_model = g_coef_meas_w_orig * coef"
    cblas_dgemv(CblasRowMajor, CblasTrans,  N_COEF, N_MEAS, 1.0, g_coef_meas_w_orig,
            N_MEAS, coef, 1, 0.0, meas_model, 1);

    // find error between meas and model
    *error = 0.0;
    for (i_meas=0; i_meas<N_MEAS; i_meas++)
    {
        *error = *error + (meas_no_coil[i_meas] -  meas_model[i_meas]) \
                * (meas_no_coil[i_meas] -  meas_model[i_meas]);
    }

    // convert current to RHS of eq
    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        source[i_grid] *= -R_MU0_DZ2[i_grid];  // BUXTON: R_MU0_DZ2=mu0 * R * d_area
    }

    //  poisson solver -> psi_plasma
    // BUXTON: calculate psi_plasma
    poisson_solver(source, flux_pls);

    // calculate coil psi on grid
    // BUXTON: "flux_total = G_GRID_COIL * coil_curr"
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_COIL, 1.0, G_GRID_COIL,
            N_COIL, coil_curr, 1, 0.0, flux_total, 1);

    // calculate vessel flux on grid
    if (N_VESS > 0)
    {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_VESS, 1.0, G_GRID_VESSEL,
                N_VESS, &coef[N_PLS], 1, 0.0, flux_vessel, 1);

        // calculate total flux = coil + plasma + vessel
        for (i_grid=0; i_grid<N_GRID; i_grid++)
        {
            flux_total[i_grid] +=  flux_pls[i_grid] + flux_vessel[i_grid];
            // BUXTON: do I need to add delta_z*d(psi)/d(psi_n)
        }
    }
    else
    {
        for (i_grid=0; i_grid<N_GRID; i_grid++)
        {
            flux_total[i_grid] +=  flux_pls[i_grid];
        }
    }

    // flux value on limiter
    lcfs_flux = find_flux_on_limiter(flux_total);

    // find x point & opt
    find_null_in_gradient_march(flux_total, opt_r, opt_z, opt_flux, &opt_n,
            xpt_r, xpt_z, xpt_flux, &xpt_n);

    // select opt
    i_opt = max_idx(opt_n, opt_flux);
    axis_flux = opt_flux[i_opt];
    axis_r = opt_r[i_opt];
    axis_z = opt_z[i_opt];

    // select xpt
    if (xpt_n > 0)
    {
        i_xpt = max_idx(xpt_n, xpt_flux);
        xpt_flux_max = xpt_flux[i_xpt];
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
    inside_lcfs(axis_r, axis_z, lcfs_r, lcfs_z, *lcfs_n, mask);

    // normalise total psi
    normalise_flux(flux_total, lcfs_flux, axis_flux, mask, flux_norm);

    // Store psi_b for later
    *flux_boundary = lcfs_flux;
  return info;
}
