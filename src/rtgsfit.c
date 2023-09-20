#include "constants.h"
#include <math.h>
#include <cblas.h>
#include "gradient.h"
#include "rtgsfit.h"
#include <lapacke.h>
#include "poisson_solver.h"
#include "find_x_point.h"
#include <stdio.h>


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

    // could use 1 - flux_norm instead of flux_norm ?????
    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        if (mask[i_grid])
        {
            basis[i_grid] = (1 - flux_norm[i_grid]) * R_GRID[i_grid];
            basis[i_grid + N_GRID] = (1 -  flux_norm[i_grid]) * INV_R_MU0[i_grid];
        }
        else
        {
            basis[i_grid] = 0.0;
            basis[i_grid + N_GRID] = 0.0;
        }        
    }
    // could use gradient from previous iteration.  should apply mask 
    gradient_z(flux_norm, &basis[2*N_GRID]);   
}



void normalise_flux(
        double* flux_norm, 
        double* flux_total, 
        double flux_lcfs, 
        double flux_axis,
        int* mask
        )
{
    double inv_flux_diff;
    int i_grid;
    inv_flux_diff = 1.0/(flux_lcfs - flux_axis);
    
    // psi norm has to be of total flux, as boundary is defined in terms of total flux !
    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        if (mask[i_grid])
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
        double* meas,
        double* coil_curr,
        double* flux_norm, 
        int* mask
        )
{
    double g_coef_meas[N_COEF*N_MEAS];
    double basis[N_COEF*N_GRID];
    int info, rank, i_grid, i_opt, i_xpt, i_meas;
    double rcond = -1.0;
    double xpt_flux_max;
    double single_vals[N_COEF];
    double source[N_GRID], meas_no_coil[N_MEAS];
    double flux_pls[N_GRID], flux_total[N_GRID];
    double lcfs_flux, axis_flux, axis_r, axis_z;
    double lcfs_r[N_LCFS_MAX], lcfs_z[N_LCFS_MAX];
    double xpt_r[N_XPT_MAX], xpt_z[N_XPT_MAX], xpt_flux[N_XPT_MAX];
    double opt_r[N_XPT_MAX], opt_z[N_XPT_MAX], opt_flux[N_XPT_MAX];
    int lcfs_n;
    int xpt_n = 0;
    int opt_n = 0;
    
    // subtract PF (vessel) contributions from measurements    
    rm_coil_from_meas(coil_curr, meas, meas_no_coil);

    // make basis    
    make_basis(flux_norm, mask, basis);

    // make meas2coeff matrix
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_COEF, N_MEAS, N_GRID, 
            1.0, basis, N_GRID, G_GRID_MEAS, N_MEAS, 0.0, g_coef_meas, N_MEAS);    

    for (i_meas=0; i_meas<N_MEAS; i_meas++)
    {
        meas_no_coil[i_meas] *= WEIGHT[i_meas];
    }

    // fit coeff or use dgelsd
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR, N_MEAS, N_COEF, 1, g_coef_meas, 
            N_MEAS, meas_no_coil, N_MEAS, single_vals, rcond, &rank);

    // apply coeff to find current
    cblas_dgemv(CblasColMajor, CblasNoTrans, N_GRID, N_COEF, 1.0, basis, 
            N_GRID, meas_no_coil, 1, 0.0, source, 1);   

    // convert current to RHS of eq   
    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        source[i_grid] *= -R_MU0_DZ2[i_grid];
    }          

    //  poisson solver -> psi_plasma
    poisson_solver(source, flux_pls);

    // calculate coil psi on grid
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_GRID, N_COIL, 1.0, G_GRID_COIL, 
            N_COIL, coil_curr, 1, 0.0, flux_total, 1);     

    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        flux_total[i_grid] += 2 * M_PI * flux_pls[i_grid];
    }          

    // find x point & opt
    find_null_in_gradient(flux_total, opt_r, opt_z, opt_flux, &opt_n, 
            xpt_r, xpt_z, xpt_flux, &xpt_n);

    // select opt
    i_opt = max_idx(opt_n, opt_flux);    
    axis_flux = opt_flux[i_opt];
    axis_r = opt_r[i_opt];
    axis_z = opt_z[i_opt];  
    
    // select xpt      
    i_xpt = max_idx(xpt_n, xpt_flux);
    xpt_flux_max = xpt_flux[i_xpt];      
    lcfs_flux = FRAC * xpt_flux_max + (1-FRAC)*axis_flux;

    // extract LCFS
    find_lcfs_rz(flux_total, lcfs_flux, lcfs_r, lcfs_z, &lcfs_n);         

    // extract inside of LCFS
    inside_lcfs(axis_r, axis_z, lcfs_r, lcfs_z, lcfs_n, mask);

    // normalise total psi                                
    normalise_flux(flux_norm, flux_total, lcfs_flux, axis_flux, mask);
}                   
   
