#include "rtgsfit.h"

#define N_GRID 2145

void rm_coil_from_meas(
        int n_meas,
        int n_coil,
        double* g_meas_coil,
        double* coil_curr,
        double* meas,
        double* meas_no_coil
        )
{
    int i_meas;
    // subtract PF (vessel) contributions from measurements    
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_meas, n_coil, 1.0, g_meas_coil, 
            n_coil, coil_curr, 1, 0.0, meas_no_coil, 1);  
    
    for (i_meas=0; i_meas<n_meas; i_meas++)
    {
        meas_no_coil[i_meas] = meas[i_meas] - meas_no_coil[i_meas];
    }
}
        
        
void make_basis(
        int n_row,
        int n_col,
        int n_grid,
        double d_row, 
        double* psi_norm,
        double* r_grid,
        double* inv_r_mu0,
        int* mask,
        double* basis
        )
{    
    int i_grid;
    
    // could use 1- psi_norm instead of psi_norm ?????
    for (i_grid=0; i_grid<n_grid; i_grid++)
    {
        if (mask[i_grid])
        {
            basis[i_grid] = (1 - psi_norm[i_grid]) * r_grid[i_grid];
            basis[i_grid + n_grid] = (1 -  psi_norm[i_grid]) * inv_r_mu0[i_grid];
        }
        else
        {
            basis[i_grid] = 0.0;
            basis[i_grid + n_grid] = 0.0;
        }
        
    }
    gradient_row(psi_norm, n_row, n_col, d_row, &basis[2*n_grid]);   
}



void normalise_psi(
        int n_grid, 
        double* psi_norm, 
        double* psi_total, 
        double psi_bound, 
        double psi_axis,
        int* mask
        )
{
    inv_psi_diff = 1.0/(psi_bound - psi_axis);
    
    // psi norm has to be of total flux, as boundary is defined in terms of total flux !
    for (i_grid=0; i_grid<n_grid; i_grid++)
    {
        if (mask[i_grid])
        {
            psi_norm[i_grid] = psi_total * inv_psi_diff;
        }
        else
        {
            psi_norm[i_grid] = 1.0;
        }
    }
}

void rtgsfit(
        int n_meas,
        double* meas,
        int n_coil,
        double* g_meas_coil,
        double* coil_curr,
        int n_row,
        int n_col,
        double* psi_norm,
        double d_row,
        double d_col,
        double* g_grid_sens,
        double* r_grid,
        double* inv_r_mu0,
        double* r_mu0_dz2,
        double* g_bound,
        double* g_grid_coeff,
        double* g_coeff_sens,
        double* lower_band,
        double* upper_band,
        int* perm_idx,
        int n_bound,
        double* g_grid_coil,
        double* r_vec,
        double* z_vec        
        )

{
    int n_grid = n_row*n_col;


    int n_coeff = 3;
    double basis[n_coeff*n_grid], g_coeff_sens[n_coeff*n_meas];
    int info, rank;
    double rcond = -1.0;
    double single_vals[n_coeff];
    double source[n_grid], meas_no_coil[n_meas];
    double psi_pls[n_grid]; psi_total[n_grid];
    int xpt_idx;
    double grad_z[n_grid], grad_r[n_grid], hess_rr[n_grid], hess_zz[n_grid];
    double hess_rz[n_grid];
    double psi_bound, psi_axis, opt_r, opt_z;
    
    
    // subtract PF (vessel) contributions from measurements    
    rm_coil_from_meas(n_meas, n_coil, g_meas_coil, coil_curr, meas, meas_no_coil);
       
    // make basis    
    make_basis(n_row, n_col, n_grid, d_row, psi_norm, r_grid, inv_r_mu0, mask, basis);
    
    // make meas2coeff matrix
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_coeff, n_meas, n_grid, 
            1.0, basis, n_grid, g_grid_sens, n_meas, 0.0, g_coeff_sens, n_meas);    

    // fit coeff
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n_meas, n_coeff, 1, g_coeff_sens, 
            n_meas, meas_no_coil, n_meas, single_vals, rcond, &rank);
            
    // apply coeff to find current
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_grid, n_coeff, 1.0, g_grid_coeff, 
            n_coeff, meas_no_coil, 1, 0.0, source, 1);   
         
    // convert current to RHS of eq   
    for (i_grid=0; i_grid<n_grid; i_grid++)
    {
        source[i_grid] *= -r_mu0_dz2;
    }          
    
    //  poisson solver -> psi_plasma
    poisson_solver(n_col, n_row, n_grid, lower_band, upper_band, source, 
        perm_idx, n_bound, g_bound, inv_r_mu0, d_row, d_col, psi_pls);
        
    // calculate coil psi on grid
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_grid, n_coil, 1.0, g_grid_coil, 
            n_coil, coil_curr, 1, 0.0, psi_total, 1);     
            
    for (i_grid=0; i_grid<n_grid; i_grid++)
    {
        psi_total[i_grid] += TWO_PI * psi_pls;
    }  
        
    // find hessian & gradients 
    gradient_row(psi, n_row, n_col, d_row, grad_z);   
    gradient_col(psi, n_row, n_col, d_col, grad_r);
    hessian_row_row(psi, n_row, n_col, d_row, hess_zz);
    hessian_col_col(psi, n_row, n_col, d_col, hess_rr);
    hessian_row_col(psi, n_row, n_col, d_row, d_col, hess_rz);
    
    // find x point & opt
    find_opt_xpt(d_col, d_row, n_row, n_col, r_vec, z_vec, psi, grad_r, grad_z, 
            hess_rr, hess_zz, hess_rz, psi_bound, psi_axis, opt_r, opt_z);
    
    // find LCFS mask
    find_lcfs_mask(d_col, d_row, n_row, n_col, r_vec, z_vec, psi_total,
            psi_bound, mask);
    
    // normalise total psi                                
    normalise_psi(n_grid, psi_norm, psi_total, psi_bound, psi_axis, mask);
}                   
   
    

    
/*void make_basis(*/
/*        int n_row,*/
/*        int n_col,*/
/*        int n_grid,*/
/*        double d_row, */
/*        double* psi_norm,*/
/*        double* r_grid,*/
/*        double* inv_r_mu0,*/
/*        int* inside_idx,*/
/*        int inside_n,*/
/*        double* basis*/
/*        )*/
/*{    */
/*    int i_grid;*/
/*    double dpsi_dz[n_grid];*/
/*    */
/*    for (i_grid=0; i_grid<3*n_grid; i_grid++)*/
/*    {*/
/*        basis[i_grid] = 0.0;*/
/*    }*/
/*    */
/*    gradient_row(psi_norm, n_row, n_col, d_row, &dpsi_dz);   */
/*    */
/*    for (i_idx=0; i_idx<inside_n; i_idx++)*/
/*    {*/
/*        basis[inside_idx[i_idx]] = (1 - psi_norm[inside_idx[i_idx]]) * r_grid[inside_idx[i_idx]];*/
/*        basis[inside_idx[i_idx] + n_grid] = (1 -  psi_norm[inside_idx[i_idx]]) * inv_r_mu0[inside_idx[i_idx]];*/
/*        basis[inside_idx[i_idx] + 2*n_grid] = dpsi_dz[inside_idx[i_idx]];*/
/*    }*/
/*}*/

    
    


