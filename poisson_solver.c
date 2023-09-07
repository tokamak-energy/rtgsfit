#include "poisson_solver.h"


/*
 * Function: hagenow_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * n_r - number of radial grid positions
 * n_z - number of vertical grid poisitons
 * n_ele - number of elements in the grid n_r*n_z
 * lower (n_ele, n_r) - Non zero subdiagonals of the lower triangular matrix
 * upper (n_ele, n_r+2) - Diagonal and superdiagonals of the upper triangle matrix
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * idx_final - array of final index of the original indexes 0:n_row-1. 
 * n_bound - number of elements on the boundary of the grid NOT REQUIRED
 * g_bound (n_bound, n_bound) - Green's matrix of boundary elements
 * inv_r_mu0 (n_bound, ) - Inverse of the major radius multipled by mu0
 * psi (n_ele, ) - Temporary variable to hold the psi matrix
 *
 * Outputs: 
 * psi_bound (n_bound, ) - psi values on the boundary
 */ 
void hagenow_bound(
        int n_r, 
        int n_z, 
        int n_ele,
        double* lower, 
        double* upper, 
        double* b_vec, 
        int* idx_final, 
        int n_bound,
        double* g_bound,
        double* inv_r_mu0,
        double d_row,
        double d_col, 
        double* psi,
        double* psi_bound
        )
{

    double dpsi_ltrb[n_bound];
    int ii;

    solve_tria(n_ele, n_r, lower, upper, b_vec, idx_final, psi);
    
    gradient_bound(psi, n_z, n_r, d_row, d_col, dpsi_ltrb);
    
    for (ii=0; ii<n_bound; ++ii)
    {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*inv_r_mu0[ii];
    }
        
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_bound, n_bound, 1.0, g_bound, 
            n_bound, dpsi_ltrb, 1, 0.0, psi_bound, 1);    
}

    
/*
 * Function: add_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * n_r - number of radial grid positions
 * n_z - number of vertical grid poisitons
 * n_ele - number of elements in the grid n_r*n_z
 * n_bound - number of elements on the boundary of the grid
 * psi_bound (n_bound, ) - psi values on the boundary
 *
 * Outputs: 
 * b_vec (n_ele, ) - Current density with flux boundary values added
 */     
void add_bound(
        int n_r, 
        int n_z, 
        int n_ele, 
        int n_bound, 
        double* psi_bound, 
        double* b_vec
        )
{
    int ii;
    
    for (ii=1; ii<n_z-1; ++ii) 
    {
        b_vec[ii*n_r] = psi_bound[ii];
    }
        
    for (ii=0; ii<n_r; ++ii) 
    {
        b_vec[ii] = psi_bound[ii+n_z];
    }
        
    for (ii=0; ii<n_z; ++ii)
    {
        b_vec[ii*n_r + n_r -1]= psi_bound[ii + n_z + n_r];
    }
          
    for (ii=0; ii<n_r; ++ii)
    {
        b_vec[(n_z-1)*n_r + ii]= psi_bound[ii + 2*n_z + n_r];
    }
        
}     



/* Function: poisson_solver
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * n_r - number of radial grid positions
 * n_z - number of vertical grid poisitons
 * n_ele - number of elements in the grid n_r*n_z
 * lower (n_ele, n_r) - Non zero subdiagonals of the lower triangular matrix
 * upper (n_ele, n_r+2) - Diagonal and superdiagonals of the upper triangle matrix
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * idx_final (n_ele, ) - array of final index of the original indexes 0:n_ele-1. 
 * n_bound - number of elements on the boundary of the grid NOT REQUIRED
 * g_bound (n_bound, n_bound) - Green's matrix of boundary elements
 * inv_r_mu0 (n_bound, ) - Inverse of the major radius x mu0 along boundary
 *
 * Outputs: 
 * out (n_ele, ) - psi values on the 2D grid
 */ 
void poisson_solver(
        int n_r, 
        int n_z, 
        int n_ele,
        double* lower, 
        double* upper, 
        double* b_vec, 
        int* idx_final, 
        int n_bound,
        double* g_bound,
        double* inv_r_mu0,
        double dz, 
        double dr, 
        double* out 
        )
{

    double psi_bound[n_bound];
    
    hagenow_bound(n_r, n_z, n_ele, lower, upper, b_vec, idx_final, n_bound,
            g_bound, inv_r_mu0, dz, dr, out, psi_bound);
            
    add_bound(n_r, n_z, n_ele, n_bound, psi_bound, b_vec);
    
    solve_tria(n_ele, n_r, lower, upper, b_vec, idx_final, out);
}  
    
    
