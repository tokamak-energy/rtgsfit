#include <cblas.h>
#include "poisson_solver.h"
#include "solve_tria.h"
#include "gradient.h"
#include "constants.h"
#include <stdio.h>
/*
 * Function: hagenow_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * lower (n_ele, N_R) - Non zero subdiagonals of the lower triangular matrix
 * upper (n_ele, N_R+2) - Diagonal and superdiagonals of the upper triangle matrix
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * idx_final - array of final index of the original indexes 0:n_row-1. 
 * N_LTRB - number of elements on the boundary of the grid NOT REQUIRED
 * G_LTRB (N_LTRB, N_LTRB) - Green's matrix of boundary elements
 * inv_r_mu0 (N_LTRB, ) - Inverse of the major radius multipled by mu0
 * psi (n_ele, ) - Temporary variable to hold the psi matrix
 *
 * Outputs: 
 * psi_bound (N_LTRB, ) - psi values on the boundary
 */ 
void hagenow_bound(     
        double* b_vec, 
        double* psi,
        double* psi_ltrb
        )
{

    double dpsi_ltrb[N_LTRB];
    int ii;

    solve_tria(b_vec, psi);
    
    gradient_bound(psi, dpsi_ltrb);
    
    for (ii=0; ii<N_LTRB; ++ii)
    {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*INV_R_LTRB_MU0[ii];
    }
        
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_LTRB, N_LTRB, 1.0, G_LTRB, 
            N_LTRB, dpsi_ltrb, 1, 0.0, psi_ltrb, 1);    
}

    
/*
 * Function: add_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * N_LTRB - number of elements on the boundary of the grid
 * psi_bound (N_LTRB, ) - psi values on the boundary
 *
 * Outputs: 
 * b_vec (n_ele, ) - Current density with flux boundary values added
 */     
void add_bound(
        double* psi_bound, 
        double* b_vec
        )
{
    int ii;
    
    // left
    for (ii=0; ii<N_Z; ++ii) 
    {
        b_vec[ii*N_R] = psi_bound[ii];
    }
        
    // top
    for (ii=1; ii<N_R-1; ++ii) 
    {
        b_vec[ii + (N_Z-1)*N_R] = psi_bound[ii+N_Z-1];
    }
    
    // right    
    for (ii=0; ii<N_Z; ++ii)
    {
        b_vec[ii*N_R + N_R -1]= psi_bound[N_Z + N_R - 3 + (N_Z - ii)];
    }
     
    //bottom     
    for (ii=1; ii<N_R-1; ++ii)
    {
        b_vec[ii]= psi_bound[2*N_Z + N_R + (N_R - ii) - 4];
    }
        
}     



/* Function: poisson_solver
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * lower (n_ele, N_R) - Non zero subdiagonals of the lower triangular matrix
 * upper (n_ele, N_R+2) - Diagonal and superdiagonals of the upper triangle matrix
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * idx_final (n_ele, ) - array of final index of the original indexes 0:n_ele-1. 
 * N_LTRB - number of elements on the boundary of the grid NOT REQUIRED
 * G_LTRB (N_LTRB, N_LTRB) - Green's matrix of boundary elements
 * inv_r_mu0 (N_LTRB, ) - Inverse of the major radius x mu0 along boundary
 *
 * Outputs: 
 * out (n_ele, ) - psi values on the 2D grid
 */ 
void poisson_solver(
        double* b_vec, 
        double* out 
        )
{

    double psi_bound[N_LTRB];

    hagenow_bound(b_vec, out, psi_bound);
            
    add_bound(psi_bound, b_vec);
    
    solve_tria(b_vec, out);
}  
    
    
