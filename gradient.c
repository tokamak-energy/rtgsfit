#include "gradient.h"    
#include "constants.h"
#include <stdio.h>

/*
 * Function: gradient_row
 * computes the gradient in the row direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * N_Z - number of rows i.e. n_grid = (N_Z * N_R)
 * N_R - number of columns i.e. n_grid = (N_Z * N_R)
 * d_row - spatial distance between rows 
 *
 * Outputs: 
 * grad_z - gradient of the 2D matrix in the row direction
 */ 
void gradient_z(
        double* arr, 
        double* grad_z
        )       
{       
 
    int i_row, i_col, idx;
    double inv_dz = 1.0/(2.0*DZ);
    
    for (i_row=0; i_row<N_Z; i_row++) 
    {
        for (i_col=0; i_col<N_R; i_col++)
        {
            
            idx = i_row*N_R + i_col;
                
            if (i_row > 0 && i_row < (N_Z - 1))
            {
                grad_z[idx] = (arr[idx + N_R] - arr[idx - N_R])*inv_dz;
            }
            else if (i_row == 0)
            {
                grad_z[idx] = (-3*arr[idx] + 4*arr[idx + N_R] - arr[idx + 2*N_R])*inv_dz; 
            }
            else
            {
                grad_z[idx] = (3*arr[idx] - 4*arr[idx - N_R] + arr[idx - 2*N_R])*inv_dz;   
            }
        }
    }
}  
          

/*
 * Function: gradient_col
 * computes the gradient in the column direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * N_Z - number of rows i.e. n_grid = (N_Z * N_R)
 * N_R - number of columns i.e. n_grid = (N_Z * N_R)
 * d_col - spatial distance between columns 
 *
 * Outputs: 
 * grad_r - gradient of the 2D matrix in the column direction
 */ 
void gradient_r(
        double* arr, 
        double* grad_r
        )
{
    
    int i_row, i_col, idx;
    double inv_dr = 1.0/(2.0*DR);
    
    for (i_row=0; i_row<N_Z; i_row++) 
    {
        for (i_col=0; i_col<N_R; i_col++)
        {
            
            idx = i_row*N_R + i_col;

            if (i_col > 0 && i_col < (N_R - 1))
            {
                grad_r[idx] = (arr[idx + 1] - arr[idx - 1])*inv_dr;
            }
            else if (i_col == 0)
            {
                grad_r[idx] = (-3*arr[idx] + 4*arr[idx + 1] - arr[idx + 2])*inv_dr;   
            }
            else
            {
                grad_r[idx] = (3*arr[idx] - 4*arr[idx - 1] + arr[idx - 2])*inv_dr;
            }
        }
    }
}
    
            
/*
 * Function: hessiaN_Z_row
 * computes the Hessian in the row direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * N_Z - number of rows i.e. n_grid = (N_Z * N_R)
 * N_R - number of columns i.e. n_grid = (N_Z * N_R)
 * d_row - spatial distance between rows 
 *
 * Outputs: 
 * hess_zz - second derivative (Hessian) of the 2D matrix in the row 
 * direction
 */         
void hessian_zz(
        double* arr, 
        double* hess_zz
        )
{
    
    int i_row, i_col, idx;
    double inv_dz2 = 1.0/(DZ*DZ);

    for (i_row=0; i_row<N_Z; i_row++) 
    {
        for (i_col=0; i_col<N_R; i_col++)
        {
            
            idx = i_row*N_R + i_col;
                
            if (i_row > 0 && i_row < (N_Z - 1))
            {
                hess_zz[idx] = (arr[idx + N_R] - 2*arr[idx] + 
                        arr[idx - N_R])*inv_dz2;
            }
            else if (i_row == 0)
            {
                hess_zz[idx] = (2*arr[idx] - 5*arr[idx + N_R] + 
                        4*arr[idx + 2*N_R] - arr[idx + 3*N_R])*inv_dz2;    
            }
            else
            {
                hess_zz[idx] = (2*arr[idx] - 5*arr[idx - N_R] + 
                        4*arr[idx - 2*N_R] - arr[idx - 3*N_R])*inv_dz2;    
            }  
        }
    }
}   

   
/*
 * Function: hessiaN_R_col
 * computes the Hessian in the col direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * N_Z - number of rows i.e. n_grid = (N_Z * N_R)
 * N_R - number of columns i.e. n_grid = (N_Z * N_R)
 * d_col - spatial distance between cols 
 *
 * Outputs: 
 * hess_rr - second derivative (Hessian) of the 2D matrix in the column
 * direction
 */         
void hessian_rr(
        double* arr, 
        double* hess_rr
        )
{
    
    int i_row, i_col, idx;
    double inv_dr2 = 1.0/(DR*DR);
    
    for (i_row=0; i_row<N_Z; i_row++) 
    {
        for (i_col=0; i_col<N_R; i_col++)
        {
            
            idx = i_row*N_R + i_col;

            if (i_col > 0 && i_col < (N_R - 1))
            {
                hess_rr[idx] = (arr[idx + 1] - 2*arr[idx] + 
                        arr[idx - 1])*inv_dr2;
            }
            else if (i_col == 0)
            {
                hess_rr[idx] = (2*arr[idx] - 5*arr[idx + 1] + 
                        4*arr[idx + 2] - arr[idx + 3])*inv_dr2;
            }
            else
            {
                hess_rr[idx] = (2*arr[idx] - 5*arr[idx - 1] + 
                        4*arr[idx - 2] - arr[idx - 3])*inv_dr2;
            }          
        }   
    }
}   


/*
 * Function: gradient_bound
 * calculates the numerator the gradient at the boundary in the inward
 * normal direction.  See Liuqe paper
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * N_Z - number of rows i.e. n_grid = (N_Z * N_R)
 * N_R - number of columns i.e. n_grid = (N_Z * N_R)
 *
 * Outputs: 
 * grad_bound - numerator of inward normal derivative at the boundary stored in 
 * the order left top right bottom
 */     
void gradient_bound(
        double* arr,
        double* grad_bound
        )
{
    
    int ii;
    double dz_dr = (DZ/DR);
    double dr_dz = (DR/DZ);
    
    // left
    for (ii=0; ii<N_Z; ++ii) 
    {
        grad_bound[ii] = dz_dr*(2*arr[ii*N_R + 1] - 0.5*arr[ii*N_R + 2]);
        // printf("%d, %f, %f \n", ii, R_GRID[ii*N_R], Z_GRID[ii*N_R]);
    }

    // top
    for (ii=1; ii<N_R-1; ++ii)
    {
        grad_bound[ii + N_Z - 1] = dr_dz*(2*arr[(N_Z - 2)*N_R + ii] - 
                0.5*arr[(N_Z - 3)*N_R + ii]);
        // printf("%d, %f, %f \n", ii + N_Z - 1,  R_GRID[(N_Z - 1)*N_R + ii], Z_GRID[(N_Z - 1)*N_R + ii]);
    }  
        
    // right  
    for (ii=0; ii<N_Z; ++ii)
    {
        grad_bound[ii + N_R + N_Z - 2] = dz_dr*(2*arr[(N_Z-ii-1)*N_R + N_R - 2] - 
                0.5*arr[(N_Z-ii-1)*N_R + N_R - 3]);
        // printf("%d, %f, %f \n", ii + N_R + N_Z - 2, R_GRID[(N_Z-ii-1)*N_R + N_R - 1], Z_GRID[(N_Z-ii-1)*N_R + N_R - 1]);
    }
    
    // bottom
    for (ii=1; ii<N_R-1; ++ii) 
    {
        grad_bound[ii + N_R + 2*N_Z - 3] = dr_dz*(2*arr[(N_R - ii - 1) + N_R] - 0.5*arr[2*N_R + (N_R - ii - 1)]);
        // printf("%d, %f, %f \n", ii + N_R + 2*N_Z - 3, R_GRID[(N_R - ii - 1)], Z_GRID[(N_R - ii - 1)]);
    }  
}  
