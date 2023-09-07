#include "gradient.h"    
        
/*
 * Function: gradient_row
 * computes the gradient in the row direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * n_row - number of rows i.e. n_grid = (n_row * n_col)
 * n_col - number of columns i.e. n_grid = (n_row * n_col)
 * d_row - spatial distance between rows 
 *
 * Outputs: 
 * grad_row - gradient of the 2D matrix in the row direction
 */ 
void gradient_row(
        double* arr, 
        int n_row,
        int n_col,
        double d_row,
        double* grad_row
        )       
{       
 
    int i_row, i_col, idx;
    double inv_d_row = 1.0/(2.0*d_row);
    
    for (i_row=0; i_row<n_row; i_row++) 
    {
        for (i_col=0; i_col<n_col; i_col++)
        {
            
            idx = i_row*n_col + i_col;
                
            if (i_row > 0 && i_row < (n_row - 1))
            {
                grad_row[idx] = (arr[idx + n_col] - arr[idx - n_col])*inv_d_row;
            }
            else if (i_row == 0)
            {
                grad_row[idx] = (-3*arr[idx] + 4*arr[idx + n_col] - arr[idx + 2*n_col])*inv_d_row; 
            }
            else
            {
                grad_row[idx] = (3*arr[idx] - 4*arr[idx - n_col] + arr[idx - 2*n_col])*inv_d_row;   
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
 * n_row - number of rows i.e. n_grid = (n_row * n_col)
 * n_col - number of columns i.e. n_grid = (n_row * n_col)
 * d_col - spatial distance between columns 
 *
 * Outputs: 
 * grad_col - gradient of the 2D matrix in the column direction
 */ 
void gradient_col(
        double* arr, 
        int n_row,
        int n_col,
        double d_col,
        double* grad_col
        )
{
    
    int i_row, i_col, idx;
    double inv_d_col = 1.0/(2.0*d_col);
    
    for (i_row=0; i_row<n_row; i_row++) 
    {
        for (i_col=0; i_col<n_col; i_col++)
        {
            
            idx = i_row*n_col + i_col;

            if (i_col > 0 && i_col < (n_col - 1))
            {
                grad_col[idx] = (arr[idx + 1] - arr[idx - 1])*inv_d_col;
            }
            else if (i_col == 0)
            {
                grad_col[idx] = (-3*arr[idx] + 4*arr[idx + 1] - arr[idx + 2])*inv_d_col;   
            }
            else
            {
                grad_col[idx] = (3*arr[idx] - 4*arr[idx - 1] + arr[idx - 2])*inv_d_col;
            }
        }
    }
}
    
            
/*
 * Function: hessian_row_row
 * computes the Hessian in the row direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * n_row - number of rows i.e. n_grid = (n_row * n_col)
 * n_col - number of columns i.e. n_grid = (n_row * n_col)
 * d_row - spatial distance between rows 
 *
 * Outputs: 
 * hess_row_row - second derivative (Hessian) of the 2D matrix in the row 
 * direction
 */         
void hessian_row_row(
        double* arr, 
        int n_row,
        int n_col,
        double d_row,
        double* hess_row_row
        )
{
    
    int i_row, i_col, idx;
    double inv_d_row2 = 1.0/(d_row*d_row);

    for (i_row=0; i_row<n_row; i_row++) 
    {
        for (i_col=0; i_col<n_col; i_col++)
        {
            
            idx = i_row*n_col + i_col;
                
            if (i_row > 0 && i_row < (n_row - 1))
            {
                hess_row_row[idx] = (arr[idx + n_col] - 2*arr[idx] + 
                        arr[idx - n_col])*inv_d_row2;
            }
            else if (i_row == 0)
            {
                hess_row_row[idx] = (2*arr[idx] - 5*arr[idx + n_col] + 
                        4*arr[idx + 2*n_col] - arr[idx + 3*n_col])*inv_d_row2;    
            }
            else
            {
                hess_row_row[idx] = (2*arr[idx] - 5*arr[idx - n_col] + 
                        4*arr[idx - 2*n_col] - arr[idx - 3*n_col])*inv_d_row2;    
            }  
        }
    }
}   

   
/*
 * Function: hessian_col_col
 * computes the Hessian in the col direction of a 2D matrix, stored in a 1D 
 * array
 *
 * Inputs:
 * arr - 2D array stored as a contiguous block (i.e. 1D array)
 * n_row - number of rows i.e. n_grid = (n_row * n_col)
 * n_col - number of columns i.e. n_grid = (n_row * n_col)
 * d_col - spatial distance between cols 
 *
 * Outputs: 
 * hess_col_col - second derivative (Hessian) of the 2D matrix in the column
 * direction
 */         
void hessian_col_col(
        double* arr, 
        int n_row,
        int n_col,
        double d_col,
        double* hess_col_col
        )
{
    
    int i_row, i_col, idx;
    double inv_d_col2 = 1.0/(d_col*d_col);
    
    for (i_row=0; i_row<n_row; i_row++) 
    {
        for (i_col=0; i_col<n_col; i_col++)
        {
            
            idx = i_row*n_col + i_col;

            if (i_col > 0 && i_col < (n_col - 1))
            {
                hess_col_col[idx] = (arr[idx + 1] - 2*arr[idx] + 
                        arr[idx - 1])*inv_d_col2;
            }
            else if (i_col == 0)
            {
                hess_col_col[idx] = (2*arr[idx] - 5*arr[idx + 1] + 
                        4*arr[idx + 2] - arr[idx + 3])*inv_d_col2;
            }
            else
            {
                hess_col_col[idx] = (2*arr[idx] - 5*arr[idx - 1] + 
                        4*arr[idx - 2] - arr[idx - 3])*inv_d_col2;
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
 * n_row - number of rows i.e. n_grid = (n_row * n_col)
 * n_col - number of columns i.e. n_grid = (n_row * n_col)
 *
 * Outputs: 
 * grad_bound - numerator of inward normal derivative at the boundary stored in 
 * the order left top right bottom
 */     
void gradient_bound(
        double* arr,
        int n_row, 
        int n_col, 
        double d_row, 
        double d_col,
        double* grad_bound
        )
{
    
    int ii;
    double d_row_col = (d_row/d_col);
    double d_col_row = (d_col/d_row);
    
    for (ii=0; ii<n_row; ++ii) 
    {
        grad_bound[ii] = d_row_col*(2*arr[ii*n_col + 1] - 0.5*arr[ii*n_col + 2]);
    }
        
    for (ii=0; ii<n_col; ++ii) 
    {
        grad_bound[ii + n_row] = d_col_row*(2*arr[ii + n_col] - 0.5*arr[2*n_col + ii]);
    }
        
    for (ii=0; ii<n_row; ++ii)
    {
        grad_bound[ii + n_col + n_row] = d_row_col*(2*arr[ii*n_col + n_col - 2] - 
                0.5*arr[ii*n_col + n_col - 3]);
    }
          
    for (ii=0; ii<n_col; ++ii)
    {
        grad_bound[ii + n_col + 2*n_row] = d_col_row*(2*arr[(n_row - 2)*n_col + ii] - 
                0.5*arr[(n_row - 3)*n_col + ii]);
    }    
}  
