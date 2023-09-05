#include "solve_tria.h"

/*
 * Function: permute
 * swaps the values of two array indexes listed in idx_a and idx_b
 *
 * Inputs:
 * arr_orig - array stored as a contiguous block of size (n_row, )
 * n_row - number of index pairs to be swapped
 * idx_final - array of final index of the original indexes 0:n_row-1. 
 *
 * Outputs: 
 * arr_swap - array with the index swapped accordingly of size (n_row, )
 */ 
void permute(
        double* arr_orig,
        int n_row, 
        int* idx_final,
        double* arr_swap
        )
{   
    int ii;
    
    for(ii=0; ii<n_row; ii++)
    {
        arr_swap[ii] = arr_orig[idx_final[ii]];
    }
}


/*
 * Function: back_sub_upper_short
 * Solves the upper triangular matrix (A) with n_rep + 1 superdiagonals by back 
 * substitution.  The matrix A has a repeating structure every n_rep rows. 
 * Solving the equation A x = b
 *
 * Inputs:
 * n_row - number of rows of the triangular matrix
 * n_rep - number of elements within each repetition block
 * upper - 2D array (A) of the upper traingular matrix (n_row, n_superdiagonal+1)
 *         == (n_row, n_rep+2) stored in a 1D contiguous block.  
 * b_vec - RHS 1D array of size (n_row, )
 *
 * Outputs: 
 * out - 1D array containing the solution (x) of length n_row
 */ 
void back_sub_upper_short( 
        int n_row, 
        int n_rep, 
        double* upper, 
        double* b_vec, 
        double* out
        )
{ 

    double sum;
    int ii, jj, idx;
    int n_col = n_rep + 2;
    
    for (ii=0; ii< n_rep; ii++)
    {
        out[ii] = b_vec[ii];
    }
        
    for (ii=n_row-n_rep; ii< n_row; ii++)
    {
        out[ii] = b_vec[ii];
    } 
    
    for (ii=(n_row-n_rep-1); ii>=n_rep; ii--)
    {
        sum = 0.0;
        for (jj=1; jj<n_col; jj++)
        {
            idx = ii*n_col + jj;
            sum = sum + upper[idx]*out[ii+jj];
        }
        out[ii] = (b_vec[ii] - sum) * upper[ii*n_col];
    }
}    
    

/*
 * Function: back_sub_lower_short
 * Solves the a lower triangular matrix with n_rep subdiagonals 
 * by back substitution. The matrix A has a repeating structure every n_rep rows. 
 * Solving the equation A x = b
 * where A is the lower triangular matrix, b is b_vec, x is out
 * 
 * Inputs:
 * n_row - number of rows of the triangular matrix
 * n_rep - number of elements within each repetition block
 * lower - 2D array (A) of the lower traingular matrix (n_row, n_subdiagonals)
 *         == (n_row, n_rep) stored in a 1D contiguous block.  
 * b_vec - RHS 1D array of size (n_row, )
 *
 * Outputs: 
 * out - 1D array containing the solution (x) of length n_row
 */    
void back_sub_lower_short( 
        int n_row, 
        int n_rep, 
        double* lower, 
        double* b_vec, 
        double* out
        )
{

    double sum;
    int ii, jj, idx;
    
    for (ii=0; ii< n_rep; ii++)
    {
        out[ii] = b_vec[ii];
    }
        
    for (ii=n_row-n_rep; ii< n_row; ii++)
    {
        out[ii] = b_vec[ii];
    }            
    
    for (ii=n_rep; ii<(n_row-n_rep); ii++)
    {
        sum = 0.0;
        for (jj=0; jj<n_rep; jj++)
        {
            idx = n_rep*ii + jj;
            sum = sum + lower[idx]*out[ii-n_rep+jj];
        }
        out[ii] = b_vec[ii] - sum;
    }
}

/*
 * Function: solve_tria
 * Solves the linear equation A x = b, by using the the LU decomposition of 
 * matrix A and solving L U x = P b.  Where P is the permutation matrix.
 * 
 * Inputs:
 * n_row - number of rows of the A, L, U, P matrices and the b vector
 * n_rep - number of elements within each repetition block
 * lower - 2D array (L) of the lower traingular matrix (n_row, n_subdiagonals)
 *         == (n_row, n_rep) stored in a 1D contiguous block.  
 * upper - 2D array (U) of the upper traingular matrix (n_row, n_superdiagonal+1)
 *         == (n_row, n_rep+2) stored in a 1D contiguous block.  
 * b_vec - RHS 1D array of size (n_row, )
 * idx_final - list of index to be swapped swap[ii] = orig[idx_final[ii]
 *
 * Outputs: 
 * out - 1D array containing the solution (x) of length n_row
 */     
void solve_tria(
        int n_row, 
        int n_rep, 
        double* lower, 
        double* upper, 
        double* b_vec, 
        int* idx_final, 
        double* out
        )
{
       
    permute(b_vec, n_row, idx_final, out);
    
    back_sub_lower_short(n_row, n_rep, lower, out, b_vec);
    
    back_sub_upper_short(n_row, n_rep, upper, b_vec, out);
    
}



