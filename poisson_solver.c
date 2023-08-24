#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>

void swap(double arr[], int n_swap, int idx_a[n_swap], int idx_b[n_swap]){
    
    int ii;
    double tmp;
    for(ii=0; ii<n_swap; ii++){
        tmp = arr[idx_a[ii]];
        arr[idx_a[ii]] = arr[idx_b[ii]];
        arr[idx_b[ii]] = tmp;        
        }
    }


void back_sub_upper_short( 
        size_t n_row, 
        size_t n_col, 
        double upper[n_row][n_col+2], 
        double b_vec[n_row], 
        double out[n_row])
        { 

    double sum;
    int ii, jj;
    
    for (ii=0; ii< n_col; ii++){
        out[ii] = b_vec[ii];
        }
        
    for (ii=n_row-n_col; ii< n_row; ii++){
        out[ii] = b_vec[ii];
        } 
    
    for (ii=(n_row-n_col-1); ii>=n_col; ii--){
        sum = 0.0;
        for (jj=1; jj<(n_col+2); jj++){
            sum = sum + upper[ii][jj]*out[ii+jj];
            }
        out[ii] = (b_vec[ii] - sum) * upper[ii][0];
        }
    }    
    
    
void back_sub_lower_short( 
        size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double b_vec[n_row], 
        double out[n_row])
        {
    /* Solves the a banded lower triangular matrix by back substitution.  
    A*x = b
    where A is the lower triangular matrix, b is b_vec, x is out
    
    Parameters
    ----------
    lower : (n_ele, bandwidth)
        2D matrix of the banded lower triangle matrix
    b_vec : (n_ele, )
        RHS of linear system to solve
    n_band :
        width of the banded matrix
    n_ele :
        length of b_vec
    out : (n_ele, )
        Solution of the linear system
        
    */   

    double sum;
    int ii, jj;
    
    for (ii=0; ii< n_col; ii++){
        out[ii] = b_vec[ii];
        }
        
    for (ii=n_row-n_col; ii< n_row; ii++){
        out[ii] = b_vec[ii];
        }            
    
    for (ii=n_col; ii<(n_row-n_col); ii++){
        sum = 0.0;
        for (jj=0; jj<n_col; jj++){
            sum = sum + lower[ii][jj]*out[ii-n_col+jj];
            }
        out[ii] = b_vec[ii] - sum;

        }
    }


void back_sub_lower_inplace( 
        size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double b_vec[n_row])
        {
    /* Solves the a banded lower triangular matrix by back substitution.  
    A*x = b
    where A is the lower triangular matrix, b is b_vec, x is out
    
    Parameters
    ----------
    lower : (n_ele, bandwidth)
        2D matrix of the banded lower triangle matrix
    b_vec : (n_ele, )
        RHS of linear system to solve
    n_band :
        width of the banded matrix
    n_ele :
        length of b_vec
    out : (n_ele, )
        Solution of the linear system
        
    */   

    double sum;
    int ii, jj;
   
    for (ii=n_col; ii<(n_row-n_col); ii++){
        sum = 0.0;
        for (jj=0; jj<n_col; jj++){
            sum = sum + lower[ii][jj]*b_vec[ii-n_col+jj];
            }
        b_vec[ii] = b_vec[ii] - sum;

        }
    }
    
void back_sub_upper_inplace( 
        size_t n_row, 
        size_t n_col, 
        double upper[n_row][n_col+2], 
        double b_vec[n_row])
    { 

    double sum;
    int ii, jj;
    
    for (ii=(n_row-n_col-1); ii>=n_col; ii--){
        sum = 0.0;
        for (jj=1; jj<(n_col+2); jj++){
            sum = sum + upper[ii][jj]*b_vec[ii+jj];
            }
        b_vec[ii] = (b_vec[ii] - sum) * upper[ii][0];
        }
    }    

void bound_grad(size_t n_r, size_t n_z, size_t n_ele, size_t n_bound, double psi[n_ele], 
        double dpsi_ltrb[n_bound])
    {
    
    int ii;
    
    for (ii=0; ii<n_z; ++ii) 
        {
        dpsi_ltrb[ii] = 2*psi[ii*n_r + 1] - 0.5*psi[ii*n_r + 2];
        }
        
    for (ii=0; ii<n_r; ++ii) 
        {
        dpsi_ltrb[ii+n_z] = 2*psi[ii + n_r] - 0.5*psi[2*n_r + ii];
        }
        
    for (ii=0; ii<n_z; ++ii)
        {
        dpsi_ltrb[ii+n_r+n_z] = 2*psi[ii*n_r + n_r - 2] - 0.5*psi[ii*n_r + n_r - 3];
        }
          
    for (ii=0; ii<n_r; ++ii)
        {
        dpsi_ltrb[ii+n_r + 2*n_z] = 2*psi[(n_z - 2)*n_r + ii] - 0.5*psi[(n_z - 3)*n_r + ii];
        }
        
    }        
        
void matrix_mul(size_t n_row, double g_bound[n_row][n_row], double dpsi_ltrb[n_row], double psi_bound[n_row])
    {
    
    int ii, jj;
    
    for (ii=0; ii<n_row; ii++)
        {
        psi_bound[ii] = g_bound[ii][0] * dpsi_ltrb[0];
        for (jj=1; jj<n_row; jj++)
            {
            psi_bound[ii] = psi_bound[ii] + g_bound[ii][jj] * dpsi_ltrb[jj];
/*            printf("%d %d %f %f \n", ii, jj, g_bound[ii][jj], dpsi_ltrb[jj]);*/
            }
        }
    }
    

    

void solve_tri(size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double upper[n_row][n_col+2], 
        double b_vec[n_row], 
        int n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        double out[n_row]
        )
    {
    
    double b_vec_out[n_row];
    
    swap(b_vec, n_swap, idx_a, idx_b);
    
    back_sub_lower_short(n_row, n_col, lower, b_vec, b_vec_out);
    
    back_sub_upper_short(n_row, n_col, upper, b_vec_out, out);
    
    }

void solve_tri_inplace(size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double upper[n_row][n_col+2], 
        double b_vec[n_row], 
        int n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap] 
        )
    {
    
    swap(b_vec, n_swap, idx_a, idx_b);
    
    back_sub_lower_inplace(n_row, n_col, lower, b_vec);
    
    back_sub_upper_inplace(n_row, n_col, upper, b_vec);
    
    }

void solve_tri_inplace_time(size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double upper[n_row][n_col+2], 
        double b_vec[n_row], 
        int n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap],
        int n_time)
    {
    
    int ii;
    clock_t start, end;
    double time_used;
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        solve_tri_inplace(n_row, n_col, lower, upper, b_vec, n_swap, idx_a, idx_b);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken per iteration: %f ms \n", time_used*1000);  
    
    }
    
    
void solve_tri_time(size_t n_row, 
        size_t n_col, 
        double lower[n_row][n_col], 
        double upper[n_row][n_col+2], 
        double b_vec[n_row], 
        int n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap],
        double out[n_row],
        int n_time)
    {
    
    int ii;
    clock_t start, end;
    double time_used;
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        solve_tri(n_row, n_col, lower, upper, b_vec, n_swap, idx_a, idx_b, out);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken per iteration: %f ms \n", time_used*1000);  
    
    }

void hagenow_bound(
        size_t n_r, 
        size_t n_z, 
        size_t n_ele,
        double lower[n_ele][n_r], 
        double upper[n_ele][n_r+2], 
        double b_vec[n_ele], 
        size_t n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        size_t n_bound,
        double g_bound[n_bound][n_bound],
        double inv_r_mu0[n_bound],
        double out[n_bound]       
        )
    {
    
    double psi[n_ele];
    double dpsi_ltrb[n_bound];
    int ii;
    
    solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, psi);
    
    bound_grad(n_r, n_z, n_ele, n_bound, psi, dpsi_ltrb);
    
    for (ii=0; ii<n_bound; ++ii)
        {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*inv_r_mu0[ii];
        }
        
    matrix_mul(n_bound, g_bound, dpsi_ltrb, out);
    
    }

void add_bound(size_t n_r, size_t n_z, size_t n_ele, size_t n_bound, 
        double psi_bound[n_bound], double b_vec[n_ele])
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
        
    b_vec[0] = 0.0;
    b_vec[n_r - 1] = 0.0;
    b_vec[n_r*(n_z-1)] = 0.0;
    b_vec[n_r*n_z - 1] = 0.0;      
    
    }
     
void gradient_z(
        size_t n_grid, 
        double psi_norm[n_grid], 
        double grad_z[n_grid],
        double dz,
        int n_r,
        int mask[n_grid])
    {
    
    int ii;
    
    for (ii=0; ii<n_r; ii++) 
        {
        if (mask[ii])
            {
            grad_z[ii] = (-3*psi_norm[ii] + 4*psi_norm[ii + n_r] - psi_norm[ii + 2*n_r])/ (2*dz);     
            }
        else
            {
            grad_z[ii] = 0.0;
            }
        }
        
    for (ii=n_r; ii<(n_grid - n_r); ii++) 
        {
        if (mask[ii]) 
            {
            grad_z[ii] = (psi_norm[ii + n_r] - psi_norm[ii - n_r]) / (2*dz);;
            }
        else 
            {
            grad_z[ii] = 0.0;
            }
        }
        
    for (ii=(n_grid - n_r); ii<n_grid; ii++)
        {
        if (mask[ii])
            {
            grad_z[ii] = (3*psi_norm[ii] - 4*psi_norm[ii + n_r] + psi_norm[ii + 2*n_r])/ (2*dz);     
            }
        else
            {
            grad_z[ii] = 0.0;
            }
        }        
    }    
    

void pprime_grid(
        size_t n_grid, 
        double r_grid[n_grid], 
        double psi_norm[n_grid], 
        int mask[n_grid], 
        double pprime[n_grid])
    {
    
    int ii;
    
    for (ii=0; ii<n_grid; ii++) 
        {
        if (mask[ii]) 
            {
            pprime[ii] = (1 - psi_norm[ii]) * r_grid[ii];
            }
        else 
            {
            pprime[ii] = 0.0;
            }
        }       
    }
    
    
void ffprime_grid(
        size_t n_grid, 
        double inv_r_mu0[n_grid], 
        double psi_norm[n_grid], 
        int mask[n_grid], 
        double ffprime[n_grid])
    {
    
    int ii;
    
    for (ii=0; ii<n_grid; ii++) 
        {
        if (mask[ii]) 
            {
            ffprime[ii] = (1 - psi_norm[ii]) * inv_r_mu0[ii];
            }
        else 
            {
            ffprime[ii] = 0.0;
            }
        }
    }
    
    
void fit_coeff(
        size_t n_grid, 
        double r_grid[n_grid], 
        double inv_r_mu0[n_grid], 
        int mask[n_grid], 
        double psi_norm[n_grid], 
        int n_r,
        double dz, 
        double basis[3*n_grid],
        int n_sens,
        double g_grid_sens[n_grid*n_sens],
        double system_mat[3*n_sens],
        int n_basis,
        double sens_val[n_sens])
    {
// 
/*    double basis[3*n_grid];*/
/*    double system_mat[n_grid]*/
/*    int ii;*/
    int info, rank;

    int n_pls = 3;
    double alpha = 1.0;
    double beta = 0.0;
    double single_vals[n_basis];
    double rcond = -1.0;
    
    printf("hi");
    pprime_grid(n_grid, r_grid, psi_norm, mask, basis);
        
    ffprime_grid(n_grid, inv_r_mu0, psi_norm, mask, &basis[n_grid]);
      
    gradient_z(n_grid, psi_norm, &basis[2*n_grid], dz, n_r, mask);      
              
    printf("%d\n", n_pls);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_pls, n_sens, n_grid, alpha, basis, n_grid, g_grid_sens, n_sens, beta, system_mat, n_sens);    
   
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n_sens, n_basis, 1, system_mat, n_sens, sens_val, n_sens, single_vals, rcond, &rank);
    printf("%d\n%d\n%f\n", info, rank, sens_val[0]);
    }
/*        */
/*    for (ii=0; ii<n_sens; ii++) */
/*        {*/
/*        for (jj=0; jj<n_pls; jj++)*/
/*            {*/
/*            system_mat[ii*(n_vess + n_pls) + n_vess +jj] = 0.0;*/
/*            for (kk=0; kk<n_grid; kk++)*/
/*                {*/
/*                system_mat[ii*(n_vess + n_pls) + n_vess +jj] += g_sens_grid[ii*n_sens + kk]*basis[jj*n_basis + kk];*/
/*                }*/
/*            }*/
/*        }*/
/*        */

    

    
void poisson_solver(
        size_t n_r, 
        size_t n_z, 
        size_t n_ele,
        double lower[n_ele][n_r], 
        double upper[n_ele][n_r+2], 
        double b_vec[n_ele], 
        size_t n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        size_t n_bound,
        double g_bound[n_bound][n_bound],
        double inv_r_mu0[n_bound],
        double out[n_bound]       
        )
    {
    
    double psi[n_ele];
    double dpsi_ltrb[n_bound];
    double psi_bound[n_bound];
    int ii;
    
    solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, psi);
    
    bound_grad(n_r, n_z, n_ele, n_bound, psi, dpsi_ltrb);
    
    for (ii=0; ii<n_bound; ++ii)
        {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*inv_r_mu0[ii];
        }
        
    matrix_mul(n_bound, g_bound, dpsi_ltrb, psi_bound);
        
    add_bound(n_r, n_z, n_ele, n_bound, psi_bound, b_vec);
    
    solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, out);
    
    }    

    
void poisson_solver2(
        size_t n_r, 
        size_t n_z, 
        size_t n_ele,
        double lower[n_ele][n_r], 
        double upper[n_ele][n_r+2], 
        double b_vec[n_ele], 
        size_t n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        size_t n_bound,
        double g_bound[n_bound*n_bound],
        double inv_r_mu0[n_bound],
        double out[n_bound]       
        )
    {
    
    double psi[n_ele];
    double dpsi_ltrb[n_bound];
    double psi_bound[n_bound];
    int ii;
    
    solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, psi);
    
    bound_grad(n_r, n_z, n_ele, n_bound, psi, dpsi_ltrb);
    
    for (ii=0; ii<n_bound; ++ii)
        {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*inv_r_mu0[ii];
        }
    
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_bound, n_bound, 1.0, g_bound, n_bound, dpsi_ltrb, 1, 0.0, psi_bound, 1);    
        
    add_bound(n_r, n_z, n_ele, n_bound, psi_bound, b_vec);
    
    solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, out);
    
    }  
    
    
void time_poisson_solver(
        size_t n_r, 
        size_t n_z, 
        size_t n_ele,
        double lower[n_ele][n_r], 
        double upper[n_ele][n_r+2], 
        double b_vec[n_ele], 
        size_t n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        size_t n_bound,
        double g_bound[n_bound][n_bound],
        double inv_r_mu0[n_bound],
        double out[n_bound],
        int n_time
        )
    {    
        
    int ii, jj;
    clock_t start, end;
    double time_used;
    double b_vec_tmp[n_ele];
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        for (jj=0; jj<n_ele; jj++)
            {
            b_vec_tmp[jj] = b_vec[jj];
            }
        poisson_solver(n_r, n_z, n_ele, lower, upper, b_vec_tmp, n_swap, idx_a, 
                idx_b, n_bound, g_bound, inv_r_mu0, out);     
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken per iteration: %f ms \n", time_used*1000);  
    
    }
    
    
void time_poisson_solver2(
        size_t n_r, 
        size_t n_z, 
        size_t n_ele,
        double lower[n_ele][n_r], 
        double upper[n_ele][n_r+2], 
        double b_vec[n_ele], 
        size_t n_swap, 
        int idx_a[n_swap], 
        int idx_b[n_swap], 
        size_t n_bound,
        double g_bound[n_bound][n_bound],
        double inv_r_mu0[n_bound],
        double out[n_bound],
        int n_time
        )
    {    
        
    int ii, jj;
    clock_t start, end;
    double time_used;
    double psi[n_ele];
    double dpsi_ltrb[n_bound];
    double psi_bound[n_bound];

    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, psi);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken solve_tri: %f ms \n", time_used*1000);  
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        bound_grad(n_r, n_z, n_ele, n_bound, psi, dpsi_ltrb);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken bound_grad: %f ms \n", time_used*1000);  
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {        
        for (jj=0; jj<n_bound; ++jj)
            {
            dpsi_ltrb[jj] = dpsi_ltrb[jj]*inv_r_mu0[jj];
            }
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken inv_r_mu0: %f ms \n", time_used*1000);  

    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {           
        matrix_mul(n_bound, g_bound, dpsi_ltrb, psi_bound);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken g_bound: %f ms \n", time_used*1000);  

    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {               
        add_bound(n_r, n_z, n_ele, n_bound, psi_bound, b_vec);
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken psi_bound: %f ms \n", time_used*1000);  

    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {           
        solve_tri(n_ele, n_r, lower, upper, b_vec, n_swap, idx_a, idx_b, out);    
        }
    end = clock();
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken solve_tri: %f ms \n", time_used*1000);  
        
    }
    
void time_cblas(int n_time, size_t n_bound, double g_bound[n_bound*n_bound], double psi_bound[n_bound], double g_bound_2d[n_bound][n_bound])
    {
    
    int ii, jj, kk;
    clock_t start, end;
    double time_used;
    double y[n_bound];
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n_bound, n_bound, 1.0, g_bound, n_bound, psi_bound, 1, 0.0, y, 1);
        }
    end = clock(); 
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken cblas: %f ms \n", time_used*1000);  
    
    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        matrix_mul(n_bound, g_bound_2d, psi_bound, y);
        }
    end = clock(); 
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken matmul: %f ms \n", time_used*1000);  

    start = clock();
    for (ii=0; ii<n_time; ii++) 
        {
        for (jj=0; jj<n_bound; jj++) 
            {
            y[jj] = 0.0;
            for (kk=0; kk < n_bound; kk++)
                {
                y[jj] += g_bound[jj*n_bound + kk] * psi_bound[kk];
                }
            }
        }
    end = clock(); 
    time_used = ((double)(end-start))/(n_time*CLOCKS_PER_SEC);
    printf("time taken matmul2: %f ms \n", time_used*1000);  
            
    }
/*void print_array(size_t n_row, size_t n_col, double array[n_row][n_col])*/
/*    {*/
/*    int ii, jj;*/
/*    for (ii = 0; ii < n_row; ii++){*/
/*        for (jj = 0; jj < n_col; jj++){*/
/*            printf("%f ", array[ii][jj]);*/
/*            }*/
/*        printf("\n");*/
/*        }*/
/*    }*/

/*void back_sub_lower( */
/*        size_t n_row, */
/*        size_t n_col, */
/*        double lower[n_row][n_row], */
/*        double b_vec[n_row], */
/*        double out[n_row])*/
/*        {*/
/*     Solves the a banded lower triangular matrix by back substitution.  */
/*    A*x = b*/
/*    where A is the lower triangular matrix, b is b_vec, x is out*/
/*    */
/*    Parameters*/
/*    ----------*/
/*    lower : (n_ele, bandwidth)*/
/*        2D matrix of the banded lower triangle matrix*/
/*    b_vec : (n_ele, )*/
/*        RHS of linear system to solve*/
/*    n_band :*/
/*        width of the banded matrix*/
/*    n_ele :*/
/*        length of b_vec*/
/*    out : (n_ele, )*/
/*        Solution of the linear system*/
/*        */
/*       */

/*    double sum;*/
/*    int ii, jj;*/
/*    */
/*    for (ii=0; ii< n_col; ii++){*/
/*        out[ii] = b_vec[ii];*/
/*        }*/
/*        */
/*    for (ii=n_row-n_col; ii< n_row; ii++){*/
/*        out[ii] = b_vec[ii];*/
/*        }            */
/*    */
/*    for (ii=n_col; ii<(n_row-n_col); ii++){*/
/*        sum = 0.0;*/
/*        for (jj=(ii-n_col); jj<ii; jj++){*/
/*            sum = sum + lower[ii][jj]*out[jj];*/
/*            }*/
/*        out[ii] = b_vec[ii] - sum;*/

/*        }*/
/*    }    */
/*void back_sub_lower2( */
/*        size_t n_row, */
/*        size_t n_col, */
/*        double** lower, */
/*        double* b_vec, */
/*        double* out)*/
/*        {*/
/*     Solves the a banded lower triangular matrix by back substitution.  */
/*    A*x = b*/
/*    where A is the lower triangular matrix, b is b_vec, x is out*/
/*    */
/*    Parameters*/
/*    ----------*/
/*    lower : (n_ele, bandwidth)*/
/*        2D matrix of the banded lower triangle matrix*/
/*    b_vec : (n_ele, )*/
/*        RHS of linear system to solve*/
/*    n_band :*/
/*        width of the banded matrix*/
/*    n_ele :*/
/*        length of b_vec*/
/*    out : (n_ele, )*/
/*        Solution of the linear system*/
/*        */
/*    */

/*    double sum;*/
/*    int ii, jj;*/
/*    */
/*    for (ii=n_col; ii<(n_row-n_col); ii++){*/
/*        sum = 0.0;*/
/*        for (jj=(ii-n_col); jj<ii; jj++){*/
/*            sum = sum + *( *(lower + ii) + jj) * *(out+ii);*/
/*            if (ii==737){*/
/*                printf("%d, %d, %f\n", ii, jj, *( *(lower + ii) + jj));*/
/*                }*/
/*            }*/
/*        *(out + ii) = *(b_vec + ii) - sum;*/

/*        }*/
/*    printf("\n");*/
/*    printf("ii=%d, jj=%d, lower[ii][jj]=%f*/
/*    }*/

// https://stackoverflow.com/questions/21043682/populate-an-array-at-compile-time-from-file      
    
/*void bound_grad(size_t n_r, size_t n_z, size_t, n_ele, size_t n_bound, double psi[n_ele], */
/*        double dpsi_ltrb[n_bound])*/
/*    {*/
/*    */
/*    int ii;*/
/*    */
/*    for (ii=0; ii<n_row; ++ii) */
/*        {*/
/*        dpsi_ltrb[ii] = 2*psi[ii][1] - 0.5*psi[ii][2];*/
/*        }*/
/*        */
/*    for (ii=0; ii<n_row; ++ii)*/
/*        {*/
/*        dpsi_ltrb[ii+n_col+n_row] = 2*psi[ii][n_col-2] - 0.5*psi[ii][n_col-3];*/
/*        }*/
/*        */
/*    for (ii=0; ii<n_col; ++ii) */
/*        {*/
/*        dpsi_ltrb[ii+n_row] = 2*psi[1][ii] - 0.5*psi[2][ii];*/
/*        }*/
/*        */
/*    for (ii=0; ii<n_col; ++ii)*/
/*        {*/
/*        dpsi_ltrb[ii+n_col+2*n_row] = 2*psi[n_col-2][ii] - 0.5*psi[n_col-3][ii];*/
/*        }*/
/*        */
/*    }        */
    
