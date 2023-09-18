#include <lapacke.h>
#include <cblas.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

int main()
{
    int n_meas = 70;
    int n_coef = 20;
    int n_grid = n_meas*n_coef;
    double arr[n_grid];
    double x_vec[n_coef];
    double b_vec[n_meas];
    int i_grid, i_coef;
    int info, rank;
    double rcond = -1.0;
    double single_vals[n_coef];  
    int n_time = 10000;
    int i_time;  
    clock_t start, end;
    
    srand(0);
    
    for (i_grid=0; i_grid<n_grid; i_grid++)
    {
        arr[i_grid] = rand();
    }
    
    for (i_coef=0; i_coef<n_coef; i_coef++)
    {
        x_vec[i_coef] = rand();
    }
 
    start = clock();
    for (i_time=0; i_time<n_time; i_time++)
    {
        cblas_dgemv(CblasColMajor, CblasNoTrans, n_meas, n_coef, 1.0, arr, 
            n_meas, x_vec, 1, 0.0, b_vec, 1);  
            
        info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n_meas, n_coef, 1, arr, 
                n_meas, b_vec, n_meas, single_vals, rcond, &rank);
    }
    end = clock();
            
    for (i_coef=0; i_coef<n_coef; i_coef++)
    {
        printf("%f, %f, %f \n", x_vec[i_coef], b_vec[i_coef], x_vec[i_coef]-b_vec[i_coef]);
    }    
    
    printf("time %f ms\n", 1000.0*((double)(end-start))/(n_time*CLOCKS_PER_SEC));
    return 0;
    
}    
    
