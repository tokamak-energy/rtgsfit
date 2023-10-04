#include <stdio.h>
#include <time.h>
#include "constants.h"
#include "rtgsfit.h"


int main()
{
    const char meas_file[] = "../data/meas.txt";
    const char flux_norm_file[] = "../data/flux_norm.txt";
    const char mask_file[] = "../data/mask.txt";
    const char curr_file[] = "../data/coil_curr.txt";  
    
    const char flux_total_file[] = "../data/psi_total.txt";  
    
    double meas[N_MEAS];
    double flux_norm[N_GRID];
    int mask[N_GRID];
    double coil_curr[N_COIL];
    FILE* fptr;
    struct timespec ts_start, ts_end;
    int n_time = 1000;
    double nanosec;
    double flux_total[N_GRID];
    double error[n_time];
    double lcfs_r[N_LCFS_MAX];
    double lcfs_z[N_LCFS_MAX];
    int lcfs_n = 0;
    double coef[N_COEF];
    
    
    int ii;
    
    fptr = fopen(meas_file, "r");
    for (ii = 0; ii<N_MEAS; ii++)
    {
        fscanf(fptr, "%lf", &meas[ii]);
        printf("meas[%d], %lf\n", ii, meas[ii]);
    }
    fclose(fptr);
    
    fptr = fopen(flux_norm_file, "r");
    for (ii = 0; ii<N_GRID; ii++)
    {
        fscanf(fptr, "%lf", &flux_norm[ii]);
        printf("flux_norm[%d], %lf\n", ii, flux_norm[ii]);
    }
    fclose(fptr);
    
    fptr = fopen(mask_file, "r");
    for (ii = 0; ii<N_GRID; ii++)
    {
        fscanf(fptr, "%d", &mask[ii]);
        printf("mask[%d], %d\n", ii, mask[ii]);
    }
    fclose(fptr);
    
    fptr = fopen(curr_file, "r");
    for (ii = 0; ii<N_COIL; ii++)
    {
        fscanf(fptr, "%lf", &coil_curr[ii]);
        printf("coil_curr[%d], %lf\n", ii, coil_curr[ii]);
    }
    fclose(fptr);    
    
    rtgsfit(meas, coil_curr, flux_norm, mask, flux_total, &error[ii], lcfs_r, 
        lcfs_z, &lcfs_n, coef);
                
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts_start);
    for (ii = 0; ii<n_time; ii++)
    {
        rtgsfit(meas, coil_curr, flux_norm, mask, flux_total, &error[ii], lcfs_r, 
                lcfs_z, &lcfs_n, coef);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts_end);
    
    nanosec = (double)(1000000000*(ts_end.tv_sec - ts_start.tv_sec) + 
            (ts_end.tv_nsec - ts_start.tv_nsec))/ ((double) n_time);
    
    for (ii=0; ii<n_time; ii++)
    {
        printf("error[%d]: %lf\n", ii, error[ii]);
    }

    printf("time: %lf ns\n", ts_end.tv_sec, ts_start.tv_sec, ts_end.tv_nsec, ts_start.tv_nsec, nanosec);
        
    return 0;
}
    
