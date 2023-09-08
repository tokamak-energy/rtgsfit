#ifndef FIND_X_POINT_H_
#define FIND_X_POINT_H_

#include <math.h>

void find_zero_on_edge(double *grad_patch, double thresh, double *cross_row, 
        double *cross_col, int *count);
        
        
void find_null_in_gradient(double dr, double dz, int n_row, int n_col,
        double *r_arr, double *z_arr, double *psi, double *grad_r,
        double *grad_z, double *hess_rr, double *hess_zz, double *hess_rz,
        double *opt_r, double *opt_z, double *opt_psi, int *i_opt,
        double *xpt_r, double *xpt_z, double *xpt_psi, int *i_xpt);
        
        
void find_lcfs(double dr, double dz, int n_row, int n_col, double *r_arr,
        double *z_arr, double *psi, double psi_bound, double thresh, 
        double *r_lcfs, double *z_lcfs, double *n_lcfs);

        
void inside_lcfs(double dr, double dz, int n_row, int n_col, double r_opt,
        double z_opt, double thresh, double *r_arr, double *z_arr,
        double *psi, double psi_bound, double *r_lcfs, double *z_lcfs,
        double *n_lcfs, int *idx, int *n_idx);                    

#endif
