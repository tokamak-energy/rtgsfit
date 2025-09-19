#ifndef FIND_X_POINT_H_
#define FIND_X_POINT_H_

#include <stdint.h>

double lin_intrp(double* psi, int32_t idx, double dist_to_null_r,
        double dist_to_null_z, double abs_dist_null_r, double abs_dist_null_z);

double lin_intrp_2(double* psi, int32_t idx, double frac_dist_null_r,
        double frac_dist_null_z);

void find_null_in_gradient(double *psi,
        double *opt_r, double *opt_z, double *opt_psi, int32_t *i_opt,
        double *xpt_r, double *xpt_z, double *xpt_psi, int32_t *i_xpt);


void find_lcfs_rz(double *psi, double psi_lcfs,
        double *r_lcfs, double *z_lcfs, int32_t *n_lcfs);


int inside_lcfs(double r_opt, double z_opt, double *r_lcfs, 
        double *z_lcfs, int32_t n_lcfs, int32_t *mask);


void find_null_in_gradient_march(double* flux, double* opt_r, double* opt_z,
        double* opt_flux, int32_t* opt_n, double* xpt_r, double* xpt_z,
        double* xpt_flux, int32_t* xpt_n);

#endif

/*int max_idx(int n_arr, double* arr);*/


/*void find_opt_xpt(double* psi, double* grad_z, double* psi_lcfs,*/
/*        double* psi_axis, double* opt_chosen_r, double* opt_chosen_z);*/
/*void find_xpoint_cross(double a_start_row, double a_end_row, double a_start_col,*/
/*        double a_end_col, double b_start_row, double b_end_row, double b_start_col,*/
/*        double b_end_col, double* col_off, double* row_off, int* flag);*/
