#ifndef RTGSFIT_H_
#define RTGSFIT_H_

#include <cblas.h>
/*#include <lapacke.h>*/
#include "gradient.h"


void rm_coil_from_meas(int n_meas, int n_coil, double* g_meas_coil,
        double* coil_curr, double* meas, double* meas_no_coil);
        
        
void make_basis(int n_row, int n_col, int n_grid, double d_row, 
        double* psi_norm, double* r_grid, double* inv_r_mu0, double* basis);


/*void rtgsfit(int n_meas, double* meas, int n_coil, double* g_meas_coil,*/
/*        double* coil_curr, int n_row, int n_col, double* psi_norm,*/
/*        double d_row, double* g_grid_sens, double* r_grid, double* inv_r_mu0,*/
/*        double* meas_no_coil);*/


#endif

