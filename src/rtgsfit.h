#ifndef RTGSFIT_H_
#define RTGSFIT_H_

int max_idx(int n_arr, double* arr);
        
        
void rm_coil_from_meas(double* coil_curr, double* meas, double* meas_no_coil);

              
void make_basis(double* psi_norm, int* mask, double* basis);


void normalise_flux(double* flux_total, double flux_lcfs, 
        double flux_axis,int* mask, double* flux_norm);
        
        
void rtgsfit(double* meas, double* coil_curr, double* flux_norm, int* mask,
        double* flux_total, double* error);

#endif

