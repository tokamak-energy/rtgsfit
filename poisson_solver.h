#ifndef POISSON_SOLVER_H_   /* Include guard */
#define POISSON_SOLVER_H_

#include <cblas.h>
#include "solve_tria.h"
#include "gradient.h"

void hagenow_bound(int n_r, int n_z, int n_ele, double* lower, double* upper, 
        double* b_vec, int* idx_final, int n_bound, double* g_bound,
        double* inv_r_mu0, double d_row, double d_col, double* psi, 
        double* psi_bound);
        
void add_bound(int n_r, int n_z, int n_ele, int n_bound, double* psi_bound, 
        double* b_vec);
        
void poisson_solver(int n_r, int n_z, int n_ele, double* lower, double* upper, 
        double* b_vec, int* idx_final, int n_bound, double* g_bound,
        double* inv_r_mu0, double dz, double dr, double* out);        

#endif 
