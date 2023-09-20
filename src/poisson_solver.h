#ifndef POISSON_SOLVER_H_   /* Include guard */
#define POISSON_SOLVER_H_

void hagenow_bound(double* b_vec, double* psi, double* psi_ltrb);
        
void add_bound(double* psi_bound, double* b_vec);
        
void poisson_solver(double* b_vec, double* out);        

#endif 
