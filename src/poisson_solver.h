#ifndef POISSON_SOLVER_H_   /* Include guard */
#define POISSON_SOLVER_H_

#define POISSON_SOLVER_METHOD_UNSET (-1)
#define POISSON_SOLVER_METHOD_LU     0
#define POISSON_SOLVER_METHOD_FAST   1

int poisson_solver_init(void);

int poisson_solver_method(void);

void hagenow_bound(double* b_vec, double* psi, double* psi_ltrb);
        
void add_bound(double* psi_bound, double* b_vec);
        
void poisson_solver(double* b_vec, double* out);        

void poisson_solver_lu(double* b_vec, double* out);

#endif 
