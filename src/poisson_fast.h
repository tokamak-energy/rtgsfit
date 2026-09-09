#ifndef POISSON_FAST_H_   /* Include guard */
#define POISSON_FAST_H_

/* Status codes returned by poisson_fast_init() and poisson_fast_status(). */
#define POISSON_FAST_OK                  0
#define POISSON_FAST_NOT_INITIALISED    -1
#define POISSON_FAST_ERR_GRID_TOO_SMALL  1
#define POISSON_FAST_ERR_ALLOC           2
#define POISSON_FAST_ERR_COEFFICIENTS    3
#define POISSON_FAST_ERR_PIVOT           4
#define POISSON_FAST_ERR_VALIDATION      5

int poisson_fast_init(void);

int poisson_fast_status(void);

double poisson_fast_max_deviation(void);

void poisson_fast_free(void);

void poisson_fast_solve(const double* b_vec, double* out);

void poisson_fast_poisson_solver(double* b_vec, double* out);

#endif
