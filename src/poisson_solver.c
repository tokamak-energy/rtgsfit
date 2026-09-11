#include <cblas.h>
#include "poisson_solver.h"
#include "poisson_fast.h"
#include "gradient.h"
#include "constants.h"
#include <stdio.h>

/* Set once poisson_solver_init() has run (successfully or not). */
static int s_init_attempted = 0;

/*
 * Function: poisson_solver_init
 * Builds the tables of the separable Poisson solver (poisson_fast.c) from
 * POISSON_A, POISSON_B and POISSON_C and validates them. Returns
 * POISSON_FAST_OK (0) on success, otherwise a nonzero status code; the solves
 * then return zero flux (see poisson_solver). Runs automatically on the first
 * solve, but calling it at start-up keeps the one-off table construction (a
 * few milliseconds) out of the first real-time cycle and lets the caller check
 * the status.
 */
int poisson_solver_init(void)
{
    s_init_attempted = 1;
    return poisson_fast_init();
}

/*
 * Function: is_ready
 * Lazily initialises and reports whether the solver is usable.
 */
static int is_ready(void)
{
    if (!s_init_attempted)
    {
        poisson_solver_init();
    }
    return poisson_fast_status() == POISSON_FAST_OK;
}

/*
 * Function: hagenow_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * N_LTRB - number of elements on the boundary of the grid NOT REQUIRED
 * G_LTRB (N_LTRB, N_LTRB) - Green's matrix of boundary elements
 * inv_r_mu0 (N_LTRB, ) - Inverse of the major radius multipled by mu0
 * psi (n_ele, ) - Temporary variable to hold the psi matrix
 *
 * Outputs: 
 * psi_bound (N_LTRB, ) - psi values on the boundary
 */ 
void hagenow_bound(     
        double* b_vec, 
        double* psi,
        double* psi_ltrb
        )
{

    double dpsi_ltrb[N_LTRB];
    int ii;

    if (!is_ready())
    {
        for (ii = 0; ii < N_GRID; ++ii) psi[ii] = 0.0;
        for (ii = 0; ii < N_LTRB; ++ii) psi_ltrb[ii] = 0.0;
        return;
    }

    poisson_fast_solve(b_vec, psi);
    
    gradient_bound(psi, dpsi_ltrb);
    
    for (ii=0; ii<N_LTRB; ++ii)
    {
        dpsi_ltrb[ii] = dpsi_ltrb[ii]*INV_R_LTRB_MU0[ii];
    }
        
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_LTRB, N_LTRB, 1.0, G_LTRB, 
            N_LTRB, dpsi_ltrb, 1, 0.0, psi_ltrb, 1);    
}

    
/*
 * Function: add_bound
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * N_LTRB - number of elements on the boundary of the grid
 * psi_bound (N_LTRB, ) - psi values on the boundary
 *
 * Outputs: 
 * b_vec (n_ele, ) - Current density with flux boundary values added
 */     
void add_bound(
        double* psi_bound, 
        double* b_vec
        )
{
    int ii;
    
    // left
    for (ii=0; ii<N_Z; ++ii) 
    {
        b_vec[ii*N_R] = psi_bound[ii];
    }
        
    // top
    for (ii=1; ii<N_R-1; ++ii) 
    {
        b_vec[ii + (N_Z-1)*N_R] = psi_bound[ii+N_Z-1];
    }
    
    // right    
    for (ii=0; ii<N_Z; ++ii)
    {
        b_vec[ii*N_R + N_R -1]= psi_bound[N_Z + N_R - 3 + (N_Z - ii)];
    }
     
    //bottom     
    for (ii=1; ii<N_R-1; ++ii)
    {
        b_vec[ii]= psi_bound[2*N_Z + N_R + (N_R - ii) - 4];
    }
        
}     



/* Function: poisson_solver
 * determines the boundary flux values of the boundary using the Hagenow method
 *
 * Inputs:
 * N_R - number of radial grid positions
 * N_Z - number of vertical grid poisitons
 * n_ele - number of elements in the grid N_R*N_Z
 * b_vec (n_ele, ) - Current density with zero values on the boundary
 * N_LTRB - number of elements on the boundary of the grid NOT REQUIRED
 * G_LTRB (N_LTRB, N_LTRB) - Green's matrix of boundary elements
 * inv_r_mu0 (N_LTRB, ) - Inverse of the major radius x mu0 along boundary
 *
 * Outputs: 
 * out (n_ele, ) - psi values on the 2D grid
 */ 
void poisson_solver(
        double* b_vec, 
        double* out 
        )
{
    static int reported = 0;

    if (!is_ready())
    {
        /* The Poisson tables could not be built from the constants: return a
         * zero plasma flux (finite, and rtgsfit() then reports ERR_NO_AXIS)
         * rather than an undefined one. poisson_solver_init() gives the status. */
        if (!reported)
        {
            fprintf(stderr, "poisson_solver: solver not initialised (status %d), returning zero flux\n",
                    poisson_fast_status());
            reported = 1;
        }
        for (int ii = 0; ii < N_GRID; ++ii) out[ii] = 0.0;
        return;
    }

    poisson_fast_poisson_solver(b_vec, out);
}  
    
    
