/*
 * Free-boundary Poisson solver for the poloidal flux produced by the plasma current.
 *
 * We solve  Delta* psi = -mu0 R j_phi  on the rectangular grid, with the free-boundary condition
 * that psi is produced by j_phi alone. The method follows von Hagenow and Lackner, as described
 * in Moret et al. 2015, section 4.2:
 *   1. Solve with psi = 0 on the edge of the grid.
 *   2. From the normal derivative of that solution on the edge, the edge Green's function gives
 *      the true psi on the edge.
 *   3. Solve again, with that psi as the Dirichlet condition on the edge.
 *
 * The 5-point finite-difference operator with identity rows on the edge is factorised once, at
 * preshot, as P A = L U. Because the factors are banded they are stored compactly:
 *   lower_band[n_grid][n_r]       the n_r sub-diagonals of L, whose diagonal is one
 *   upper_band[n_grid][n_r + 2]   column 0 is 1 / diag(U); columns 1..n_r+1 are super-diagonals
 *   permutation[n_grid]           row i of P b is b[permutation[i]]
 */
#ifndef RTGSFIT_POISSON_SOLVER_H
#define RTGSFIT_POISSON_SOLVER_H

#include <stdint.h>

#include "grid.h"

struct poisson_operator {
    struct rtgsfit_grid grid;
    const double *lower_band;     /* [n_grid * n_r] */
    const double *upper_band;     /* [n_grid * (n_r + 2)] */
    const int32_t *permutation;   /* [n_grid] */
    const double *edge_inv_mu0_r; /* [n_edge] 1 / (mu0 R) on the edge */
    const double *edge_g_edge;    /* [n_edge * n_edge] psi on the edge per unit edge current */
};

/* Number of doubles of workspace needed by `poisson_solve_free_boundary()`. */
static inline int32_t poisson_work_size(const struct rtgsfit_grid *grid)
{
    return 2 * grid_n_grid(grid) + 2 * grid_n_edge(grid);
}

/*
 * Solve A psi = rhs with the banded LU factors. This is the fixed-boundary solve: the edge rows
 * of A are identity rows, so psi on the edge equals rhs on the edge.
 *
 * rhs[n_rhs], psi[n_psi] and work[n_work] are grid arrays (length n_grid).
 */
int32_t poisson_solve_fixed_boundary(const struct poisson_operator *op,
                                     const double rhs[],
                                     int32_t n_rhs,
                                     double psi[],
                                     int32_t n_psi,
                                     double work[],
                                     int32_t n_work);

/*
 * Solve the free-boundary problem.
 *
 * rhs[n_rhs] is -mu0 R j_phi d_z^2 on the grid; its edge values are ignored (treated as zero).
 * psi[n_psi] receives the plasma psi on the grid. work[n_work] must hold at least
 * `poisson_work_size()` doubles.
 */
int32_t poisson_solve_free_boundary(const struct poisson_operator *op,
                                    const double rhs[],
                                    int32_t n_rhs,
                                    double psi[],
                                    int32_t n_psi,
                                    double work[],
                                    int32_t n_work);

#endif /* RTGSFIT_POISSON_SOLVER_H */
