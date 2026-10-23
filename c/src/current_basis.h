/*
 * Plasma current basis functions.
 *
 * The toroidal current density is modelled as  j_phi(R, Z) = sum_k x_k basis_k(R, Z)  with
 *   basis_0 = (1 - psi_norm) R               p' term
 *   basis_1 = (1 - psi_norm) / (mu0 R)       FF' term
 *   basis_2 = d(psi_norm)/dZ                 rigid vertical shift of the current
 * inside the plasma core (mask == 1 and mask_limiter == 1), and zero elsewhere.
 * This is the transpose of T_yg in eq. 61 of Moret et al. 2015, without the d_r d_z factor.
 */
#ifndef RTGSFIT_CURRENT_BASIS_H
#define RTGSFIT_CURRENT_BASIS_H

#include <stdint.h>

#include "grid.h"

enum { CURRENT_BASIS_N_DOF = 3 };

/*
 * basis[n_basis] is row-major [CURRENT_BASIS_N_DOF][n_grid]. Every other array is a grid array.
 */
int32_t current_basis_compute(const struct rtgsfit_grid *grid,
                              const double psi_norm[],
                              int32_t n_psi_norm,
                              const int32_t mask[],
                              int32_t n_mask,
                              const int32_t mask_limiter[],
                              int32_t n_mask_limiter,
                              const double r_flat[],
                              int32_t n_r_flat,
                              const double inv_mu0_r[],
                              int32_t n_inv_mu0_r,
                              double basis[],
                              int32_t n_basis);

#endif /* RTGSFIT_CURRENT_BASIS_H */
