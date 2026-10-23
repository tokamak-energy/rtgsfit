#include "current_basis.h"

#include "finite_difference.h"

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
                              int32_t n_basis)
{
    const int32_t status = grid_check(grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    const int32_t n_grid = grid_n_grid(grid);
    RTGSFIT_CHECK_SIZE(n_psi_norm, n_grid);
    RTGSFIT_CHECK_SIZE(n_mask, n_grid);
    RTGSFIT_CHECK_SIZE(n_mask_limiter, n_grid);
    RTGSFIT_CHECK_SIZE(n_r_flat, n_grid);
    RTGSFIT_CHECK_SIZE(n_inv_mu0_r, n_grid);
    RTGSFIT_CHECK_SIZE(n_basis, CURRENT_BASIS_N_DOF * n_grid);

    double *basis_p_prime = basis;
    double *basis_ff_prime = basis + n_grid;
    double *basis_shift = basis + 2 * n_grid;

    const int32_t err = finite_difference_d_dz(grid, psi_norm, n_grid, basis_shift, n_grid);
    if (err != RTGSFIT_OK) {
        return err;
    }
    for (int32_t i = 0; i < n_grid; i++) {
        if (mask[i] && mask_limiter[i]) {
            basis_p_prime[i] = (1.0 - psi_norm[i]) * r_flat[i];
            basis_ff_prime[i] = (1.0 - psi_norm[i]) * inv_mu0_r[i];
        } else {
            basis_p_prime[i] = 0.0;
            basis_ff_prime[i] = 0.0;
            basis_shift[i] = 0.0;
        }
    }
    return RTGSFIT_OK;
}
