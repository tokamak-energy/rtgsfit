#include "poisson_solver.h"

#include <cblas.h>
#include <string.h>

#include "finite_difference.h"

/* Forward substitution with the unit lower-triangular factor L. */
static void solve_lower(const struct poisson_operator *op, const double b[], double x[])
{
    const int32_t n_r = op->grid.n_r;
    const int32_t n_grid = grid_n_grid(&op->grid);

    /* The first and last n_r rows (the bottom and top edges) are identity rows. */
    for (int32_t i = 0; i < n_r; i++) {
        x[i] = b[i];
        x[n_grid - n_r + i] = b[n_grid - n_r + i];
    }
    for (int32_t i = n_r; i < n_grid - n_r; i++) {
        const double *band = &op->lower_band[(size_t)i * n_r];
        double sum = 0.0;
        for (int32_t j = 0; j < n_r; j++) {
            sum += band[j] * x[i - n_r + j];
        }
        x[i] = b[i] - sum;
    }
}

/* Back substitution with the upper-triangular factor U. */
static void solve_upper(const struct poisson_operator *op, const double b[], double x[])
{
    const int32_t n_r = op->grid.n_r;
    const int32_t n_band = n_r + 2;
    const int32_t n_grid = grid_n_grid(&op->grid);

    for (int32_t i = 0; i < n_r; i++) {
        x[i] = b[i];
        x[n_grid - n_r + i] = b[n_grid - n_r + i];
    }
    for (int32_t i = n_grid - n_r - 1; i >= n_r; i--) {
        const double *band = &op->upper_band[(size_t)i * n_band];
        /* The last super-diagonal of the last interior row would point one past the grid. */
        const int32_t j_end = (i + n_band <= n_grid) ? n_band : n_grid - i;
        double sum = 0.0;
        for (int32_t j = 1; j < j_end; j++) {
            sum += band[j] * x[i + j];
        }
        x[i] = (b[i] - sum) * band[0];
    }
}

int32_t poisson_solve_fixed_boundary(const struct poisson_operator *op,
                                     const double rhs[],
                                     int32_t n_rhs,
                                     double psi[],
                                     int32_t n_psi,
                                     double work[],
                                     int32_t n_work)
{
    const int32_t status = grid_check(&op->grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    const int32_t n_grid = grid_n_grid(&op->grid);
    RTGSFIT_CHECK_SIZE(n_rhs, n_grid);
    RTGSFIT_CHECK_SIZE(n_psi, n_grid);
    RTGSFIT_CHECK_MIN_SIZE(n_work, n_grid);

    for (int32_t i = 0; i < n_grid; i++) {
        psi[i] = rhs[op->permutation[i]];
    }
    solve_lower(op, psi, work);
    solve_upper(op, work, psi);
    return RTGSFIT_OK;
}

/*
 * Visit the edge of the grid in edge order (see `finite_difference_edge_d_dn()`), and
 * copy `edge_values[i_edge]` into `grid_values` at that point, or zero it if `edge_values` is NULL.
 */
static void set_edge(const struct rtgsfit_grid *grid,
                     const double edge_values[],
                     double grid_values[])
{
    const int32_t n_r = grid->n_r;
    const int32_t n_z = grid->n_z;
    int32_t i_edge = 0;

    for (int32_t i_z = 0; i_z < n_z; i_z++) {
        grid_values[grid_index(grid, 0, i_z)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_r = 1; i_r < n_r - 1; i_r++) {
        grid_values[grid_index(grid, i_r, n_z - 1)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_z = n_z - 1; i_z >= 0; i_z--) {
        grid_values[grid_index(grid, n_r - 1, i_z)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_r = n_r - 2; i_r >= 1; i_r--) {
        grid_values[grid_index(grid, i_r, 0)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
}

int32_t poisson_solve_free_boundary(const struct poisson_operator *op,
                                    const double rhs[],
                                    int32_t n_rhs,
                                    double psi[],
                                    int32_t n_psi,
                                    double work[],
                                    int32_t n_work)
{
    const int32_t status = grid_check(&op->grid);
    if (status != RTGSFIT_OK) {
        return status;
    }
    const int32_t n_grid = grid_n_grid(&op->grid);
    const int32_t n_edge = grid_n_edge(&op->grid);
    RTGSFIT_CHECK_SIZE(n_rhs, n_grid);
    RTGSFIT_CHECK_SIZE(n_psi, n_grid);
    RTGSFIT_CHECK_MIN_SIZE(n_work, poisson_work_size(&op->grid));

    double *rhs_work = work;
    double *lu_work = work + n_grid;
    double *edge_d_psi_dn = work + 2 * n_grid;
    double *edge_psi = edge_d_psi_dn + n_edge;

    /* 1. Fixed-boundary solve with psi = 0 on the edge. */
    memcpy(rhs_work, rhs, sizeof(double) * (size_t)n_grid);
    set_edge(&op->grid, NULL, rhs_work);
    int32_t err = poisson_solve_fixed_boundary(op, rhs_work, n_grid, psi, n_grid, lu_work, n_grid);
    if (err != RTGSFIT_OK) {
        return err;
    }

    /*
     * 2. The fixed-boundary solution corresponds to a surface current on the edge proportional to
     *    d(psi)/dn / (mu0 R). Removing that current gives psi on the edge via the edge Green's
     *    function (signs and the d_r, d_z factors are folded into EDGE_G_EDGE at preshot).
     */
    err = finite_difference_edge_d_dn(&op->grid, psi, n_grid, edge_d_psi_dn, n_edge);
    if (err != RTGSFIT_OK) {
        return err;
    }
    for (int32_t i = 0; i < n_edge; i++) {
        edge_d_psi_dn[i] *= op->edge_inv_mu0_r[i];
    }
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_edge, n_edge, 1.0, op->edge_g_edge, n_edge,
                edge_d_psi_dn, 1, 0.0, edge_psi, 1);

    /* 3. Fixed-boundary solve with the free-boundary psi on the edge. */
    set_edge(&op->grid, edge_psi, rhs_work);
    return poisson_solve_fixed_boundary(op, rhs_work, n_grid, psi, n_grid, lu_work, n_grid);
}
