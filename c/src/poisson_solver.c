#include "poisson_solver.h"

#include <cblas.h>
#include <string.h>

#include "finite_difference.h"

static int32_t operator_check(const struct poisson_operator *op)
{
    if (!op)
    {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    const int32_t status = grid_check(&op->grid);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    const int32_t n_grid = grid_n_grid(&op->grid);
    const int32_t n_edge = grid_n_edge(&op->grid);
    RTGSFIT_CHECK_SIZE(op->n_lower_band, (int64_t)n_grid * op->grid.n_r);
    RTGSFIT_CHECK_SIZE(op->n_upper_band, (int64_t)n_grid * (op->grid.n_r + 2));
    RTGSFIT_CHECK_SIZE(op->n_permutation, n_grid);
    RTGSFIT_CHECK_SIZE(op->n_edge_inv_mu0_r, n_edge);
    RTGSFIT_CHECK_SIZE(op->n_edge_g_edge, (int64_t)n_edge * n_edge);
    RTGSFIT_CHECK_ARRAY(op->lower_band, op->n_lower_band);
    RTGSFIT_CHECK_ARRAY(op->upper_band, op->n_upper_band);
    RTGSFIT_CHECK_ARRAY(op->permutation, op->n_permutation);
    RTGSFIT_CHECK_ARRAY(op->edge_inv_mu0_r, op->n_edge_inv_mu0_r);
    RTGSFIT_CHECK_ARRAY(op->edge_g_edge, op->n_edge_g_edge);
    return RTGSFIT_OK;
}

/* Forward substitution with the unit lower-triangular factor L. */
static int32_t solve_lower(const struct poisson_operator *op,
                           const double b[],
                           const int32_t n_b,
                           double x[],
                           const int32_t n_x)
{
    const int32_t n_r = op->grid.n_r;
    const int32_t n_grid = grid_n_grid(&op->grid);
    RTGSFIT_CHECK_SIZE(n_b, n_grid);
    RTGSFIT_CHECK_SIZE(n_x, n_grid);
    RTGSFIT_CHECK_ARRAY(b, n_b);
    RTGSFIT_CHECK_ARRAY(x, n_x);

    /* The first and last n_r rows (the bottom and top edges) are identity rows. */
    for (int32_t i = 0; i < n_r; i++)
    {
        x[i] = b[i];
        x[n_grid - n_r + i] = b[n_grid - n_r + i];
    }
    for (int32_t i = n_r; i < n_grid - n_r; i++)
    {
        const double *band = &op->lower_band[(size_t)i * n_r];
        double sum = 0.0;
        for (int32_t j = 0; j < n_r; j++)
        {
            sum += band[j] * x[i - n_r + j];
        }
        x[i] = b[i] - sum;
    }
    return RTGSFIT_OK;
}

/* Back substitution with the upper-triangular factor U. */
static int32_t solve_upper(const struct poisson_operator *op,
                           const double b[],
                           const int32_t n_b,
                           double x[],
                           const int32_t n_x)
{
    const int32_t n_r = op->grid.n_r;
    const int32_t n_band = n_r + 2;
    const int32_t n_grid = grid_n_grid(&op->grid);
    RTGSFIT_CHECK_SIZE(n_b, n_grid);
    RTGSFIT_CHECK_SIZE(n_x, n_grid);
    RTGSFIT_CHECK_ARRAY(b, n_b);
    RTGSFIT_CHECK_ARRAY(x, n_x);

    for (int32_t i = 0; i < n_r; i++)
    {
        x[i] = b[i];
        x[n_grid - n_r + i] = b[n_grid - n_r + i];
    }
    for (int32_t i = n_grid - n_r - 1; i >= n_r; i--)
    {
        const double *band = &op->upper_band[(size_t)i * n_band];
        /* The last super-diagonal of the last interior row would point one past the grid. */
        const int32_t j_end = (i + n_band <= n_grid) ? n_band : n_grid - i;
        double sum = 0.0;
        for (int32_t j = 1; j < j_end; j++)
        {
            sum += band[j] * x[i + j];
        }
        x[i] = (b[i] - sum) * band[0];
    }
    return RTGSFIT_OK;
}

int32_t poisson_solve_fixed_boundary(const struct poisson_operator *op,
                                     const double rhs[],
                                     const int32_t n_rhs,
                                     double psi[],
                                     const int32_t n_psi,
                                     double work[],
                                     const int32_t n_work)
{
    int32_t status = operator_check(op);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    const int32_t n_grid = grid_n_grid(&op->grid);
    RTGSFIT_CHECK_SIZE(n_rhs, n_grid);
    RTGSFIT_CHECK_SIZE(n_psi, n_grid);
    RTGSFIT_CHECK_MIN_SIZE(n_work, n_grid);
    RTGSFIT_CHECK_ARRAY(rhs, n_rhs);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(work, n_work);

    for (int32_t i = 0; i < n_grid; i++)
    {
        if (op->permutation[i] < 0 || op->permutation[i] >= n_rhs)
        {
            return RTGSFIT_ERR_ARGUMENT_SIZE;
        }
    }
    for (int32_t i = 0; i < n_grid; i++)
    {
        psi[i] = rhs[op->permutation[i]];
    }
    status = solve_lower(op, psi, n_psi, work, n_grid);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    return solve_upper(op, work, n_grid, psi, n_psi);
}

/*
 * Visit the edge of the grid in edge order (see `finite_difference_edge_d_dn()`), and
 * copy `edge_values[i_edge]` into `grid_values` at that point, or zero it if `edge_values` is NULL.
 */
static int32_t set_edge(const struct rtgsfit_grid *grid,
                        const double edge_values[],
                        const int32_t n_edge_values,
                        double grid_values[],
                        const int32_t n_grid_values)
{
    RTGSFIT_CHECK_SIZE(n_grid_values, grid_n_grid(grid));
    RTGSFIT_CHECK_SIZE(n_edge_values, edge_values ? grid_n_edge(grid) : 0);
    RTGSFIT_CHECK_ARRAY(grid_values, n_grid_values);
    const int32_t n_r = grid->n_r;
    const int32_t n_z = grid->n_z;
    int32_t i_edge = 0;

    for (int32_t i_z = 0; i_z < n_z; i_z++)
    {
        grid_values[grid_index(grid, 0, i_z)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_r = 1; i_r < n_r - 1; i_r++)
    {
        grid_values[grid_index(grid, i_r, n_z - 1)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_z = n_z - 1; i_z >= 0; i_z--)
    {
        grid_values[grid_index(grid, n_r - 1, i_z)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    for (int32_t i_r = n_r - 2; i_r >= 1; i_r--)
    {
        grid_values[grid_index(grid, i_r, 0)] = edge_values ? edge_values[i_edge] : 0.0;
        i_edge++;
    }
    return RTGSFIT_OK;
}

int32_t poisson_solve_free_boundary(const struct poisson_operator *op,
                                    const double rhs[],
                                    const int32_t n_rhs,
                                    double psi[],
                                    const int32_t n_psi,
                                    double work[],
                                    const int32_t n_work)
{
    const int32_t status = operator_check(op);
    if (status != RTGSFIT_OK)
    {
        return status;
    }
    const int32_t n_grid = grid_n_grid(&op->grid);
    const int32_t n_edge = grid_n_edge(&op->grid);
    RTGSFIT_CHECK_SIZE(n_rhs, n_grid);
    RTGSFIT_CHECK_SIZE(n_psi, n_grid);
    RTGSFIT_CHECK_MIN_SIZE(n_work, poisson_work_size(&op->grid));
    RTGSFIT_CHECK_ARRAY(rhs, n_rhs);
    RTGSFIT_CHECK_ARRAY(psi, n_psi);
    RTGSFIT_CHECK_ARRAY(work, n_work);

    double *rhs_work = work;
    double *lu_work = work + n_grid;
    double *edge_d_psi_dn = work + 2 * n_grid;
    double *edge_psi = edge_d_psi_dn + n_edge;

    /* 1. Fixed-boundary solve with psi = 0 on the edge. */
    memcpy(rhs_work, rhs, sizeof(double) * (size_t)n_grid);
    int32_t err = set_edge(&op->grid, NULL, 0, rhs_work, n_grid);
    if (err != RTGSFIT_OK)
    {
        return err;
    }
    err = poisson_solve_fixed_boundary(op, rhs_work, n_grid, psi, n_grid, lu_work, n_grid);
    if (err != RTGSFIT_OK)
    {
        return err;
    }

    /*
     * 2. The fixed-boundary solution corresponds to a surface current on the edge proportional to
     *    d(psi)/dn / (mu0 R). Removing that current gives psi on the edge via the edge Green's
     *    function (signs and the d_r, d_z factors are folded into EDGE_G_EDGE at preshot).
     */
    err = finite_difference_edge_d_dn(&op->grid, psi, n_grid, edge_d_psi_dn, n_edge);
    if (err != RTGSFIT_OK)
    {
        return err;
    }
    for (int32_t i = 0; i < n_edge; i++)
    {
        edge_d_psi_dn[i] *= op->edge_inv_mu0_r[i];
    }
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n_edge, n_edge, 1.0, op->edge_g_edge, n_edge,
                edge_d_psi_dn, 1, 0.0, edge_psi, 1);

    /* 3. Fixed-boundary solve with the free-boundary psi on the edge. */
    err = set_edge(&op->grid, edge_psi, n_edge, rhs_work, n_grid);
    if (err != RTGSFIT_OK)
    {
        return err;
    }
    return poisson_solve_fixed_boundary(op, rhs_work, n_grid, psi, n_grid, lu_work, n_grid);
}
