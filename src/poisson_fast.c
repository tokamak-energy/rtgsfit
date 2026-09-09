/*
 * File: poisson_fast.c
 * --------------------
 * Direct solver for the RT-GSFit Poisson (Grad-Shafranov) system.
 *
 * The operator is the LIUQE eq. (45) five-point stencil on the uniform (R, Z)
 * grid, supplied in constants.c as POISSON_A, POISSON_B and POISSON_C, with
 * identity rows on the grid boundary (Dirichlet data is passed in the boundary
 * entries of the right-hand side b):
 *
 *   x[i-1][j] + x[i+1][j] + POISSON_B[j] x[i][j-1] + POISSON_A[j] x[i][j+1] - POISSON_C[j] x[i][j] = b[i][j]
 *
 * for interior points (i = 1..N_Z-2, j = 1..N_R-2). The vertical coupling is
 * the same on every row and the other coefficients depend on the column j
 * only, so the vertical direction can be diagonalised with a discrete sine
 * transform (DST-I with period N = N_Z-1), after which every vertical mode k is
 * an independent tridiagonal system in the radial direction:
 *
 *   b_j X[k][j-1] + (2 cos(pi k / N) - c_j) X[k][j] + a_j X[k][j+1] = B[k][j]
 *
 * This is the classic "Fourier plus tridiagonal" direct method for separable
 * elliptic problems (Hockney 1965; Jardin, Computational Methods in Plasma
 * Physics, 2010, sec. 3.3). The transform is a dense matrix product (two
 * cblas_dgemm calls per direction, split by mode parity), so any grid size
 * works and the cost is O(N_Z^2 N_R) in small cache-resident products.
 *
 * At initialisation the coefficients are checked to be finite, every mode's
 * tridiagonal is checked to have safe pivots, and the solver is checked
 * against the stencil residual and against a plain two-solve Hagenow flow on
 * test right-hand sides. On failure poisson_fast_init() returns a nonzero
 * status and poisson_solver() returns zero flux.
 *
 * The Hagenow boundary calculation needs the first (zero-boundary) solution only
 * on the two rows and two columns next to the grid boundary, and the second
 * solve differs from the first only through the boundary values. The Hagenow
 * flow therefore does one forward transform, evaluates the first solution just
 * where gradient_bound() reads it, transforms the boundary-only correction
 * analytically, and does one inverse transform for the final flux.
 *
 * Storage layout: the sine transform is split into odd and even modes so that
 * the (m x m) transform matrix becomes two (m/2 x m/2) blocks (S[k][N-i] =
 * (-1)^(k+1) S[k][i]); transformed data are stored with the odd modes first,
 * then the even modes ("split mode order"). All work buffers are allocated once
 * in poisson_fast_init(); the real-time path does not allocate.
 */
#include "poisson_fast.h"
#include "constants.h"
#include "gradient.h"
#include "poisson_solver.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Relative tolerance for the tridiagonal pivots. */
#define PIVOT_RTOL 1.0e-10
/* Relative tolerance for the residual and self-consistency checks at initialisation. */
#define VALIDATION_RTOL 1.0e-9
/* Rows and columns next to the boundary that gradient_bound() reads. */
#define N_EDGE 2

static int s_status = POISSON_FAST_NOT_INITIALISED;
static double s_max_dev = -1.0;

/* Grid geometry: m interior rows, n interior columns, transform period N. */
static int s_m, s_n, s_N, s_h, s_hp, s_m_odd, s_m_even;
static int s_n_edge_rows, s_edge_rows[2 * N_EDGE];
static int s_n_edge_cols, s_edge_cols[2 * N_EDGE];

/* Operator: vertical coupling e, and a_j, b_j, d_j = -c_j for interior columns. */
static double s_e;
static double *s_a, *s_b, *s_d;      /* n each, indexed by interior column */

/* Sine transform blocks and per-mode tridiagonal factors (split mode order). */
static double *s_s_odd;              /* m_odd x hp */
static double *s_s_even;             /* m_even x h */
static double *s_s_edge_rows;        /* n_edge_rows x m, scaled inverse transform rows */
static double *s_sin_k;              /* m, sin(pi k / N) in split mode order */
static double *s_mult;               /* m x n */
static double *s_inv_w;              /* m x n */

/* Work buffers. */
static double *s_f;                  /* m x n interior right-hand side */
static double *s_fe;                 /* hp x n */
static double *s_fo;                 /* h x n */
static double *s_fhat;               /* m x n, split mode order */
static double *s_fhat2;              /* m x n, split mode order */
static double *s_ge;                 /* hp x n */
static double *s_go;                 /* h x n */
static double *s_psi1;               /* N_GRID, first Hagenow solution (partial) */
static double *s_bnd;                /* N_GRID, boundary values from add_bound */
static double *s_col_in;             /* m x 2 N_EDGE */
static double *s_col_ge;             /* hp x 2 N_EDGE */
static double *s_col_go;             /* h x 2 N_EDGE */
static double *s_row_out;            /* 2 N_EDGE x n */
static double *s_vec_e;              /* hp */
static double *s_vec_o;              /* h */
static double *s_vec_hat;            /* m */

static int is_boundary(int i, int j)
{
    return (i == 0 || i == N_Z - 1 || j == 0 || j == N_R - 1);
}

static int mode_of_row(int kk)
{
    /* Split mode order: rows 0..m_odd-1 hold k = 1, 3, 5, ...; the rest hold k = 2, 4, ... */
    return (kk < s_m_odd) ? (2 * kk + 1) : (2 * (kk - s_m_odd) + 2);
}

/* Distinct indices among the first and last N_EDGE of 0..count-1. */
static int edge_indices(int count, int* idx)
{
    int n_idx = 0;
    for (int q = 0; q < 2 * N_EDGE; q++)
    {
        const int cand = (q < N_EDGE) ? q : (count - 1 - (q - N_EDGE));
        int dup = 0;
        if (cand < 0 || cand >= count) continue;
        for (int t = 0; t < n_idx; t++)
        {
            if (idx[t] == cand) dup = 1;
        }
        if (!dup) idx[n_idx++] = cand;
    }
    return n_idx;
}

void poisson_fast_free(void)
{
    free(s_a); free(s_b); free(s_d);
    free(s_s_odd); free(s_s_even); free(s_s_edge_rows); free(s_sin_k); free(s_mult); free(s_inv_w);
    free(s_f); free(s_fe); free(s_fo); free(s_fhat); free(s_fhat2); free(s_ge); free(s_go);
    free(s_psi1); free(s_bnd); free(s_col_in); free(s_col_ge); free(s_col_go); free(s_row_out);
    free(s_vec_e); free(s_vec_o); free(s_vec_hat);
    s_a = s_b = s_d = NULL;
    s_s_odd = s_s_even = s_s_edge_rows = s_sin_k = s_mult = s_inv_w = NULL;
    s_f = s_fe = s_fo = s_fhat = s_fhat2 = s_ge = s_go = NULL;
    s_psi1 = s_bnd = s_col_in = s_col_ge = s_col_go = s_row_out = NULL;
    s_vec_e = s_vec_o = s_vec_hat = NULL;
    s_status = POISSON_FAST_NOT_INITIALISED;
    s_max_dev = -1.0;
}

/*
 * Function: load_operator
 * Copies the interior stencil coefficients from constants.c and checks that
 * they are finite. Returns a status code.
 */
static int load_operator(void)
{
    s_e = 1.0;
    for (int jj = 0; jj < s_n; jj++)
    {
        const int j = jj + 1;
        s_a[jj] = POISSON_A[j];
        s_b[jj] = POISSON_B[j];
        s_d[jj] = -POISSON_C[j];
        if (!isfinite(s_a[jj]) || !isfinite(s_b[jj]) || !isfinite(s_d[jj]))
        {
            return POISSON_FAST_ERR_COEFFICIENTS;
        }
    }
    return POISSON_FAST_OK;
}

/*
 * Function: build_tables
 * Builds the sine transform blocks and the tridiagonal factors for every mode.
 */
static int build_tables(void)
{
    const double pi_over_n = M_PI / (double)s_N;
    const double scale = 2.0 / (double)s_N;

    for (int kk = 0; kk < s_m_odd; kk++)
    {
        const int k = 2 * kk + 1;
        for (int ii = 0; ii < s_hp; ii++)
        {
            s_s_odd[kk * s_hp + ii] = sin(pi_over_n * (double)k * (double)(ii + 1));
        }
    }
    for (int kk = 0; kk < s_m_even; kk++)
    {
        const int k = 2 * kk + 2;
        for (int ii = 0; ii < s_h; ii++)
        {
            s_s_even[kk * s_h + ii] = sin(pi_over_n * (double)k * (double)(ii + 1));
        }
    }
    for (int q = 0; q < s_n_edge_rows; q++)
    {
        const int i = s_edge_rows[q] + 1;
        for (int kk = 0; kk < s_m; kk++)
        {
            s_s_edge_rows[q * s_m + kk] = scale * sin(pi_over_n * (double)i * (double)mode_of_row(kk));
        }
    }
    for (int kk = 0; kk < s_m; kk++)
    {
        s_sin_k[kk] = sin(pi_over_n * (double)mode_of_row(kk));
    }

    for (int kk = 0; kk < s_m; kk++)
    {
        const int k = mode_of_row(kk);
        const double lam = 2.0 * s_e * cos(pi_over_n * (double)k);
        double* mult = &s_mult[kk * s_n];
        double* inv_w = &s_inv_w[kk * s_n];
        double w = s_d[0] + lam;
        double pivot_scale = fabs(s_d[0]) + fabs(lam) + fabs(s_a[0]) + fabs(s_b[0]);

        if (fabs(w) <= PIVOT_RTOL * pivot_scale) return POISSON_FAST_ERR_PIVOT;
        mult[0] = 0.0;
        inv_w[0] = 1.0 / w;
        for (int jj = 1; jj < s_n; jj++)
        {
            mult[jj] = s_b[jj] * inv_w[jj - 1];
            w = s_d[jj] + lam - mult[jj] * s_a[jj - 1];
            pivot_scale = fabs(s_d[jj]) + fabs(lam) + fabs(s_a[jj]) + fabs(s_b[jj]);
            if (fabs(w) <= PIVOT_RTOL * pivot_scale) return POISSON_FAST_ERR_PIVOT;
            inv_w[jj] = 1.0 / w;
        }
    }
    return POISSON_FAST_OK;
}

/*
 * Function: dst_forward
 * Sine transform of the m x n interior array s_f along the vertical direction.
 * Output s_fhat in split mode order.
 */
static void dst_forward(void)
{
    const int n = s_n;

    for (int ii = 0; ii < s_h; ii++)
    {
        const double* f_lo = &s_f[ii * n];
        const double* f_hi = &s_f[(s_m - 1 - ii) * n];
        double* fe = &s_fe[ii * n];
        double* fo = &s_fo[ii * n];
        for (int jj = 0; jj < n; jj++)
        {
            fe[jj] = f_lo[jj] + f_hi[jj];
            fo[jj] = f_lo[jj] - f_hi[jj];
        }
    }
    if (s_hp > s_h)
    {
        memcpy(&s_fe[s_h * n], &s_f[s_h * n], (size_t)n * sizeof(double));
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, s_m_odd, n, s_hp,
                1.0, s_s_odd, s_hp, s_fe, n, 0.0, s_fhat, n);
    if (s_m_even > 0)
    {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, s_m_even, n, s_h,
                    1.0, s_s_even, s_h, s_fo, n, 0.0, &s_fhat[s_m_odd * n], n);
    }
}

/*
 * Function: dst_forward_vector
 * Sine transform of one vertical vector v (length m) into vhat (length m,
 * split mode order).
 */
static void dst_forward_vector(const double* v, double* vhat)
{
    for (int ii = 0; ii < s_h; ii++)
    {
        s_vec_e[ii] = v[ii] + v[s_m - 1 - ii];
        s_vec_o[ii] = v[ii] - v[s_m - 1 - ii];
    }
    if (s_hp > s_h)
    {
        s_vec_e[s_h] = v[s_h];
    }
    cblas_dgemv(CblasRowMajor, CblasNoTrans, s_m_odd, s_hp, 1.0, s_s_odd, s_hp, s_vec_e, 1, 0.0, vhat, 1);
    if (s_m_even > 0)
    {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, s_m_even, s_h, 1.0, s_s_even, s_h, s_vec_o, 1, 0.0,
                    &vhat[s_m_odd], 1);
    }
}

/*
 * Function: dst_inverse
 * Inverse sine transform of s_fhat (split mode order), written into the
 * interior of the N_Z x N_R output array.
 */
static void dst_inverse(double* out)
{
    const int n = s_n;
    const double scale = 2.0 / (double)s_N;

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, s_hp, n, s_m_odd,
                1.0, s_s_odd, s_hp, s_fhat, n, 0.0, s_ge, n);
    if (s_m_even > 0)
    {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, s_h, n, s_m_even,
                    1.0, s_s_even, s_h, &s_fhat[s_m_odd * n], n, 0.0, s_go, n);
    }

    for (int ii = 0; ii < s_h; ii++)
    {
        const double* ge = &s_ge[ii * n];
        const double* go = &s_go[ii * n];
        double* x_lo = &out[(ii + 1) * N_R + 1];
        double* x_hi = &out[(s_m - ii) * N_R + 1];
        for (int jj = 0; jj < n; jj++)
        {
            x_lo[jj] = scale * (ge[jj] + go[jj]);
            x_hi[jj] = scale * (ge[jj] - go[jj]);
        }
    }
    if (s_hp > s_h)
    {
        const double* ge = &s_ge[s_h * n];
        double* x_mid = &out[(s_h + 1) * N_R + 1];
        for (int jj = 0; jj < n; jj++)
        {
            x_mid[jj] = scale * ge[jj];
        }
    }
}

/*
 * Function: evaluate_edges
 * Evaluates the inverse transform of s_fhat only on the interior rows and
 * columns next to the grid boundary (the ones gradient_bound() reads), writing
 * them into psi.
 */
static void evaluate_edges(double* psi)
{
    const int n = s_n, m = s_m;
    const int nc = s_n_edge_cols;
    const double scale = 2.0 / (double)s_N;

    /* Rows next to the top and bottom boundary. */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, s_n_edge_rows, n, m,
                1.0, s_s_edge_rows, m, s_fhat, n, 0.0, s_row_out, n);
    for (int q = 0; q < s_n_edge_rows; q++)
    {
        memcpy(&psi[(s_edge_rows[q] + 1) * N_R + 1], &s_row_out[q * n], (size_t)n * sizeof(double));
    }

    /* Columns next to the left and right boundary. */
    for (int kk = 0; kk < m; kk++)
    {
        for (int q = 0; q < nc; q++)
        {
            s_col_in[kk * nc + q] = s_fhat[kk * n + s_edge_cols[q]];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, s_hp, nc, s_m_odd,
                1.0, s_s_odd, s_hp, s_col_in, nc, 0.0, s_col_ge, nc);
    if (s_m_even > 0)
    {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, s_h, nc, s_m_even,
                    1.0, s_s_even, s_h, &s_col_in[s_m_odd * nc], nc, 0.0, s_col_go, nc);
    }
    for (int ii = 0; ii < s_h; ii++)
    {
        for (int q = 0; q < nc; q++)
        {
            const int j = s_edge_cols[q] + 1;
            psi[(ii + 1) * N_R + j] = scale * (s_col_ge[ii * nc + q] + s_col_go[ii * nc + q]);
            psi[(m - ii) * N_R + j] = scale * (s_col_ge[ii * nc + q] - s_col_go[ii * nc + q]);
        }
    }
    if (s_hp > s_h)
    {
        for (int q = 0; q < nc; q++)
        {
            psi[(s_h + 1) * N_R + s_edge_cols[q] + 1] = scale * s_col_ge[s_h * nc + q];
        }
    }
}

/*
 * Function: solve_modes
 * Solves the tridiagonal system of every vertical mode in place in fhat
 * (m x n, split mode order). Four modes are advanced together so their
 * recurrences overlap.
 */
static void solve_modes(double* fhat)
{
    const int n = s_n;

    for (int kk0 = 0; kk0 < s_m; kk0 += 4)
    {
        const int nb = (s_m - kk0 < 4) ? (s_m - kk0) : 4;
        double* f[4];
        const double* mult[4];
        const double* inv_w[4];

        for (int q = 0; q < 4; q++)
        {
            const int kk = (q < nb) ? (kk0 + q) : kk0;
            f[q] = &fhat[kk * n];
            mult[q] = &s_mult[kk * n];
            inv_w[q] = &s_inv_w[kk * n];
        }

        if (nb == 4)
        {
            for (int jj = 1; jj < n; jj++)
            {
                f[0][jj] -= mult[0][jj] * f[0][jj - 1];
                f[1][jj] -= mult[1][jj] * f[1][jj - 1];
                f[2][jj] -= mult[2][jj] * f[2][jj - 1];
                f[3][jj] -= mult[3][jj] * f[3][jj - 1];
            }
            for (int q = 0; q < 4; q++)
            {
                f[q][n - 1] *= inv_w[q][n - 1];
            }
            for (int jj = n - 2; jj >= 0; jj--)
            {
                const double a = s_a[jj];
                f[0][jj] = (f[0][jj] - a * f[0][jj + 1]) * inv_w[0][jj];
                f[1][jj] = (f[1][jj] - a * f[1][jj + 1]) * inv_w[1][jj];
                f[2][jj] = (f[2][jj] - a * f[2][jj + 1]) * inv_w[2][jj];
                f[3][jj] = (f[3][jj] - a * f[3][jj + 1]) * inv_w[3][jj];
            }
        }
        else
        {
            for (int q = 0; q < nb; q++)
            {
                for (int jj = 1; jj < n; jj++)
                {
                    f[q][jj] -= mult[q][jj] * f[q][jj - 1];
                }
                f[q][n - 1] *= inv_w[q][n - 1];
                for (int jj = n - 2; jj >= 0; jj--)
                {
                    f[q][jj] = (f[q][jj] - s_a[jj] * f[q][jj + 1]) * inv_w[q][jj];
                }
            }
        }
    }
}

/*
 * Function: set_boundary_values
 * out[p] = b_vec[p] on the grid boundary (identity rows).
 */
static void set_boundary_values(const double* b_vec, double* out)
{
    for (int j = 0; j < N_R; j++)
    {
        out[j] = b_vec[j];
        out[(N_Z - 1) * N_R + j] = b_vec[(N_Z - 1) * N_R + j];
    }
    for (int i = 1; i < N_Z - 1; i++)
    {
        out[i * N_R] = b_vec[i * N_R];
        out[i * N_R + N_R - 1] = b_vec[i * N_R + N_R - 1];
    }
}

/*
 * Function: build_interior_rhs
 * s_f = interior of b_vec with the boundary contributions (boundary values
 * taken from x) moved to the right-hand side.
 */
static void build_interior_rhs(const double* b_vec, const double* x)
{
    const int m = s_m, n = s_n;

    for (int ii = 0; ii < m; ii++)
    {
        memcpy(&s_f[ii * n], &b_vec[(ii + 1) * N_R + 1], (size_t)n * sizeof(double));
    }
    for (int jj = 0; jj < n; jj++)
    {
        s_f[jj] -= s_e * x[jj + 1];
        s_f[(m - 1) * n + jj] -= s_e * x[(N_Z - 1) * N_R + jj + 1];
    }
    for (int ii = 0; ii < m; ii++)
    {
        s_f[ii * n] -= s_b[0] * x[(ii + 1) * N_R];
        s_f[ii * n + n - 1] -= s_a[n - 1] * x[(ii + 1) * N_R + N_R - 1];
    }
}

/*
 * Function: poisson_fast_solve
 * Solves A x = b_vec: the boundary entries of b_vec are Dirichlet data and the
 * interior entries are the right-hand side. All N_GRID entries of out are
 * written.
 */
void poisson_fast_solve(const double* b_vec, double* out)
{
    set_boundary_values(b_vec, out);
    build_interior_rhs(b_vec, out);
    dst_forward();
    solve_modes(s_fhat);
    dst_inverse(out);
}

/*
 * Function: poisson_fast_poisson_solver
 * Fast version of poisson_solver(): Hagenow boundary calculation followed by
 * the solve with the computed boundary flux. Same inputs, outputs and side
 * effects as poisson_solver() (b_vec receives the boundary flux values).
 */
void poisson_fast_poisson_solver(double* b_vec, double* out)
{
    const int m = s_m, n = s_n;
    double dpsi_ltrb[N_LTRB];
    double psi_ltrb[N_LTRB];

    /* First solve (zero interior boundary flux), evaluated only where needed. */
    set_boundary_values(b_vec, s_psi1);
    build_interior_rhs(b_vec, s_psi1);
    dst_forward();
    solve_modes(s_fhat);
    evaluate_edges(s_psi1);

    /* Boundary flux from the normal derivative (as hagenow_bound). */
    gradient_bound(s_psi1, dpsi_ltrb);
    for (int ii = 0; ii < N_LTRB; ++ii)
    {
        dpsi_ltrb[ii] = dpsi_ltrb[ii] * INV_R_LTRB_MU0[ii];
    }
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N_LTRB, N_LTRB, 1.0, G_LTRB,
                N_LTRB, dpsi_ltrb, 1, 0.0, psi_ltrb, 1);
    add_bound(psi_ltrb, b_vec);
    set_boundary_values(b_vec, s_bnd);

    /* Second solve differs from the first only through the boundary values:
     * its transformed right-hand side is the transform of the boundary
     * correction, which is nonzero on the outermost interior ring only. */
    for (int j = 0; j < N_R; j++)
    {
        s_bnd[j] -= s_psi1[j];
        s_bnd[(N_Z - 1) * N_R + j] -= s_psi1[(N_Z - 1) * N_R + j];
    }
    for (int i = 1; i < N_Z - 1; i++)
    {
        s_bnd[i * N_R] -= s_psi1[i * N_R];
        s_bnd[i * N_R + N_R - 1] -= s_psi1[i * N_R + N_R - 1];
    }
    for (int kk = 0; kk < m; kk++)
    {
        /* Bottom row i = 1 has S[k][1] = sin(pi k / N); top row i = m has (-1)^(k+1) sin(pi k / N). */
        const double sign = (kk < s_m_odd) ? 1.0 : -1.0;
        const double c = -s_e * s_sin_k[kk];
        double* f2 = &s_fhat2[kk * n];
        for (int jj = 0; jj < n; jj++)
        {
            f2[jj] = c * (s_bnd[jj + 1] + sign * s_bnd[(N_Z - 1) * N_R + jj + 1]);
        }
    }
    for (int ii = 0; ii < m; ii++)
    {
        s_col_in[ii] = s_bnd[(ii + 1) * N_R];
    }
    dst_forward_vector(s_col_in, s_vec_hat);
    for (int kk = 0; kk < m; kk++)
    {
        s_fhat2[kk * n] -= s_b[0] * s_vec_hat[kk];
    }
    for (int ii = 0; ii < m; ii++)
    {
        s_col_in[ii] = s_bnd[(ii + 1) * N_R + N_R - 1];
    }
    dst_forward_vector(s_col_in, s_vec_hat);
    for (int kk = 0; kk < m; kk++)
    {
        s_fhat2[kk * n + n - 1] -= s_a[n - 1] * s_vec_hat[kk];
    }
    solve_modes(s_fhat2);
    for (int q = 0; q < m * n; q++)
    {
        s_fhat[q] += s_fhat2[q];
    }

    set_boundary_values(b_vec, out);
    dst_inverse(out);
}

/*
 * Function: relative_deviation
 * max |x - y| / max |x| over N_GRID entries.
 */
static double relative_deviation(const double* x, const double* y)
{
    double max_x = 0.0, max_diff = 0.0;
    for (int p = 0; p < N_GRID; p++)
    {
        if (fabs(x[p]) > max_x) max_x = fabs(x[p]);
        if (fabs(x[p] - y[p]) > max_diff) max_diff = fabs(x[p] - y[p]);
    }
    return (max_x > 0.0) ? max_diff / max_x : max_diff;
}

/*
 * Function: stencil_residual
 * max |A x - b| over the grid relative to max |b| (identity rows on the
 * boundary, the stencil in the interior).
 */
static double stencil_residual(const double* x, const double* b)
{
    double max_r = 0.0, max_b = 0.0;
    for (int p = 0; p < N_GRID; p++)
    {
        if (fabs(b[p]) > max_b) max_b = fabs(b[p]);
    }
    for (int i = 0; i < N_Z; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            const int p = i * N_R + j;
            double r;
            if (is_boundary(i, j))
            {
                r = x[p] - b[p];
            }
            else
            {
                r = s_e * (x[p - N_R] + x[p + N_R]) + s_b[j - 1] * x[p - 1] + s_a[j - 1] * x[p + 1]
                  + s_d[j - 1] * x[p] - b[p];
            }
            if (fabs(r) > max_r) max_r = fabs(r);
        }
    }
    return (max_b > 0.0) ? max_r / max_b : max_r;
}

/*
 * Function: validate
 * Checks the solver against the stencil residual on two right-hand sides
 * (zero and nonzero boundary data), and the shortcut Hagenow flow against a
 * plain two-solve Hagenow flow. Records the largest relative deviation.
 */
static int validate(void)
{
    double* b = malloc((size_t)N_GRID * sizeof(double));
    double* b2 = malloc((size_t)N_GRID * sizeof(double));
    double* x_ref = malloc((size_t)N_GRID * sizeof(double));
    double* x_fast = malloc((size_t)N_GRID * sizeof(double));
    double dpsi_ltrb[N_LTRB];
    double psi_ltrb[N_LTRB];
    int status = POISSON_FAST_OK;

    if (b == NULL || b2 == NULL || x_ref == NULL || x_fast == NULL)
    {
        free(b); free(b2); free(x_ref); free(x_fast);
        return POISSON_FAST_ERR_ALLOC;
    }

    s_max_dev = 0.0;
    for (int trial = 0; trial < 3; trial++)
    {
        for (int i = 0; i < N_Z; i++)
        {
            for (int j = 0; j < N_R; j++)
            {
                const int p = i * N_R + j;
                if (trial == 1)
                {
                    b[p] = cos(0.37 * (double)p) + 0.1 * sin(1.3 * (double)i) * cos(0.7 * (double)j);
                }
                else
                {
                    b[p] = is_boundary(i, j) ? 0.0
                         : sin(3.1 * (double)i / (double)N_Z) * cos(2.3 * (double)j / (double)N_R) + 0.5;
                }
            }
        }
        double dev;
        if (trial < 2)
        {
            poisson_fast_solve(b, x_fast);
            dev = stencil_residual(x_fast, b);
        }
        else
        {
            /* Plain Hagenow flow: two full solves. */
            memcpy(b2, b, (size_t)N_GRID * sizeof(double));
            poisson_fast_solve(b2, x_ref);
            gradient_bound(x_ref, dpsi_ltrb);
            for (int ii = 0; ii < N_LTRB; ++ii)
            {
                dpsi_ltrb[ii] = dpsi_ltrb[ii] * INV_R_LTRB_MU0[ii];
            }
            cblas_dgemv(CblasRowMajor, CblasNoTrans, N_LTRB, N_LTRB, 1.0, G_LTRB,
                        N_LTRB, dpsi_ltrb, 1, 0.0, psi_ltrb, 1);
            add_bound(psi_ltrb, b2);
            poisson_fast_solve(b2, x_ref);
            poisson_fast_poisson_solver(b, x_fast);
            dev = relative_deviation(x_ref, x_fast);
            const double dev_b = relative_deviation(b2, b);
            if (dev_b > dev) dev = dev_b;
        }
        if (dev > s_max_dev) s_max_dev = dev;
    }
    if (!(s_max_dev <= VALIDATION_RTOL)) status = POISSON_FAST_ERR_VALIDATION;

    free(b); free(b2); free(x_ref); free(x_fast);
    return status;
}

/*
 * Function: poisson_fast_init
 * Loads the stencil from constants.c, builds the transform tables and
 * validates the solver. Returns POISSON_FAST_OK when the solver can be used,
 * otherwise a nonzero status code. Safe to call more than once; the real-time
 * path never calls malloc.
 */
int poisson_fast_init(void)
{
    poisson_fast_free();

    if (N_Z < 4 || N_R < 3)
    {
        s_status = POISSON_FAST_ERR_GRID_TOO_SMALL;
        return s_status;
    }

    s_m = N_Z - 2;
    s_n = N_R - 2;
    s_N = N_Z - 1;
    s_h = s_m / 2;
    s_hp = s_h + (s_m % 2);
    s_m_odd = (s_m + 1) / 2;
    s_m_even = s_m / 2;
    s_n_edge_rows = edge_indices(s_m, s_edge_rows);
    s_n_edge_cols = edge_indices(s_n, s_edge_cols);

    const size_t h1 = (size_t)(s_h > 0 ? s_h : 1);
    const size_t me1 = (size_t)(s_m_even > 0 ? s_m_even : 1);
    s_a = calloc((size_t)s_n, sizeof(double));
    s_b = calloc((size_t)s_n, sizeof(double));
    s_d = calloc((size_t)s_n, sizeof(double));
    s_s_odd = calloc((size_t)s_m_odd * (size_t)s_hp, sizeof(double));
    s_s_even = calloc(me1 * h1, sizeof(double));
    s_s_edge_rows = calloc((size_t)s_n_edge_rows * (size_t)s_m, sizeof(double));
    s_sin_k = calloc((size_t)s_m, sizeof(double));
    s_mult = calloc((size_t)s_m * (size_t)s_n, sizeof(double));
    s_inv_w = calloc((size_t)s_m * (size_t)s_n, sizeof(double));
    s_f = calloc((size_t)s_m * (size_t)s_n, sizeof(double));
    s_fe = calloc((size_t)s_hp * (size_t)s_n, sizeof(double));
    s_fo = calloc(h1 * (size_t)s_n, sizeof(double));
    s_fhat = calloc((size_t)s_m * (size_t)s_n, sizeof(double));
    s_fhat2 = calloc((size_t)s_m * (size_t)s_n, sizeof(double));
    s_ge = calloc((size_t)s_hp * (size_t)s_n, sizeof(double));
    s_go = calloc(h1 * (size_t)s_n, sizeof(double));
    s_psi1 = calloc((size_t)N_GRID, sizeof(double));
    s_bnd = calloc((size_t)N_GRID, sizeof(double));
    s_col_in = calloc((size_t)s_m * (size_t)(2 * N_EDGE), sizeof(double));
    s_col_ge = calloc((size_t)s_hp * (size_t)(2 * N_EDGE), sizeof(double));
    s_col_go = calloc(h1 * (size_t)(2 * N_EDGE), sizeof(double));
    s_row_out = calloc((size_t)(2 * N_EDGE) * (size_t)s_n, sizeof(double));
    s_vec_e = calloc((size_t)s_hp, sizeof(double));
    s_vec_o = calloc(h1, sizeof(double));
    s_vec_hat = calloc((size_t)s_m, sizeof(double));
    if (s_a == NULL || s_b == NULL || s_d == NULL || s_s_odd == NULL ||
        s_s_even == NULL || s_s_edge_rows == NULL || s_sin_k == NULL || s_mult == NULL ||
        s_inv_w == NULL || s_f == NULL || s_fe == NULL || s_fo == NULL || s_fhat == NULL ||
        s_fhat2 == NULL || s_ge == NULL || s_go == NULL || s_psi1 == NULL || s_bnd == NULL ||
        s_col_in == NULL || s_col_ge == NULL || s_col_go == NULL || s_row_out == NULL ||
        s_vec_e == NULL || s_vec_o == NULL || s_vec_hat == NULL)
    {
        poisson_fast_free();
        s_status = POISSON_FAST_ERR_ALLOC;
        return s_status;
    }

    s_status = load_operator();
    if (s_status == POISSON_FAST_OK) s_status = build_tables();
    if (s_status == POISSON_FAST_OK) s_status = validate();

    if (s_status != POISSON_FAST_OK)
    {
        const int keep_status = s_status;
        const double keep_dev = s_max_dev;
        poisson_fast_free();
        s_status = keep_status;
        s_max_dev = keep_dev;
        fprintf(stderr, "poisson_fast_init: Poisson solver could not be initialised (status %d)\n",
                s_status);
    }
    return s_status;
}

int poisson_fast_status(void)
{
    return s_status;
}

double poisson_fast_max_deviation(void)
{
    return s_max_dev;
}
