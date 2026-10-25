/*
 * rtgsfit_bench.c : benchmark and state-dump driver for rtgsfit().
 *
 * Phase 1 ("converge"): call rtgsfit() n_iter times from the supplied initial
 *   state, exactly like tests/rtgsfit_verify_analytic/replay_rtgsfit.py, and
 *   dump per-iteration scalars plus flux_total so two builds can be compared.
 * Phase 2 ("steady"): from the converged state, call rtgsfit() n_rep times,
 *   re-copying the in/out arrays before every call (same protocol as
 *   tests/rtgsfit_vs_gsfit/tests/test_timing_rtgsfit.py), timing each call
 *   with CLOCK_MONOTONIC_RAW and, when built with -DENABLE_RT_TIMING, reading
 *   the library's per-section CPU timers.
 *
 * Usage: rtgsfit_bench <input_dir> <n_iter> <n_rep> <out_dir>
 * input_dir must contain meas_pcs.txt (N_SENS_PCS values), coil_curr.txt
 * (N_COIL), flux_norm.txt (N_GRID) and mask.txt (N_GRID), e.g. as written by
 * tests/rtgsfit_verify_analytic/src/rtgsfit_verify_analytic/export_bench_inputs.py.
 * Build with tests/makefile_bench after building src/ with the same constants.c.
 * OpenBLAS is pinned to one thread inside the program.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <cblas.h>
#include "constants.h"
#include "rtgsfit.h"

#define T_NTIMERS 18
static const char *section_names[T_NTIMERS] = {
    "meas_prep", "basis", "meas_matrix", "weighting", "ls_fit", "source",
    "model_meas", "chi2", "poisson", "coil_flux", "vessel_flux", "xpts_and_axis",
    "xpts_sort", "limiter", "lcfs", "inside", "normalise", "total"};

static double now_us(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

static void read_doubles(const char *dir, const char *name, int n, double *out)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    for (int i = 0; i < n; i++)
        if (fscanf(f, "%lf", &out[i]) != 1) { fprintf(stderr, "short read %s at %d\n", path, i); exit(2); }
    fclose(f);
}

static void read_ints(const char *dir, const char *name, int n, int32_t *out)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    for (int i = 0; i < n; i++)
        if (fscanf(f, "%d", &out[i]) != 1) { fprintf(stderr, "short read %s at %d\n", path, i); exit(2); }
    fclose(f);
}

typedef struct {
    double *meas_pcs, *coil_curr, *flux_norm, *flux_total, *lcfs_r, *lcfs_z, *coef, *meas_model;
    double *xpt_r, *xpt_z, *xpt_flux;
    int32_t *mask;
    double chi_sq_err, flux_boundary, plasma_current, r_mag_axis, z_mag_axis, mag_axis_flux;
    double r_cur_centroid, z_cur_centroid;
    int32_t lcfs_n, lcfs_err_code, xpt_n, xpt_diverted;
    int lapack_dgelss_info;
} state_t;

static void state_alloc(state_t *s)
{
    s->meas_pcs = calloc(N_SENS_PCS, sizeof(double));
    s->coil_curr = calloc(N_COIL, sizeof(double));
    s->flux_norm = calloc(N_GRID, sizeof(double));
    s->flux_total = calloc(N_GRID, sizeof(double));
    s->lcfs_r = calloc(N_LCFS_MAX, sizeof(double));
    s->lcfs_z = calloc(N_LCFS_MAX, sizeof(double));
    s->coef = calloc(N_COEF, sizeof(double));
    s->meas_model = calloc(N_MEAS, sizeof(double));
    s->xpt_r = calloc(N_XPT_MAX, sizeof(double));
    s->xpt_z = calloc(N_XPT_MAX, sizeof(double));
    s->xpt_flux = calloc(N_XPT_MAX, sizeof(double));
    s->mask = calloc(N_GRID, sizeof(int32_t));
}

static void state_free(state_t *s)
{
    free(s->meas_pcs); free(s->coil_curr); free(s->flux_norm); free(s->flux_total);
    free(s->lcfs_r); free(s->lcfs_z); free(s->coef); free(s->meas_model);
    free(s->xpt_r); free(s->xpt_z); free(s->xpt_flux); free(s->mask);
}

static void call(state_t *s)
{
    rtgsfit(s->meas_pcs, s->coil_curr, s->flux_norm, s->mask, s->flux_total,
            &s->chi_sq_err, s->lcfs_r, s->lcfs_z, &s->lcfs_n, s->coef,
            &s->flux_boundary, &s->plasma_current, &s->lcfs_err_code,
            &s->lapack_dgelss_info, s->meas_model, N_MEAS, &s->r_mag_axis,
            &s->z_mag_axis, &s->mag_axis_flux, &s->r_cur_centroid,
            &s->z_cur_centroid, s->xpt_r, s->xpt_z, s->xpt_flux, N_XPT_MAX,
            &s->xpt_n, &s->xpt_diverted);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void print_stats(FILE *f, const char *label, double *v, int n)
{
    double *s = malloc(n * sizeof(double));
    memcpy(s, v, n * sizeof(double));
    qsort(s, n, sizeof(double), cmp_double);
    double sum = 0, sq = 0;
    for (int i = 0; i < n; i++) { sum += s[i]; sq += s[i] * s[i]; }
    double mean = sum / n, sd = sqrt(fmax(sq / n - mean * mean, 0.0));
    #define PCT(p) s[(int)fmin(n - 1, floor((p) * (n - 1) + 0.5))]
    fprintf(f, "%-14s %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f\n", label,
            mean, PCT(0.5), sd, s[0], PCT(0.95), PCT(0.99), PCT(0.999), s[n - 1]);
    free(s);
}

int main(int argc, char **argv)
{
    if (argc < 5) { fprintf(stderr, "usage: %s workload_dir n_iter n_rep out_dir\n", argv[0]); return 1; }
    const char *wdir = argv[1];
    int n_iter = atoi(argv[2]);
    int n_rep = atoi(argv[3]);
    const char *odir = argv[4];
    char path[4096];

    openblas_set_num_threads(1);
    printf("grid %dx%d (N_GRID %d) N_MEAS %d N_COEF %d N_COIL %d N_VESS %d N_LTRB %d openblas_threads %d\n",
           N_R, N_Z, N_GRID, N_MEAS, N_COEF, N_COIL, N_VESS, N_LTRB, openblas_get_num_threads());

    state_t s;
    state_alloc(&s);
    read_doubles(wdir, "meas_pcs.txt", N_SENS_PCS, s.meas_pcs);
    read_doubles(wdir, "coil_curr.txt", N_COIL, s.coil_curr);
    read_doubles(wdir, "flux_norm.txt", N_GRID, s.flux_norm);
    read_ints(wdir, "mask.txt", N_GRID, s.mask);

    /* ---------------- Phase 1: convergence from the initial state ---------------- */
    snprintf(path, sizeof path, "%s/converge.csv", odir);
    FILE *fc = fopen(path, "w");
    snprintf(path, sizeof path, "%s/converge_flux_total.bin", odir);
    FILE *fb = fopen(path, "wb");
    if (!fc || !fb) { fprintf(stderr, "cannot write to %s\n", odir); return 2; }
    fprintf(fc, "iter,wall_us,err_code,dgelss_info,chi_sq_err,plasma_current,r_mag_axis,z_mag_axis,mag_axis_flux,flux_boundary,xpt_n,xpt_diverted,mask_sum,r_cur_centroid,z_cur_centroid\n");
    for (int it = 0; it < n_iter; it++) {
        double t0 = now_us();
        call(&s);
        double t1 = now_us();
        long mask_sum = 0;
        for (int i = 0; i < N_GRID; i++) mask_sum += s.mask[i];
        fprintf(fc, "%d,%.3f,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%ld,%.17g,%.17g\n", it + 1, t1 - t0,
                s.lcfs_err_code, s.lapack_dgelss_info, s.chi_sq_err, s.plasma_current, s.r_mag_axis,
                s.z_mag_axis, s.mag_axis_flux, s.flux_boundary, s.xpt_n, s.xpt_diverted, mask_sum,
                s.r_cur_centroid, s.z_cur_centroid);
        fwrite(s.flux_total, sizeof(double), N_GRID, fb);
    }
    fclose(fc); fclose(fb);
    /* final state dump: flux_norm, mask, coef, meas_model, xpt arrays */
    snprintf(path, sizeof path, "%s/final_state.bin", odir);
    FILE *ff = fopen(path, "wb");
    fwrite(s.flux_norm, sizeof(double), N_GRID, ff);
    fwrite(s.flux_total, sizeof(double), N_GRID, ff);
    fwrite(s.mask, sizeof(int32_t), N_GRID, ff);
    fwrite(s.coef, sizeof(double), N_COEF, ff);
    fwrite(s.meas_model, sizeof(double), N_MEAS, ff);
    fwrite(s.xpt_r, sizeof(double), N_XPT_MAX, ff);
    fwrite(s.xpt_z, sizeof(double), N_XPT_MAX, ff);
    fwrite(s.xpt_flux, sizeof(double), N_XPT_MAX, ff);
    fclose(ff);
    printf("converge: %d iterations, final err_code %d dgelss_info %d chi2 %.6g Ip %.6g axis (%.6g, %.6g) psi_ax %.6g psi_b %.6g xpt_n %d\n",
           n_iter, s.lcfs_err_code, s.lapack_dgelss_info, s.chi_sq_err, s.plasma_current, s.r_mag_axis,
           s.z_mag_axis, s.mag_axis_flux, s.flux_boundary, s.xpt_n);

    if (n_rep <= 0)
    {
        state_free(&s);
        return 0;
    }

    /* ---------------- Phase 2: steady-state repeated calls ---------------- */
    double *flux_norm_conv = malloc(N_GRID * sizeof(double));
    int32_t *mask_conv = malloc(N_GRID * sizeof(int32_t));
    double *meas_fixed = malloc(N_SENS_PCS * sizeof(double));
    double *coil_fixed = malloc(N_COIL * sizeof(double));
    memcpy(flux_norm_conv, s.flux_norm, N_GRID * sizeof(double));
    memcpy(mask_conv, s.mask, N_GRID * sizeof(int32_t));
    memcpy(meas_fixed, s.meas_pcs, N_SENS_PCS * sizeof(double));
    memcpy(coil_fixed, s.coil_curr, N_COIL * sizeof(double));

    double *wall = malloc(n_rep * sizeof(double));
    double *sect = malloc((size_t)n_rep * T_NTIMERS * sizeof(double));
    double *ref_flux = malloc(N_GRID * sizeof(double));
    double *ref_coef = malloc(N_COEF * sizeof(double));
    int ref_err = 0;
    long n_mismatch = 0, n_err = 0;

    for (int r = 0; r < n_rep; r++) {
        memcpy(s.flux_norm, flux_norm_conv, N_GRID * sizeof(double));
        memcpy(s.mask, mask_conv, N_GRID * sizeof(int32_t));
        memcpy(s.meas_pcs, meas_fixed, N_SENS_PCS * sizeof(double));
        memcpy(s.coil_curr, coil_fixed, N_COIL * sizeof(double));
#ifdef ENABLE_RT_TIMING
        rtgsfit_timing_reset();
#endif
        double t0 = now_us();
        call(&s);
        double t1 = now_us();
        wall[r] = t1 - t0;
#ifdef ENABLE_RT_TIMING
        rtgsfit_timing_get(&sect[(size_t)r * T_NTIMERS], T_NTIMERS);
#endif
        if (r == 0) {
            memcpy(ref_flux, s.flux_total, N_GRID * sizeof(double));
            memcpy(ref_coef, s.coef, N_COEF * sizeof(double));
            ref_err = s.lcfs_err_code;
        } else {
            if (memcmp(ref_flux, s.flux_total, N_GRID * sizeof(double)) != 0 ||
                memcmp(ref_coef, s.coef, N_COEF * sizeof(double)) != 0) n_mismatch++;
        }
        if (s.lcfs_err_code != ref_err) n_err++;
    }

    snprintf(path, sizeof path, "%s/steady.csv", odir);
    FILE *fs = fopen(path, "w");
    fprintf(fs, "rep,wall_us");
    for (int i = 0; i < T_NTIMERS; i++) fprintf(fs, ",%s_us", section_names[i]);
    fprintf(fs, "\n");
    for (int r = 0; r < n_rep; r++) {
        fprintf(fs, "%d,%.3f", r, wall[r]);
        for (int i = 0; i < T_NTIMERS; i++) fprintf(fs, ",%.3f", sect[(size_t)r * T_NTIMERS + i]);
        fprintf(fs, "\n");
    }
    fclose(fs);

    printf("steady: %d calls, err_code %d, bitwise-identical outputs across calls: %s (%ld mismatches), err-code changes %ld\n",
           n_rep, ref_err, n_mismatch == 0 ? "yes" : "NO", n_mismatch, n_err);
    printf("%-14s %9s %9s %9s %9s %9s %9s %9s %9s  (us)\n", "section", "mean", "median", "std", "min", "p95", "p99", "p99.9", "max");
    print_stats(stdout, "wall_total", wall, n_rep);
#ifdef ENABLE_RT_TIMING
    double *col = malloc(n_rep * sizeof(double));
    for (int i = 0; i < T_NTIMERS; i++) {
        for (int r = 0; r < n_rep; r++) col[r] = sect[(size_t)r * T_NTIMERS + i];
        print_stats(stdout, section_names[i], col, n_rep);
    }
    free(col);
#endif
    free(flux_norm_conv); free(mask_conv); free(meas_fixed); free(coil_fixed);
    free(wall); free(sect); free(ref_flux); free(ref_coef);
    state_free(&s);
    return 0;
}
