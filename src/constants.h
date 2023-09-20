#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/* N_ROW : number of rows */
extern const int N_Z;

/* N_COL : number of rows */
extern const int N_R;

/* N_GRID : number of grid points (N_COL * N_ROW) */
extern const int N_GRID;

/* N_MEAS : number of measurements */
extern const int N_MEAS;

/* N_COIL : number of PF coils */
extern const int N_COIL;

/* N_BOUND : number of boundary points + 4  CHECK - should this be like this?*/
extern const int N_LTRB;

/* N_LCFS_MAX : max storage allocated for find LCFS coordinates */ 
extern const int N_LCFS_MAX;

/* N_XPT_MAX : max storage allocated for find xpts/opts */ 
extern const int N_XPT_MAX;

/* N_COEF : number of basis coefficients */ 
extern const int N_COEF;

/* FRAC : relative weight between flux_lcfs and flux_axis */
extern const double FRAC;

/* THRESH : numerical precision threshold */
extern const double THRESH;

/* D_ROW : difference in Z_VEC */
extern const double DZ;

/* D_COL : difference in R_VEC */
extern const double DR;

/* R_VEC : R values (length N_COL) */
extern const double R_VEC[]; 

/* Z_VEC : Z values (length N_ROW) */
extern const double Z_VEC[];

/* WEIGHTS: weight values (length N_MEAS) */
extern const double WEIGHT[];

/* N_LIMIT : number of limit points */
extern const int N_LIMIT;

/* N_INTRP : number of interpolation points */
extern const int N_INTRP;

/* LIMIT_IDX : limiter index (N_LIMIT, N_INTRP) */
extern const int LIMIT_IDX[];

/* LIMIT_WEIGHT : limiter interpolation weights (N_LIMIT, N_INTRP) */
extern const double LIMIT_WEIGHT[];

/* MASK_LIM: mask of the points inside the vessel */
extern const int MASK_LIM[];

/* R_GRID : R values (length N_GRID) */
extern const double R_GRID[];

/* Z_GRID : Z values (length N_GRID) */
extern const double Z_GRID[];

/* INV_R_MU0 : 1/(mu0 * R_GRID) */
extern const double INV_R_MU0[];

/* INV_R_MU0 : 1/(mu0 * R_GRID) for boundary (LTRB)*/
extern const double INV_R_LTRB_MU0[];

/* R_MU0_DZ2 : mu0 * R_GRID * D_ROW * D_ROW */
extern const double R_MU0_DZ2[];

/* G_GRID_MEAS : greens matrix (N_GRID, N_MEAS) */
extern const double G_GRID_MEAS[];

/* G_BOUND : greens matrix (N_BOUND, N_BOUND) */
extern const double G_LTRB[];

/* G_MEAS_COIL : greens matrix (N_MEAS, N_COIL) */
extern const double G_MEAS_COIL[];

/* G_GRID_COIL : greens matrix (N_GRID, N_COIL) */
extern const double G_GRID_COIL[];

/* LOWER_BAND : lower triangular matrix (N_GRID, N_COL) CHECK!*/
extern const double LOWER_BAND[];

/* UPPER_BAND : upper traingular matrix (N_GRID, N_COL+2) CHECK!*/
extern const double UPPER_BAND[];

/* PERM_IDX : idx that are to be swapped (N_GRID, ) */
extern const int PERM_IDX[];
    
#endif

