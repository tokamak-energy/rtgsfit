
/* N_ROW : number of rows */
const int N_Z = 0;

/* N_COL : number of rows */
const int N_R = 0;

/* N_GRID : number of grid points (N_COL * N_ROW) */
const int N_GRID = 0;

/* N_MEAS : number of measurements */
const int N_MEAS = 0;

/* N_COIL : number of PF coils */
const int N_COIL = 0;

/* N_BOUND : number of boundary points + 4 */
const int N_LTRB = 0;

/* N_LCFS_MAX : max storage allocated for find LCFS coordinates */ 
const int N_LCFS_MAX = 0;

/* N_XPT_MAX : max storage allocated for find xpts/opts */ 
const int N_XPT_MAX = 0;

const int N_COEF = 0;

/* FRAC : relative weight between flux_lcfs and flux_axis */
const double FRAC = 0.0;

/* THRESH : numerical precision threshold */
const double THRESH = 0.0;

/* D_ROW : difference in Z_VEC */
const double DZ = 0.0;

/* D_COL : difference in R_VEC */
const double DR = 0.0;

/* R_VEC : R values (length N_COL) */
const double R_VEC[] = {}; 

/* Z_VEC : Z values (length N_ROW) */
const double Z_VEC[] = {};

/* WEIGHT: weight values (length N_MEAS) */
const double WEIGHT[] = {};

/* N_LIMIT : number of limit points */
const int N_LIMIT = 0;

/* N_INTRP : number of interpolation points */
const int N_INTRP = 0;

/* LIMIT_IDX : limiter index (N_LIMIT, N_INTRP) */
const int LIMIT_IDX[] = {};

/* LIMIT_WEIGHT : limiter interpolation weights (N_LIMIT, N_INTRP) */
const double LIMIT_WEIGHT[] = {};

/* MASK_LIM: mask of the points inside the vessel */
const int MASK_LIM[] = {};

/* R_GRID : R values (length N_GRID) */
const double R_GRID[] = {};

/* Z_GRID : Z values (length N_GRID) */
const double Z_GRID[] = {};

/* INV_R_MU0 : 1/(mu0 * R_GRID) */
const double INV_R_MU0[] = {};

/* INV_R_MU0 : 1/(mu0 * R_GRID) for boundary (LTRB)*/
const double INV_R_LTRB_MU0[] = {};

/* R_MU0_DZ2 : mu0 * R_GRID * D_ROW * D_ROW */
const double R_MU0_DZ2[] = {};

/* G_GRID_MEAS : greens matrix (N_GRID, N_MEAS) */
const double G_GRID_MEAS[] = {};

/* G_MEAS_COIL : greems matrix (N_MEAS, N_COIL) */
const double G_MEAS_COIL[] = {};

/* G_BOUND : greens matrix (N_BOUND, N_BOUND) */
const double G_LTRB[] = {};

/* G_GRID_COIL : greens matrix (N_GRID, N_COIL) */
const double G_GRID_COIL[] = {};

/* LOWER_BAND : lower triangular matrix (N_GRID, N_COL) CHECK!*/
const double LOWER_BAND[] = {};

/* UPPER_BAND : upper traingular matrix (N_GRID, N_COL+2) CHECK!*/
const double UPPER_BAND[] = {};

/* PERM_IDX : idx that are to be swapped (N_GRID, ) */
const int PERM_IDX[] = {};    



