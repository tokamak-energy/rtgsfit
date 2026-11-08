/* Array sizes. They are enum constants (not const int variables) so that they
   are integer constant expressions and can be used as fixed array dimensions.
   This file is generated from constants_sizes.h.in at the same time as
   constants.c. */
#ifndef CONSTANTS_SIZES_H_
#define CONSTANTS_SIZES_H_

/* N_ROW : number of rows */
enum { N_Z = 65 };

/* N_COL : number of rows */
enum { N_R = 33 };

/* N_GRID : number of grid points (N_COL * N_ROW) */
enum { N_GRID = 2145 };

/* N_MEAS : number of measurements */
enum { N_MEAS = 72 };

/* N_FLUX_LOOPS : number of flux loop measurements */
enum { N_FLUX_LOOPS = 26 };

/* N_BP_PROBES : number of bp-probe measurements */
enum { N_BP_PROBES = 21 };

/* N_ROGOWSKI_COILS : number of Rogowski coil measurements */
enum { N_ROGOWSKI_COILS = 9 };

/* N_COIL : number of PF coils */
enum { N_COIL = 8 };

/* N_LTRB : number of boundary points */
enum { N_LTRB = 192 };

/* N_LCFS_MAX : max storage allocated for find LCFS coordinates */ 
enum { N_LCFS_MAX = 500 };

/* N_XPT_MAX : max storage allocated for find xpts/opts */ 
enum { N_XPT_MAX = 50 };

/* N_COEF : number of basis coefficients */ 
enum { N_COEF = 25 };

/* N_PLS : number of basis coefficients */ 
enum { N_PLS = 3 };

/* N_VESS : number of basis coefficients */ 
enum { N_VESS = 22 };

/* N_SENS_PCS : number of sensors read from PCS (i.e. the "good" sensors) */ 
enum { N_SENS_PCS = 56 };

/* N_LIMIT : number of limit points */
enum { N_LIMIT = 110 };

/* N_INTRP : number of interpolation points */
enum { N_INTRP = 4 };

#endif
