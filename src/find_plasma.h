// find_plasma.h
#ifndef FIND_PLASMA_H
#define FIND_PLASMA_H

#include <stdint.h>

void find_nulls(
    double* flux,
    double* opt_r,
    double* opt_z,
    double* opt_flux,
    int32_t* opt_n,
    double* xpt_r,
    double* xpt_z,
    double* xpt_flux,
    int32_t* xpt_n
);

#endif // FIND_PLASMA_H
