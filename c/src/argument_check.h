/*
 * Argument-size checks used at the top of every function that takes arrays.
 *
 * Every array argument is passed together with its length, and the function checks the length
 * before touching the array. A mismatch is a programming error (the library and the caller were
 * built from different preshots, for example), and is reported as RTGSFIT_ERR_ARGUMENT_SIZE.
 */
#ifndef RTGSFIT_ARGUMENT_CHECK_H
#define RTGSFIT_ARGUMENT_CHECK_H

#include <stdint.h>
#include <stddef.h>

#include "rtgsfit.h"

/* Return RTGSFIT_ERR_ARGUMENT_SIZE from the calling function if `actual != expected`. */
#define RTGSFIT_CHECK_SIZE(actual, expected)                                                       \
    do                                                                                             \
    {                                                                                              \
        if ((int64_t)(actual) < 0 || (int64_t)(actual) != (int64_t)(expected))                     \
        {                                                                                          \
            return RTGSFIT_ERR_ARGUMENT_SIZE;                                                      \
        }                                                                                          \
    } while (0)

/* Return RTGSFIT_ERR_ARGUMENT_SIZE from the calling function if `actual < minimum`. */
#define RTGSFIT_CHECK_MIN_SIZE(actual, minimum)                                                    \
    do                                                                                             \
    {                                                                                              \
        if ((int64_t)(actual) < 0 || (int64_t)(actual) < (int64_t)(minimum))                       \
        {                                                                                          \
            return RTGSFIT_ERR_ARGUMENT_SIZE;                                                      \
        }                                                                                          \
    } while (0)

/* Empty arrays may be NULL; nonempty arrays must have a valid pointer. */
#define RTGSFIT_CHECK_ARRAY(array, length)                                                         \
    do                                                                                             \
    {                                                                                              \
        if ((length) < 0 || ((length) > 0 && (array) == NULL))                                     \
        {                                                                                          \
            return RTGSFIT_ERR_ARGUMENT_SIZE;                                                      \
        }                                                                                          \
    } while (0)

#endif /* RTGSFIT_ARGUMENT_CHECK_H */
