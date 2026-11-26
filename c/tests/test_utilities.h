/*
 * Minimal test helpers for the C unit tests (no external test framework).
 */
#ifndef RTGSFIT_TEST_UTILITIES_H
#define RTGSFIT_TEST_UTILITIES_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "argument_check.h"

static int32_t test_failures = 0;

#define CHECK(condition)                                                                           \
    do                                                                                             \
    {                                                                                              \
        if (!(condition))                                                                          \
        {                                                                                          \
            fprintf(stderr, "FAILED %s:%d: %s\n", __FILE__, __LINE__, #condition);                 \
            test_failures++;                                                                       \
        }                                                                                          \
    } while (0)

#define CHECK_CLOSE(actual, expected, tolerance)                                                   \
    do                                                                                             \
    {                                                                                              \
        const double check_actual_ = (actual);                                                     \
        const double check_expected_ = (expected);                                                 \
        if (!(fabs(check_actual_ - check_expected_) <= (tolerance)))                               \
        {                                                                                          \
            fprintf(stderr, "FAILED %s:%d: %s = %.17g, expected %.17g\n", __FILE__, __LINE__,      \
                    #actual, check_actual_, check_expected_);                                      \
            test_failures++;                                                                       \
        }                                                                                          \
    } while (0)

#define RUN_TEST(test_function)                                                                    \
    do                                                                                             \
    {                                                                                              \
        const int32_t failures_before_ = test_failures;                                            \
        test_function();                                                                           \
        printf("%s %s\n", test_failures == failures_before_ ? "PASS" : "FAIL", #test_function);    \
    } while (0)

/* Read exactly n comma- or newline-separated numbers from `path`. Returns 0 on success. */
static inline int32_t read_csv(const char *path, double values[], int32_t n)
{
    RTGSFIT_CHECK_MIN_SIZE(n, 1);
    RTGSFIT_CHECK_ARRAY(values, n);
    FILE *file = fopen(path, "r");
    if (!file)
    {
        fprintf(stderr, "could not open %s\n", path);
        return -1;
    }
    int32_t count = 0;
    double value = 0.0;
    while (count < n && fscanf(file, " %lf ,", &value) == 1)
    {
        values[count++] = value;
    }
    const int32_t extra = fscanf(file, " %lf", &value);
    fclose(file);
    if (count != n || extra == 1)
    {
        fprintf(stderr, "%s: expected %d values\n", path, (int)n);
        return -1;
    }
    return 0;
}

static inline int32_t read_csv_int(const char *path, int32_t values[], int32_t n)
{
    RTGSFIT_CHECK_MIN_SIZE(n, 1);
    RTGSFIT_CHECK_ARRAY(values, n);
    double *buffer = malloc(sizeof(double) * (size_t)n);
    if (!buffer)
    {
        return -1;
    }
    const int32_t status = read_csv(path, buffer, n);
    for (int32_t i = 0; status == 0 && i < n; i++)
    {
        values[i] = (int32_t)buffer[i];
    }
    free(buffer);
    return status;
}

#endif /* RTGSFIT_TEST_UTILITIES_H */
