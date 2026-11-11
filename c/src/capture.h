/*
 * NetCDF capture writer shared by replay and the PCS.
 *
 * The caller buffers variables during the run, then passes their names, dimensions and contiguous
 * data to rtgsfit_capture_write() once the run has finished. Strings are arrays of UTF-8 strings.
 */
#ifndef RTGSFIT_CAPTURE_H
#define RTGSFIT_CAPTURE_H

#include <stddef.h>
#include <stdint.h>

enum rtgsfit_capture_type
{
    RTGSFIT_CAPTURE_FLOAT64 = 1,
    RTGSFIT_CAPTURE_INT32 = 2,
    RTGSFIT_CAPTURE_STRING = 3
};

struct rtgsfit_capture_dimension
{
    const char *name;
    size_t size;
};

struct rtgsfit_capture_variable
{
    const char *group;
    const char *name;
    enum rtgsfit_capture_type type;
    size_t n_dimensions;
    const char *const *dimensions;
    const void *values;
    size_t n_values;
};

/*
 * Write an rtgsfit-capture-2 NetCDF-4 file atomically. The writer sources compiled preshot nodes
 * from librtgsfit and accepts the other 14 preshot nodes in `variables`. It includes both the
 * solver CRC and a CRC over the complete embedded preshot. Returns 0 on success and -1 on failure;
 * call rtgsfit_capture_last_error().
 */
int32_t rtgsfit_capture_write(const char *path,
                              const char *settings_json,
                              const char *preshot_source,
                              const char *code_version,
                              int32_t timing_enabled,
                              const struct rtgsfit_capture_dimension dimensions[],
                              size_t n_dimensions,
                              const struct rtgsfit_capture_variable variables[],
                              size_t n_variables);

const char *rtgsfit_capture_last_error(void);

#endif
