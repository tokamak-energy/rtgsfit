#define _POSIX_C_SOURCE 200809L

#include "capture.h"

#include <errno.h>
#include <fcntl.h>
#include <netcdf.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "preshot.h"
#include "rtgsfit.h"

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "capture checksums use little-endian numeric bytes"
#endif

#define CAPTURE_SCHEMA     "rtgsfit-capture-2"
#define CAPTURE_ERROR_SIZE 512

static char capture_error[CAPTURE_ERROR_SIZE];
static atomic_uint_fast64_t temporary_counter = 0;

static const char *const size_names[RTGSFIT_SIZE_COUNT] = {
    "N_R",          "N_Z",         "N_GRID",     "N_EDGE",     "N_LIMITER",    "N_COIL",
    "N_MEAS",       "N_FLUX_LOOP", "N_BP_PROBE", "N_ROGOWSKI", "N_CONSTRAINT", "N_REG",
    "N_FIT_TARGET", "N_PLASMA",    "N_PASSIVE",  "N_FIT_DOF",  "N_XPT_MAX",
};

static const char *const dim_n_r[] = {"N_R"};
static const char *const dim_n_z[] = {"N_Z"};
static const char *const dim_n_grid[] = {"N_GRID"};
static const char *const dim_n_edge[] = {"N_EDGE"};
static const char *const dim_n_limiter[] = {"N_LIMITER"};
static const char *const dim_n_fit_target[] = {"N_FIT_TARGET"};
static const char *const dim_grid_coil[] = {"N_GRID", "N_COIL"};
static const char *const dim_grid_passive[] = {"N_GRID", "N_PASSIVE"};
static const char *const dim_edge_edge[] = {"N_EDGE", "N_EDGE"};
static const char *const dim_limiter_interp[] = {"N_LIMITER", "N_INTERP"};
static const char *const dim_constraint_coil[] = {"N_CONSTRAINT", "N_COIL"};
static const char *const dim_constraint_meas[] = {"N_CONSTRAINT", "N_MEAS"};
static const char *const dim_fit_target_grid[] = {"N_FIT_TARGET", "N_GRID"};
static const char *const dim_fit_target_passive[] = {"N_FIT_TARGET", "N_PASSIVE"};
static const char *const dim_grid_r[] = {"N_GRID", "N_R"};
static const char *const dim_grid_upper[] = {"N_GRID", "N_R_PLUS_2"};

#define CAPTURE_NODE(group_, name_, type_, ndims_, dims_, values_)                                 \
    {group_, name_, type_, ndims_, dims_, values_, 0}
#define CAPTURE_SCALAR(group_, name_, type_, values_)                                              \
    CAPTURE_NODE(group_, name_, type_, 0, NULL, values_)

static const struct rtgsfit_capture_variable compiled_nodes[] = {
    CAPTURE_NODE("PRESHOT/GRID", "R", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_r, GRID_R),
    CAPTURE_NODE("PRESHOT/GRID", "Z", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_z, GRID_Z),
    CAPTURE_SCALAR("PRESHOT/GRID", "D_R", RTGSFIT_CAPTURE_FLOAT64, &GRID_D_R),
    CAPTURE_SCALAR("PRESHOT/GRID", "D_Z", RTGSFIT_CAPTURE_FLOAT64, &GRID_D_Z),
    CAPTURE_NODE("PRESHOT/GRID", "R_FLAT", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_grid, GRID_R_FLAT),
    CAPTURE_NODE("PRESHOT/GRID", "Z_FLAT", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_grid, GRID_Z_FLAT),
    CAPTURE_NODE("PRESHOT/GRID", "INV_R", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_grid, GRID_INV_R),
    CAPTURE_NODE(
        "PRESHOT/GRID", "MASK_LIMITER", RTGSFIT_CAPTURE_INT32, 1, dim_n_grid, GRID_MASK_LIMITER),
    CAPTURE_NODE("PRESHOT/GRID", "G_COIL", RTGSFIT_CAPTURE_FLOAT64, 2, dim_grid_coil, GRID_G_COIL),
    CAPTURE_NODE(
        "PRESHOT/GRID", "G_PASSIVE", RTGSFIT_CAPTURE_FLOAT64, 2, dim_grid_passive, GRID_G_PASSIVE),
    CAPTURE_NODE("PRESHOT/EDGE", "R", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_edge, EDGE_R),
    CAPTURE_NODE("PRESHOT/EDGE", "Z", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_edge, EDGE_Z),
    CAPTURE_NODE(
        "PRESHOT/EDGE", "INV_MU0_R", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_edge, EDGE_INV_MU0_R),
    CAPTURE_NODE("PRESHOT/EDGE", "G_EDGE", RTGSFIT_CAPTURE_FLOAT64, 2, dim_edge_edge, EDGE_G_EDGE),
    CAPTURE_NODE("PRESHOT/LIMITER", "R", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_limiter, LIMITER_R),
    CAPTURE_NODE("PRESHOT/LIMITER", "Z", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_limiter, LIMITER_Z),
    CAPTURE_NODE("PRESHOT/LIMITER",
                 "GRID_INDEX",
                 RTGSFIT_CAPTURE_INT32,
                 2,
                 dim_limiter_interp,
                 LIMITER_GRID_INDEX),
    CAPTURE_NODE("PRESHOT/LIMITER",
                 "GRID_WEIGHT",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_limiter_interp,
                 LIMITER_GRID_WEIGHT),
    CAPTURE_NODE(
        "PRESHOT/POISSON", "RHS_SCALE", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_grid, POISSON_RHS_SCALE),
    CAPTURE_NODE("PRESHOT/POISSON",
                 "LOWER_BAND",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_grid_r,
                 POISSON_LOWER_BAND),
    CAPTURE_NODE("PRESHOT/POISSON",
                 "UPPER_BAND",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_grid_upper,
                 POISSON_UPPER_BAND),
    CAPTURE_NODE("PRESHOT/POISSON",
                 "PERMUTATION",
                 RTGSFIT_CAPTURE_INT32,
                 1,
                 dim_n_grid,
                 POISSON_PERMUTATION),
    CAPTURE_NODE("PRESHOT/CONSTRAINT",
                 "G_COIL",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_constraint_coil,
                 CONSTRAINT_G_COIL),
    CAPTURE_NODE("PRESHOT/CONSTRAINT",
                 "REPLACEMENT",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_constraint_meas,
                 CONSTRAINT_REPLACEMENT),
    CAPTURE_NODE("PRESHOT/FIT", "WEIGHT", RTGSFIT_CAPTURE_FLOAT64, 1, dim_n_fit_target, FIT_WEIGHT),
    CAPTURE_NODE(
        "PRESHOT/FIT", "G_GRID", RTGSFIT_CAPTURE_FLOAT64, 2, dim_fit_target_grid, FIT_G_GRID),
    CAPTURE_NODE("PRESHOT/FIT",
                 "G_PASSIVE",
                 RTGSFIT_CAPTURE_FLOAT64,
                 2,
                 dim_fit_target_passive,
                 FIT_G_PASSIVE),
    CAPTURE_SCALAR(
        "PRESHOT/SETTINGS", "XPT_FRACTION", RTGSFIT_CAPTURE_FLOAT64, &SETTINGS_XPT_FRACTION),
    CAPTURE_SCALAR("PRESHOT/SETTINGS", "TOLERANCE", RTGSFIT_CAPTURE_FLOAT64, &SETTINGS_TOLERANCE),
    CAPTURE_SCALAR("PRESHOT/SETTINGS", "IP_MIN", RTGSFIT_CAPTURE_FLOAT64, &SETTINGS_IP_MIN),
};

static const char *const extra_preshot_nodes[] = {
    "PRESHOT/CONSTRAINT:NAMES", "PRESHOT/FIT:TARGET_NAMES", "PRESHOT/FIT:DOF_NAMES",
    "PRESHOT/FIT:DOF_UNITS",    "PRESHOT/COIL:NAMES",       "PRESHOT/COIL:SIGNALS",
    "PRESHOT/COIL:MATRIX",      "PRESHOT/MEAS:NAMES",       "PRESHOT/INITIAL:PSI_NORM",
    "PRESHOT/INITIAL:MASK",     "PRESHOT/INFO:SETTINGS",    "PRESHOT/INFO:MU0",
    "PRESHOT/INFO:PSI_TO_WB",   "PRESHOT/INFO:SOURCE",
};

static const struct
{
    const char *group;
    const char *name;
} preshot_node_order[] = {
    {"PRESHOT/GRID", "R"},
    {"PRESHOT/GRID", "Z"},
    {"PRESHOT/GRID", "D_R"},
    {"PRESHOT/GRID", "D_Z"},
    {"PRESHOT/GRID", "R_FLAT"},
    {"PRESHOT/GRID", "Z_FLAT"},
    {"PRESHOT/GRID", "INV_R"},
    {"PRESHOT/GRID", "MASK_LIMITER"},
    {"PRESHOT/GRID", "G_COIL"},
    {"PRESHOT/GRID", "G_PASSIVE"},
    {"PRESHOT/EDGE", "R"},
    {"PRESHOT/EDGE", "Z"},
    {"PRESHOT/EDGE", "INV_MU0_R"},
    {"PRESHOT/EDGE", "G_EDGE"},
    {"PRESHOT/LIMITER", "R"},
    {"PRESHOT/LIMITER", "Z"},
    {"PRESHOT/LIMITER", "GRID_INDEX"},
    {"PRESHOT/LIMITER", "GRID_WEIGHT"},
    {"PRESHOT/POISSON", "RHS_SCALE"},
    {"PRESHOT/POISSON", "LOWER_BAND"},
    {"PRESHOT/POISSON", "UPPER_BAND"},
    {"PRESHOT/POISSON", "PERMUTATION"},
    {"PRESHOT/CONSTRAINT", "NAMES"},
    {"PRESHOT/CONSTRAINT", "G_COIL"},
    {"PRESHOT/CONSTRAINT", "REPLACEMENT"},
    {"PRESHOT/FIT", "TARGET_NAMES"},
    {"PRESHOT/FIT", "DOF_NAMES"},
    {"PRESHOT/FIT", "DOF_UNITS"},
    {"PRESHOT/FIT", "WEIGHT"},
    {"PRESHOT/FIT", "G_GRID"},
    {"PRESHOT/FIT", "G_PASSIVE"},
    {"PRESHOT/SETTINGS", "XPT_FRACTION"},
    {"PRESHOT/SETTINGS", "TOLERANCE"},
    {"PRESHOT/SETTINGS", "IP_MIN"},
    {"PRESHOT/COIL", "NAMES"},
    {"PRESHOT/COIL", "SIGNALS"},
    {"PRESHOT/COIL", "MATRIX"},
    {"PRESHOT/MEAS", "NAMES"},
    {"PRESHOT/INITIAL", "PSI_NORM"},
    {"PRESHOT/INITIAL", "MASK"},
    {"PRESHOT/INFO", "SETTINGS"},
    {"PRESHOT/INFO", "MU0"},
    {"PRESHOT/INFO", "PSI_TO_WB"},
    {"PRESHOT/INFO", "SOURCE"},
};

enum
{
    N_COMPILED_NODES = sizeof(compiled_nodes) / sizeof(compiled_nodes[0]),
    N_EXTRA_PRESHOT_NODES = sizeof(extra_preshot_nodes) / sizeof(extra_preshot_nodes[0]),
    N_PRESHOT_NODES = sizeof(preshot_node_order) / sizeof(preshot_node_order[0]),
};
_Static_assert(N_PRESHOT_NODES + RTGSFIT_SIZE_COUNT == 61,
               "capture preshot order must cover the complete data dictionary");

#undef CAPTURE_SCALAR
#undef CAPTURE_NODE

static uint32_t crc32_update(uint32_t crc, const void *data, const size_t size)
{
    const unsigned char *bytes = data;
    crc = ~crc;
    for (size_t i = 0; i < size; i++)
    {
        crc ^= bytes[i];
        for (int bit = 0; bit < 8; bit++)
        {
            crc = (crc & 1u) ? 0xEDB88320u ^ (crc >> 1) : crc >> 1;
        }
    }
    return ~crc;
}

static int find_variable(const struct rtgsfit_capture_variable variables[],
                         const size_t n_variables,
                         const char *group,
                         const char *name)
{
    for (size_t i = 0; i < n_variables; i++)
    {
        if (variables[i].group != NULL && variables[i].name != NULL &&
            strcmp(variables[i].group, group) == 0 && strcmp(variables[i].name, name) == 0)
        {
            return (int)i;
        }
    }
    return -1;
}

static uint32_t checksum_variable(uint32_t checksum,
                                  const struct rtgsfit_capture_variable *variable,
                                  const char *key)
{
    const char *type_name = variable->type == RTGSFIT_CAPTURE_FLOAT64
                                ? "float64"
                                : (variable->type == RTGSFIT_CAPTURE_INT32 ? "int32" : "str");
    checksum = crc32_update(checksum, key, strlen(key) + 1);
    checksum = crc32_update(checksum, type_name, strlen(type_name) + 1);
    if (variable->type == RTGSFIT_CAPTURE_FLOAT64)
    {
        return crc32_update(checksum, variable->values, variable->n_values * sizeof(double));
    }
    if (variable->type == RTGSFIT_CAPTURE_INT32)
    {
        return crc32_update(checksum, variable->values, variable->n_values * sizeof(int32_t));
    }
    const char *const *values = variable->values;
    for (size_t i = 0; i < variable->n_values; i++)
    {
        const uint64_t length = strlen(values[i]);
        checksum = crc32_update(checksum, &length, sizeof(length));
        checksum = crc32_update(checksum, values[i], (size_t)length);
    }
    return checksum;
}

static uint32_t checksum_preshot(const struct rtgsfit_capture_variable variables[],
                                 const size_t n_variables)
{
    uint32_t checksum = 0;
    char key[128];
    for (size_t i = 0; i < RTGSFIT_SIZE_COUNT; i++)
    {
        const int index = find_variable(variables, n_variables, "PRESHOT/SIZES", size_names[i]);
        if (index < 0)
        {
            return 0;
        }
        snprintf(key, sizeof(key), "PRESHOT.SIZES.%s", size_names[i]);
        checksum = checksum_variable(checksum, &variables[index], key);
    }
    for (size_t i = 0; i < N_PRESHOT_NODES; i++)
    {
        const int index = find_variable(variables, n_variables, preshot_node_order[i].group,
                                        preshot_node_order[i].name);
        if (index < 0)
        {
            return 0;
        }
        snprintf(key, sizeof(key), "%s.%s", preshot_node_order[i].group,
                 preshot_node_order[i].name);
        for (char *character = key; *character != '\0'; character++)
        {
            if (*character == '/')
            {
                *character = '.';
            }
        }
        checksum = checksum_variable(checksum, &variables[index], key);
    }
    return checksum;
}

const char *rtgsfit_capture_last_error(void)
{
    return capture_error;
}

static int32_t set_error(const char *message)
{
    snprintf(capture_error, sizeof(capture_error), "%s", message);
    return -1;
}

static int32_t set_netcdf_error(const char *operation, const int status)
{
    snprintf(capture_error, sizeof(capture_error), "%s: %s", operation, nc_strerror(status));
    return -1;
}

static int valid_group(const char *group)
{
    return strcmp(group, "INPUT") == 0 || strcmp(group, "OUTPUT") == 0 ||
           (strncmp(group, "PRESHOT/", 8) == 0 && group[8] != '\0' &&
            strchr(group + 8, '/') == NULL);
}

static int get_variable_group(const int root_groups[],
                              const struct rtgsfit_capture_variable *variable,
                              int *group_id)
{
    if (strcmp(variable->group, "INPUT") == 0)
    {
        *group_id = root_groups[1];
        return NC_NOERR;
    }
    if (strcmp(variable->group, "OUTPUT") == 0)
    {
        *group_id = root_groups[2];
        return NC_NOERR;
    }
    const char *name = variable->group + 8;
    int status = nc_inq_ncid(root_groups[0], name, group_id);
    if (status == NC_ENOGRP)
    {
        status = nc_def_grp(root_groups[0], name, group_id);
    }
    return status;
}

static int find_dimension(const struct rtgsfit_capture_dimension dimensions[],
                          const size_t n_dimensions,
                          const int root,
                          const char *name,
                          int *dimension_id,
                          size_t *size)
{
    for (size_t i = 0; i < n_dimensions; i++)
    {
        if (strcmp(dimensions[i].name, name) == 0)
        {
            *size = dimensions[i].size;
            return nc_inq_dimid(root, name, dimension_id);
        }
    }
    return NC_EBADDIM;
}

static int validate_variable(const struct rtgsfit_capture_variable *variable,
                             const struct rtgsfit_capture_dimension dimensions[],
                             const size_t n_dimensions)
{
    if (variable->group == NULL || variable->group[0] == '\0' || variable->name == NULL ||
        variable->name[0] == '\0' || !valid_group(variable->group) || variable->n_dimensions > 8 ||
        (variable->n_dimensions != 0 && variable->dimensions == NULL))
    {
        return set_error("invalid capture variable descriptor");
    }
    if (variable->type != RTGSFIT_CAPTURE_FLOAT64 && variable->type != RTGSFIT_CAPTURE_INT32 &&
        variable->type != RTGSFIT_CAPTURE_STRING)
    {
        snprintf(capture_error, sizeof(capture_error), "unsupported type for %s/%s",
                 variable->group, variable->name);
        return -1;
    }
    size_t expected_values = 1;
    for (size_t i = 0; i < variable->n_dimensions; i++)
    {
        size_t dimension_size = 0;
        int found = 0;
        if (variable->dimensions[i] == NULL)
        {
            return set_error("capture variable is missing a dimension name");
        }
        for (size_t j = 0; j < n_dimensions; j++)
        {
            if (strcmp(variable->dimensions[i], dimensions[j].name) == 0)
            {
                dimension_size = dimensions[j].size;
                found = 1;
                break;
            }
        }
        if (!found)
        {
            snprintf(capture_error, sizeof(capture_error), "unknown dimension %s for %s/%s",
                     variable->dimensions[i], variable->group, variable->name);
            return -1;
        }
        if (dimension_size != 0 && expected_values > SIZE_MAX / dimension_size)
        {
            snprintf(capture_error, sizeof(capture_error), "element count overflows for %s/%s",
                     variable->group, variable->name);
            return -1;
        }
        expected_values *= dimension_size;
    }
    if (expected_values != variable->n_values ||
        (variable->values == NULL && variable->n_values != 0))
    {
        snprintf(capture_error, sizeof(capture_error),
                 "element count or data pointer is invalid for %s/%s", variable->group,
                 variable->name);
        return -1;
    }
    if (variable->type == RTGSFIT_CAPTURE_STRING)
    {
        const char *const *values = variable->values;
        for (size_t i = 0; i < variable->n_values; i++)
        {
            if (values[i] == NULL)
            {
                snprintf(capture_error, sizeof(capture_error), "null string in %s/%s",
                         variable->group, variable->name);
                return -1;
            }
        }
    }
    return 0;
}

static int make_temporary_path(const char *path, char **temporary)
{
    const size_t length = strlen(path) + 64;
    char *candidate = malloc(length);
    if (candidate == NULL)
    {
        return set_error("out of memory creating temporary capture path");
    }

    for (int attempt = 0; attempt < 100; attempt++)
    {
        const uint64_t counter = atomic_fetch_add(&temporary_counter, 1);
        snprintf(candidate, length, "%s.tmp.%ld.%llu", path, (long)getpid(),
                 (unsigned long long)counter);
        const int descriptor = open(candidate, O_CREAT | O_EXCL | O_WRONLY, 0600);
        if (descriptor >= 0)
        {
            close(descriptor);
            *temporary = candidate;
            return 0;
        }
        if (errno != EEXIST)
        {
            snprintf(capture_error, sizeof(capture_error), "cannot create %s: %s", candidate,
                     strerror(errno));
            free(candidate);
            return -1;
        }
    }
    free(candidate);
    return set_error("could not create a unique temporary capture file");
}

int32_t rtgsfit_capture_write(const char *path,
                              const char *settings_json,
                              const char *preshot_source,
                              const char *code_version,
                              const int32_t timing_enabled,
                              const struct rtgsfit_capture_dimension dimensions[],
                              const size_t n_dimensions,
                              const struct rtgsfit_capture_variable variables[],
                              size_t n_variables)
{
    capture_error[0] = '\0';
    if (path == NULL || path[0] == '\0' || settings_json == NULL || preshot_source == NULL ||
        code_version == NULL || dimensions == NULL || n_dimensions == 0 || variables == NULL ||
        n_variables == 0 || (timing_enabled != 0 && timing_enabled != 1))
    {
        return set_error("invalid capture arguments");
    }
    for (size_t i = 0; i < n_dimensions; i++)
    {
        if (dimensions[i].name == NULL || dimensions[i].name[0] == '\0')
        {
            return set_error("capture dimension has no name");
        }
    }

    int status;
    int root = -1;
    int groups[3] = {-1, -1, -1};
    char *temporary = NULL;
    struct rtgsfit_capture_variable *combined_variables = NULL;
    const char *group_names[3] = {"PRESHOT", "INPUT", "OUTPUT"};
    int32_t compiled_sizes[RTGSFIT_SIZE_COUNT];
    struct rtgsfit_capture_variable size_variables[RTGSFIT_SIZE_COUNT];

    if (rtgsfit_get_sizes(compiled_sizes, RTGSFIT_SIZE_COUNT) != 0)
    {
        return set_error("cannot get the compiled solver dimensions");
    }
    for (size_t size_index = 0; size_index < RTGSFIT_SIZE_COUNT; size_index++)
    {
        int found_size = 0;
        for (size_t j = 0; j < n_dimensions; j++)
        {
            if (strcmp(dimensions[j].name, size_names[size_index]) == 0)
            {
                if (dimensions[j].size != (size_t)compiled_sizes[size_index])
                {
                    snprintf(capture_error, sizeof(capture_error),
                             "dimension %s does not match the compiled solver",
                             size_names[size_index]);
                    return -1;
                }
                found_size = 1;
                break;
            }
        }
        if (!found_size)
        {
            snprintf(capture_error, sizeof(capture_error), "capture is missing dimension %s",
                     size_names[size_index]);
            return -1;
        }
    }

    const struct
    {
        const char *name;
        size_t expected;
    } derived_dimensions[] = {
        {"N_INTERP", PLASMA_BOUNDARY_N_INTERP},
        {"N_R_PLUS_2", (size_t)compiled_sizes[RTGSFIT_SIZE_N_R] + 2},
    };
    for (size_t derived_index = 0;
         derived_index < sizeof(derived_dimensions) / sizeof(derived_dimensions[0]);
         derived_index++)
    {
        int found_derived = 0;
        for (size_t j = 0; j < n_dimensions; j++)
        {
            if (strcmp(dimensions[j].name, derived_dimensions[derived_index].name) == 0)
            {
                if (dimensions[j].size != derived_dimensions[derived_index].expected)
                {
                    snprintf(capture_error, sizeof(capture_error), "dimension %s is inconsistent",
                             derived_dimensions[derived_index].name);
                    return -1;
                }
                found_derived = 1;
                break;
            }
        }
        if (!found_derived)
        {
            snprintf(capture_error, sizeof(capture_error), "capture is missing dimension %s",
                     derived_dimensions[derived_index].name);
            return -1;
        }
    }

    size_t extra_count = 0;
    for (size_t i = 0; i < n_variables; i++)
    {
        if (variables[i].group != NULL && strncmp(variables[i].group, "PRESHOT/", 8) == 0)
        {
            extra_count++;
            for (size_t j = 0; variables[i].name != NULL && j < N_COMPILED_NODES; j++)
            {
                if (strcmp(variables[i].group, compiled_nodes[j].group) == 0 &&
                    strcmp(variables[i].name, compiled_nodes[j].name) == 0)
                {
                    snprintf(capture_error, sizeof(capture_error),
                             "compiled preshot node %s/%s must not be supplied", variables[i].group,
                             variables[i].name);
                    return -1;
                }
            }
        }
    }
    if (extra_count != N_EXTRA_PRESHOT_NODES)
    {
        snprintf(capture_error, sizeof(capture_error),
                 "capture needs all %d non-compiled preshot nodes", N_EXTRA_PRESHOT_NODES);
        return -1;
    }
    for (size_t i = 0; i < N_EXTRA_PRESHOT_NODES; i++)
    {
        const char *separator = strchr(extra_preshot_nodes[i], ':');
        int found = 0;
        if (separator == NULL)
        {
            return set_error("invalid internal preshot node table");
        }
        const size_t group_length = (size_t)(separator - extra_preshot_nodes[i]);
        for (size_t j = 0; j < n_variables; j++)
        {
            if (variables[j].group != NULL && variables[j].name != NULL &&
                strlen(variables[j].group) == group_length &&
                strncmp(variables[j].group, extra_preshot_nodes[i], group_length) == 0 &&
                strcmp(variables[j].name, separator + 1) == 0)
            {
                found = 1;
                break;
            }
        }
        if (!found)
        {
            snprintf(capture_error, sizeof(capture_error), "capture is missing preshot node %s",
                     extra_preshot_nodes[i]);
            return -1;
        }
    }

    const size_t total_variables = n_variables + RTGSFIT_SIZE_COUNT + N_COMPILED_NODES;
    combined_variables = calloc(total_variables, sizeof(*combined_variables));
    if (combined_variables == NULL)
    {
        return set_error("out of memory preparing compiled preshot variables");
    }
    for (size_t i = 0; i < RTGSFIT_SIZE_COUNT; i++)
    {
        size_variables[i] = (struct rtgsfit_capture_variable){
            .group = "PRESHOT/SIZES",
            .name = size_names[i],
            .type = RTGSFIT_CAPTURE_INT32,
            .n_dimensions = 0,
            .dimensions = NULL,
            .values = &compiled_sizes[i],
            .n_values = 1,
        };
        combined_variables[i] = size_variables[i];
    }
    memcpy(combined_variables + RTGSFIT_SIZE_COUNT, compiled_nodes, sizeof(compiled_nodes));
    memcpy(combined_variables + RTGSFIT_SIZE_COUNT + N_COMPILED_NODES, variables,
           n_variables * sizeof(*variables));
    for (size_t i = 0; i < RTGSFIT_SIZE_COUNT + N_COMPILED_NODES; i++)
    {
        struct rtgsfit_capture_variable *variable = &combined_variables[i];
        size_t count = 1;
        for (size_t j = 0; j < variable->n_dimensions; j++)
        {
            size_t dimension_size = 0;
            for (size_t k = 0; k < n_dimensions; k++)
            {
                if (strcmp(variable->dimensions[j], dimensions[k].name) == 0)
                {
                    dimension_size = dimensions[k].size;
                    break;
                }
            }
            count *= dimension_size;
        }
        variable->n_values = count;
    }
    variables = combined_variables;
    n_variables = total_variables;
    for (size_t i = 0; i < n_variables; i++)
    {
        if (validate_variable(&variables[i], dimensions, n_dimensions) != 0)
        {
            free(combined_variables);
            return -1;
        }
    }
    const uint32_t full_preshot_checksum = checksum_preshot(variables, n_variables);

    if (make_temporary_path(path, &temporary) != 0)
    {
        free(combined_variables);
        return -1;
    }
    status = nc_create(temporary, NC_NETCDF4 | NC_CLOBBER, &root);
    if (status != NC_NOERR)
    {
        set_netcdf_error("cannot create NetCDF-4 capture", status);
        goto failed;
    }

    status = nc_put_att_text(root, NC_GLOBAL, "schema", strlen(CAPTURE_SCHEMA), CAPTURE_SCHEMA);
    if (status == NC_NOERR)
    {
        status = nc_put_att_text(root, NC_GLOBAL, "flux_convention", strlen("native"), "native");
    }
    if (status == NC_NOERR)
    {
        status =
            nc_put_att_text(root, NC_GLOBAL, "settings_json", strlen(settings_json), settings_json);
    }
    if (status == NC_NOERR)
    {
        status = nc_put_att_text(root, NC_GLOBAL, "preshot_source", strlen(preshot_source),
                                 preshot_source);
    }
    if (status == NC_NOERR)
    {
        status =
            nc_put_att_text(root, NC_GLOBAL, "code_version", strlen(code_version), code_version);
    }
    if (status == NC_NOERR)
    {
        const uint32_t checksum = rtgsfit_get_preshot_checksum();
        status = nc_put_att_uint(root, NC_GLOBAL, "solver_checksum", NC_UINT, 1, &checksum);
    }
    if (status == NC_NOERR)
    {
        status = nc_put_att_uint(root, NC_GLOBAL, "preshot_checksum", NC_UINT, 1,
                                 &full_preshot_checksum);
    }
    if (status == NC_NOERR)
    {
        status = nc_put_att_int(root, NC_GLOBAL, "timing_enabled", NC_INT, 1, &timing_enabled);
    }
    if (status != NC_NOERR)
    {
        set_netcdf_error("cannot write capture metadata", status);
        goto failed;
    }

    for (size_t i = 0; i < n_dimensions; i++)
    {
        if (dimensions[i].name == NULL || dimensions[i].name[0] == '\0')
        {
            set_error("capture dimension has no name");
            goto failed;
        }
        for (size_t j = 0; j < i; j++)
        {
            if (strcmp(dimensions[i].name, dimensions[j].name) == 0)
            {
                snprintf(capture_error, sizeof(capture_error), "duplicate capture dimension %s",
                         dimensions[i].name);
                goto failed;
            }
        }
        const size_t length = dimensions[i].size == 0 ? NC_UNLIMITED : dimensions[i].size;
        status = nc_def_dim(root, dimensions[i].name, length, NULL);
        if (status != NC_NOERR)
        {
            set_netcdf_error("cannot define capture dimension", status);
            goto failed;
        }
    }

    for (size_t i = 0; i < 3; i++)
    {
        status = nc_def_grp(root, group_names[i], &groups[i]);
        if (status != NC_NOERR)
        {
            set_netcdf_error("cannot define capture group", status);
            goto failed;
        }
    }

    int *variable_ids = calloc(n_variables, sizeof(*variable_ids));
    int **dimension_ids = calloc(n_variables, sizeof(*dimension_ids));
    int *variable_groups = calloc(n_variables, sizeof(*variable_groups));
    if (variable_ids == NULL || dimension_ids == NULL || variable_groups == NULL)
    {
        free(variable_ids);
        free(dimension_ids);
        free(variable_groups);
        set_error("out of memory defining capture variables");
        goto failed;
    }

    for (size_t i = 0; i < n_variables; i++)
    {
        const struct rtgsfit_capture_variable *variable = &variables[i];
        if (variable->group == NULL || variable->name == NULL ||
            (variable->values == NULL && variable->n_values != 0) ||
            !valid_group(variable->group) || variable->n_dimensions > 8 ||
            (variable->n_dimensions != 0 && variable->dimensions == NULL))
        {
            snprintf(capture_error, sizeof(capture_error), "invalid capture variable at index %zu",
                     i);
            goto variable_failed;
        }
        status = get_variable_group(groups, variable, &variable_groups[i]);
        if (status != NC_NOERR)
        {
            set_netcdf_error("cannot define capture subgroup", status);
            goto variable_failed;
        }
        nc_type type;
        switch (variable->type)
        {
            case RTGSFIT_CAPTURE_FLOAT64:
                type = NC_DOUBLE;
                break;
            case RTGSFIT_CAPTURE_INT32:
                type = NC_INT;
                break;
            case RTGSFIT_CAPTURE_STRING:
                type = NC_STRING;
                break;
            default:
                snprintf(capture_error, sizeof(capture_error), "unsupported type for %s/%s",
                         variable->group, variable->name);
                goto variable_failed;
        }

        dimension_ids[i] = calloc(variable->n_dimensions == 0 ? 1 : variable->n_dimensions,
                                  sizeof(*dimension_ids[i]));
        if (dimension_ids[i] == NULL)
        {
            set_error("out of memory defining capture variable dimensions");
            goto variable_failed;
        }
        size_t expected_values = 1;
        for (size_t j = 0; j < variable->n_dimensions; j++)
        {
            size_t dimension_size = 0;
            if (variable->dimensions[j] == NULL)
            {
                snprintf(capture_error, sizeof(capture_error), "missing dimension name for %s/%s",
                         variable->group, variable->name);
                goto variable_failed;
            }
            status = find_dimension(dimensions, n_dimensions, root, variable->dimensions[j],
                                    &dimension_ids[i][j], &dimension_size);
            if (status != NC_NOERR)
            {
                snprintf(capture_error, sizeof(capture_error), "unknown dimension %s for %s/%s",
                         variable->dimensions[j], variable->group, variable->name);
                goto variable_failed;
            }
            if (dimension_size != 0 && expected_values > SIZE_MAX / dimension_size)
            {
                snprintf(capture_error, sizeof(capture_error), "element count overflows for %s/%s",
                         variable->group, variable->name);
                goto variable_failed;
            }
            expected_values *= dimension_size;
        }
        if (expected_values != variable->n_values)
        {
            snprintf(capture_error, sizeof(capture_error),
                     "element count does not match dimensions for %s/%s", variable->group,
                     variable->name);
            goto variable_failed;
        }
        status = nc_def_var(variable_groups[i], variable->name, type, (int)variable->n_dimensions,
                            dimension_ids[i], &variable_ids[i]);
        if (status != NC_NOERR)
        {
            set_netcdf_error("cannot define capture variable", status);
            goto variable_failed;
        }
    }

    status = nc_enddef(root);
    if (status != NC_NOERR)
    {
        set_netcdf_error("cannot finish capture definition", status);
        goto variable_failed;
    }

    for (size_t i = 0; i < n_variables; i++)
    {
        const struct rtgsfit_capture_variable *variable = &variables[i];
        if (variable->n_values == 0)
        {
            continue;
        }
        if (variable->type == RTGSFIT_CAPTURE_FLOAT64)
        {
            status = nc_put_var_double(variable_groups[i], variable_ids[i], variable->values);
        }
        else if (variable->type == RTGSFIT_CAPTURE_INT32)
        {
            status = nc_put_var_int(variable_groups[i], variable_ids[i], variable->values);
        }
        else
        {
            status = nc_put_var_string(variable_groups[i], variable_ids[i],
                                       (const char **)variable->values);
        }
        if (status != NC_NOERR)
        {
            snprintf(capture_error, sizeof(capture_error), "cannot write %s/%s: %s",
                     variable->group, variable->name, nc_strerror(status));
            goto variable_failed;
        }
    }

    for (size_t i = 0; i < n_variables; i++)
    {
        free(dimension_ids[i]);
    }
    free(dimension_ids);
    free(variable_ids);
    free(variable_groups);
    status = nc_close(root);
    root = -1;
    if (status != NC_NOERR)
    {
        set_netcdf_error("cannot close capture", status);
        goto failed;
    }
    if (rename(temporary, path) != 0)
    {
        snprintf(capture_error, sizeof(capture_error), "cannot replace %s: %s", path,
                 strerror(errno));
        goto failed;
    }
    free(temporary);
    free(combined_variables);
    return 0;

variable_failed:
    for (size_t i = 0; i < n_variables; i++)
    {
        free(dimension_ids[i]);
    }
    free(dimension_ids);
    free(variable_ids);
    free(variable_groups);
failed:
    if (root >= 0)
    {
        nc_close(root);
    }
    if (temporary != NULL)
    {
        unlink(temporary);
        free(temporary);
    }
    free(combined_variables);
    return -1;
}
