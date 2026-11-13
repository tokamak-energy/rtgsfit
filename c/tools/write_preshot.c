/*
 * write_preshot: writes the preshot as C source, which c/Makefile compiles into librtgsfit.
 *
 *   write_preshot_mdsplus --mdsplus <server> <pulse> <run name> <output directory>
 *   write_preshot --file <preshot file> <output directory>
 *   write_preshot --list
 *
 * `--mdsplus` reads `\RTGSFIT::TOP.<run name>.PRESHOT` of <pulse> from an MDSplus server. This is
 * how the PCS builds RT-GSFit, with no Python. `--file` reads a preshot file written by Python
 * (`rtgsfit.database_writers.c_source.write_preshot_file`), for preshots that are not in MDSplus,
 * such as the analytic and FreeGSNKE ones. `--list` prints the node table, which
 * tests/test_database_io.py checks against `rtgsfit.data_dictionary.C_NODES`.
 *
 * Two files are written into the output directory (normally c/src):
 *   preshot_sizes.h  every size as an enum constant, so that arrays have constant sizes (#97)
 *   preshot_data.c   every array declared in preshot.h, row-major and flattened, with 17
 *                    significant digits so that each double is reproduced exactly, and
 *                    PRESHOT_CHECKSUM, the CRC-32 of all the values
 * The whole preshot is read before either file is written, and a file is replaced only when its
 * content changes, so that make recompiles only what changed.
 *
 * The preshot file is "RTGSFIT PRESHOT 1\n" and "SOURCE <INFO:SOURCE>\n", then for each node, in
 * the order of SIZES and NODES below, "<GROUP>:<NAME> <float64|int32> <count>\n" followed by the
 * values, raw and little-endian, and finally "END\n".
 */
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "plasma_boundary.h"
#include "rtgsfit.h"

#ifdef RTGSFIT_MDSPLUS
#include <mdslib.h>
#endif

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "the preshot file holds little-endian values"
#endif

/* The sizes, in the order of rtgsfit.data_dictionary.SIZE_NAMES and enum rtgsfit_size_index. */
static const struct
{
    const char *name;
    const char *description;
} SIZES[] = {
    {"N_R", "Number of grid points in R"},
    {"N_Z", "Number of grid points in Z"},
    {"N_GRID", "Number of grid points, N_R * N_Z"},
    {"N_EDGE", "Number of grid points on the edge of the grid, 2 N_R + 2 N_Z - 4"},
    {"N_LIMITER", "Number of limiter points"},
    {"N_COIL", "Number of PF coil circuits"},
    {"N_MEAS", "Number of PCS measurements, the length of meas"},
    {"N_FLUX_LOOP", "Number of flux loop constraints"},
    {"N_BP_PROBE", "Number of Bp probe constraints"},
    {"N_ROGOWSKI", "Number of Rogowski coil constraints"},
    {"N_CONSTRAINT", "Number of constraints, N_FLUX_LOOP + N_BP_PROBE + N_ROGOWSKI"},
    {"N_REG", "Number of regularisation rows of the least-squares problem"},
    {"N_FIT_TARGET", "Number of rows of the least-squares problem, N_CONSTRAINT + N_REG"},
    {"N_PLASMA", "Number of plasma degrees of freedom"},
    {"N_PASSIVE", "Number of passive-structure degrees of freedom"},
    {"N_FIT_DOF", "Number of unknowns of the least-squares problem, N_PLASMA + N_PASSIVE"},
    {"N_XPT_MAX", "Maximum number of x-points that are stored"},
};

enum
{
    N_SIZES = sizeof(SIZES) / sizeof(SIZES[0])
};
_Static_assert((int)N_SIZES == (int)RTGSFIT_SIZE_COUNT, "SIZES must match enum rtgsfit_size_index");

enum type
{
    FLOAT64,
    INT32
};

/*
 * The arrays and scalars declared in preshot.h, in the order of rtgsfit.data_dictionary.C_NODES.
 * The C name is GROUP_NAME. `dims` are the names of the dimensions; a scalar has none.
 */
static const struct node
{
    const char *group;
    const char *name;
    enum type type;
    const char *dims[2];
} NODES[] = {
    {"GRID", "R", FLOAT64, {"N_R"}},
    {"GRID", "Z", FLOAT64, {"N_Z"}},
    {"GRID", "D_R", FLOAT64, {NULL}},
    {"GRID", "D_Z", FLOAT64, {NULL}},
    {"GRID", "R_FLAT", FLOAT64, {"N_GRID"}},
    {"GRID", "Z_FLAT", FLOAT64, {"N_GRID"}},
    {"GRID", "INV_R", FLOAT64, {"N_GRID"}},
    {"GRID", "MASK_LIMITER", INT32, {"N_GRID"}},
    {"GRID", "G_COIL", FLOAT64, {"N_GRID", "N_COIL"}},
    {"GRID", "G_PASSIVE", FLOAT64, {"N_GRID", "N_PASSIVE"}},
    {"EDGE", "R", FLOAT64, {"N_EDGE"}},
    {"EDGE", "Z", FLOAT64, {"N_EDGE"}},
    {"EDGE", "INV_MU0_R", FLOAT64, {"N_EDGE"}},
    {"EDGE", "G_EDGE", FLOAT64, {"N_EDGE", "N_EDGE"}},
    {"LIMITER", "R", FLOAT64, {"N_LIMITER"}},
    {"LIMITER", "Z", FLOAT64, {"N_LIMITER"}},
    {"LIMITER", "GRID_INDEX", INT32, {"N_LIMITER", "N_INTERP"}},
    {"LIMITER", "GRID_WEIGHT", FLOAT64, {"N_LIMITER", "N_INTERP"}},
    {"POISSON", "RHS_SCALE", FLOAT64, {"N_GRID"}},
    {"POISSON", "LOWER_BAND", FLOAT64, {"N_GRID", "N_R"}},
    {"POISSON", "UPPER_BAND", FLOAT64, {"N_GRID", "N_R_PLUS_2"}},
    {"POISSON", "PERMUTATION", INT32, {"N_GRID"}},
    {"CONSTRAINT", "G_COIL", FLOAT64, {"N_CONSTRAINT", "N_COIL"}},
    {"CONSTRAINT", "REPLACEMENT", FLOAT64, {"N_CONSTRAINT", "N_MEAS"}},
    {"FIT", "WEIGHT", FLOAT64, {"N_FIT_TARGET"}},
    {"FIT", "G_GRID", FLOAT64, {"N_FIT_TARGET", "N_GRID"}},
    {"FIT", "G_PASSIVE", FLOAT64, {"N_FIT_TARGET", "N_PASSIVE"}},
    {"SETTINGS", "XPT_FRACTION", FLOAT64, {NULL}},
    {"SETTINGS", "TOLERANCE", FLOAT64, {NULL}},
    {"SETTINGS", "IP_MIN", FLOAT64, {NULL}},
};

enum
{
    N_NODES = sizeof(NODES) / sizeof(NODES[0]),
    VALUES_PER_LINE = 4,
    MAX_TEXT = 512
};

/* preshot.h pads only these dimensions to one element when they are zero (PRESHOT_AT_LEAST_ONE). */
static const char *const OPTIONAL_DIMS[] = {"N_PASSIVE"};

static const char *const TYPE_NAMES[] = {[FLOAT64] = "float64", [INT32] = "int32"};
static const char *const C_TYPES[] = {[FLOAT64] = "double", [INT32] = "int32_t"};
static const size_t TYPE_BYTES[] = {[FLOAT64] = sizeof(double), [INT32] = sizeof(int32_t)};

struct preshot
{
    char source[MAX_TEXT]; /* INFO:SOURCE, for a comment in preshot_data.c */
    int32_t sizes[N_SIZES];
    void *values[N_NODES]; /* padded as in preshot.h; the padding is zeros */
    uint32_t checksum;     /* CRC-32 of the values of every size and node, as read */
};

/* Where the preshot is read from: a file, or (run_name != NULL) the connected MDSplus server. */
struct source
{
    FILE *file;
    const char *path;
    const char *run_name;
};

_Noreturn static void fail(const char *format, ...) __attribute__((format(printf, 1, 2)));

static void fail(const char *format, ...)
{
    va_list arguments;
    va_start(arguments, format);
    fputs("write_preshot: ", stderr);
    vfprintf(stderr, format, arguments);
    fputc('\n', stderr);
    va_end(arguments);
    exit(EXIT_FAILURE);
}

/* The CRC-32 of zlib (Python's `zlib.crc32`), continued from `crc`. */
static uint32_t crc32_update(uint32_t crc, const void *data, const size_t n_bytes)
{
    static uint32_t table[256];
    if (table[1] == 0)
    {
        for (uint32_t i = 0; i < 256; i++)
        {
            uint32_t c = i;
            for (int bit = 0; bit < 8; bit++)
            {
                c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            }
            table[i] = c;
        }
    }
    const unsigned char *bytes = data;
    crc = ~crc;
    for (size_t i = 0; i < n_bytes; i++)
    {
        crc = table[(crc ^ bytes[i]) & 0xFFu] ^ (crc >> 8);
    }
    return ~crc;
}

/* The value of a dimension, given the sizes. N_INTERP and N_R_PLUS_2 are derived, not stored. */
static int64_t dimension(const char *name, const int32_t sizes[])
{
    if (strcmp(name, "N_INTERP") == 0)
    {
        return PLASMA_BOUNDARY_N_INTERP;
    }
    if (strcmp(name, "N_R_PLUS_2") == 0)
    {
        return dimension("N_R", sizes) + 2;
    }
    for (int i = 0; i < N_SIZES; i++)
    {
        if (strcmp(name, SIZES[i].name) == 0)
        {
            return sizes[i];
        }
    }
    fail("unknown dimension %s", name);
}

static int is_optional(const char *dim)
{
    for (size_t i = 0; i < sizeof(OPTIONAL_DIMS) / sizeof(OPTIONAL_DIMS[0]); i++)
    {
        if (strcmp(dim, OPTIONAL_DIMS[i]) == 0)
        {
            return 1;
        }
    }
    return 0;
}

/* The number of values of `node`; with `padded`, as declared in preshot.h (zero dims become 1). */
static int64_t count_values(const struct node *node, const int32_t sizes[], const int padded)
{
    int64_t count = 1;
    for (int k = 0; k < 2 && node->dims[k] != NULL; k++)
    {
        const int64_t n = dimension(node->dims[k], sizes);
        if (n == 0 && !is_optional(node->dims[k]))
        {
            fail("%s:%s has %s = 0, which preshot.h does not allow", node->group, node->name,
                 node->dims[k]);
        }
        count *= (padded && n == 0) ? 1 : n;
    }
    return count;
}

static void file_read_node(struct source *source,
                           const char *group,
                           const char *name,
                           const enum type type,
                           const int64_t count,
                           void *values)
{
    char line[MAX_TEXT], found_name[MAX_TEXT], found_type[16];
    long long found_count = -1;
    if (fgets(line, sizeof(line), source->file) == NULL ||
        sscanf(line, "%511s %15s %lld", found_name, found_type, &found_count) != 3)
    {
        fail("%s: cannot read the header of %s:%s", source->path, group, name);
    }
    char expected_name[MAX_TEXT];
    snprintf(expected_name, sizeof(expected_name), "%s:%s", group, name);
    if (strcmp(found_name, expected_name) != 0 || strcmp(found_type, TYPE_NAMES[type]) != 0 ||
        found_count != count)
    {
        fail("%s: expected %s %s %" PRId64 ", found %s %s %lld", source->path, expected_name,
             TYPE_NAMES[type], count, found_name, found_type, found_count);
    }
    if (count > 0 && fread(values, TYPE_BYTES[type], (size_t)count, source->file) != (size_t)count)
    {
        fail("%s: the file ends inside %s", source->path, expected_name);
    }
}

#ifdef RTGSFIT_MDSPLUS
static void mdsplus_connect(char *server, const int32_t pulse)
{
    if (MdsConnect(server) == -1)
    {
        fail("cannot connect to the MDSplus server %s", server);
    }
    char tree[] = "RTGSFIT";
    int shot = pulse;
    const int status = MdsOpen(tree, &shot);
    if (!(status & 1))
    {
        fail("cannot open the RTGSFIT tree of pulse %d on %s (MDSplus status %d)", shot, server,
             status);
    }
}

static void mdsplus_read_node(struct source *source,
                              const char *group,
                              const char *name,
                              const enum type type,
                              const int64_t count,
                              void *values)
{
    char path[MAX_TEXT], expression[MAX_TEXT + 16];
    snprintf(path, sizeof(path), "\\RTGSFIT::TOP.%s.PRESHOT.%s:%s", source->run_name, group, name);
    int dtype_long = DTYPE_LONG, dtype_double = DTYPE_DOUBLE, null = 0, length = 0;
    int stored = -1;
    snprintf(expression, sizeof(expression), "SIZE(%s)", path);
    int descriptor = descr(&dtype_long, &stored, &null);
    int status = MdsValue(expression, &descriptor, &null, &length);
    if (!(status & 1))
    {
        fail("cannot read %s (MDSplus status %d)", path, status);
    }
    if (stored != count)
    {
        fail("%s has %d values, expected %" PRId64, path, stored, count);
    }
    if (count == 0)
    {
        return;
    }
    int n = (int)count;
    /* Convert on the server, so that the values arrive as IEEE doubles and 32-bit integers. */
    snprintf(expression, sizeof(expression), type == FLOAT64 ? "FT_FLOAT(%s)" : "LONG(%s)", path);
    descriptor = descr(type == FLOAT64 ? &dtype_double : &dtype_long, values, &n, &null);
    status = MdsValue(expression, &descriptor, &null, &length);
    if (!(status & 1) || length != n)
    {
        fail("cannot read %s (MDSplus status %d, %d of %d values)", path, status, length, n);
    }
}

static void mdsplus_read_source(const char *run_name, char source[], const int n_source)
{
    char expression[MAX_TEXT];
    snprintf(expression, sizeof(expression), "\\RTGSFIT::TOP.%s.PRESHOT.INFO:SOURCE", run_name);
    int dtype_cstring = DTYPE_CSTRING, null = 0, length = 0, size = n_source - 1;
    memset(source, 0, (size_t)n_source);
    int descriptor = descr(&dtype_cstring, source, &null, &size);
    const int status = MdsValue(expression, &descriptor, &null, &length);
    if (!(status & 1))
    {
        fail("cannot read %s (MDSplus status %d)", expression, status);
    }
    source[length < n_source ? length : n_source - 1] = '\0';
}
#endif

static void read_node(struct source *source,
                      const char *group,
                      const char *name,
                      const enum type type,
                      const int64_t count,
                      void *values)
{
#ifdef RTGSFIT_MDSPLUS
    if (source->run_name != NULL)
    {
        mdsplus_read_node(source, group, name, type, count, values);
        return;
    }
#endif
    file_read_node(source, group, name, type, count, values);
}

/* Read and check every size and node from `source`. */
static void read_preshot(struct source *source, struct preshot *preshot)
{
    preshot->checksum = 0;
    for (int i = 0; i < N_SIZES; i++)
    {
        read_node(source, "SIZES", SIZES[i].name, INT32, 1, &preshot->sizes[i]);
        if (preshot->sizes[i] < 0)
        {
            fail("SIZES:%s = %" PRId32 " is negative", SIZES[i].name, preshot->sizes[i]);
        }
        preshot->checksum = crc32_update(preshot->checksum, &preshot->sizes[i], sizeof(int32_t));
    }
    for (int i = 0; i < N_NODES; i++)
    {
        const struct node *node = &NODES[i];
        const int64_t count = count_values(node, preshot->sizes, 0);
        preshot->values[i] =
            calloc((size_t)count_values(node, preshot->sizes, 1), TYPE_BYTES[node->type]);
        if (preshot->values[i] == NULL)
        {
            fail("out of memory for %s:%s", node->group, node->name);
        }
        read_node(source, node->group, node->name, node->type, count, preshot->values[i]);
        if (node->type == FLOAT64)
        {
            const double *values = preshot->values[i];
            for (int64_t k = 0; k < count; k++)
            {
                if (!isfinite(values[k]))
                {
                    fail("%s:%s[%" PRId64 "] is not finite", node->group, node->name, k);
                }
            }
        }
        preshot->checksum = crc32_update(preshot->checksum, preshot->values[i],
                                         (size_t)count * TYPE_BYTES[node->type]);
    }
}

static void read_preshot_file(const char *path, struct preshot *preshot)
{
    struct source source = {.file = fopen(path, "rb"), .path = path, .run_name = NULL};
    if (source.file == NULL)
    {
        fail("cannot open %s: %s", path, strerror(errno));
    }
    char line[MAX_TEXT];
    if (fgets(line, sizeof(line), source.file) == NULL || strcmp(line, "RTGSFIT PRESHOT 1\n") != 0)
    {
        fail("%s is not an RT-GSFit preshot file (version 1)", path);
    }
    if (fgets(line, sizeof(line), source.file) == NULL || strncmp(line, "SOURCE ", 7) != 0)
    {
        fail("%s: the second line must be SOURCE <INFO:SOURCE>", path);
    }
    char *newline = strchr(line, '\n');
    if (newline == NULL)
    {
        fail("%s: SOURCE line is too long", path);
    }
    *newline = '\0';
    const size_t source_length = strlen(line + 7);
    if (source_length >= sizeof(preshot->source))
    {
        fail("%s: SOURCE value is too long", path);
    }
    memcpy(preshot->source, line + 7, source_length + 1);
    read_preshot(&source, preshot);
    if (fgets(line, sizeof(line), source.file) == NULL || strcmp(line, "END\n") != 0 ||
        fgetc(source.file) != EOF)
    {
        fail("%s: expected END after the last node", path);
    }
    fclose(source.file);
}

/* A comment cannot contain its own terminator, so a star-slash in the text gets a space. */
static void write_comment_text(FILE *out, const char *text)
{
    for (const char *c = text; *c != '\0'; c++)
    {
        fputc(*c == '\n' ? ' ' : *c, out);
        if (c[0] == '*' && c[1] == '/')
        {
            fputc(' ', out);
        }
    }
}

static void write_value(FILE *out, const enum type type, const void *values, const int64_t k)
{
    if (type == INT32)
    {
        fprintf(out, "%" PRId32, ((const int32_t *)values)[k]);
        return;
    }
    char text[32];
    snprintf(text, sizeof(text), "%.17g", ((const double *)values)[k]);
    /* Keep it a floating-point literal: 1 is written 1.0 */
    fprintf(out, "%s%s", text, strpbrk(text, ".e") == NULL ? ".0" : "");
}

static void write_sizes_header(FILE *out, const struct preshot *preshot)
{
    int width = 0;
    for (int i = 0; i < N_SIZES; i++)
    {
        const int length = (int)strlen(SIZES[i].name);
        width = length > width ? length : width;
    }
    fputs("/* Generated by c/tools/write_preshot.c. Do not edit. */\n"
          "/* clang-format off */\n"
          "#ifndef RTGSFIT_PRESHOT_SIZES_H\n"
          "#define RTGSFIT_PRESHOT_SIZES_H\n"
          "\n"
          "enum {\n",
          out);
    for (int i = 0; i < N_SIZES; i++)
    {
        fprintf(out, "    %-*s = %8" PRId32 ", /* %s */\n", width, SIZES[i].name, preshot->sizes[i],
                SIZES[i].description);
    }
    fputs("};\n\n#endif /* RTGSFIT_PRESHOT_SIZES_H */\n", out);
}

static void write_data_source(FILE *out, const struct preshot *preshot)
{
    fputs("/* Generated by c/tools/write_preshot.c. Do not edit. */\n"
          "/* clang-format off */\n"
          "/* ",
          out);
    write_comment_text(out, preshot->source);
    fprintf(out,
            " */\n#include \"preshot.h\"\n\nconst uint32_t PRESHOT_CHECKSUM = 0x%08" PRIX32 "u;\n",
            preshot->checksum);
    for (int i = 0; i < N_NODES; i++)
    {
        const struct node *node = &NODES[i];
        const char *c_type = C_TYPES[node->type];
        if (node->dims[0] == NULL)
        {
            fprintf(out, "\nconst %s %s_%s = ", c_type, node->group, node->name);
            write_value(out, node->type, preshot->values[i], 0);
            fputs(";\n", out);
            continue;
        }
        const int64_t count = count_values(node, preshot->sizes, 1);
        fprintf(out, "\nconst %s %s_%s[%" PRId64 "] = {\n    ", c_type, node->group, node->name,
                count);
        for (int64_t k = 0; k < count; k++)
        {
            write_value(out, node->type, preshot->values[i], k);
            if (k + 1 < count)
            {
                fputs((k + 1) % VALUES_PER_LINE == 0 ? ",\n    " : ", ", out);
            }
        }
        fputs("\n};\n", out);
    }
}

static int files_equal(const char *path_a, const char *path_b)
{
    FILE *a = fopen(path_a, "rb");
    FILE *b = fopen(path_b, "rb");
    int equal = a != NULL && b != NULL;
    while (equal)
    {
        static char buffer_a[1 << 16], buffer_b[1 << 16];
        const size_t n_a = fread(buffer_a, 1, sizeof(buffer_a), a);
        const size_t n_b = fread(buffer_b, 1, sizeof(buffer_b), b);
        equal = n_a == n_b && memcmp(buffer_a, buffer_b, n_a) == 0;
        if (n_a == 0)
        {
            break;
        }
    }
    if (a != NULL)
    {
        fclose(a);
    }
    if (b != NULL)
    {
        fclose(b);
    }
    return equal;
}

/* Write `directory`/`name` with `write`, replacing the file only if its content changed. */
static void write_if_changed(const char *directory,
                             const char *name,
                             void (*write)(FILE *, const struct preshot *),
                             const struct preshot *preshot)
{
    const size_t directory_length = strlen(directory);
    const size_t name_length = strlen(name);
    char *path = malloc(directory_length + name_length + 2);
    char *temporary = malloc(directory_length + name_length + 8);
    if (path == NULL || temporary == NULL)
    {
        free(path);
        free(temporary);
        fail("out of memory creating output path");
    }
    snprintf(path, directory_length + name_length + 2, "%s/%s", directory, name);
    snprintf(temporary, directory_length + name_length + 8, "%s/.%s.tmp", directory, name);
    FILE *out = fopen(temporary, "w");
    if (out == NULL)
    {
        fail("cannot write %s: %s", temporary, strerror(errno));
    }
    write(out, preshot);
    if (ferror(out) || fclose(out) != 0)
    {
        fail("cannot write %s", temporary);
    }
    if (files_equal(temporary, path))
    {
        remove(temporary);
        printf("%s is unchanged\n", path);
    }
    else if (rename(temporary, path) != 0)
    {
        fail("cannot replace %s: %s", path, strerror(errno));
    }
    else
    {
        printf("Wrote %s\n", path);
    }
    free(path);
    free(temporary);
}

static void list_nodes(void)
{
    for (int i = 0; i < N_SIZES; i++)
    {
        printf("SIZES:%s int32\n", SIZES[i].name);
    }
    for (int i = 0; i < N_NODES; i++)
    {
        printf("%s:%s %s", NODES[i].group, NODES[i].name, TYPE_NAMES[NODES[i].type]);
        for (int k = 0; k < 2 && NODES[i].dims[k] != NULL; k++)
        {
            printf(" %s", NODES[i].dims[k]);
        }
        printf("\n");
    }
}

_Noreturn static void usage(void)
{
    fail("usage:\n"
#ifdef RTGSFIT_MDSPLUS
         "  write_preshot_mdsplus --mdsplus <server> <pulse> <run name> <output directory>\n"
#endif
         "  write_preshot --file <preshot file> <output directory>\n"
         "  write_preshot --list");
}

int main(int argc, char *argv[])
{
    static struct preshot preshot;
    const char *output_dir = NULL;
    if (argc == 2 && strcmp(argv[1], "--list") == 0)
    {
        list_nodes();
        return EXIT_SUCCESS;
    }
    if (argc == 4 && strcmp(argv[1], "--file") == 0)
    {
        read_preshot_file(argv[2], &preshot);
        output_dir = argv[3];
    }
    else if (argc == 6 && strcmp(argv[1], "--mdsplus") == 0)
    {
#ifdef RTGSFIT_MDSPLUS
        char *end = NULL;
        const long pulse = strtol(argv[3], &end, 10);
        if (*argv[3] == '\0' || *end != '\0' || pulse <= 0 || pulse > INT32_MAX)
        {
            fail("the pulse must be a positive integer, not %s", argv[3]);
        }
        mdsplus_connect(argv[2], (int32_t)pulse);
        struct source source = {.file = NULL, .path = NULL, .run_name = argv[4]};
        mdsplus_read_source(source.run_name, preshot.source, (int)sizeof(preshot.source));
        read_preshot(&source, &preshot);
        MdsDisconnect();
        output_dir = argv[5];
#else
        fail("this write_preshot was built without MDSplus; use write_preshot_mdsplus");
#endif
    }
    else
    {
        usage();
    }
    write_if_changed(output_dir, "preshot_sizes.h", write_sizes_header, &preshot);
    write_if_changed(output_dir, "preshot_data.c", write_data_source, &preshot);
    return EXIT_SUCCESS;
}
