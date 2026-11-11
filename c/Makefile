# RT-GSFit C library.
#
#   make SHOT=<pulse> RUN_NAME=<run>   read the preshot from MDSplus, then build (no Python needed)
#   make PRESHOT_FILE=<path>           read the preshot from a file written by Python, then build
#   make                               build with the preshot already generated in src/
#   make test                          build and run the C unit tests (no preshot needed)
#   make clean                         remove build/ and the generated preshot
#
# Options: MDSPLUS_SERVER=smaug, MDSPLUS_DIR=/usr/local/mdsplus, DEBUG=1 (-O0, gprof),
# TIMING=1 (per-section timers), ARCH_FLAGS=...
#
# The preshot is compiled into the library: build/write_preshot (c/tools/write_preshot.c) writes
# src/preshot_sizes.h and src/preshot_data.c, which are not version-controlled, and rewrites them
# only when they change. The output is build/librtgsfit.so. There is one active preshot at a time,
# so run one build at a time.

CC ?= gcc
DEBUG ?= 0
TIMING ?= 0
ARCH_FLAGS ?= -march=native
MDSPLUS_SERVER ?= smaug
MDSPLUS_DIR ?= /usr/local/mdsplus

C_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SRC_DIR := $(C_DIR)/src
TEST_DIR := $(C_DIR)/tests
TOOLS_DIR := $(C_DIR)/tools
BUILD_DIR := $(C_DIR)/build
OBJ_DIR := $(BUILD_DIR)/obj

BLAS_CFLAGS ?= -I/usr/include/openblas
BLAS_LIBS ?= -lopenblas

WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wvla -Wstrict-prototypes -Warray-bounds
CFLAGS_COMMON := -std=c11 -fPIC -pthread -g -fno-omit-frame-pointer $(WARNINGS) $(BLAS_CFLAGS)

ifeq ($(DEBUG), 1)
  CFLAGS := $(CFLAGS_COMMON) -O0 -pg -DDEBUG
  LDFLAGS := -shared -pthread -g -pg
else
  CFLAGS := $(CFLAGS_COMMON) -O3 $(ARCH_FLAGS)
  LDFLAGS := -shared -pthread -g
endif
ifeq ($(TIMING), 1)
  CFLAGS += -DRTGSFIT_TIMING
endif

# Modules that do not depend on the preshot.
MODULES := finite_difference current_basis least_squares poisson_solver magnetic_nulls measurement_processing \
           plasma_boundary timing
MODULE_OBJS := $(patsubst %,$(OBJ_DIR)/%.o,$(MODULES))
CAPTURE_OBJ := $(OBJ_DIR)/capture.o
# Only rtgsfit.c and preshot_data.c depend on the generated sizes.
HEADERS := $(filter-out $(SRC_DIR)/preshot_sizes.h,$(wildcard $(SRC_DIR)/*.h))
GENERATED := $(SRC_DIR)/preshot_sizes.h $(SRC_DIR)/preshot_data.c

LIBRARY := $(BUILD_DIR)/librtgsfit.so
GENERATOR := $(BUILD_DIR)/write_preshot
GENERATOR_MDSPLUS := $(BUILD_DIR)/write_preshot_mdsplus
TOOL_CFLAGS := -std=c11 -O2 -g $(WARNINGS) -I$(SRC_DIR)

# Every object depends on FLAGS_STAMP, which is rewritten when the flags change, or the CPU that
# -march=native means, so that a change of options or of machine rebuilds everything.
FLAGS_STAMP := $(OBJ_DIR)/flags
BUILD_FLAGS := $(CC) $(CFLAGS) $(LDFLAGS) $(BLAS_LIBS) \
               $(if $(findstring native,$(ARCH_FLAGS)),$(shell $(CC) $(ARCH_FLAGS) -Q --help=target 2>/dev/null | cksum))

.PHONY: all preshot tools tools-mdsplus test clean FORCE

all: $(LIBRARY)

tools: $(GENERATOR)

tools-mdsplus: $(GENERATOR_MDSPLUS)

# --------------------------------------------------------------------------------------------
# The preshot. `preshot` always runs the generator; make then checks again whether the generated
# files changed, and recompiles only if they did.
# --------------------------------------------------------------------------------------------
ifneq ($(and $(SHOT),$(PRESHOT_FILE)),)
  $(error give SHOT or PRESHOT_FILE, not both)
endif

ifdef SHOT
  ifndef RUN_NAME
    $(error RUN_NAME is needed with SHOT, e.g. make SHOT=99000230 RUN_NAME=TEST03)
  endif
preshot: $(GENERATOR_MDSPLUS)
	$< --mdsplus $(MDSPLUS_SERVER) $(SHOT) $(RUN_NAME) $(SRC_DIR)
$(GENERATED): preshot
	@:
else ifdef PRESHOT_FILE
preshot: $(GENERATOR)
	$< --file $(PRESHOT_FILE) $(SRC_DIR)
$(GENERATED): preshot
	@:
else
preshot:
	$(error give SHOT=<pulse> RUN_NAME=<run> or PRESHOT_FILE=<path>)
$(GENERATED):
	@echo "No preshot in $(SRC_DIR). Run 'make SHOT=<pulse> RUN_NAME=<run>' to read one from" >&2
	@echo "MDSplus, or build from Python (rtgsfit.library.build_library)." >&2
	@exit 1
endif

$(GENERATOR): $(TOOLS_DIR)/write_preshot.c $(HEADERS) | $(BUILD_DIR)
	$(CC) $(TOOL_CFLAGS) $< -o $@

$(GENERATOR_MDSPLUS): $(TOOLS_DIR)/write_preshot.c $(HEADERS) | $(BUILD_DIR)
	$(CC) $(TOOL_CFLAGS) -DRTGSFIT_MDSPLUS -isystem $(MDSPLUS_DIR)/include $< -o $@ \
	    -L$(MDSPLUS_DIR)/lib -Wl,-rpath,$(MDSPLUS_DIR)/lib -lMdsLib

# --------------------------------------------------------------------------------------------
# The library
# --------------------------------------------------------------------------------------------
$(BUILD_DIR) $(OBJ_DIR) $(BUILD_DIR)/tests:
	mkdir -p $@

$(FLAGS_STAMP): FORCE | $(OBJ_DIR)
	@echo '$(BUILD_FLAGS)' | cmp -s - $@ || echo '$(BUILD_FLAGS)' > $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) $(FLAGS_STAMP) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/rtgsfit.o: $(SRC_DIR)/rtgsfit.c $(HEADERS) $(SRC_DIR)/preshot_sizes.h $(FLAGS_STAMP) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(CAPTURE_OBJ): $(SRC_DIR)/capture.c $(SRC_DIR)/capture.h $(SRC_DIR)/rtgsfit.h $(HEADERS) \
                $(SRC_DIR)/preshot_sizes.h $(FLAGS_STAMP) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# The data is only read, so -O0 loses nothing and compiles the large arrays much faster.
$(OBJ_DIR)/preshot_data.o: $(SRC_DIR)/preshot_data.c $(SRC_DIR)/preshot_sizes.h $(HEADERS) $(FLAGS_STAMP) \
                           | $(OBJ_DIR)
	$(CC) $(CFLAGS) -O0 -c $< -o $@

$(LIBRARY): $(OBJ_DIR)/rtgsfit.o $(OBJ_DIR)/preshot_data.o $(MODULE_OBJS) $(CAPTURE_OBJ) $(FLAGS_STAMP)
	$(CC) $(LDFLAGS) -Wl,--no-undefined -o $@ $(filter %.o,$^) -lm $(BLAS_LIBS) -lnetcdf

# --------------------------------------------------------------------------------------------
# Unit tests: each tests/test_<name>.c is linked against the preshot-independent modules.
# --------------------------------------------------------------------------------------------
TESTS := $(patsubst $(TEST_DIR)/%.c,$(BUILD_DIR)/tests/%,$(wildcard $(TEST_DIR)/test_*.c))

$(BUILD_DIR)/tests/%: $(TEST_DIR)/%.c $(TEST_DIR)/test_utilities.h $(HEADERS) $(MODULE_OBJS) $(FLAGS_STAMP) \
                      | $(BUILD_DIR)/tests
	$(CC) $(CFLAGS) -I$(SRC_DIR) -o $@ $< $(MODULE_OBJS) -lm $(BLAS_LIBS)

test: $(TESTS)
	@set -e; for t in $(TESTS); do echo "== $$t"; (cd $(TEST_DIR) && $$t); done

clean:
	rm -rf $(BUILD_DIR) $(GENERATED)

FORCE:
