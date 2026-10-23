# RT-GSFit C library.
#
#   make                                 build $(BUILD_DIR)/librtgsfit.so from $(PRESHOT_DIR)
#   make SHOT=<pulse> RUN_NAME=<RUNxx>   first write the preshot from MDSplus (TE network only)
#   make test                            build and run the C unit tests (no preshot needed)
#   make clean
#
# Options: DEBUG=1 (-O0, gprof), TIMING=1 (per-section timers), ARCH_FLAGS=..., PYTHON=...
#
# The preshot directory must contain `preshot_sizes.h` and `preshot_data.c`, written by
# `rtgsfit.database_writers.c_source` (for example via `python -m rtgsfit preshot ...`).

CC ?= gcc
PYTHON ?= python3
DEBUG ?= 0
TIMING ?= 0
ARCH_FLAGS ?= -march=native

C_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SRC_DIR := $(C_DIR)/src
TEST_DIR := $(C_DIR)/tests
PRESHOT_DIR ?= $(C_DIR)/preshot
BUILD_DIR ?= $(C_DIR)/build
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
MODULES := finite_difference current_basis least_squares poisson_solver magnetic_nulls \
           plasma_boundary timing
MODULE_OBJS := $(patsubst %,$(OBJ_DIR)/%.o,$(MODULES))
HEADERS := $(wildcard $(SRC_DIR)/*.h)

LIBRARY := $(BUILD_DIR)/librtgsfit.so

.PHONY: all preshot test clean

all: $(LIBRARY)

ifdef SHOT
$(PRESHOT_DIR)/preshot_sizes.h $(PRESHOT_DIR)/preshot_data.c: preshot
preshot:
	$(PYTHON) -m rtgsfit write-c-preshot --pulse $(SHOT) --run-name $(RUN_NAME) \
	    --output-dir $(PRESHOT_DIR)
else
$(PRESHOT_DIR)/preshot_sizes.h $(PRESHOT_DIR)/preshot_data.c:
	$(error No preshot in $(PRESHOT_DIR). Either run 'make SHOT=<pulse> RUN_NAME=<RUNxx>' \
	    (TE network) or write one with 'python -m rtgsfit write-c-preshot ...')
endif

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/rtgsfit.o: $(SRC_DIR)/rtgsfit.c $(HEADERS) $(PRESHOT_DIR)/preshot_sizes.h | $(OBJ_DIR)
	$(CC) $(CFLAGS) -I$(PRESHOT_DIR) -c $< -o $@

$(OBJ_DIR)/preshot_data.o: $(PRESHOT_DIR)/preshot_data.c $(HEADERS) \
                           $(PRESHOT_DIR)/preshot_sizes.h | $(OBJ_DIR)
	$(CC) $(CFLAGS) -O0 -I$(SRC_DIR) -I$(PRESHOT_DIR) -c $< -o $@

$(LIBRARY): $(OBJ_DIR)/rtgsfit.o $(OBJ_DIR)/preshot_data.o $(MODULE_OBJS)
	$(CC) $(LDFLAGS) -Wl,--no-undefined -o $@ $^ -lm $(BLAS_LIBS)

# --------------------------------------------------------------------------------------------
# Unit tests: each tests/test_<name>.c is linked against the preshot-independent modules.
# --------------------------------------------------------------------------------------------
TESTS := $(patsubst $(TEST_DIR)/%.c,$(BUILD_DIR)/tests/%,$(wildcard $(TEST_DIR)/test_*.c))

$(BUILD_DIR)/tests:
	mkdir -p $@

$(BUILD_DIR)/tests/%: $(TEST_DIR)/%.c $(TEST_DIR)/test_utilities.h $(MODULE_OBJS) \
                      | $(BUILD_DIR)/tests
	$(CC) $(CFLAGS) -I$(SRC_DIR) -o $@ $< $(MODULE_OBJS) -lm $(BLAS_LIBS)

test: $(TESTS)
	@set -e; for t in $(TESTS); do echo "== $$t"; (cd $(TEST_DIR) && $$t); done

clean:
	rm -rf $(BUILD_DIR)
