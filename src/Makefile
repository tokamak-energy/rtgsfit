# --------------------------------------------------
# Required runtime parameters
# --------------------------------------------------
ifndef SHOT
  ifndef RUN_NAME
    ifneq ($(MAKECMDGOALS),clean)
      $(error SHOT and RUN_NAME are not defined. Please provide values using \
              'make SHOT=<integer> RUN_NAME=<string>')
    endif
  endif
endif

# --------------------------------------------------
# Compiler
# --------------------------------------------------
CC = gcc

# --------------------------------------------------
# Build mode
# --------------------------------------------------
ifndef DEBUG
  DEBUG = 0
endif

# --------------------------------------------------
# Timing (off by default)
# --------------------------------------------------
ifndef TIMING
  TIMING = 0
endif

$(info DEBUG is set to $(DEBUG))
$(info TIMING is set to $(TIMING))

# --------------------------------------------------
# Common flags (good for perf)
# --------------------------------------------------
COMMON_CFLAGS = \
  -fPIC \
  -Wall \
  -Warray-bounds \
  -std=gnu11 \
  -g \
  -fno-omit-frame-pointer \
  -pthread

COMMON_LDFLAGS = -shared -pthread

# --------------------------------------------------
# Mode-specific flags
# --------------------------------------------------
ifeq ($(DEBUG), 1)
  # DEBUG + gprof (instrumented, slow)
  CFLAGS  = $(COMMON_CFLAGS) -O0 -g -pg -DDEBUG
  LDFLAGS = -shared -g -pg
else
  # RELEASE + perf (real performance)
  CFLAGS  = $(COMMON_CFLAGS) -O3 -march=native \
            -g -fno-omit-frame-pointer
  LDFLAGS = -shared -g
endif

# --------------------------------------------------
# Optional runtime timing
# --------------------------------------------------
ifeq ($(TIMING), 1)
  CFLAGS += -DENABLE_RT_TIMING
endif

# --------------------------------------------------
# Objects and libraries
# --------------------------------------------------
OBJS = \
  constants.o \
  solve_tria.o \
  gradient.o \
  poisson_solver.o \
  find_x_point.o \
  rtgsfit.o \
  find_plasma.o

LIB     = ../lib
TEST_DIR = ../tests
SHARED  = $(patsubst %.o,$(LIB)/lib%.so,$(OBJS))

# --------------------------------------------------
# Default target
# --------------------------------------------------
all: constants.c $(SHARED)

# --------------------------------------------------
# Shared libraries
# --------------------------------------------------
$(LIB)/libfind_x_point.so: find_x_point.o gradient.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm

$(LIB)/libfind_plasma.so: find_plasma.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/libpoisson_solver.so: poisson_solver.o gradient.o solve_tria.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lopenblas

$(LIB)/libsolve_tria.so: solve_tria.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 

$(LIB)/libgradient.so: gradient.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 

$(LIB)/libconstants.so: constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 

$(LIB)/librtgsfit.so: \
  rtgsfit.o gradient.o constants.o poisson_solver.o find_x_point.o solve_tria.o find_plasma.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -lopenblas

# --------------------------------------------------
# Object compilation
# --------------------------------------------------
$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) -I/usr/include/openblas -c $< -o $@

# --------------------------------------------------
# Generated source
# --------------------------------------------------
constants.c: constants.h
	python ../utility/mdsplus_const_to_file.py SHOT=$(SHOT) RUN_NAME=$(RUN_NAME)

# --------------------------------------------------
# Tests
# --------------------------------------------------
test_find_plasma: $(TEST_DIR)/test_find_plasma.c find_plasma.o constants.o
	$(CC) $(CFLAGS) -o test_find_plasma \
        $(TEST_DIR)/test_find_plasma.c find_plasma.o constants.o -lm

# --------------------------------------------------
# Cleanup
# --------------------------------------------------
clean:
	rm -f $(SHARED) $(OBJS) constants.c *.pyc test_*