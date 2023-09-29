# Check if FILENAME is defined as a command line argument
ifndef DATAFILE
    ifneq ($(MAKECMDGOALS),clean)
    $(error DATAFILE is not defined. Please provide a value using 'make DATAFILE=<your_file_name>')
    endif
endif

CC = gcc
DEBUG = 1

ifeq ($(DEBUG), 1)
	CFLAGS = -fPIC -Wall -Warray-bounds -pg -std=gnu11
	LDFLAGS = -pg -shared
else
	CFLAGS = -fPIC -Wall -Warray-bounds -std=gnu11 -march=native -Ofast -funroll-loops -frename-registers -flto
	LDFLAGS = -shared -Wl,-Ofast -pthread
endif

OBJS = constants.o solve_tria.o gradient.o poisson_solver.o find_x_point.o rtgsfit.o
LIB = ../lib
SHARED = $(patsubst %.o,$(LIB)/lib%.so,$(OBJS))
all: constants.c $(SHARED)

$(LIB)/libfind_x_point.so: find_x_point.o gradient.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm

$(LIB)/libpoisson_solver.so: poisson_solver.o gradient.o solve_tria.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,-Wl,--no-undefined -o $@ $^ -lopenblas

$(LIB)/libsolve_tria.so: solve_tria.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/libgradient.so: gradient.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/libconstants.so: constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/librtgsfit.so: rtgsfit.o gradient.o constants.o poisson_solver.o find_x_point.o solve_tria.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -lopenblas

$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

rtgsfit.o: rtgsfit.c
	$(CC) $(CFLAGS) -o $@  -c $< -I/usr/include/lapacke

constants.c: constants.h
	python ../utility/const_to_file.py $(DATAFILE)

clean:
	rm -f $(SHARED) $(OBJS)  constants.c *.pyc
