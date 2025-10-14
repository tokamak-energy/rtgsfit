# Check if FILENAME is defined as a command line argument
ifndef SHOT
	ifndef RUN_NAME
		ifneq ($(MAKECMDGOALS),clean)
		$(error SHOT and RUN_NAME are not defined. Please provide values using 'make SHOT=<integer> RUN_NAME=<string>')
		endif
	endif
endif

CC = gcc
ifndef DEBUG
	DEBUG = 0
endif
# print the value of DEBUG
$(info DEBUG is set to $(DEBUG))

ifeq ($(DEBUG), 1)
	CFLAGS = -fPIC -Wall -Warray-bounds -pg -std=gnu11 -Ofast -pthread -g3 -O0 -DDEBUG
	LDFLAGS = -pg -shared -g3
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
	$(CC) $(LDFLAGS) -Wl,-soname,-Wl,--no-undefined -o $@ $^ -L/usr/lib/x86_64-linux-gnu/openblas

$(LIB)/libsolve_tria.so: solve_tria.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/libgradient.so: gradient.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/libconstants.so: constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

$(LIB)/librtgsfit.so: rtgsfit.o gradient.o constants.o poisson_solver.o find_x_point.o solve_tria.o
	$(CC) $(LDFLAGS) $(CFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -L/usr/lib/x86_64-linux-gnu/openblas

$(OBJS): %.o : %.c
	#$(CC) $(CFLAGS) -I/usr/include/openblas -c $< -o $@
	$(CC) $(CFLAGS) -I/usr/include/x86_64-linux-gnu -c $< -o $@

rtgsfit.o: rtgsfit.c
	#$(CC) $(CFLAGS) -I/usr/include/openblas -o $@  -c $< 
	$(CC) $(CFLAGS) -I/usr/include/x86_64-linux-gnu -o $@  -c $< 

constants.c: constants.h
	python ../utility/mdsplus_const_to_file.py SHOT=$(SHOT) RUN_NAME=$(RUN_NAME)

clean:
	rm -f $(SHARED) $(OBJS)  constants.c *.pyc
