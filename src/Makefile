# Check if FILENAME is defined as a command line argument
#ifndef DATAFILE
#$(error DATAFILE is not defined. Please provide a value using 'make DATAFILE=<your_file_name>')
#endif
CC = gcc
DEBUG = 0
CFLAGS = -fPIC -Wall -Warray-bounds -pg -std=gnu11
LDFLAGS = -pg -shared
OBJS = constants.o solve_tria.o gradient.o poisson_solver.o find_x_point.o rtgsfit.o 
SHARED = $(patsubst %.o,%.so,$(OBJS))

all: constants.c $(SHARED)

find_x_point.so: find_x_point.o gradient.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm
	
poisson_solver.so: poisson_solver.o gradient.o solve_tria.o constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,poisson_solver.so -Wl,--no-undefined -o \
	poisson_solver.so poisson_solver.o solve_tria.o gradient.o constants.o -lopenblas

solve_tria.so: solve_tria.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 
	
gradient.so: gradient.o  constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 
	
constants.so: constants.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

rtgsfit.so: rtgsfit.o gradient.o constants.o poisson_solver.o find_x_point.o solve_tria.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -lopenblas
	
$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@	

rtgsfit.o: rtgsfit.c
	$(CC) $(CFLAGS) -o $@  -c $< -I/usr/include/lapacke 
	
constants.c: constants_template.c $(DATAFILE)
	python ../utility/const_to_file.py $(DATAFILE)
	
clean:
	rm -f $(SHARED) $(OBJS)  constants.c *.pyc

