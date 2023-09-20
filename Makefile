
all:
    



#CC = gcc
#DEBUG = 0
#CFLAGS = -fPIC -Wall -Warray-bounds -pg -std=gnu11
#LDFLAGS = -pg -shared


##SRCS = $(wildcard *.c)
##HEADERS = $(wildcard *.h)
##OBJS = $(patsubst %.c,%.o,$(SRCS))
##SHARED = $(patsubst %.o,%.so,$(HEADERS))

#OBJS = constants.o solve_tria.o gradient.o poisson_solver.o find_x_point.o rtgsfit.o dgelss.o
#SHARED = $(patsubst %.o,%.so,$(OBJS))
#CTYPES = $(addprefix c_, $(patsubst %.so,%.py,$(SHARED)))

#all: $(CTYPES)

#$(CTYPES): c_%.py : %.so %.h
#	ctypesgen -o  $@ -l $^

#find_x_point.so: find_x_point.o gradient.o constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm
#	
#poisson_solver.so: poisson_solver.o gradient.o solve_tria.o constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,poisson_solver.so -Wl,--no-undefined -o \
#	poisson_solver.so poisson_solver.o solve_tria.o gradient.o constants.o -lopenblas

#solve_tria.so: solve_tria.o  constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 
#	
#gradient.so: gradient.o  constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ 
#	
#constants.so: constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^

#$(OBJS): %.o : %.c
#	$(CC) $(CFLAGS) -c $< -o $@

#rtgsfit.so: rtgsfit.o gradient.o constants.o poisson_solver.o find_x_point.o solve_tria.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -lopenblas

#rtgsfit.o: rtgsfit.c
#	$(CC) $(CFLAGS) -o $@  -c $< -I/usr/include/lapacke 
#	
#dgelss.so: dgelss.o  constants.o
#	$(CC) $(LDFLAGS) -Wl,-soname,$@ -Wl,--no-undefined -o $@ $^ -lm -lopenblas

#dgelss.o: dgelss.c
#	$(CC) $(CFLAGS) -o $@  -c $< -I/usr/include/lapacke 
#	
#c_cblas.py: 
#	ctypesgen -o c_cblas.py -llibcblas.so /usr/include/cblas.h
#	
#clean:
#	rm -f $(SHARED) $(OBJS) $(CTYPES) constants.c *.pyc



#all: $(SHARED)

#c_poisson_solver.py: poisson_solver.so
#	ctypesgen -o c_poisson_solver.py -l poisson_solver.so poisson_solver.h

#$(SHARED): %.so : %.o
#	$(CC) $(LDFLAGS) $^ -o $@  
	
#poisson_solver.so: poisson_solver.o gradient.so solve_tria.so
#	$(CC) $(LDFLAGS) -Wl,-soname,poisson_solver.so -Wl,--no-undefined -o poisson_solver.so poisson_solver.o solve_tria.so gradient.so -lopenblas
