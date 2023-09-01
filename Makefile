CC = gcc
DEBUG = 0
CFLAGS = -fPIC -Wall -Warray-bounds -pg
LDFLAGS = -pg -shared

SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c,%.o,$(SRCS))
SHARED = $(patsubst %.o,%.so,$(OBJS))
CTYPES = $(addprefix c_, $(patsubst %.so,%.py,$(SHARED)))

all: $(CTYPES)

$(CTYPES): c_%.py : %.so %.h
	ctypesgen -o  $@ -l $^

$(SHARED): %.so : %.o
	$(CC) $(LDFLAGS) $^ -o $@  
	
$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(SHARED) $(OBJS) $(CTYPES)

   
