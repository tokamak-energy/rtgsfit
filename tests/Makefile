SHARED := $(wildcard ../lib/*.so)
CTYPES = $(patsubst ../lib/%,c_%, $(patsubst %.so,%.py,$(SHARED)))

all: $(CTYPES)

$(CTYPES): c_lib%.py : ../lib/lib%.so ../src/%.h
	ctypesgen -o  $@ -l $^

clean:
	rm c_*.py
