SHARED := $(wildcard ../src/*.so)
CTYPES = $(patsubst ../src/%,c_%, $(patsubst %.so,%.py,$(SHARED)))

all: $(CTYPES)

$(CTYPES): c_%.py : ../src/%.so ../src/%.h
	ctypesgen -o  $@ -l $^

clean:
	rm c_*.py
