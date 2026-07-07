# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CC ?= cc
CXX ?= c++
OPT ?= -O2 -g
CPPFLAGS += -I../../src/include

MISC_LIB := ../../src/misc.a

# linear_leaf links the full engine (bartcore headers pull in the R print API
# through external/io.h), so it needs R's flags and all three static libraries.
R_CPPFLAGS := $(shell R CMD config --cppflags)
R_LDFLAGS := $(shell R CMD config --ldflags)
ENGINE_LIBS := ../../src/misc.a ../../src/external.a ../../src/rc.a

all: bench

bench: bench.c partition_u8.c partition_u8.h $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ bench.c partition_u8.c $(MISC_LIB) -lpthread -lm

sparse: sparse.c $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ sparse.c $(MISC_LIB) -lpthread -lm

linear_leaf: linear_leaf.cpp $(ENGINE_LIBS)
	$(CXX) $(OPT) -std=gnu++20 -I../../src -I../../src/include $(R_CPPFLAGS) -o $@ linear_leaf.cpp $(ENGINE_LIBS) $(R_LDFLAGS) -lpthread -lm

$(MISC_LIB):
	@echo "error: $(MISC_LIB) not found; run 'R CMD INSTALL .' from the package root first" >&2; exit 1

clean:
	rm -f bench sparse linear_leaf

.PHONY: all clean
