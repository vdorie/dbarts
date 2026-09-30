# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CC ?= cc
CXX ?= c++
OPT ?= -O2 -g
CPPFLAGS += -I../../src/include

MISC_LIB := ../../src/misc.a

all: bench

bench: bench.c partition_u8.c partition_u8.h $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ bench.c partition_u8.c $(MISC_LIB) -lpthread -lm

grow_from_root: grow_from_root.c $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ grow_from_root.c $(MISC_LIB) -lpthread -lm

sparse: sparse.c $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ sparse.c $(MISC_LIB) -lpthread -lm

monotone_count: monotone_count.cpp
	$(CXX) $(OPT) -std=c++17 -o $@ monotone_count.cpp

monotone_ratio: monotone_ratio.cpp
	$(CXX) $(OPT) -std=c++17 -o $@ monotone_ratio.cpp

$(MISC_LIB):
	@echo "error: $(MISC_LIB) not found; run 'R CMD INSTALL .' from the package root first" >&2; exit 1

clean:
	rm -f bench sparse grow_from_root monotone_count monotone_ratio

.PHONY: all clean
