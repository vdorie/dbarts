# Requires an in-place package build first (R CMD INSTALL . leaves src/misc.a).

CC ?= cc
OPT ?= -O2 -g
CPPFLAGS += -I../../src/include

MISC_LIB := ../../src/misc.a

all: bench

bench: bench.c partition_u8.c partition_u8.h $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ bench.c partition_u8.c $(MISC_LIB) -lpthread -lm

sparse: sparse.c $(MISC_LIB)
	$(CC) $(OPT) $(CPPFLAGS) -o $@ sparse.c $(MISC_LIB) -lpthread -lm

$(MISC_LIB):
	@echo "error: $(MISC_LIB) not found; run 'R CMD INSTALL .' from the package root first" >&2; exit 1

clean:
	rm -f bench sparse

.PHONY: all clean
