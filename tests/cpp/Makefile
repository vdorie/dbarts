# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CXX ?= c++
OPT ?= -O2 -g
# external/io.h (ext_printf) pulls in R's print headers
R_CPPFLAGS := $(shell R CMD config --cppflags)
CPPFLAGS += -I../../src/include -I../../src $(R_CPPFLAGS)
CXXFLAGS += -std=gnu++20 -Wall -Wextra

LIBS := ../../src/misc.a ../../src/external.a ../../src/rc.a
# external.a's rng object references the R API even for non-native
# algorithms; satisfy it with libR (never called at runtime here).
R_LDFLAGS := $(shell R CMD config --ldflags)

all: test_bartcore

test_bartcore: test_bartcore.cpp ../../src/bartcore/*.hpp $(LIBS)
	$(CXX) $(OPT) $(CXXFLAGS) $(CPPFLAGS) -o $@ test_bartcore.cpp $(LIBS) $(R_LDFLAGS) -lpthread -lm

run: test_bartcore
	./test_bartcore

clean:
	rm -f test_bartcore

.PHONY: all run clean
