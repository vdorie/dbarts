# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CXX ?= c++
OPT ?= -O2 -g
# external/io.h (ext_printf) pulls in R's print headers
R_CPPFLAGS := $(shell R CMD config --cppflags)
# -I../../inst/include: the shipped header, which the bridge's shared
# internals include for the flat C API's per-draw struct
CPPFLAGS += -I../../src/include -I../../src -I../../inst/include $(R_CPPFLAGS)
# Scaffolding builds append their own defines here rather than through
# CPPFLAGS, which a command-line assignment would replace wholesale:
#   make EXTRA_CPPFLAGS=-DBARTCORE_MOVE_CENSUS
CPPFLAGS += $(EXTRA_CPPFLAGS)
# -Werror=switch is a gate, not a style: the engine's switches over
# ResponseFamily carry no default arm, so a new family that is not
# classified at every one of them must break this build rather than
# inherit gaussian's behaviour. The package build only warns (CRAN
# forbids -Werror there), which is why the hard failure lives here.
CXXFLAGS += -std=gnu++20 -Wall -Wextra -Werror=switch \
            -Werror=return-type -MMD -MP

LIBS := ../../src/misc.a ../../src/external.a ../../src/rc.a
# external.a's rng object references the R API even for non-native
# algorithms; satisfy it with libR (never called at runtime here).
R_LDFLAGS := $(shell R CMD config --ldflags)

# One translation unit per suite plus the shared assert/fixture helpers and
# main()'s registration; -MMD -MP below tracks each one's engine-header
# dependencies so only the affected objects rebuild on a header touch.
SOURCES := main.cpp assert.cpp common.cpp test_data.cpp test_tree.cpp \
           test_scan.cpp test_grow.cpp test_moves.cpp test_interaction.cpp \
           test_model.cpp test_sampler.cpp test_shape.cpp test_facade.cpp \
           test_state.cpp test_capi.cpp \
           test_ensemble.cpp test_fuzz.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS := $(SOURCES:.cpp=.d)

# The shipped header's prototype view must compile as plain C with NO R include
# path (a consumer that includes it before any R header, and the plain-C
# specification stage it is being kept ready for), and as C++ for the LinkingTo
# consumers that are C++. Both are syntax-only: nothing is linked, and nothing
# here may add -I$(R_CPPFLAGS), which is the whole point of the target.
HEADER := ../../inst/include/dbarts/dbarts.h
CC ?= cc

all: header-compiles test_bartcore

# -Wno-unused-function: the header's one static inline helper
# (dbarts_dense_predictor_source) has no caller when the header is compiled
# alone, which is the only way this target compiles it.
header-compiles:
	$(CC) -std=c99 -pedantic -Wall -Wno-unused-function -x c -fsyntax-only $(HEADER)
	$(CXX) -std=c++11 -Wno-unused-function -x c++ -fsyntax-only $(HEADER)

test_bartcore: $(OBJECTS) $(LIBS)
	$(CXX) $(OPT) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(OBJECTS) $(LIBS) $(R_LDFLAGS) -lpthread -lm

%.o: %.cpp
	$(CXX) $(OPT) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

-include $(DEPS)

run: test_bartcore
	./test_bartcore

clean:
	rm -f test_bartcore $(OBJECTS) $(DEPS)
	rm -rf test_bartcore.dSYM

.PHONY: all header-compiles run clean
