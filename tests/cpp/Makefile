# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CXX ?= c++
OPT ?= -O2 -g
# external/io.h (ext_printf) pulls in R's print headers
R_CPPFLAGS := $(shell R CMD config --cppflags)
CPPFLAGS += -I../../src/include -I../../src $(R_CPPFLAGS)
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
           test_state.cpp \
           test_ensemble.cpp test_fuzz.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS := $(SOURCES:.cpp=.d)

all: test_bartcore

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

.PHONY: all run clean
