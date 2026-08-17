# Requires an in-place package build first (R CMD INSTALL . leaves src/*.a).

CXX ?= c++
OPT ?= -O2 -g
# external/io.h (ext_printf) pulls in R's print headers
R_CPPFLAGS := $(shell R CMD config --cppflags)
CPPFLAGS += -I../../src/include -I../../src $(R_CPPFLAGS)
CXXFLAGS += -std=gnu++20 -Wall -Wextra -MMD -MP

LIBS := ../../src/misc.a ../../src/external.a ../../src/rc.a
# external.a's rng object references the R API even for non-native
# algorithms; satisfy it with libR (never called at runtime here).
R_LDFLAGS := $(shell R CMD config --ldflags)

# One translation unit per suite plus the shared assert/fixture helpers and
# main()'s registration; -MMD -MP below tracks each one's engine-header
# dependencies so only the affected objects rebuild on a header touch.
SOURCES := main.cpp assert.cpp common.cpp test_data.cpp test_tree.cpp \
           test_scan.cpp test_grow.cpp test_moves.cpp test_interaction.cpp \
           test_model.cpp test_sampler.cpp test_shape.cpp test_state.cpp \
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

.PHONY: all run clean
