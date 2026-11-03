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
# inherit gaussian's behavior. The package build only warns (CRAN
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
           test_ensemble.cpp test_monotone.cpp test_fuzz.cpp test_nmath.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS := $(SOURCES:.cpp=.d)

# The shipped header's prototype view must compile as plain C with NO R include
# path (a consumer that includes it before any R header, and the plain-C
# specification stage it is being kept ready for), and as C++ for the LinkingTo
# consumers that are C++. Both are syntax-only: nothing is linked, and nothing
# here may add -I$(R_CPPFLAGS), which is the whole point of the target.
HEADER := ../../inst/include/dbarts/dbarts.h
CC ?= cc

# --- The standalone build (DBARTS_HOST_STANDALONE): no R header, no libR ---
#
# The recipe a binding for another language copies: the support libraries and
# the copy of R's math (standalone/nmath), configured by the checked-in headers
# in standalone/config since configure reads R CMD config, compiled into
# build-standalone/ and never into src/, where a later R CMD INSTALL would link
# a stray object into dbarts.so. The C compiler and flags default to R's own:
# only those make the copy R's code bit for bit (contraction of a * b + c, for
# one) and the kernels the R build's; a build without R sets RCC and RCFLAGS.
RCC ?= $(or $(shell R CMD config CC 2>/dev/null),cc)
RCFLAGS ?= $(shell R CMD config CFLAGS 2>/dev/null)
AR ?= ar

STANDALONE := ../../standalone
SA_BUILD := build-standalone
# standalone/config ahead of src/include, so its misc/types.h wins; a source's
# own "config.h" resolves next to the source first, so in a tree configure has
# run in, misc and external read configure's copy of the same values
SA_CPPFLAGS := -DDBARTS_HOST_STANDALONE -DNDEBUG -I$(STANDALONE)/config -I../../src/include

NMATH_DIR := $(STANDALONE)/nmath
NMATH_SOURCES := $(wildcard $(NMATH_DIR)/*.c)
NMATH_OBJECTS := $(patsubst $(NMATH_DIR)/%.c,$(SA_BUILD)/nmath/%.o,$(NMATH_SOURCES))
NMATH_LIB := $(SA_BUILD)/nmath/libnmath.a
NMATH_CPPFLAGS := -DMATHLIB_STANDALONE -include errno.h -include $(NMATH_DIR)/prefix.h \
                  -I$(NMATH_DIR)/include

all: header-compiles test_bartcore r-math-symbols standalone-math

# -Wno-unused-function: the header's one static inline helper
# (dbarts_dense_predictor_source) has no caller when the header is compiled
# alone, which is the only way this target compiles it.
header-compiles:
	$(CC) -std=c99 -pedantic -Wall -Wno-unused-function -x c -fsyntax-only $(HEADER)
	$(CXX) -std=c++11 -Wno-unused-function -x c++ -fsyntax-only $(HEADER)

# The copy of R's math rides along so the nmath suite can compare it with
# libR's on its grid; its prefix keeps the two apart.
test_bartcore: $(OBJECTS) $(LIBS) $(NMATH_LIB)
	$(CXX) $(OPT) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(OBJECTS) $(LIBS) $(NMATH_LIB) $(R_LDFLAGS) -lpthread -lm

%.o: %.cpp
	$(CXX) $(OPT) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

-include $(DEPS)

$(SA_BUILD)/nmath $(SA_BUILD)/misc $(SA_BUILD)/external:
	mkdir -p $@

$(SA_BUILD)/nmath/%.o: $(NMATH_DIR)/%.c $(NMATH_DIR)/prefix.h $(NMATH_DIR)/nmath.h \
                       $(NMATH_DIR)/dpq.h | $(SA_BUILD)/nmath
	$(RCC) $(RCFLAGS) $(NMATH_CPPFLAGS) -c -o $@ $<

# Every symbol the copy defines must carry the prefix, or a program that also
# links R or libRmath would see two definitions.
$(NMATH_LIB): $(NMATH_OBJECTS)
	@unprefixed=$$(nm -g $(NMATH_OBJECTS) | awk 'NF == 3 && $$2 != "U" { print $$3 }' | \
	  grep -v -E '^_?dbarts_nmath_'); \
	if [ -n "$$unprefixed" ]; then \
	  echo "error: the copy of R's math defines unprefixed symbols:" $$unprefixed >&2; exit 1; \
	fi
	rm -f $@
	$(AR) -crs $@ $(NMATH_OBJECTS)

# misc through its own Makefile, read from here with vpath at src/misc, and
# the per-file instruction-set flags configure passes on this architecture
SA_ARCH := $(shell uname -m)
ifeq ($(SA_ARCH),x86_64)
  SA_SIMD_FLAGS := SSE2_FLAG=-msse2 SSE4_1_FLAG=-msse4.1 AVX_FLAG=-mavx AVX2_FLAG=-mavx2
endif
ifneq (,$(filter arm64 aarch64,$(SA_ARCH)))
  SA_SIMD_FLAGS := NEON_FLAG="-mfloat-abi=softfp -mfpu=neon"
endif
SA_MISC_LIB := $(SA_BUILD)/misc/libmisc.a

$(SA_MISC_LIB): FORCE | $(SA_BUILD)/misc
	$(MAKE) -C $(SA_BUILD)/misc -f $(CURDIR)/standalone-misc.mk \
	  -f $(abspath ../../src/misc/Makefile) library \
	  MISC_SRC=$(abspath ../../src/misc) INCLUDE_DIR=$(abspath ../../src/include) \
	  MISC_LIB=libmisc.a CC="$(RCC)" ALL_CFLAGS="$(RCFLAGS)" AR="$(AR)" ARFLAGS=-crs \
	  ALL_CPPFLAGS="-DNDEBUG -I$(abspath $(STANDALONE)/config/misc) \
	    -I$(abspath $(STANDALONE)/config) -I$(abspath ../../src/include)" \
	  $(SA_SIMD_FLAGS)

# external without its R-only parts: io_standalone.c in place of io.c
SA_EXTERNAL_SOURCES := ../../src/external/randomBase.c ../../src/external/randomNorm.c \
                       ../../src/external/random.c $(STANDALONE)/io_standalone.c
SA_EXTERNAL_OBJECTS := $(addprefix $(SA_BUILD)/external/,$(notdir $(SA_EXTERNAL_SOURCES:.c=.o)))
SA_EXTERNAL_LIB := $(SA_BUILD)/external/libexternal.a
SA_EXTERNAL_CPPFLAGS := $(SA_CPPFLAGS) -I$(STANDALONE)/config/external

$(SA_BUILD)/external/%.o: ../../src/external/%.c | $(SA_BUILD)/external
	$(RCC) $(RCFLAGS) $(SA_EXTERNAL_CPPFLAGS) -MMD -MP -c -o $@ $<

$(SA_BUILD)/external/%.o: $(STANDALONE)/%.c | $(SA_BUILD)/external
	$(RCC) $(RCFLAGS) $(SA_EXTERNAL_CPPFLAGS) -MMD -MP -c -o $@ $<

$(SA_EXTERNAL_LIB): $(SA_EXTERNAL_OBJECTS)
	rm -f $@
	$(AR) -crs $@ $(SA_EXTERNAL_OBJECTS)

-include $(SA_EXTERNAL_OBJECTS:.o=.d) $(SA_BUILD)/standalone_math.d

# compiled and linked in two steps, so no .dSYM directory is left on macOS
$(SA_BUILD)/standalone_math.o: standalone_math.c | $(SA_BUILD)/external
	$(RCC) $(RCFLAGS) $(SA_CPPFLAGS) -MMD -MP -c -o $@ $<

$(SA_BUILD)/standalone_math: $(SA_BUILD)/standalone_math.o $(SA_EXTERNAL_LIB) $(SA_MISC_LIB) \
                             $(NMATH_LIB)
	$(RCC) $(RCFLAGS) -o $@ $(SA_BUILD)/standalone_math.o $(SA_EXTERNAL_LIB) $(SA_MISC_LIB) \
	  $(NMATH_LIB) -lpthread -lm

# Builds, runs, and checks the binary references nothing of R.
standalone-math: $(SA_BUILD)/standalone_math
	./$(SA_BUILD)/standalone_math
	@rsyms=$$(nm -u $(SA_BUILD)/standalone_math | \
	  grep -E '^ *U? *_?(Rf_|R_|Rprintf|REprintf|unif_rand|norm_rand)'); \
	if [ -n "$$rsyms" ]; then \
	  echo "error: the standalone binary references R:" $$rsyms >&2; exit 1; \
	fi
	@echo "standalone-math: no R symbol in $(SA_BUILD)/standalone_math"

# The R package's dbarts.so (left in src/ by the in-place build) calls R's own
# math, never the copy: each of the six names the engine and the support
# library use is R's, and no dbarts_nmath_ symbol is in the library. Draws
# cannot tell the two apart where the copy is exact, so this reads the symbols.
R_MATH_NAMES := Rf_pnorm5 Rf_qnorm5 Rf_dnorm4 Rf_dt Rf_qt Rf_qchisq
r-math-symbols:
	@so=../../src/dbarts.so; \
	if [ ! -f $$so ]; then echo "error: $$so not found; run 'R CMD INSTALL .' first" >&2; exit 1; fi; \
	undefined=$$(nm -u $$so); \
	for name in $(R_MATH_NAMES); do \
	  echo "$$undefined" | grep -q -E "^ *U? *_?$$name\$$" || \
	    { echo "error: $$so does not call R's $$name" >&2; exit 1; }; \
	done; \
	if nm $$so | grep -q dbarts_nmath_; then \
	  echo "error: $$so carries the standalone copy of R's math" >&2; exit 1; \
	fi; \
	echo "r-math-symbols: $$so calls R's $(R_MATH_NAMES)"

FORCE:

run: test_bartcore
	./test_bartcore

clean:
	rm -f test_bartcore $(OBJECTS) $(DEPS)
	rm -rf test_bartcore.dSYM
	rm -f build-standalone/nmath/*.o build-standalone/nmath/libnmath.a
	rm -f build-standalone/misc/*.o build-standalone/misc/libmisc.a
	rm -f build-standalone/external/*.o build-standalone/external/*.d
	rm -f build-standalone/external/libexternal.a
	rm -f build-standalone/standalone_math build-standalone/standalone_math.o
	rm -f build-standalone/standalone_math.d
	-rmdir build-standalone/nmath build-standalone/misc build-standalone/external
	-rmdir build-standalone

.PHONY: all header-compiles run clean standalone-math r-math-symbols FORCE
