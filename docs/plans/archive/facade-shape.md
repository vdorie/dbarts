# facade-shape

status: LANDED 40082c7 2026-08-06, all CI green incl. sanitizers
agent: opus
rng: neutral
budget: ~800 lines (facade.hpp, ~143 call sites, new tests/cpp TU)

## Goal

SamplerBase's nullary count/capability accessors collapse into one
SamplerShape POD returned by a single virtual, filled on demand: a
capability costs one field + one fill line, not declaration +
forward + override. Bitwise-neutral.

## Context

- SamplerBase (facade.hpp): every nullary count/capability accessor
  collapses, hasVarianceForest through savedTreeCapacity, 21 fields.
  Excluded: numVarianceTrees, a dead virtual (the bridge sizes from
  options.numVarianceTrees), deleted outright; currentSampleNum and
  fitScale (state, not shape); rng, data, latents, sigma,
  sumOfSquaredResiduals (live state); parameterized accessors.
- Callers (census 2026-08-05): 116 reads in R_interface_bartcore.cpp
  plus 27 in C_interface.cpp (the dbarts.h implementation, what
  stan4bart/bartCause exercise); a few tests/cpp SamplerBase sites.
  benchmarks/kernels/linear_leaf.cpp compiles unchanged (stale-binary
  caveat applies). Surfaced by the 2026-08-05 facade review.

## Decision

Build the collapse, or keep the per-feature triple? Recommend build:
the surface grows with every arc, each member a three-site edit, and
a swapped field is caught by the shape-vs-impl assertion instead of a
silently mis-sized buffer. Critique record (2026-08-05, refuting): NO
BLOCKER - no mutate-then-read entry point, no hot path (fill 5.74 ns),
prototype passed tests/cpp, tinytest, and the equivalence trio
bitwise. Six amendments folded in here.

## Constraints

- FILLED ON DEMAND: shape() returns by value, assembled at call time
  in SamplerFacade from impl_. No cached SamplerShape member anywhere
  - a forgotten refresh after a future mutator is a silent bug.
- SamplerShape is trivially destructible, enforced by a static_assert
  beside the definition: bridge entry points Rf_error/longjmp past
  destructors, so an owning field would leak on every error path.
- The typed layer keeps its individual methods; Sampler<L> gains a
  constexpr usesFunctionLeaves (= L::hasFunctionParams) so the test
  oracle is uniform across all 21 fields.
- X-macro forwarding DECLINED (recorded in TODO); leafCovariateColumns
  stays a borrowed pointer, today's lifetime. dbarts.h is untouched
  (no bartcore type in it, only its impl TU is swept): no ABI
  checklist.
- Caller discipline: fetch shape() at the top of each entry point;
  shared refusal helpers refill internally (~6 ns).
- Struct-layout changes corrupt stale .o exactly like vtable changes
  (verified: a stale C_interface.o segfaults): R CMD INSTALL
  --preclean and tests/cpp make clean are mandatory.

## Steps

1. Define SamplerShape + the static_assert in facade.hpp (accessor
   doc comments migrate to the fields); add virtual shape() const,
   implemented once in SamplerFacade; add usesFunctionLeaves to
   Sampler<L>.
2. Delete the 21 virtuals, their forwards, and dead numVarianceTrees;
   the compiler enumerates every caller.
3. Sweep R_interface_bartcore.cpp and C_interface.cpp (and any TU the
   compiler flags): one fetch per entry point; helpers refill.
4. New tests/cpp/test_shape.cpp asserting every shape field against
   the typed impl via facade.impl(), over construction paths held BY
   CONCRETE FACADE TYPE (constant gaussian, BCF, multinomial at
   minimum; factory dispatch out of scope). Register it: Makefile
   SOURCES, common.hpp, main.cpp.

## Verification

- cd tests/cpp && make clean && make && ./test_bartcore: all pass.
  Local ASAN build (detect_container_overflow=0) clean; CI sanitizer
  watched to green before landing more on top.
- R CMD INSTALL --preclean .; full tinytest (3509, test-capi.R
  included); equivalence trio bitwise via the dedicated harnesses.

## Landing

- 40082c7 (609+/245-): SamplerShape (21 fields) + shape() in
  facade.hpp; the 21 virtuals, their forwards, and dead
  numVarianceTrees deleted; 116 sites swept in
  R_interface_bartcore.cpp, 27 in C_interface.cpp; registered
  tests/cpp/test_shape.cpp oracles all 21 fields against a new const
  impl() overload over SEVEN construction paths (constant gaussian,
  linear, GP, ordinal, heteroscedastic, BCF, multinomial - beyond the
  plan's minimum so every non-default field value is exercised).
- Gates independently re-run at review on top of the implementer's
  battery: tests/cpp clean-build 180 ok incl. the shape suite; local
  ASAN/UBSAN run clean; INSTALL --preclean zero compiler warnings;
  tinytest 3509/0; equivalence trio bitwise; all six CI legs green
  incl. sanitizers.
- Notes: the joint per-observation entry fetches one shape per
  sampler (inherent to its multi-sampler loop, cold); the ASAN -O1
  leg caught a non-inline static const in the new test that -O2 hid -
  keep the ASAN build in the tests/cpp gate list.
