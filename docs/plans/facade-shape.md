# facade-shape

agent: opus
rng: neutral
budget: ~700 lines (facade.hpp -22 virtuals/+1 POD, ~119 bridge call
  sites, one new tests/cpp TU)

## Goal

SamplerBase's nullary count/capability accessors collapse into one
SamplerShape POD returned by a single virtual, filled on demand: a
future capability costs one POD field plus one fill line instead of
declaration + facade forward + impl override. Bitwise-neutral.

## Context

- SamplerBase (facade.hpp): the 22 collapsing members are the nullary
  count/capability accessors, hasVarianceForest through
  savedTreeCapacity (numVarianceTrees, numLeafCovariates,
  leafCovariateColumns, usesFunctionLeaves, family, numGroups,
  kIsSampled, usesDart, numChains, numThreads, numForests,
  supportsResponseMutation, numReportedLocations,
  numVariableCountForests, numCutpoints, testFitsAreDefined, numTrees,
  numObservations, numPredictors, numTestObservations).
- Staying virtual: currentSampleNum and fitScale (state, not shape);
  rng, data, latents, sigma, sumOfSquaredResiduals (live state or
  computation); every parameterized accessor.
- Callers: ~119 reads in src/R_interface_bartcore.cpp (census
  2026-08-05); a few tests/cpp SamplerBase-interface sites. Surfaced
  by the 2026-08-05 facade review (TODO facade-shape).

## Decision

Build the collapse, or keep the per-feature triple? Recommend build:
the surface is 22 members and still growing (supportsResponseMutation
with FIX-B; numReportedLocations, numVariableCountForests, numCutpoints
with the family arcs), each a three-site edit today; after, each is one
field plus one fill line, and a swapped field is caught by the
tests/cpp shape-vs-impl assertion instead of by a silently mis-sized
bridge buffer. Evidence that would change it: the step 3 audit finding
pervasive mutate-then-read entry points.

## Constraints

- FILLED ON DEMAND: shape() returns by value, assembled at call time
  in SamplerFacade from impl_ methods. No cached SamplerShape member
  anywhere - a forgotten cache refresh after a future mutator is a
  silent runtime bug; the callers size buffers, they are not hot.
- The typed sampler layer keeps its individual methods; only the
  type-erased virtuals collapse.
- X-macro forwarding generation DECLINED (recorded in TODO).
  leafCovariateColumns stays a borrowed pointer field with today's
  lifetime. No dbarts.h touch (src-internal only): no ABI checklist.
- Caller discipline: fetch shape() at the top of each bridge entry
  point; re-fetch after any in-function mutator (step 3 lists these).
- Struct-layout changes corrupt stale .o exactly like vtable changes:
  R CMD INSTALL --preclean and tests/cpp make clean are mandatory.

## Steps

1. Define SamplerShape in facade.hpp (accessor doc comments migrate to
   the fields); add virtual SamplerShape shape() const to SamplerBase,
   implemented once in SamplerFacade from impl_.
2. Delete the 22 virtuals and their SamplerFacade forwards; the
   compiler now enumerates every caller.
3. Sweep src/R_interface_bartcore.cpp (and any TU the compiler flags,
   tests/cpp included): one shape fetch per entry point; audit
   mutate-then-read.
4. New tests/cpp/test_shape.cpp: for each construction path the
   component tests already exercise, assert every shape field equals
   the typed impl's accessor.

## Verification

- cd tests/cpp && make clean && make && ./test_bartcore: all pass.
  Local ASAN build (detect_container_overflow=0) clean; CI sanitizer
  watched to green before landing more on top.
- R CMD INSTALL --preclean .; full tinytest (3484); equivalence trio
  bitwise via the dedicated harnesses (equivalence-7903855,
  bcf-equivalence-99205ee, multinomial-equivalence-ec2a3d0).
