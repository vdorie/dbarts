# mutation-fuzzing

agent: opus
rng: neutral (tests only)
budget: ~400 lines in tests/cpp

## Goal

The transactional mutation surface survives randomized adversarial use:
a property test drives long random sequences of mutations and checks
the invariants after every step, so rollback and re-routing bugs
surface before users find them.

## Context

- The surface: setPredictor (full/column/forced), updatePredictor,
  setCell, per-observation sessions (commit and abandon),
  setData (resizing), setOffset/setResponse/setWeights/setSigma,
  setCutPoints, setTestPredictors, getState/setState round trips -
  the package's most state-heavy code (src/bartcore/sampler.hpp,
  data.hpp).
- Existing mutation tests are example-based; the audit found the
  surface is a faithful classic port whose failure modes are
  state-consistency ones - exactly what property testing finds.

## Constraints

- Deterministic: seeded ext_rng stream drives the op sequence; a
  failure prints the seed and the op trace for replay (no
  Date/clock anywhere).
- Invariants after every op: every tree well-formed (existing
  helpers); node observation counts sum to n; totalFits equals the
  tree-order sum of fits (to tolerance); rejected transactions leave
  a bitwise-identical getState; accepted ones keep run() functional
  (one sweep, finite draws).
- Budget-bounded: ~2000 ops across mixed configurations (families,
  categorical/sparse/MIA columns, linear leaves) within ~30s; runs in
  tests/cpp (no R churn), wired into cpp-tests CI.
- Out of scope: the R5 layer (thin; exercised by tinytest), fuzzing
  the C API marshaling.

## Steps

1. Op vocabulary + weights table; a generator producing valid-but-
   nasty inputs (columns that force rollback, quantile-infeasible
   cuts, resizing setData, degenerate constants, all-same categorical
   columns).
2. Invariant checkers reusing the component-test helpers; the
   before/after state capture for rejected ops.
3. Runner with seed/trace replay; N seeds in CI, more locally via an
   argument.
4. Triage anything it finds into separate fix items (do not fix
   in-band beyond trivial).

## Verification

- The fuzzer passes clean over >= 20 seeds locally, 3 in CI.
- Deliberately breaking a rollback path (locally, reverted) makes it
  fail with a usable trace - prove the harness can catch.

## Landing note (2026-07-07, 44cee98)

Landed: testMutationFuzzer in tests/cpp - 8 configurations x 13
weighted ops with invariants after every op (well-formed trees with
exact range nesting and occupancy, totalFits resum, rejected
transactions bitwise-stable via fingerprint, run finite, semantic
state round-trip); replayable seed + op trace on failure; main takes
a seed count (CI runs 3, default 20, clean at 100; whole binary
~2.4s). Sabotage-proofed: a removed rollback restore fails with a
usable trace. Deviations: setCell is store-level, not on the Sampler
surface, so not fuzzed; quantile-infeasible cuts exercised via
invalid categorical codes. FOUND two latent state-serialization
edges, triaged to fuzz-state-roundtrip (not fixed in-band): forced
re-cut of a near-constant column installs a non-ascending uniform
grid, and a category-mutated live mask fails the canonical-gauge
check on flatten - both make setState refuse the sampler's own state
(cleanly; nothing corrupts). Gates: full C++ suite + fuzzer 20 seeds,
diff confined to tests/cpp + cpp-tests.yaml, tinytest All ok 2468.
