# tests-cpp-split

agent: sonnet
rng: neutral (test infrastructure)
budget: file moves + ~60 lines of Makefile

## Goal

tests/cpp compiles incrementally and runs subsets: per-area translation
units with dependency tracking and a name filter, replacing the
6075-line single-TU monolith when its compile time starts hurting.

## Context

- tests/cpp/test_bartcore.cpp: 78 test functions, one TU, hand-rolled
  check/checkNear; Makefile already tracks ../../src/bartcore/*.hpp as
  prerequisites but rebuilds everything on any header touch.
- Trigger condition: this item is deliberately dormant until edit-
  compile-test latency on the engine headers becomes an irritant
  (VD's call); splitting early buys little.

## Constraints

- Zero-dependency stays (no test framework); the micro-assert helpers
  move to a shared header.
- CI (cpp-tests.yaml) keeps a single binary entry point.
- Out of scope: rshim.cpp (equivalence harness shim, parked here;
  leave it).

## Steps

1. Split by area: test_{data,tree,moves,model,sampler,state}.cpp + a
   main.cpp registering suites; shared assert header.
2. Makefile: -MMD -MP dependency files; a filter argument
   (./test_bartcore sampler) selecting suites by prefix.
3. cpp-tests.yaml unchanged except the build producing the same binary
   name.

## Verification

- All 78 tests still run and pass (count them in the runner output).
- Touching one engine header recompiles only dependent TUs.
