# gp-cache-test-flake

agent: opus
rng: neutral (test-only)
budget: ~40 lines

## Goal

tests/cpp's testGPLeafKernelCache passes in every build: the
"warm cache matches cold clone after member re-route" check stops
depending on allocation layout.

## Context

- Found 2026-07-07 while landing linear-leaf-reuse: the check fails
  deterministically in some builds at clean HEAD (6/6 on one main
  build, 0/N on the implementer's clean build of the same code) and
  is unrelated to the linear cache (GP path only).
- Known mechanism (docs and memory both record it): SIMD reductions
  split allocation-dependently, so two separately allocated but
  logically identical accumulations can differ in the last bit. The
  test compares a warm-cached GP fit bitwise against a freshly
  cloned cold one - different allocations, so the comparison is
  layout roulette.

## Constraints

- Keep the contract strong: the cache must still be shown to serve
  values bitwise-equal to a scan WITHIN one object, where buffers
  are shared and the comparison is sound (the linear-leaf cache test
  in test_model.cpp shows the shape: warm vs cold on the SAME leaf
  after an explicit invalidate).
- No engine changes; tests/cpp/test_model.cpp only.

## Steps

1. Rework the re-route check: compare warm lookup vs post-invalidate
   rescan on the same leaf object instead of vs a cloned instance
   (or force the clone's buffers through the same allocation
   pattern if that is simpler and still deterministic).
2. Verify the reworked test fails when the cache is deliberately
   poisoned (sanity that it still tests something).

## Verification

- tests/cpp passes repeatedly across a fresh rebuild (make clean,
  rebuild, run several times).
