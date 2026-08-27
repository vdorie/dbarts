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

## Status (2026-07-07)

Landed (ee4d58c). The re-route check now runs on a single
GPGaussianLeaf: warm the kernel cache over members [0, 50), re-route
the node onto the disjoint [100, 150), and compare the memcmp-driven
rebuild against a rescan taken after regatherTrainingCovariates drops
the cache - shared buffers make the comparison bitwise-deterministic.
Root cause confirmed: the GP leaf math is fully scalar; the flake was
the two separately allocated samplers' misc SIMD fit reductions
splitting differently, exactly the known layout-roulette mechanism.
Poison check ran: a 1e-6 perturbation on the cache serve path made the
reworked test fail, then was reverted (step 2 satisfied). The
sampler-level designated-mutation check (warm vs cold clone) was left
as-is per scope; if it ever flakes by the same mechanism, rework it
the same way. Gates: baseline reproduced the failure 5/5 at HEAD;
after the fix, the model suite passed 5x plus 3x on a clean rebuild
(implementer) and 3x plus a full run on the main tree (reviewer).
tinytest/equivalence skipped: the diff is tests/cpp-only, outside the
installed package. Review fix: the function-head comment still
described the deleted cold-clone cycle; rewritten in the landing.
