# post-mutation-assertions

agent: sonnet
rng: neutral (tests only)
budget: ~150 lines of tests

## Goal

The older linear/gp mutation tests assert something falsifiable after
mutating: a mutated sampler's continued run agrees statistically with a
from-scratch fit on the post-mutation data, instead of is.finite.

## Context

- Near-tautologies live in the older linear/gp mutation tests
  (TODO note, 2026-07-06); mutate-then-serialize covers the bitwise
  side already.
- The comparison must be statistical, not identical: setPredictor
  re-quantizes differently from a fresh build - the MIA setPredictor
  test has the working pattern to copy.

## Constraints

- Statistical bounds must be robust across seeds (the flat-hyperprior
  lesson from the flip: knife-edge bounds get replaced, not tightened).
  Calibrate any bound by running >= 20 seeds locally and leaving 4-sigma
  headroom.
- Runtime budget: keep added test time under ~10s per file.
- Out of scope: engine changes; new mutation surface.

## Steps

1. Inventory is.finite-only assertions in the linear/gp mutation test
   files; list them in the review report.
2. For each: mutate, run, compare posterior mean predictions against a
   fresh sampler built on the mutated data (same seed protocol as the
   MIA test), with a calibrated tolerance.
3. Keep one cheap structural assertion per mutation as a smoke layer
   (tree well-formedness via the existing helpers, not is.finite).

## Verification

- New assertions pass across 20 local seeds (documented in the test
  comment only if the bound needs justification).
- Full tinytest runtime increase measured and reported.
