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

## Landing note (2026-07-07, fb7bceb)

Landed: the one is.finite-only mutation assertion in each of
test-linear-leaves.R and test-gp-leaves.R replaced by a structural
getTrees well-formedness check plus an RMSE comparison of the
mutated-continued fit against a fresh fit on the same mutated data;
bounds calibrated over 60/30 seeds with mean+4sd headroom, basis
recorded in test comments. Sabotage demonstration: a misaligned
(shuffled) predictor fails the new assertion cleanly - a stale-scale
sabotage does NOT, because linear leaves absorb a scale mutation into
beta, which is why the RMSE form was needed. Gates: suite twice at
2470/0, diff confined to inst/tinytest, lint zero. ~67 lines.
