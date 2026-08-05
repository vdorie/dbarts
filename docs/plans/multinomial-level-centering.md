# multinomial-level-centering

status: OPEN (memo first; found 2026-08-04 by sbc-family-tiers'
        raw-f_ik arm, docs/plans/sbc-family-tiers.md Results)
agent: opus
rng: posterior-changing IF the fix ships (multinomial channels only;
     every other family untouched)
budget: memo, then ~30-80 engine lines if the memo says fix

## Goal

A verdict on `MultinomialForestCombiner::afterCombine`'s level-centering
draw: either the exact per-forest-level full conditional replaces the
current approximation and the raw-f_ik SBC arm passes, or the
approximation is proven inert for the stationary distribution and
documented as such. Today it is neither: the draw's precision sums
invV * n over observations, its own comment concedes the exact
conditional would carry a leaf-count correction (observations sharing a
leaf are not independent prior draws, over-counting by roughly
n / #leaves), and SBC now supplies evidence - all three ranked f_ik
cells are U-shaped (chisqP 0.000) at both chain-length points and do
NOT shrink at 3x spacing while their ACF is clean by lag 9, so mixing
does not explain them. The identified softmax probabilities PASS at
both points, so user-facing output is not known to be affected - but
the tree update conditions on the level through the centered fits, so
inertness is unproven, not established.

## Context

- src/bartcore/combiner.hpp `MultinomialForestCombiner::afterCombine`
  (the draw and its own caveat comment).
- docs/plans/sbc-family-tiers.md Results (the finding, the ladder
  table, the pre-registration); docs/design/multinomial.md
  (level-centering rationale).
- benchmarks/R/multinomial-exact.R checks posterior MEANS at a tiny
  config; SBC is the whole-posterior check that caught this.

## Decision (memo first)

Derive the exact per-forest-level full conditional under the
constant-leaf prior (the level shifts every leaf of forest k equally:
the conditional precision should count leaves, not observations, plus
the likelihood term). Quantify the discrepancy at the SBC config, and
answer: does the approximation distort only the unidentified level's
distribution (then: fix cheap, or document), or does it feed back into
tree structure via the centered C_ik (then: fix, full gates)? Evidence
that ends the item without a code change: a proof of inertness for
every identified functional plus a documented-approximation note in
the design doc.

## Constraints

- Any draw change re-records multinomial-equivalence and runs
  benchmarks/R/multinomial-exact.R plus the SBC multinomial arm
  (raw f_ik must PASS at the matrix band); all other families bitwise.
- No public-surface change.

## Steps

1. Memo: exact conditional derivation, magnitude at the SBC config,
   inertness analysis; record here; VD signs off on fix vs document.
2. If fix: implement, re-run the SBC multinomial arm (f_ik and softmax
   both), re-record multinomial-equivalence, exact gate, design note.

## Verification

- Memo reviewed against docs/design/multinomial.md's centering
  derivation; if fix ships: full gate battery + the SBC arm
  `Rscript benchmarks/R/sbc.R multinom 200 150 30` with every
  functional PASS.
