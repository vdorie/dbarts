# multinomial-level-centering

status: LANDED ec2a3d0 2026-08-05 (see Landing below; fix implemented
        under VD's proceed-at-discretion grant per the Memo's unanimous
        verdict - two blind derivations + adversarial critique with a
        counterfactual SBC run)
        (found 2026-08-04 by sbc-family-tiers' raw-f_ik arm,
        docs/plans/sbc-family-tiers.md Results)
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

## Memo (2026-08-05)

Process: two independent blind derivations (opus), reconciled by the
orchestrator; then an adversarial critique (opus) instructed to
refute, which instead validated the fix end to end on a patched
scratch build. The repo was untouched throughout.

THE DEFECT. The prior lives on LEAF VALUES (per-leaf sd
s = nodeScale/(k sqrt(m)), 0.272 at the SBC config), not on total
fits (tau = 1.924). The implemented draw treats the n*K total fits as
independent N(0, tau^2) prior draws. For the tree-0 absorption the
code actually performs, the exact conditional is
prec = sum_k L_k0/s_k^2, mean = -(sum_k S_k0/s_k^2)/prec (L_k0, S_k0
tree-0 occupied-leaf count and value sum). At the SBC config the
PRECISION error nearly cancels (1.21-1.23x, a coincidence: n/m = 3 vs
mean occupied leaves/tree 2.46). The real defect is the MEAN:
-num/prec is the precision-weighted grand fit mean, so the move
deterministically subtracts the WHOLE level every sweep, leaving it
pinned at N(0, tau^2/(nK)) - sd 0.091 (predicted 0.0907, measured
0.089-0.097) against a correct stationary spread 0.874-0.917. That
mean belongs to no valid absorption, and the critique closed the
loopholes: it is not rescued as MH, overrelaxation, PX-DA, or by the
composite sweep (stationary sd(level) 0.091 vs 0.879/0.910 under two
independent exact schemes on the same build).

EVIDENCE FIT. Per-cell posterior shrinkage rho = 0.690-0.701 predicts
ecdfDiff 0.118 +/- 0.022 against the recorded 0.111/0.114/0.117; the
level is redrawn iid every sweep (ACF ~ 0), which is why 3x spacing
cannot shrink the U; identified content moves < 0.4%, which is why
the softmax arm passes at both chain lengths. The 20-first/29-last
bin asymmetry is 1-2 bin-sd at R=200 - no second defect. (memoA's
pooled-rho underprediction of 0.098 was its estimator, refuted by the
critique's per-cell measurement.)

INERTNESS: not provable, leak second order. The PG psi = f - C and
every tree t != 0 are exactly shift-equivariant; the one leak is tree
0's constant-leaf conditional, non-equivariant under the leaf prior
(the birth log-ratio picks up -4.90 c^2 at the SBC config; tree-0
mean leaf value fluctuates sd 0.26 vs 0.02 for trees 1..49). Effect
on identified functionals: < 0.4%, zero-mean in the level.

THE FIX (unanimous): uniform absorption - add c/m_k to every occupied
leaf of every tree - with its exact conditional
prec = sum_k L_k/(m_k^2 s_k^2), num = sum_k S_k/(m_k s_k^2) (L_k, S_k
over ALL trees' occupied leaves, the fillBottom pattern BCF's
afterCombine already uses), c = -num/prec + N(0,1)/sqrt(prec). One
normal draw as today, so no stream reshuffle beyond the changed
values. ~25-30 lines in afterCombine. Besides exactness it is the
better mixing device: level ACF 0.325/0.035 (lag 1/5) vs 0.965/0.849
for the exact-tree-0 variant, and at the intercept-only exact-gate
configuration it reduces to an exact independence sampler drawing the
level from its true marginal N(0, tau^2/K).

COUNTERFACTUAL VALIDATION (critique, patched build, real SBC arm,
same seeds): the implemented-draw build reproduces the recorded run
to the last digit (f cells 0.1107/0.1144/0.1170, chisqP 0.000, the
f.1.1 histogram 20...29); the fixed build clears all three raw cells
at even the per-functional 5% band (0.0635/0.0733/0.0356, chisqP
0.735/0.229/0.383) with softmax p unchanged-PASS both ways;
multinomial-exact.R passes both ways with identical gaps.

BLOCKING CORRECTIONS for the implementer (critique):
1. tests/cpp/test_sampler.cpp ~2052-2113 must gain real single-node
   trees and a sized muByTree - a bounds guard alone makes the move a
   no-op, the softmax-invariance check vacuous, and shifts the shared
   rng stream for every later testMultinomial block.
2. The accumulation loop needs a non-const Forest& (it writes
   tree.bottomScratch); bound the tree loop by
   min(numTrees, trees.size(), muByTree.size()).
3. Draw-changing beyond the sweep: bartcore_createMultinomialCounts
   shares the combiner (grouped-count channels change) and
   growForestFromRoot (chain.hpp ~1268) calls afterCombine (XBART
   warm start changes) - re-record ALL multinomial equivalence paths;
   every other family stays bitwise.
4. The total = sum-of-tree-fits invariant is preserved TO ROUNDING
   (m * fl(c/m) != c), as today - the comment must not claim exact.
5. docs/design/multinomial.md's "exact Gaussian full conditional ...
   (f_ik ~ N(0, tau^2))" sentence is the defect in prose; rewrite as
   the leaf-space statement whichever way VD forks.
Non-blocking, worth a comment: sampleNodeParametersFromPrior draws
empty bottom nodes too (inert today, contradicts the support argument
in spirit); the blanket muByTree shift is safe only because
multinomial chains are ConstantGaussianLeaf-only.

Full memos, critique, and scripts: session scratchpad (disposable);
this section is the durable record.

## Landing (2026-08-05, ec2a3d0)

Implemented per the Memo with all five blocking corrections honored:
uniform absorption (+c/m_k on every occupied leaf of every tree, the
fillBottom pattern), exact leaf-space conditional, one normal draw with
the pre-draw early return; non-const Forest& with the tree loop bounded
by min(numTrees, trees.size(), muByTree.size()) (numTreesOf helper);
the tests/cpp K=3 combiner fixture gained real single-node trees plus a
non-vacuity check (the shift demonstrably moves every fit and leaf by
the same c); the invariant comment says "to rounding"; the design doc's
level-centering passage rewritten in leaf space, including the second
stale sentence (the calibration section's claim that the centering
conditional reads tau - it reads the per-leaf sd).

Gates, run by the implementer and re-run independently from a fresh
--preclean install before landing: tests/cpp from clean all pass; ASAN/
UBSAN tests/cpp leg clean (implementer); tinytest 3474/0 (no snapshot
regeneration needed - the multinomial suites assert invariants, not
draws); equivalence-7903855 27/27 identical draws; bcf-equivalence-
99205ee bitwise on every channel; multinomial-equivalence-8c2b5fc
MISMATCHES on all 3 scenarios x all 5 channels (tree structure moved
too - expected for a draw change feeding the centered fits), re-recorded
as multinomial-equivalence-ec2a3d0.rds, self-reproducible;
multinomial-exact.R full PASS (max gaps 0.0003-0.0008 vs tols
0.008-0.015); SBC acceptance `sbc.R multinom 200 150 30` (549 s) every
functional PASS at band 0.1282 - the three raw f cells at
0.0688/0.0824/0.0675 (chisqP 0.657/0.031/0.066) against the recorded
failures 0.111/0.114/0.117 (chisqP 0.000); softmax probabilities PASS
unchanged. CI sanitizer watched to green post-push per the sanitizer
discipline.
