# perturb: a same-variable cut move

Status: PROPOSED, 2026-09-07.

A fourth tree kernel that keeps a node's split variable and displaces only its cut, by a small fixed number of grid positions.
[4.2 A same-variable cut move ("perturb") - the first tree-space candidate](tree-mixing-proposals.md#42-a-same-variable-cut-move-perturb---the-first-tree-space-candidate)
ranked it first among tree-space candidates on a first-principles argument and weak external evidence; Stage 0 has since run and its
cut probe prices the move directly on this sampler, in section 6.1's 2026-09-07 addendum
([6.1 Stage 0 - the move census (pilot; no kill criterion)](tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)).
Every census number below is from it.

**Premise, not reopened here.** The swap move was removed pre-release (VD, 2026-09-07), its slice landing first and taking the one
bundled baseline re-record
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)). The mixture is therefore `birth_death 0.6,
change 0.4, birth 0.5`, `perturb` the fourth named element of `proposal.probs` and the third STRUCTURAL probability, and Stage 2's
control arm that mixture. Every count below is written against the post-removal tree: the surface work is the removal's inverse, so
each site the removal emptied is one a `perturbProbability` refills.

## 1. What change does, and why a displacement is a different move

[`changeMove`](../../src/bartcore/moves.hpp) picks uniformly among interior nodes, redraws the split variable from the prior
([`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp)), then draws a cut uniformly over the descendant-valid set with
the skeleton below held fixed. A pure cut move happens only when the redraw lands back on the incumbent, about `1/p_avail` per
proposal ([2. What the sampler can and cannot do today](tree-mixing-proposals.md#2-what-the-sampler-can-and-cannot-do-today)).

The census measured what that costs. Change accepts 4.07 percent of its proposals at the default cell, 1.61 at low noise, 5.37 at
wide, 2.80 at the causal-forest cell, and its rejections are not close calls: the median LOG-LIKELIHOOD difference among REJECTED
scored change proposals - what [`rejectionTable`](../../benchmarks/R/move-census.R) reports, likelihood only, conditional on
rejection - is -56.61 at the default cell and -135.19 at low noise, 0.78 and 0.27 percent of them within one log unit of zero. That
answers the fork 6.1 pre-registered against the temperature family: the bottleneck is proposal accuracy, not scale. The same run
probed the displacement - at each interior node a change proposal visited, hold the variable, move the cut +/-1, 2, 4, 8 positions
on a snapshot, take the FULL MH log ratio (subtree-below prior plus resolved likelihood), restore:

    cell       |1|    |2|    |4|    |8|   median full log ratio at |1|
    default   38.34  23.76  12.46  7.69                          -1.83
    lownoise  26.43  14.75   6.88  3.83                          -4.44
    wide      47.50  30.74  17.40 10.20                          -1.04
    bcf       34.08  21.02  10.91  7.17                          -2.40

The two are not the same statistic - the probe's is unconditional and carries the prior term, change's is conditional on rejection
and carries none - but they are comparable in the direction that matters, since change rejects 96 percent of its scored proposals so
its unconditional median sits near its rejected one, and an O(1) prior term cannot close a gap in the tens. One position sits in a
workable 26 to 48 percent band in every cell, two only at `wide`, so the displacement buys back roughly 50 to 130 log units: the
random-walk Metropolis tuning argument made concrete, step size being the only free parameter any dbarts structural move has.

Two findings are not in the move's favour. Change acceptance RISES from depth 0 to depth 1 in every cell - 2.91 to 6.83 percent at
the default cell, and on to 7.55 and 9.21 at depths 2 and 3 there, though `wide` falls back to 2.01 at depth 3 on 149 proposals -
the opposite of what
[3.3 A high split cannot be changed once the tree is deep (CODE-DERIVED HYPOTHESIS)](tree-mixing-proposals.md#33-a-high-split-cannot-be-changed-once-the-tree-is-deep-code-derived-hypothesis)
predicts, so a cut move cannot be sold on unfreezing high nodes; and `wide`, where two positions also works, is where change ALREADY
does best, so a width rule keyed to it would be keyed to the easy cell. What survives is
[3.2 Structure freezes when the noise level is low (ESTABLISHED)](tree-mixing-proposals.md#32-structure-freezes-when-the-noise-level-is-low-established)
and, indirectly,
[3.5 Tree size moves by a random walk (CODE-DERIVED HYPOTHESIS)](tree-mixing-proposals.md#35-tree-size-moves-by-a-random-walk-code-derived-hypothesis).

## 2. The move

**Eligible nodes: the interior nodes whose rule is on an ordinal column. Nothing else.** The selected set must be a tree function a
cut displacement cannot move, its reciprocal having to cancel between `T` and `T'`; that one does, a displacement changing no node's
variable and no node's shape. **A, draw among ALL interior nodes and no-op on a categorical one**: simplest, but wastes a proposal
per categorical node, most of them in a mixed design. **B, filter to the ordinal interior nodes**: one pass over `fillNotBottom` and
a `splitsBySubset` lookup each, about five lines. **RECOMMEND B.**

**A width filter is forbidden, and not for cost: it is not invariant, so it biases the sampler.** Counterexample: the root splits
`x1` at `c = 2` and its left child splits `x1` at index 0 with no further `x1` below. The root's window is `[lo, hi] = [1, ...]`, so
`c' = 1` is legal, while the left child's own interval `[0, c-1]` is width 2 at `c = 2` and width 1 at `c' = 1`. The selected set's
size moves with the proposal and `1/|S|` stops cancelling; the same construction runs upward through a same-variable ancestor, whose
`Walker::minMax` reads this node. A degenerate interval is a no-op inside the kernel, never a node the selector skips.

A categorical rule has no cut to displace; its analogue, a single-category flip in the direction mask, needs its own validity walk
([`categoricalSubtreeIsValid`](../../src/bartcore/moves.hpp)) and correction, and is not proposed here - so perturb is INERT on an
all-categorical design, the null control
[6.3 Stage 2 - benefit, with matched exposure](tree-mixing-proposals.md#63-stage-2---benefit-with-matched-exposure) requires.

### 2.1 The proposal, and the boundary

Let `[lo, hi]` be [`findGoodOrdinalRules`](../../src/bartcore/moves.hpp) at the node on its own variable - the interval keeping
every same-variable descendant satisfiable - and `c` the current split index. Define `W(c) = {j in [lo, hi] : 0 < |j - c| <= w}` and
draw `c'` uniformly from it. Then

    |W(c)| = min(hi, c + w) - max(lo, c - w)
    logProposalCorrection = log|W(c)| - log|W(c')|

`|W(c)| >= 1` whenever `hi > lo`, the current rule always lying inside its own valid interval; `hi == lo` gives an empty window and
a no-op. The correction is EXACT and needs no new machinery: [`Tree::splitInterval`](../../src/bartcore/tree.hpp) and
`findGoodOrdinalRules` both ignore the node's OWN rule and read only ancestors and descendants, neither of which a cut displacement
touches, so `[lo, hi]` is identical on `T` and `T'` and the reverse count is taken on the unmodified tree. At `w = 1` the window is
`{c-1, c+1}` clipped, `|W|` is 1 or 2, and the correction is bounded by `log 2`. Equalling or crossing an ancestor's or descendant's
cut on the same variable is IMPOSSIBLE, not handled: `splitInterval` sets the bounds one index inside every ancestor's, and
`findGoodOrdinalRules` sets `lower = leftMax + 1`, `upper = rightMin - 1` against the whole subtree.

Three alternatives to A, the clip-and-correct scheme above. **B, reflect at the ends:** folded about `lo - 0.5` and `hi + 0.5` it
really is symmetric with correction 1, but a step off the end returns the node to itself, so at `w = 1` it spends half its boundary
proposals on self-transitions that A converts into real moves. **C, propose on `{-w..w}\{0}` and no-op outside `[lo, hi]`:**
symmetric on its support, correction identically 1, fewest lines, but it throws proposals away exactly where the interval is narrow,
a cut pinned by a descendant and a nudge worth most. **D, uniform over the whole valid interval:** correction identically 1 and no
window, but this is the change move restricted to one variable, rejections at -56.61. **RECOMMEND A**: it wastes no proposal, its
counts come from a function the move already calls, and at `w = 1` the correction cannot exceed `log 2`. Its cost is an arithmetic
path firing only at the interval ends, which is why section 4's gate must place cuts there. The missing-value direction bit
([`Rule`](../../src/bartcore/tree.hpp)) carries over unchanged, contributing `log 2` to both sides of the prior ratio.

### 2.2 Acceptance, the veto, and the grid

`changeMove`'s shape, with one simplification: the node's own prior factors cancel EXACTLY rather than against a proposal density.
`splitVariableLogProbability` reads ancestors only, and
[`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp) is `-log|SI|` over the ancestor-constrained interval,
which an unchanged variable leaves fixed - so `changeMove`'s whole per-side `|Valid|/|SI|` machinery collapses to the window ratio
alone. So

    alpha = exp( B(T') - B(T) + logL(T') - logL(T) + logProposalCorrection )

with `B` the [`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp) of the subtree STRICTLY BELOW the node, which
`changeMove` already computes. It is not zero: a moved cut moves each descendant's own `splitInterval` and can move a variable in or
out of availability below. Three of `changeMove`'s checks drop - no mask pool (an ordinal rule allocates no words), no interaction
walk (`tree.interactionSubtreeIsValid` tests co-occurrence and order of split VARIABLES, which a displacement never moves), no
stranding check (`[lo, hi]` strands none) - and [`maintainMonotoneLeafStore`](../../src/bartcore/chain.hpp) switches only on birth
and death, so the new `StepType` falls through as change does. The snapshot is
[`SubtreeSnapshot`](../../src/bartcore/tree.hpp) from [`MoveScratch`](../../src/bartcore/moves.hpp), node CONTENTS for a fixed id
set, all a shape-preserving move needs.

`[lo, hi]` guarantees logical satisfiability, never occupancy: a displaced cut can empty a descendant leaf, raising the proposal's
rank in [`resolveVetoRank`](../../src/bartcore/moves.hpp) and driving `alpha` to 0 whenever the CURRENT branch is the better ranked
one ([Which move paths can create an empty leaf](empty-leaf-veto.md#which-move-paths-can-create-an-empty-leaf)); from an
already-vetoed state the ordering runs the other way and a rank-improving proposal is accepted outright. At `w = 1` only one cut
bin's rows cross, so the exposure is small, and it is inside the probe's band, which folded `resolveVetoRank`'s `-Inf` into its log
ratio. PERTURB'S OWN veto share is the unrecorded `vetoed.pct` column of section 7; the 0.07 to 0.25 percent of a cell's rejections
6.1 reports is birth, death, change and the then-live swap's, quoted here only as an order of magnitude.

`w` counts GRID POSITIONS, and the default `useQuantiles = FALSE` lays `n.cuts` uniform cuts over the observed range whatever a
column's distinct-value count ([`fillCutsOverRange`](../../src/bartcore/data.hpp)). On a coarse or discrete column ADJACENT SPLIT
INDICES THEN INDUCE THE SAME PARTITION, and a `w = 1` perturb between two of them re-routes zero rows: likelihood difference exactly
0, prior difference small, `alpha` at or near 1. The failure is a NULL MOVE accepted for nothing, not merely a smaller physical
step, and section 5 cell 3 is where it bites. The census cells do not exhibit it - Friedman, continuous, n = 5000, 100 cuts, ~50
rows per bin - so the band is uncontaminated but does not transfer to a coarse design; a mid-chain grid rebuild moves the physical
step again ([`mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)). The move does NOT read
[`scanOrdinalCuts`](../../src/bartcore/scan.hpp); an informed perturb needs the reverse window's scan too, a door priced with
section 4.5's. The variance forest runs the identical kernel ([`sweepVarianceForest`](../../src/bartcore/chain.hpp)), so perturb
reaches it free.

### 2.3 Window width

**A, a compile-time constant `w = 1`**, the only value in the census band in all four cells; cost, no runtime grid, so a width arm
needs a private build. **B, a fraction of the node's interval** (OpenBT ships 10 percent, Pratola fixes 85, no comparison anywhere);
cost, the census measured POSITIONS, so a fraction reproduces none of its numbers, and 10 percent of the FULL default grid is about
ten positions, where acceptance is 3.83 to 10.20 percent - though an interior node's interval is narrower than the grid, so the
identification is loose. **C, a user knob** (`perturb.width`, a new formal and slot); cost, a knob no user can set from evidence,
plus its Rd, validity and refusals. **RECOMMEND A**, any width arm run on a private `-D` build as the census itself ran
(`R_MAKEVARS_USER` appending to `CPPFLAGS`, a private library); C stays additive and pre-release costs nothing if a confirmatory
width ever proves cell-dependent.

## 3. The mixture, and the surface

`proposal.probs` gains `perturb` beside `birth_death`, `change` and `birth`.

**Where the share comes from. A, from change** - `birth_death 0.6, change 0.4 - d, perturb d`, proposal count fixed. **B, from
birth_death**: costs the only dimension-changing moves, already at 10 to 13 percent acceptance. **C, an extra pass on top**
(OpenBT's shape): [`metropolisJumpForTree`](../../src/bartcore/moves.hpp) runs exactly once per tree per sweep in
[`Chain`](../../src/bartcore/chain.hpp)'s loop, so a second pass is a second call site - an engine change with its own surface and
slice. **RECOMMEND A**; C is dropped, and section 5 drops the arms that depended on it.

**Default share: 0**, bitwise-neutral (section 6), which satisfies
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
absolute null-control gate by construction.

Surface. SIX R spellings of the mixture: [`defaultProposalProbs`](../../R/model.R); the [`dbarts`](../../R/dbarts.R), `bart2`
(R/bart.R) and [`dbartsSpec`](../../R/spec.R) formals; and TWO literals inside the monotone branch - the comparison default and the
birth/death-only rewrite (["'monotone' forces birth/death-only proposals"](../../R/spec.R)), the first of which does NOT read
`defaultProposalProbs`, so leaving it stale makes the refusal compare against a vector that no longer exists. Both `all.equal`
branches, monotone and treatment-forest, live in [`resolveSamplerSpec`](../../R/spec.R), which `dbarts()`, `bart2()` and
`dbartsSpec()` all route through, so a defect there fires from every entry point; only the treatment-forest branch reads
`defaultProposalProbs`. [`dbartsModel`](../../R/A_class.R) gains a `p.perturb` slot, prototype and sum-to-one validity. In C++, a
`perturbProbability` takes the eight sites the removal emptied in [`parseModel`](../../src/R_interface_bartcore.cpp)'s file (parsed
struct, read, sum check, creation printout, options copy, two-forest refusal, forest spec, multinomial parameters), the ten across
[`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp), and the three in
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp) - `Forest` being what
[`MoveContext`](../../src/bartcore/moves.hpp) is built from, so without it the kernel is unreachable and without the two specs it is
silently zero in BCF and multinomial fits. Four Rd files: ["proposal.probs"](../../man/dbarts.Rd),
["proposal.probs"](../../man/bart2.Rd), man/dbartsSpec.Rd and ["proposalprobs"](../../man/bart.Rd) (`bart()` itself takes `NULL`, so
the Rd alone). Plus tests/cpp/test_moves.cpp for the interval-invariance assertion, inst/NEWS.Rd, and five tinytest files: the
pinned default vector in test-argument-surface.R, the two-forest refusal's own literal in test-bcf-creation.R, the monotone slot
reads in test-spec.R and test-monotone.R, and test-sum-to-one-tolerance.R below. **At least twenty files, not eleven.** The mixture
rides none of the stored state - `storeState` writes forests, sigma, scale, latents, DART, RNG, glue and the digests, and no
proposal probability among them - so `stateFormatVersion` does not move.

**Three traps.** First, the `all.equal` comparison: widening `defaultProbs` makes `proposal.probs[names(defaultProbs)]` return `NA`
for any caller vector lacking `perturb` - verified in R - so the refusal fires SPURIOUSLY on every caller passing the documented
default, which man/dbarts.Rd, man/bart2.Rd and inst/tinytest/test-argument-surface.R all spell out and which stan4bart forwards
verbatim from `bart_args`. (inst/tinytest/test-monotone.R spells a NON-default vector, deliberately, to trigger the refusal.)
Compare the RESOLVED slots, or fill the missing name first. Second, the engine-side twin, src/R_interface_bartcore.cpp's two-forest
refusal, hard-codes the values and is staled by both the removal and the new name; not a leak, sum-to-one making a nonzero `perturb`
beside two defaults unrepresentable, but it must be restated or it refuses the new default. Third, the one-NA fill, which breaks in
BOTH directions. Post-removal `initialize` selects `c("birth_death", "change")` and perturb makes that THREE names. A caller naming
exactly ONE presents two `NA`s, the branch does not fire, and they reach the slots and fail `setValidity`. A caller naming TWO gets
`perturb` filled with the residual, which is worse because it is silent: **inst/tinytest/test-sum-to-one-tolerance.R is the gate
that must keep failing**, and its `makeModel(1e-7)` names `birth_death` and `change` and expects a sum error - under a naive
widening `perturb` becomes the single `NA`, is filled with `1e-7`, the sum is exact and the expected error disappears. Default
`perturb` to 0 when absent BEFORE the fill, and both that test and today's behaviour are preserved.

**The flat C header does not move.** `dbarts_sampler_create` takes the model as a `SEXP`, so the mixture never crosses the C ABI and
[`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) stays as it is: no `LinkingTo` consumer
recompiles. stan4bart on bartcore uses `formals(dbarts::dbartsSpec)` only for the allowed name set and forwards values from
`bart_args`, so an unnamed caller inherits the new default; bartCause on dbarts-1.0 names no proposal argument at all.

## 4. Correctness: `perturb-balance.R`

[6.2 Stage 1 - correctness (perturb-balance.R, new)](tree-mixing-proposals.md#62-stage-1---correctness-perturb-balancer-new) asks
for a per-kernel exact gate on the WITHIN-VARIABLE cut distribution, the one quantity no shipped gate reads: `change-balance.R`
computes it in [`cutReport`](../../benchmarks/R/change-balance.R) and only PRINTS it, its own pass/fail statistic being a
root-split-VARIABLE marginal at `|z| < 4` on a different problem, which this gate neither borrows nor supersedes.

**The target.** The active-row mask is a PRECISION channel, so under an all-zeros mask
[`Tree::leafVetoRank`](../../src/bartcore/tree.hpp) returns 1 for a leaf holding rows and 2 for one holding none: the rank-0 set is
empty, the likelihood difference is exactly 0 on both branches, and the kernel is reversible with respect to the CGM prior TRUNCATED
to MEMBER-OCCUPIED trees and renormalized - [`Tree::bottomNodesAreOccupied`](../../src/bartcore/tree.hpp)'s predicate, NOT
`bd-balance.R`'s `isAdmissible`, which tests positive weight and is empty under this mask.

**The design makes that truncation vacuous.** Two ordinal columns of 6 and 4 distinct values as a FULL FACTORIAL, at least one row
per cell of the 24, and `useQuantiles = TRUE`, which is what gives 5 and 3 cuts
([`finishQuantileGrid`](../../src/bartcore/data.hpp): `inducedNumCuts = numUnique - 1`; the default FALSE would lay 100 uniform cuts
over 6 values and make empty leaves common). Quantile cuts fall between consecutive distinct values, so every leaf of every
reachable tree holds a cell and the target is the plain CGM prior. Every chain starts at a bare root, the initializer taking a stump
rather than a prior draw under an all-zero composed vector, so the run needs a burn-in it would not otherwise.

**Three statistics; the space is not enumerable.** With 5 and 3 cuts there are 33,610,060,775 distinct trees, so no chi-square over
trees. (1) The root's (variable, cut) marginal plus the stump, NINE states, closed form `P(grow) x P(v) x 1/|SI_v|`: 0.095 for each
of x1's five cuts, 0.158333 for each of x2's three, 0.05 for the stump. (2) The leaf count, from a dynamic program whose state is
(remaining x1 cuts, remaining x2 cuts, DEPTH) - `growthProbability` is `base/(1 + depth)^power`, so remaining cuts alone does not
define the recursion; it gives 0.0500, 0.5523, 0.2796, 0.0913, 0.0219 for one to five leaves over a support of 1..24. (3) The (root
cut, left-child cut) joint on the SAME variable, thirteen states, closed form - the descendant-valid interval and the clipped
window, which nothing gates today.

**One gate.** States of prior mass below 0.004 are dropped by this pre-stated rule, their counts being degenerate at any feasible
run length: statistic 1 keeps all nine, statistic 2 one to five leaves plus a pooled `>= 6` bin (mass 0.00493), statistic 3 the six
states at root cut 1 and 2, dropping the seven at 3 and 4. Family size `m = 21`; batch-means z per state, Holm at family alpha 0.05,
thresholds running from `|z| = 1.96` (least strict) to `|z| = 3.04` (strictest, on the most extreme test). Run length 4 chains x
250,000 kept draws at `n.thin = 20`, batch means over 500 batches per chain: even at an autocorrelation time of 20 that is ~50,000
effective draws against the ~630 statistic 1 needs for poison (i) at `|z| = 3.04`, and ~245 expected counts in the smallest retained
state. The arm's mixture is `birth_death 0.10, change 0.10, perturb 0.80` - change retained because it moves the root's VARIABLE
directly (a perturb-only chain would have to pass through a stump, 5 percent of the prior mass) and birth/death because statistic 2
moves through nothing else, both at the smallest share keeping their own statistic non-degenerate while leaving perturb dominant on
statistics 1 and 3. `changeMove` is separately gated ([The gate](change-move-balance.md#the-gate)), so borrowing it proves nothing
about it.

**Poisons.** (i) Drop `logProposalCorrection`. The uncorrected chain is then reversible for `p(c)` proportional to `pi(c)|W(c)|`, so
moves OUT of an end are over-accepted and moves INTO one under-accepted and the boundary cuts starve. The effect is computable in
advance: on the two-leaf shape, 0.55 of the prior mass, `pi(c)` is uniform over x1's five cuts, so the poisoned conditional law is
`(1,2,2,2,1)/8` against a true 0.2 each - the end states fall to 0.125, a 37 percent relative shift, diluted by the arm's change and
birth/death shares. (ii) A one-sided `+w` window, which drives the cut to `hi` and shifts far more. Both must fail statistic 1. A
third check is a tests/cpp assertion rather than a poison: `findGoodOrdinalRules` must return the same pair before and after any
in-interval rule is installed at the node, the invariance the reverse count rests on.

**The confirmation arm.** The prior-only arm never exercises the likelihood term or the rank-1/rank-2 boundary. **A, ship it
alone**, cheap, the likelihood path being verbatim `changeMove`'s and already gated; **B, add an exact-posterior arm** on the same
grid with positive weights and NO mask - the mask is what makes the likelihood constant, so the arms cannot share a configuration -
scored against `change-balance.R`'s region dynamic program, 150 to 200 lines plus a calibration match. **RECOMMEND B**: a defect in
the descendant-valid interval that bites only when a leaf's occupancy changes is invisible to a run whose leaves are all occupied by
construction, and 6.2 asks for the exact posterior by name. An occupancy arm on `bd-balance.R`'s veto pattern is a further door.

## 5. Benefit, pre-registered

**Two arms**, matched seeds, paired: **A**, the post-removal control `birth_death 0.6, change 0.4, perturb 0`; **B**, `perturb d`
taken from change. 6.3's arms C and D are DROPPED - C is an extra `metropolisJumpForTree` pass, an engine change with no surface
here and no slice, and D exists only to price C. The alternative, pricing them as a fifth slice, is not recommended: it doubles the
engine work to measure a dosage no default would ship at.

**Dosage grid: `d` in `{0.04, 0.10, 0.16}`**, one move per tree per sweep making `d` exactly attempts per tree per sweep, so it
spans the survey's band. `d = 0.40` is EXCLUDED because it leaves change at 0, making arm B birth/death-only in every respect the
two null controls measure, so both would fail by construction rather than on the estimator
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function) records birth/death-only failing
that control outright). **Width: `w = 1` only.** The selection rule below already resolves the width from Stage 0's frozen table, so
a `w = 2` build is exploratory and is not priced in the confirmatory run.

**The selection rule, fixed before Stage 1.** The confirmatory `(w, d)` maximizes `d x accept(w)`, expected accepted perturbs per
tree per sweep at the low-noise cell, subject to `accept(w) >= 0.25` and `d <= 0.16`. Stage 0's table resolves it now: `accept(1) =
0.2643` passes, `accept(2) = 0.1475` does not, so `w = 1, d = 0.16`.

**Four cells.** These replace 6.3's four survey cells, which predate the battery.

1. **C1's independent design at 75 trees, the PRIMARY, on Trig+poly.**
   [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial) establishes that the pre-registered
   correlated arm is not the published setting and that this one is where the coverage deficit reproduces with headroom - 0.822
   against a nominal 0.95, the correlated 200-tree arm at 0.957 with none. The statistic is C1's OWN 95 percent pointwise coverage
   of true `f`; Single index is reported beside it and is not gated. A HARM is
   [6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
   -0.010 core margin. 2. **P1, section 13's rung** (n = 2000, m = 200, `sigma = 0.25`), the known-positive control: it must return
   90 percent coverage near 0.71 or no verdict is valid, which is the rung 6.4's absolute gate names. Not Pratola's n = 5000,
   `sigma^2 = 0.1` rung, whose 0.71 this is not. No win is claimed here -
   [6.3 The pathologies](benchmark-surfaces.md#63-the-pathologies) records three shipped proposal mixtures indistinguishable on it,
   and perturb is a proposal-mixture change. 3. **P2's duplicate-column null**, must-not-degrade, carrying two facts. Arm A is
   ALREADY degraded on it: 10.1's no-swap arm returns 70.8 mean switches per chain but a minimum of 0, 5 of 40 chains parked at an
   x3 root, between-chain sd 0.149 against the default's 0.051. And the cell IS section 2.2's null-move hazard -
   `surfacesDuplicateColumnNull` draws a four-value grid and the harness calls `bart2` at `n.cuts = 100L, useQuantiles = FALSE`, so
   100 uniform cuts sit over four values and about 98 percent of `w = 1` displacements re-route zero rows. **Alternative i, run the
   nulls at `useQuantiles = TRUE`**, four values giving three cuts so perturb is a real move; cost, the cell leaves 10.1's recorded
   configuration and arm A's level must be re-measured, about 3 unit fits. **Alternative ii, keep the shipped grid**, in which case
   what the cell measures is the cost of spending change's share on accepted null moves - a real shipped-configuration cost, not a
   test of mixing. **RECOMMEND i for the gate, ii beside it as a labelled arm**, three unit fits buying that cost as a measurement
   rather than an assumption. Threshold for both: arm B's mean switches within arm A's own seed range, no more than 10 of 40 chains
   parked against arm A's 5, pooled p(root on x1) within Monte Carlo error of 0.5. 4. **The all-categorical null**, where perturb is
   inert by construction - but arm B is not the control there, since it cuts change from 0.40 to 0.24 and discards 16 percent of
   every tree's proposals. The cell therefore runs against its OWN control, `birth_death 0.6, change 0.24` with 0.16 discarded, and
   every metric must sit within Monte Carlo error of it. Read against arm A instead, a failure would void the estimator family for
   arm B's construction rather than the estimator.

**What arm B must produce, in arithmetic.** At `w = 1, d = 0.16`, against the POST-REMOVAL control - the census's per-move rates
recomposed at the shipping shares, `0.6 x 5.145% + 0.4 x 1.61%` at low noise, not at the mixture that carried swap - arm A sits at
0.0373 accepted structural moves per tree per sweep and arm B at `0.0373 - 0.16 x 0.0161 + 0.16 x 0.9961 x 0.2643 = 0.0769`, 2.1x;
at the default cell 0.0771 to 0.1303, 1.7x. The required coverage effect comes from C1's OWN twenty-seed spread, not 6.4's imported
0.010: [`surfacesRange`](../../benchmarks/R/surfaces/surfaces-common.R) prints mean and min-max, so 10.4's independent-75 Trig+poly
range of 0.076 implies sd 0.0203 (`E[range] = 3.735 sd`), a per-arm SE of 0.0045 on the 20-seed mean, and a paired-difference SE
bounded by `sqrt(2) x 0.0045 = 0.0064`, attained only at zero seed correlation. **A WIN is therefore +0.026 absolute**, four times
that bound, moving Trig+poly 0.822 to 0.848 against the 0.100 that 75 to 200 trees buys on the same arm. What the chain lacks:
`accept(1)` and the 0.9961 factor are low-noise FRIEDMAN numbers applied to a C1 target the census never ran, and nothing anywhere
connects a 2.1x accepted-move rate to a coverage move. The design does not predict the bar will be met, and says so here so Stage 2
can fail honestly.

**Kill criterion.** This is the design's own, derived from
[6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered) with every departure stated. **KILL
if, at `w = 1, d = 0.16`, arm B does not improve cell 1's 95 percent coverage over arm A by more than +0.026, over at least 20
matched pairs, with a mandatory fresh-seed re-run of any flagged cell before a flag counts.** Four departures. (a) The cell is C1,
not the low-noise cell 6.4 names, because P1 could not separate three shipped mixtures and this design claims no win there; quoting
6.4 verbatim would fire the kill by construction. (b) 6.4's second conjunct, arm C against arm D, is dropped with arm C, which makes
the kill a single condition and therefore STRICTER than 6.4's disjunction. (c) "four times the measured per-replicate standard
error" is read as four times the twenty-PAIR standard error; four times the per-replicate sd would be 0.081, which no local move
could reach and which 6.4's own 0.010 scaling shows was not the intent. (d) 6.4's "KILL the default question independently" clause
needs plateau prediction error in the noise-heavy or large-n stratum; no cell here is either and none measures it, so that clause
belongs to slice 4, not to slice 3.

The asymmetry stands: passing justifies shipping the move opt-in at weight 0. Flipping the default needs the grow-from-root harm
battery, which is not in `benchmarks/` and must be reconstructed
([5. Verdict and consequences](grow-from-root-default.md#5-verdict-and-consequences)).

[6.5 Cost, honestly](tree-mixing-proposals.md#65-cost-honestly) repriced. Stage 0's instrumentation and driver are spent, both
landed. The kernel is 120 to 160 lines including 6.5's separate +15 for the window ratio, but across at least twenty files rather
than six. `perturb-balance.R` is 400 to 500 with the confirmation arm. The Stage 2 harness stays at 6.5's ~400; what 6.5 did not
price is the private-library build any width arm needs, which the confirmatory run avoids by fixing `w = 1`. Compute is unchanged,
on the order of a day.

## 6. RNG and baselines

After the removal the dispatch is `if (u < bd) ... else change`, so the perturb branch goes in SECOND, at threshold `birthOrDeath +
perturb`, with `changeMove` remaining the `else`. That is what makes `perturb = 0.0` bitwise-neutral: `bd + 0.0` is exactly `bd` in
IEEE for any finite `bd`, so the new test is the old one and control reaches `changeMove` at the same stream position. Giving change
the threshold and perturb the `else` would instead rest neutrality on the probabilities summing to exactly 1.0 and the uniform never
returning 1.0 - a coincidence, not a construction. The move-type SELECTION consumes one uniform per tree per sweep whichever branch
it takes ([`metropolisJumpForTree`](../../src/bartcore/moves.hpp); the kernels draw more), so at weight 0 no stream moves and
`benchmarks/R/equivalence.R` stays green by construction, which is what makes the correctness gate runnable before anything changes
for users.

The stream shift belongs to the swap removal, which lands first and takes the one bundled re-record of every baseline and every
hardcoded tinytest snapshot. A nonzero perturb default would be a SECOND shift, so if one is ever adopted it lands inside that same
re-record rather than paying its own. Stage 2's arm A is therefore the mixture that will ship: measuring perturb against a mixture
the release does not contain would price it against a kernel no user will run, which is also why section 5's arithmetic recomposes
the census's per-move rates at the post-removal shares.

## 7. What the census does not settle

- **No proposal correction was priced.** The probe's log ratio carries no window term: exact at an unclipped node, where the
  correction is 0, and off by at most `log 2` at a clipped one at `w = 1`, in a direction it does not report. - **The probe clamped,
  this move clips, and it weights records not nodes.** A clamped step gives the recommended window's targets at `w = 1`, not at `w
  >= 2`. `cutProbe` writes one record per distinct displacement, so an interior node contributes two `|1|` records and a boundary
  node one, and `cutTable`'s per-record mean UNDER-WEIGHTS the boundary about 2x - exactly where the omitted correction applies;
  degenerate-interval nodes leave the denominator entirely, so the 0.9961 eligibility factor misses them. Magnitudes 3, 5, 6 and 7
  come from narrow intervals only, 139 to 1140 probes against ~28000 per power of two, so acceptance at a narrow interval is
  unmeasured. - **The probe took no move**, being a one-step acceptance on the chain the SHIPPED mixture produced; a chain running
  perturb visits different trees, so its realized rate is not this one. Nor is there a depth breakdown, so whether perturb's
  acceptance rises with depth as change's does is unknown. - **Acceptance is not mixing.** Stage 0 measured no coverage, effective
  sample size, inclusion or error, and had no kill criterion. That is section 5's job. - **The veto's share of the probes WAS
  separated and not recorded.** [`cutTable`](../../benchmarks/R/move-census.R) computes `vetoed.pct` per magnitude; the 6.1
  addendum's table omits the column. Re-summarizing the existing census files records it, and slice 2 should not start without it. -
  **One grid, one tree count, one chain, no categorical and no coarse cell.** Default `n.cuts`, 75 trees (50 in the causal cell's
  treatment forest), 200 burn plus 500 sampled sweeps, continuous columns throughout - so 6.3's all-categorical null has no pilot,
  and section 2.2's null-move hazard has none either, though section 5 cell 3 is where it would be piloted.

## 8. Slices

0. **The swap removal**, a blocking prerequisite for every slice below and the source of every count in this document. DONE: it is
   recorded and landed ([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)) and it took the bundled
   re-record.
1. **The kernel, at default weight 0.** `perturbMove`, a shortened [`changeMove`](../../src/bartcore/moves.hpp), the dispatch branch
   and `StepType` enumerator, `MoveContext`; ten `chain.hpp` and three `combiner.hpp` sites; eight bridge sites including the
   two-forest refusal; six R spellings, the new slot and validity, the guarded one-NA fill, `resolveSamplerSpec`'s two `all.equal`
   comparisons; four Rd files; five tinytest files, tests/cpp/test_moves.cpp and inst/NEWS.Rd - at least twenty files. Roughly 150
   lines of code and 120 of tests: bitwise neutrality, the spurious-refusal regression under `monotone` and a treatment forest,
   test-sum-to-one-tolerance.R still failing where it fails today, a perturb-dominant run that changes cuts and never a variable or
   a shape, and the interval-invariance assertion.
2. **`perturb-balance.R`.** The prior-only arm on the full factorial, the exact-posterior confirmation arm, both poisons. Roughly
   400 to 500 lines; not startable before slice 1, or before the `vetoed.pct` column is recorded.
3. **The Stage 2 harness and run.** Two arms, the four named cells at `w = 1, d = 0.16`, twenty matched pairs. Roughly 400 lines
   plus a day of compute; the verdict is recorded here.
4. **The default share, if slice 3 passes**, carrying 6.4's second kill clause, which needs the reconstructed grow-from-root harm
   battery and lands inside the swap-removal re-record rather than after it.
