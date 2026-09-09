# perturb: a same-variable cut move

Status: PROPOSED, 2026-09-07; AMENDED 2026-09-07 (slices sized, the benefit stage re-primaried on minimum ESS); SLICE 1 LANDED 2026-09-07 (the kernel at weight zero, ab49f83a); SLICE 2 LANDED 2026-09-07 (perturb-balance.R, 30472110); SLICE 3 RUN 2026-09-07: KILL at w = 1, d = 0.16 (d73fb4e0).

Amended by [pure-c-header](../plans/pure-c-header.md#pure-c-header): the flat C header creates no sampler and
no longer declares the predictor, test-data, weight, active-row, per-forest, state,
tree-extraction or augmentation entries - each is a method on the R sampler object the
handle is now read from. The `retired:` cites below name constructs that are gone; what
this record says about the R and engine sides still holds.

A fourth tree kernel that keeps a node's split variable and displaces only its cut, by a small fixed number of grid positions.
[4.2 A same-variable cut move ("perturb") - the first tree-space candidate](tree-mixing-proposals.md#42-a-same-variable-cut-move-perturb---the-first-tree-space-candidate)
ranked it first among tree-space candidates on a first-principles argument and weak external evidence; Stage 0 has since run and its
cut probe prices the move directly on this sampler
([6.1 Stage 0 - the move census (pilot; no kill criterion)](tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)).
Every census number below is from 6.1's SECOND 2026-09-07 addendum, the re-run at the two-move kernel - which IS the shipped kernel,
swap dispatching at a probability of zero - and not from the three-move run beside it. The one exception is marked where it stands,
in section 7: the clipped-magnitude probe counts, which only the three-move addendum reports.

**Premise, not reopened here.** The swap slice landed first and took the one bundled baseline re-record
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)). It dropped swap out of the DEFAULT and then
put the move back at that default of zero
([9. Reversal: the move returns at default zero](swap-removal.md#9-reversal-the-move-returns-at-default-zero)), so the mixture is
`birth_death 0.6, swap 0, change 0.4, birth 0.5`, `perturb` the FIFTH named element of `proposal.probs` and the FOURTH structural
probability, and Stage 2's control arm that mixture. Every count in section 3 is a count against that restored tree: a
`perturbProbability` does not refill a `swapProbability` site but sits beside it, so each is a count of NEW sites and each dispatch,
initializer and literal named gains one element rather than replacing one.

## 1. What change does, and why a displacement is a different move

[`changeMove`](../../src/bartcore/moves.hpp) picks uniformly among interior nodes, redraws the split variable from the prior
([`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp)), then draws a cut uniformly over the descendant-valid set with
the skeleton below held fixed. A pure cut move happens only when the redraw lands back on the incumbent, about `1/p_avail` per
proposal ([2. What the sampler can and cannot do today](tree-mixing-proposals.md#2-what-the-sampler-can-and-cannot-do-today)).

The census measured what that costs. Change accepts 3.77 percent of its proposals at the default cell, 1.66 at low noise, 6.06 at
wide, 2.88 at the causal-forest cell, and its rejections are not close calls: the median LOG-LIKELIHOOD difference among REJECTED
scored change proposals - what [`rejectionTable`](../../benchmarks/R/move-census.R) reports, likelihood only, conditional on
rejection - is -62.34 at the default cell and -143.45 at low noise, 0.87 and 0.25 percent of them within one log unit of zero. That
answers the fork 6.1 pre-registered against the temperature family: the bottleneck is proposal accuracy, not scale. The same run
probed the displacement - at each interior node a change proposal visited, hold the variable, move the cut +/-1, 2, 4, 8 positions
on a snapshot, take the FULL MH log ratio (subtree-below prior plus resolved likelihood), restore:

    cell       |1|    |2|    |4|    |8|   median log ratio at |1|
    default   40.65  24.06  11.94  6.81                     -1.53
    lownoise  27.15  13.73   5.63  3.07                     -4.83
    wide      51.38  34.82  20.77 12.91                     -0.76
    bcf       33.34  21.09  12.31  8.09                     -2.46

The two are not the same statistic - the probe's is unconditional and carries the prior term, change's is conditional on rejection
and carries none - but they are comparable in the direction that matters, since change rejects 96.1 percent of its scored proposals
at the default cell and 93.5 to 98.3 percent across the four, so its unconditional median sits near its rejected one, and an O(1)
prior term cannot close a gap in the tens. One position sits in a workable 27 to 51 percent band in every cell, two only at
`wide`, so the displacement buys back roughly 60 to 140 log units: the
random-walk Metropolis tuning argument made concrete, step size being the only free parameter any dbarts structural move has.

Two findings are not in the move's favour. Change acceptance RISES from depth 0 to depth 1 in every cell - 2.46 to 6.84 percent at
the default cell, and on to 8.44 at depth 2 there, though it falls back to 7.09 at depth 3 on 268 proposals - the opposite of
what
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
all-categorical design. That is what
[6.3 Stage 2 - benefit, with matched exposure](tree-mixing-proposals.md#63-stage-2---benefit-with-matched-exposure) asks a null
control to exploit; section 5 explains why the arm that takes its share from change cannot be that control, and what replaces it.

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
window, but this is the change move restricted to one variable, rejections at -62.34. **RECOMMEND A**: it wastes no proposal, its
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
and death, so the new `StepType` falls through as swap and change do. The snapshot is
[`SubtreeSnapshot`](../../src/bartcore/tree.hpp) from [`MoveScratch`](../../src/bartcore/moves.hpp), node CONTENTS for a fixed id
set, all a shape-preserving move needs.

`[lo, hi]` guarantees logical satisfiability, never occupancy: a displaced cut can empty a descendant leaf, raising the proposal's
rank in [`resolveVetoRank`](../../src/bartcore/moves.hpp) and driving `alpha` to 0 whenever the CURRENT branch is the better ranked
one ([Which move paths can create an empty leaf](empty-leaf-veto.md#which-move-paths-can-create-an-empty-leaf)); from an
already-vetoed state the ordering runs the other way and a rank-improving proposal is accepted outright. At `w = 1` only one cut
bin's rows cross, so the exposure is small, and it is inside the probe's band, which folded `resolveVetoRank`'s `-Inf` into its log
ratio. PERTURB'S OWN veto share is the unrecorded `vetoed.pct` column of section 7; the 0.13 to 0.27 percent of a cell's rejections
6.1's two-move re-run reports is birth, death and change's, quoted here only as an order of magnitude.

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
ten positions, where acceptance is 3.07 to 12.91 percent - though an interior node's interval is narrower than the grid, so the
identification is loose. **C, a user knob** (`perturb.width`, a new formal and slot); cost, a knob no user can set from evidence,
plus its Rd, validity and refusals. **RECOMMEND A**, any width arm run on a private `-D` build as the census itself ran
(`R_MAKEVARS_USER` appending to `CPPFLAGS`, a private library); C stays additive and pre-release costs nothing if a confirmatory
width ever proves cell-dependent.

## 3. The mixture, and the surface

`proposal.probs` gains `perturb` beside `birth_death`, `swap`, `change` and `birth`: the FIFTH name and the FOURTH structural
probability.

**Where the share comes from. A, from change** - `birth_death 0.6, swap 0, change 0.4 - d, perturb d`, proposal count fixed. **B,
from birth_death**: costs the only dimension-changing moves, already at 10 to 13 percent acceptance. **C, an extra pass on top**
(OpenBT's shape): [`metropolisJumpForTree`](../../src/bartcore/moves.hpp) runs exactly once per tree per sweep in
[`Chain`](../../src/bartcore/chain.hpp)'s loop, so a second pass is a second call site - an engine change with its own surface and
slice, priced and declined in section 5. **RECOMMEND A.**

**Default share: 0**, bitwise-neutral (section 6), which satisfies
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
absolute null-control gate by construction.

### 3.1 Every site, enumerated

The removal record lists the analogous sites for swap ([2. What is removed](swap-removal.md#2-what-is-removed),
[3. The surface](swap-removal.md#3-the-surface)) and the restore put every one of them back, so the enumeration below is that list
read as an addition.

**R, five files, six default vectors.** [`defaultProposalProbs`](../../R/model.R); the [`dbarts`](../../R/dbarts.R), `bart2`
(R/bart.R) and [`dbartsSpec`](../../R/spec.R) formals; and TWO literals inside the monotone branch of
[`resolveSamplerSpec`](../../R/spec.R) - the comparison default and the birth/death-only rewrite
(["'monotone' forces birth/death-only proposals"](../../R/spec.R)). The first of those does NOT read `defaultProposalProbs`, so
leaving it stale makes the refusal compare against a vector that no longer exists. Both `all.equal` branches, monotone and
treatment-forest, live in `resolveSamplerSpec`, which `dbarts()`, `bart2()` and `dbartsSpec()` all route through, so a defect there
fires from every entry point; only the treatment-forest branch reads `defaultProposalProbs`, which is what widens its comparison
under the caller's feet - trap 1 below.
[`dbartsModel`](../../R/A_class.R) gains a `p.perturb` slot, a prototype of 0 and a fourth term in validity's sum. Four further
literals carry the structural name set inside R/model.R's initializer - the subset, the `names(probs) <-` assignment, the all-NA
fallback subset and the slot write - plus the fill rule below.

**How a fifth name enters the fill.** The rule landed at 6934e487 runs over three names: one unnamed structural element takes the
residual, two unnamed with `swap` among them still resolve because swap takes its zero and the other the residual
(["two unnamed, one of them swap"](../../inst/tinytest/test-proposal-probs.R)), and `swap` named alone is an error because the
birth/death-versus-change split is undetermined. A fourth structural name cannot simply join that set, and the way it fails is
SILENT. `c(birth_death = 0.5, change = 0.4)` resolves today to `swap 0.1`. Widen the subset to four names and that same call
presents TWO unnamed elements, so `sum(unnamed) == 2L` still zeroes `swap` and the residual then falls to the one `NA` left:
`perturb 0.1`, a vector summing to 1 that `setValidity` accepts and no test reads - a nonzero perturb default nobody asked for. The
variant that zeroes BOTH zero-default names instead lands at `swap 0, perturb 0` and a sum of 0.9, which `setValidity` does catch;
that one is loud, and it is not the hazard. **RECOMMEND: `perturb` resolves BEFORE the
three-name fill and never enters it.** An unnamed `perturb` is 0 - its default is a number, not a share - the residual is then taken
against `1 - perturb` rather than 1, and the existing three-name rule runs verbatim. Every resolution above is preserved element for
element, and the only widening is the error: `c(perturb = 0.16)` alone leaves the same split undetermined as `swap` alone, so that
message names both zero-default moves. The alternative, folding `perturb` into the fill's name set and zeroing the surplus unnamed
elements, is rejected for the `c(birth_death = 0.5, change = 0.4)` regression in either of its forms.

**C++: the kernel's own file, and TWENTY-ONE probability sites outside it.** In src/bartcore/moves.hpp: `perturbMove` beside
[`swapMove`](../../src/bartcore/moves.hpp), a `perturbProbability` on [`MoveContext`](../../src/bartcore/moves.hpp) - which carries
birth/death, swap and birth today and not change, so a fourth field really does have to reach it - a fifth enumerator on
[`StepType`](../../src/bartcore/moves.hpp), the dispatch branch (section 6), the census hooks the other kernels carry - a
[`BARTCORE_CENSUS_NOOP`](../../src/bartcore/moves.hpp) record at each early return, of which `changeMove` has four and `swapMove`
two, plus one shape record and one proposal record on the scored path - and the census legend that names the moves. Then TEN in
[`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp): three struct
fields, four copies into a forest (creation, `setModel`, `buildSpecifiedForest`, `buildMultinomialForest`), the variance forest's
copy and the two `MoveContext` initializers. THREE in
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp) - `Forest` being what `MoveContext` is
built from, so without it the kernel is unreachable, and without the two specs it is silently zero in BCF and multinomial fits.
EIGHT in [`parseModel`, `refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp): the parsed struct's field, the
`p.perturb` slot read, the sum check's fourth term, the creation printout - one `ext_printf`, format string and argument list
together - the `SamplerOptions` copy, the two-forest refusal's hard-coded mixture, the forest spec copy and `setModel`'s
`ModelParameters` copy.

**tests/cpp, two files.** FIFTEEN positional `MoveContext` initializers - twelve in tests/cpp/test_moves.cpp, three in
tests/cpp/test_interaction.cpp - each of which must gain an element. The probability block sits ahead of `const double* weights`, so
an initializer left short binds that pointer to a `double`: a hard compile error, not a silent zero, which is what makes this count
safe to trust. The five loose probability assignments (three in tests/cpp/test_model.cpp, two in tests/cpp/test_sampler.cpp) need no
twin, a new field with a zero default resolving itself. tests/cpp/test_moves.cpp also carries slice 1's interval-invariance
assertion.

**Four Rd files**: ["proposal.probs"](../../man/dbarts.Rd), ["proposal.probs"](../../man/bart.Rd) and man/dbartsSpec.Rd - usage
line and argument text in each - and ["proposalprobs"](../../man/bartBT.Rd), where `bartBT()` takes `NULL` so the Rd alone moves.
**SIX tinytest files**: ["a caller-supplied three-move mixture"](../../inst/tinytest/test-proposal-probs.R), the surface's own test
and the one that pins the fill; the pinned default vector in test-argument-surface.R; the two-forest refusal's literal in
test-bcf-creation.R; the slot reads in test-spec.R and test-monotone.R; and test-sum-to-one-tolerance.R below. **Plus inst/NEWS.Rd.**

**Twenty-two files**: four C++, five R, four Rd, six tinytest, two tests/cpp and NEWS. Roughly 150 lines of kernel including
section 2.1's window ratio, about 60 more spread across the twenty-one surface sites, and 120 of tests. The benchmark harnesses are
not shipped surface and mostly do not move: the eleven one-tree exact-posterior gates that pass an explicit mixture name three
structural probabilities and keep working, `perturb` filling 0 under the rule above; section 4's arm is a new mixture, and section
5's arms are new files.

**Three traps.** First, the `all.equal` comparison: widening `defaultProbs` makes `proposal.probs[names(defaultProbs)]` return `NA`
for any caller vector lacking `perturb` - verified in R - so the refusal fires SPURIOUSLY on every caller passing the documented
default, which man/dbarts.Rd, man/bart2.Rd and inst/tinytest/test-argument-surface.R all spell out and which stan4bart forwards
verbatim from `bart_args`. (inst/tinytest/test-monotone.R spells a NON-default vector, deliberately, to trigger the refusal.)
Compare the RESOLVED slots, or fill the missing name first. Second, the engine-side twin,
[`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp)'s two-forest refusal, hard-codes the mixture and is
staled by the new name - but staled only in what it says, not in what it does. It tests the four probabilities it names, so the new
default passes it exactly as today, `perturb 0` being untested; and sum-to-one makes a nonzero `perturb` beside three unchanged
defaults unrepresentable, so a nonzero one always trips a term already there. It neither leaks nor refuses the new default.
Restating it is completeness - the refusal should name the probability it now omits - not correctness. Third, the fill, which
breaks silently rather than loudly: **inst/tinytest/test-sum-to-one-tolerance.R is the gate that must keep failing**, and its `makeModel(1e-7)` names all THREE
structural probabilities and expects a sum error - under a naive widening `perturb` becomes the single `NA`, is filled with `1e-7`,
the sum is exact and the expected error disappears. Resolving `perturb` to 0 ahead of the fill preserves it.

**The flat C header does not move.** `dbarts_sampler_create` takes the model as a `SEXP`, so the mixture never crosses the C ABI and
retired: [`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) stays as it is: no `LinkingTo` consumer
recompiles. Nor does the stored state: `storeState` writes forests, sigma, scale, latents, DART, RNG, glue and the digests, and no
proposal probability among them, so `stateFormatVersion` does not move. stan4bart on bartcore uses `formals(dbarts::dbartsSpec)`
only for the allowed name set and forwards values from `bart_args`, so an unnamed caller inherits the new default; bartCause on
dbarts-1.0 names no proposal argument at all.

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
define the recursion; it gives 0.0500, 0.5523, 0.2796, 0.0913, 0.0219 for one to five leaves over a support of 1..24. Both closed
forms and the program draw the split VARIABLE uniformly among the variables that still have a cut at the node, which is what
[`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp) does with no `splitProbabilities` set: it counts the available
variables with `collectAvailableVariables` and picks one of that count. A variable exhausted below its ancestors drops out of the
draw rather than wasting it, and the vector above reproduces only under that convention: drawing uniformly over both columns
regardless of availability moves every state from three leaves up, which is where the two conventions first differ. (3) The (root
cut, left-child cut) joint on the SAME variable, thirteen states, closed form - the descendant-valid interval and the clipped
window, which nothing gates today.

**One gate.** States of prior mass below 0.004 are dropped by this pre-stated rule, their counts being degenerate at any feasible
run length: statistic 1 keeps all nine, statistic 2 one to five leaves plus a pooled `>= 6` bin (mass 0.00493), statistic 3 the six
states at root cut 1 and 2, dropping the seven at 3 and 4. Family size `m = 21`; batch-means z per state, Holm at family alpha 0.05,
thresholds running from `|z| = 1.96` (least strict) to `|z| = 3.04` (strictest, on the most extreme test). Run length 4 chains x
250,000 kept draws at `n.thin = 20`, batch means over 500 batches per chain: even at an autocorrelation time of 20 that is ~50,000
effective draws against the ~630 statistic 1 needs for poison (i) at `|z| = 3.04`, and ~245 expected counts in the smallest retained
state. That 630 is the UNDILUTED floor - it is what a 37 percent shift needs - and the realized shift is smaller, poison (i)'s 37
percent being the perturb-only limit and the arm's change and birth/death shares pulling the chain back toward the true law, so the
margin above the floor is deliberate rather than slack. The length is affordable because the likelihood is never evaluated under the
mask: a sweep costs the move machinery alone, on one tree over 24 rows. The arm's mixture is
`birth_death 0.10, change 0.10, perturb 0.80` - change retained because it moves the root's VARIABLE
directly (a perturb-only chain would have to pass through a stump, 5 percent of the prior mass) and birth/death because statistic 2
moves through nothing else, both at the smallest share keeping their own statistic non-degenerate while leaving perturb dominant on
statistics 1 and 3. `changeMove` is separately gated ([The gate](change-move-balance.md#the-gate)), so borrowing it proves nothing
about it.

**The z, and the burn-in.** A state's indicator is averaged within a chain over 500 batches of 500 consecutive kept draws; `s_c` is
the standard deviation of chain `c`'s 500 batch means, so that chain's null variance is `s_c^2 / 500`, batching rather than a
binomial formula being what absorbs the within-chain autocorrelation
([`batchMeanSE`](../../benchmarks/R/change-balance.R) is the same estimator at one chain and 400 batches). The four chains are
independent, so they pool as the mean of the four estimates against `sqrt(sum_c s_c^2 / 500) / 4`, and `z` is the pooled deviation
from the closed-form or dynamic-program mass over that. **Burn-in: 20,000 sweeps per chain**, 1,000 kept draws, discarded before
the batching. It is sized off the LEAF COUNT, the slowest statistic in the family: birth/death holds 0.10 of the one proposal a tree
gets per sweep, so a dimension change is offered about once in ten sweeps and the leaf count is the only statistic that has to
random-walk out of the stump the initializer leaves, while statistics 1 and 3 move on every accepted perturb. 20,000 sweeps is
roughly 2,000 dimension proposals against a support of 1..24, and it costs 0.4 percent of the chain. The script does not assume it:
it reads the first lag at which each statistic's kept-draw autocorrelation falls under 0.1, the ladder
[`firstUnder`](../../benchmarks/R/sbc.R) uses, and refuses to score a run whose burn-in is under fifty of those lags for any
statistic, doubling and re-running instead.

**Poisons.** (i) Drop `logProposalCorrection`. The uncorrected chain is then reversible for `p(c)` proportional to `pi(c)|W(c)|`, so
moves OUT of an end are over-accepted and moves INTO one under-accepted and the boundary cuts starve. The effect is computable in
advance: on the two-leaf shape, 0.55 of the prior mass, `pi(c)` is uniform over x1's five cuts, so the poisoned conditional law is
`(1,2,2,2,1)/8` against a true 0.2 each - the end states fall to 0.125, a 37 percent relative shift, diluted by the arm's change and
birth/death shares. (ii) A one-sided `+w` window. `|W|` is then 1 everywhere but at `hi`, where it is 0, so on the same two-leaf shape every proposal
is `c -> c+1` at `alpha = 1` and the cut is absorbed at the top: x1's five cuts go from 0.2 each to a law whose lower four states
hold only the mass change and birth/death re-seed them with, an order of magnitude past poison (i)'s 37 percent. Both must fail
statistic 1. A
third check is a tests/cpp assertion rather than a poison: `findGoodOrdinalRules` must return the same pair before and after any
in-interval rule is installed at the node, the invariance the reverse count rests on.

**The confirmation arm.** The prior-only arm never exercises the likelihood term or the rank-1/rank-2 boundary. **A, ship it
alone**, cheap, the likelihood path being verbatim `changeMove`'s and already gated; **B, add an exact-posterior arm** on the same
grid with positive weights and NO mask - the mask is what makes the likelihood constant, so the arms cannot share a configuration -
scored against `change-balance.R`'s region dynamic program, 150 to 200 lines plus a calibration match. **RECOMMEND B**: a defect in
the descendant-valid interval that bites only when a leaf's occupancy changes is invisible to a run whose leaves are all occupied by
construction, and 6.2 asks for the exact posterior by name. An occupancy arm on `bd-balance.R`'s veto pattern is a further door.

## 5. Benefit, pre-registered

**Two prerequisites, both outside this design.** First, the battery's absolute gate: `benchmarks/R/surfaces/P1-friedman.R` must
return 90 percent coverage near 0.71 in the control arm before any verdict here is valid
([6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically));
it now exists and reads 0.725 held-out at the shipped default
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)),
so the gate is in force.
Second, the chain configuration of the primary cell, which is now settled and is what the rest of this section is built on.

### 5.1 The chain configuration, and what it makes the primary statistic

[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s 0.822 is a SINGLE-chain reading at
the paper's 1000 burn-in and 2500 kept. Three chain-configuration arms were run on 2026-09-07 over the same twenty seeds, the same
independent design and the same 75 trees, varying nothing else. Trig+poly; min ESS is over C1's own 25 fixed points, summed over
chains, and the between-chain column is the across-chain sd of a chain's posterior mean of `f` over the pooled posterior sd, median
over those points - near 0 the chains agree, near 1 each sits in its own place and pooling is what widens the interval.

    chains                 95% coverage        min ESS summed  per-chain min ESS  between-chain
    4 x 500                0.961(0.945-0.977)  15(8-31)        2(1-2)             0.78
    4 x 2500               0.959(0.937-0.975)  18(9-46)        2(1-2)             0.65
    1 x 25000              0.902(0.858-0.928)  2(1-3)          2(1-3)             -
    1 x 25000, first 2500  0.818(0.770-0.866)  -               -                  -

The fourth row is the third one's own fit read at its first 2500 draws, the single-chain reading at THIS kernel; 10.4's recorded
single-chain arm reads 0.822 with a summed and per-chain minimum ESS of 2(1-4), but it was measured at swap 0.1 and is not a row of
this comparison. There is no wall column: 10.4's host carried a 1-minute load of 8 to 15 throughout, so the script's times carry no
timing claim there and none is made here.

Three readings, all of which bear on this design. **The shipped four-chain default reaches nominal coverage at 75 trees**, and a
single chain does not, even at ten times the paper's length; the same fit read at its first 2500 draws returns 0.818, so the gain is
chain length and not a different fit. **So coverage is not the statistic this benefit stage can win on**: arm A sits at
0.961 against a nominal 0.95 with no headroom left. **And the mixing symptom survives the pooling intact.** A chain's own minimum
ESS is 2 whether it is 500 draws long or 25000, and five times the length moves the summed figure only from 15 to 18. That is a
statistic counting how many separate places the chains find, which is what a proposal mixture acts on, and it reads 15 out of 2000
kept draws. The He-Hahn coverage deficit is read as a MIXING symptom and mixing is the lever (VD, 2026-09-07), so
**the primary is the minimum ESS over C1's 25 fixed points, summed over four chains, on the independent design at 75 trees,
Trig+poly**, and the benefit cell runs the SHIPPED four-chain configuration, `n.chains = 4` at 500 burn-in and 500 kept - a fifth of
the compute of the 1000 + 2500 configuration it matches, and the one that reads the symptom. Single
index is reported beside it and is not gated. The tree-count arm stays in the battery as a comparison only; the default is parked.

**Summed, against 10.4's own sentence.** 10.4 closes "the measure a kernel change has to move on this cell is the per-chain ESS",
and this design gates the summed figure instead. The two are the same measurement at two scales: the chain count is FIXED at four
by this pre-registration, the same four in arm A and arm B, so the summed figure is the per-chain ESS added over a constant, and
adding chains cannot move it because no arm may add one. What differs is resolution. The per-chain column reads 2(1-2) in every
four-chain arm - one integer, a range of one, no spread a paired difference can be taken over - while the summed column reads
15(8-31) and 18(9-46), a spread the +8 bar below is derived from and can resolve. So the per-chain figure is what 10.4's sentence
warns is not to be inflated by chain count, and holding the count fixed is how that warning is honoured; the summed figure is the
readable form of it. The per-chain minimum is reported beside the primary in every cell, and a summed gain that leaves it at 2 is
reported as such.

The three arms are `independent75pool4`, `independent75pool4long` and `independent75long` in
[`arms`](../../benchmarks/R/surfaces/C1-he-hahn.R), landed at fef6dca6 and their readouts recorded at 4cf7b94b; slice 3's harness
reuses the first as the benefit cell rather than adding one.

**The bar the cell can resolve.** [`surfacesRange`](../../benchmarks/R/surfaces/surfaces-common.R) prints mean and min-max over the
twenty seeds, so 15(8-31) is a range of 23. At `E[range] = 3.735 sd` for twenty draws that is sd 6.2, a per-arm SE of 1.4 on the
twenty-seed mean, and a paired-difference SE bounded above by `sqrt(2) x 1.4 = 2.0`, attained only at zero seed correlation and
strictly smaller under matched seeds. **Four times that bound is +7.8, so the bar is +8 summed minimum ESS**, 15 to 23, a ratio of
1.5. Two caveats, both stated rather than buried: the range-to-sd conversion assumes normality and minimum ESS is right-skewed, so
6.2 overstates the spread and the bar is conservative in the direction that costs the move; and
[6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered)'s "four times the measured
per-replicate standard error" is read here as four times the twenty-PAIR standard error, which is the scale every margin
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)
tabulates is stated at, four times a per-replicate sd being 25 - larger than the statistic's whole observed range.

### 5.2 Arms, dosage, and cells

**Two arms**, matched seeds, paired: **A**, the shipped mixture `birth_death 0.6, swap 0, change 0.4, perturb 0`; **B**, `perturb d`
taken from change. 6.3's arms C and D are DROPPED, and priced first so the drop is a decision rather than an omission. Arm C is an
OpenBT-style extra perturb pass, and it is an engine change, not a mixture setting: [`Chain`](../../src/bartcore/chain.hpp) calls
[`metropolisJumpForTree`](../../src/bartcore/moves.hpp) exactly once per tree per sweep, so a second pass is a second call site
carrying a `perturb.passes` control with its own formal, slot, validity, Rd and bitwise-neutrality argument at zero passes - a fifth
slice of perhaps 60 to 80 engine lines and a surface tour of its own - and arm D, which exists only to spend C's extra
compute on plain sweeps, doubles the Stage 2 run. **RECOMMEND dropping both.** They measure a dosage no default would ship at, and
6.4's kill is a CONJUNCTION of two failures - B fails to beat A AND C fails to beat D - so surviving it needs only one of the two
to succeed. Dropping C leaves a single failure condition, which the move must clear on its own: STRICTER than the criterion it
departs from, not looser.

**Dosage grid: `d` in `{0.04, 0.10, 0.16}`**, one move per tree per sweep making `d` exactly attempts per tree per sweep, so it
spans the survey's band. `d = 0.40` is EXCLUDED because it leaves change at 0, and change is the only move that supplies variable
switching at one tree: cell 3's must-not-degrade gate would then fail outright rather than on the estimator
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function) records birth/death-only at 0
switches in 40 of 40 chains). **Width: `w = 1` only.** The selection rule below resolves the width from Stage 0's frozen table, so a
`w = 2` build is exploratory and is not priced in the confirmatory run.

**The selection rule, fixed before Stage 1.** The confirmatory `(w, d)` maximizes `d x accept(w)`, expected accepted perturbs per
tree per sweep at the low-noise cell, subject to `accept(w) >= 0.25` and `d <= 0.16`. Stage 0's table resolves it now: `accept(1) =
0.2715` passes, `accept(2) = 0.1373` does not, so `w = 1, d = 0.16`. The rule is written with the table in hand and clears `w = 1`
by 0.0215, so it selects rather than discovers; that is the point of freezing it here, and it is why the confirmatory run needs one
build and not two.

**Three cells.** These replace 6.3's four survey cells, which predate the battery.

1. **C1's independent design at 75 trees, the PRIMARY, on Trig+poly**, four chains at 500 + 500, twenty matched pairs. The primary
   is 5.1's summed minimum ESS at the +8 bar. Four secondaries, each at
   [6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
   own margin, all must-not-degrade: 95 percent coverage at -0.010 absolute (arm A at 0.961), held-out RMSE at a ratio above 1.02,
   interval length reported, and wall time per sweep at a ratio above 1.05 - the last of which is what stops an ESS win bought
   purely with compute, and is why the primary is stated in draws rather than per second.
2. **P1, the prerequisite absolute gate, not a win cell.** Section 13's rung, n = 2000, m = 200, `sigma = 0.25`: its 90 percent
   coverage must return near 0.71 in the control arm or no verdict from any cell is valid. Not Pratola's n = 5000, `sigma^2 = 0.1`
   rung, whose published 53 percent this is not. No win is claimed here -
   [6.3 The pathologies](benchmark-surfaces.md#63-the-pathologies) records three proposal mixtures indistinguishable on it,
   and perturb is a proposal-mixture change.
3. **P2's duplicate-column cell, must-not-degrade, and it is honestly compromised in two ways.** Arm A is ALREADY degraded on it:
   10.1's no-swap arm - numerically the shipped mixture - returns 70.8 mean switches per chain but a MINIMUM of 0, 5 of 40 chains
   parked at an x3 root, between-chain sd 0.149 against the swap-carrying 0.051. And the cell IS section 2.2's null-move hazard:
   [`surfacesDuplicateColumnNull`](../../benchmarks/R/surfaces/surfaces-common.R) draws a four-value grid and the harness calls
   `bart2` at bart2's defaults, naming neither `n.cuts` nor `useQuantiles`, so the 100 uniform cuts of `n.cuts = 100L` and
   `useQuantiles = FALSE` sit over four values, only two of the 99 adjacent index pairs straddle a value, and about 98 percent
   of `w = 1` displacements re-route zero rows. **Alternative i, run the cell at
   `useQuantiles = TRUE`**, four values giving three cuts so perturb is a real move; cost, the cell leaves 10.1's recorded
   configuration and arm A's level must be re-measured, about 3 unit fits. **Alternative ii, keep the shipped grid**, in which case
   what the cell measures is the cost of spending change's share on accepted null moves - a real shipped-configuration cost, not a
   test of mixing. **RECOMMEND i for the gate, ii beside it as a labelled arm.** The threshold cannot be stated on the minimum,
   which arm A already reads 0: arm B's mean switches per chain within arm A's own seed range, chains parked no more than 10 of 40
   against arm A's 5 (a two-standard-error allowance at that binomial), and pooled p(root on x1) within Monte Carlo error of 0.5.

**The controls, replacing 6.3's all-categorical null.** 6.3 asks for a cell where the move provably cannot act and voids the
estimator family if the arms differ there. Arm B cannot pass such a cell, and the reason is where its share comes from: perturb is
inert on an all-categorical design (section 2), but arm B still cuts change from 0.40 to 0.24 and discards 16 percent of every
tree's proposals, so any difference is arm B's construction and not the estimator's. Taking the share from birth_death instead only
moves which move is diluted. The null is therefore taken in two pieces the move can pass. **(a) The bitwise null**, which is 6.4's
own first absolute gate and is strictly stronger than a statistical one: at `perturb 0` arm B must be BITWISE identical to arm A
under benchmarks/R/equivalence.R. Slice 1 already runs it. **(b) A sham comparison**, arm A against arm A at twenty FRESH sampler
seeds on cell 1, whose paired minimum-ESS difference must sit inside the +8 bar - the only control that calibrates the bar against
the harness rather than against a formula, at the cost of twenty more fits. The all-categorical cell stays, RELABELLED a dilution
arm rather than a null: arm B against `birth_death 0.6, swap 0, change 0.24` with 0.16 discarded, measuring what a discarded
proposal share costs where perturb cannot act. It reports; it gates nothing and voids nothing.

### 5.3 What arm B must produce, and the kill

**In arithmetic.** The two-move re-run measures arm A directly - it IS arm A's kernel - so nothing has to be recomposed: its pooled
acceptance is 0.0388 accepted structural moves per tree per sweep at low noise and 0.0790 at the default cell, one proposal per tree
per sweep making the pooled rate exactly that. At `w = 1, d = 0.16` arm B is
`0.0388 - 0.16 x 0.0166 + 0.16 x 0.9927 x 0.2715 = 0.0793` at low noise, 2.0x, and
`0.0790 - 0.16 x 0.0377 + 0.16 x 0.9754 x 0.4065 = 0.1364` at the default cell, 1.7x - the forgone change term at its measured
acceptance, then the perturb term at `d` times the eligible share (one minus change's no-op rate) times `accept(1)`. What the chain
lacks is worse for a minimum-ESS primary than it was for coverage: every factor above is a FRIEDMAN number, at low noise and at the
default cell, applied to a C1 target the census never ran, and nothing anywhere connects a 2.0x accepted-move rate to a move in the
minimum ESS over 25 fixed points. The design does not predict the bar will be met, and says so here so slice 3 can fail honestly.

**Kill criterion.** This is the design's own, derived from
[6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered) with every departure stated. **KILL
if, at `w = 1, d = 0.16` on the shipped four-chain configuration, arm B does not improve cell 1's summed minimum ESS over arm A by
more than +8, over at least 20 matched pairs, with wall time per sweep no worse than 6.4's 1.05 ratio, and with a mandatory
fresh-seed re-run of any flagged cell before a flag counts.** Five departures. (a) The statistic is minimum ESS, not coverage,
the He-Hahn deficit being read as a mixing symptom (VD, 2026-09-07); at the shipped chain configuration coverage is 0.961 and has
no headroom to win in. (b) The cell is C1, not the low-noise cell 6.4 names, because P1 could not separate three shipped mixtures and
this design claims no win there; quoting 6.4 verbatim would fire the kill by construction. (c) 6.4's second conjunct, arm C against
arm D, is dropped with arm C, so the kill fires on B's failure alone rather than on B's and C's together, and the move can no longer
survive on the conjunct it is not running: STRICTER than 6.4. (d) "four times the measured per-replicate standard error" is read as
four times the twenty-PAIR standard error, derived in 5.1 from the cell's own recorded spread rather than imported from
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
0.010, which is benchmark-surfaces' coverage margin and not this statistic's. (e) 6.4's "KILL the default question
independently" clause needs plateau prediction error in the noise-heavy or large-n stratum; no cell here is either and none measures it, so that clause belongs to
slice 4, not to slice 3.

The asymmetry stands: passing justifies shipping the move opt-in at weight 0. Flipping the default needs the grow-from-root harm
battery, which is not in benchmarks/ and must be reconstructed
([5. Verdict and consequences](grow-from-root-default.md#5-verdict-and-consequences)).

[6.5 Cost, honestly](tree-mixing-proposals.md#65-cost-honestly) repriced. Stage 0's instrumentation and driver are spent, both
landed. The kernel is 135 to 175 lines - 6.5's 120 to 160 plus its separate +15 for the window ratio - across twenty-two files
rather than six. `perturb-balance.R` is 400 to 500 with the confirmation arm, the lower half of 6.5's 400 to 600, narrowed because
the prior-only arm needs no likelihood oracle of its own. The Stage 2 harness stays at 6.5's ~400, and now carries the
sham arm and the P2 quantile-grid arm as well; the four-chain C1 arm it reuses is already landed. What 6.5 did not price is the private-library build any width
arm needs, which the confirmatory run avoids by fixing `w = 1`. Compute is unchanged, on the order of a day: cell 1 is 20 matched
pairs of four-chain fits at 500 + 500, the sham arm is 20 more, and cell 3's re-measurement of arm A is 3. No per-fit time is
quoted; 10.4's host was loaded throughout its run and its times carry no claim.

**Pilot measurement (2026-09-07).** Arm B has now been measured on cell 1 at the pre-registered configuration - `w = 1, d = 0.16`,
the shipped four chains at 500 + 500, Trig+poly - over twenty matched pairs, as a PILOT ahead of slice 3 and not its confirmatory
run ([10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)). Summed minimum ESS moves +0.1
+/- 8.1 against control (t 0.06), against the +8 bar above; RMSE ratio is 1.020 on Trig+poly and 1.030 on Single index. On these
numbers the kill criterion above would fire.

**Verdict (2026-09-07).** What the pilot owed has now run. The sham arm - control (b), `independent75pool4sham`, the control
against itself at sampler seeds offset 1000 on the same twenty data seeds - reads summed minimum ESS 14.82 against 12.51, a
paired difference of -2.3 +/- 9.7 (8 of 20, t -1.07, paired SE 2.17), inside the +8 bar; the bar is four times that paired SE,
which is 8.7, so +8 is marginally optimistic rather than wrong. The P1 control rung re-ran at 60 fits and reads default
0.725 (0.682-0.760) held-out, birthdeath 0.709, swap 0.728, identical to
[10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07),
so the absolute gate stands. The mandatory fresh-seed re-run
([6.1 The rule, stated operationally](benchmark-surfaces.md#61-the-rule-stated-operationally)), seeds 21 to 40, 80 fits:
Trig+poly, the primary, reads summed minimum ESS 16.46 against 14.60, -1.9 +/- 6.5 (9 of 20, t -1.27), per-chain 1.55 against
1.61 (+0.06, t 0.92), coverage -0.004 +/- 0.009 (t -2.00), RMSE ratio 1.042 (t 3.02), held-out RMSE ratio 1.032 (t 2.36), wall
ratio 1.041 - it reproduces the pilot's null and sharpens it negative (the pilot read +0.1 +/- 8.1). Every clause of the kill
criterion is satisfied: twenty matched pairs, twice; wall per sweep at 1.041 sits inside 6.4's 1.05 bar, so the criterion turns
on ESS alone; the sham control passes inside the bar; the fresh-seed re-run reproduces the flagged null; the P1 gate is in
force. Arm B does not improve cell 1's summed minimum ESS by more than +8 at either block, and held-out RMSE on the primary
regresses to 1.032, past
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
1.02 margin with the one-sided bound excluding the null. **KILL, at `w = 1, d = 0.16`.** Residue, recorded not as survival:
Single index, ungated, gains +9.5 summed and +0.37 per chain on the fresh block (19 of 20, t 5.35 and t 6.47), with its own
1.027 RMSE regression (t 6.72) - a gain the pilot's own Single index reading did not show (+1.9, t 0.74 there). Slice 4, a
nonzero default share, is closed by this kill; the kernel stays in the tree at default weight 0 as an opt-in whose measured
benefit on the pre-registered cell is nil. Whether it is removed before release is the maintainer's call and is not decided
here.

## 6. RNG and baselines

The restored dispatch is `if (u < bd) ... else if (u < bd + swap) ... else change`
([`metropolisJumpForTree`](../../src/bartcore/moves.hpp)), so the perturb branch goes in THIRD, at threshold `birthOrDeath + swap +
perturb`, with `changeMove` remaining the `else`. That is what makes `perturb = 0.0` bitwise-neutral: `(bd + swap) + 0.0` is exactly
`bd + swap` in IEEE for any finite sum, so the new test IS the second test, it fails wherever that one failed, and control reaches
`changeMove` at the same stream position. It is the identity the restore already used in the other direction
([9. Reversal: the move returns at default zero](swap-removal.md#9-reversal-the-move-returns-at-default-zero)). Giving change the
threshold and perturb the `else` would instead rest neutrality on the probabilities summing to exactly 1.0 and the uniform never
returning 1.0 - a coincidence, not a construction. The move-type SELECTION consumes one uniform per tree per sweep whichever branch
it takes (the kernels draw more), so at weight 0 no stream moves and benchmarks/R/equivalence.R stays green by construction, which
is what makes the correctness gate runnable before anything changes for users.

**No baseline is re-recorded.** The removal took the one bundled re-record of every baseline and every hardcoded tinytest snapshot,
and the restore, resting on the same identity, took none. Slice 1 takes none either. A nonzero perturb default would be a stream
shift of its own and is slice 4's to pay for. Stage 2's arm A is therefore the mixture that will ship: measuring perturb against a
mixture the release does not contain would price it against a kernel no user will run, which is also why section 5's arithmetic
reads the census's two-move re-run rather than the three-move figures beside it.

## 7. What the census does not settle

- **No proposal correction was priced.** The probe's log ratio carries no window term: exact at an unclipped node, where the
  correction is 0, and off by at most `log 2` at a clipped one at `w = 1`, in a direction it does not report. - **The probe clamped,
  this move clips, and it weights records not nodes.** A clamped step gives the recommended window's targets at `w = 1`, not at `w
  >= 2`. `cutProbe` writes one record per distinct displacement, so an interior node contributes two `|1|` records and a boundary
  node one, and `cutTable`'s per-record mean UNDER-WEIGHTS the boundary about 2x - exactly where the omitted correction applies;
  degenerate-interval nodes leave the denominator entirely, so section 5's eligibility factor misses them. Magnitudes 3, 5, 6 and 7
  come from narrow intervals only - 139 to 1140 probes against ~28000 per power of two, the one figure here taken from the
  THREE-move addendum, the two-move re-run tabulating no count off the powers of two - so acceptance at a narrow interval is
  unmeasured either way. - **The probe took no move**, being a one-step acceptance on the chain the SHIPPED mixture produced; a
  chain running perturb visits different trees, so its realized rate is not this one. Nor is there a depth breakdown, so whether perturb's
  acceptance rises with depth as change's does is unknown. - **Acceptance is not mixing.** Stage 0 measured no coverage, effective
  sample size, inclusion or error, and had no kill criterion. That is section 5's job, and section 5.1's chain-configuration run is
  what turned its primary from coverage into minimum ESS. - **The veto's share of the probes WAS
  separated and not recorded.** [`cutTable`](../../benchmarks/R/move-census.R) computes `vetoed.pct` per magnitude; the 6.1
  addendum's table omits the column. Re-summarizing the existing census files records it. It bears on SLICE 3, not slice 2: it says
  what share of the `accept(1) = 0.2715` the selection rule clears its threshold on is veto rather than likelihood, and slice 2's
  own veto exposure is fixed by occupancy alone under the all-zero mask, which the design makes vacuous by construction. -
  **One grid, one tree count, one chain, no categorical and no coarse cell.** Default `n.cuts`, 75 trees (50 in the causal cell's
  treatment forest), 200 burn plus 500 sampled sweeps, continuous columns throughout - so section 5's dilution arm has no
  pilot, and section 2.2's null-move hazard has none either, though section 5 cell 3 is now where it is piloted, and the price of
  not piloting it first is that the cell's gate runs on a grid 10.1 never measured arm A on.

## 8. Slices

**Prerequisites and order.**

1. **The swap restore, LANDED** (3ecd7f47, recorded at
   [9. Reversal: the move returns at default zero](swap-removal.md#9-reversal-the-move-returns-at-default-zero)), together with the
   structural fill rule it needed (6934e487). Section 3's enumeration counts against that tree and nothing below is startable
   against another.
2. **benchmarks/R/surfaces/P1-friedman.R, the battery's absolute gate, MET** (0.725 held-out in the control arm,
   [10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)).
   Slice 3 re-runs that rung as its own control arm at the time it runs.
3. **The C1 chain configuration, SETTLED** (section 5.1): the benefit cell runs the shipped four chains at 500 + 500 and reads
   summed minimum ESS. This fixes slice 3's primary statistic and its bar.
4. **The census's `vetoed.pct` column**, re-summarized from the existing census files (section 7). It blocks slice 3: the selection
   rule reads `accept(1)` against a 0.25 threshold and this says how much of that acceptance is the veto. Slice 2 does not wait on
   it, its own veto exposure being fixed by occupancy under the all-zero mask.

Then, in order:

1. **The kernel, at default weight 0.** `perturbMove`, a shortened [`changeMove`](../../src/bartcore/moves.hpp), the dispatch branch
   and `StepType` enumerator, `MoveContext`'s fourth probability; ten `chain.hpp` and three `combiner.hpp` sites; eight bridge sites
   including the two-forest refusal; six R default vectors, the new slot, prototype and validity term, `perturb` resolved ahead of
   the three-name fill, `resolveSamplerSpec`'s two `all.equal` comparisons; fifteen positional `MoveContext` initializers in
   tests/cpp; four Rd files; six tinytest files and inst/NEWS.Rd - twenty-two files (section 3.1). Roughly 150 lines of code and 120
   of tests: bitwise neutrality, the spurious-refusal regression under `monotone` and a treatment forest,
   test-sum-to-one-tolerance.R still failing where it fails today, the fill preserving every resolution
   test-proposal-probs.R pins, a perturb-dominant run that changes cuts and never a variable or a shape, and the
   interval-invariance assertion.

   **Landed** (ab49f83a, 2026-09-07). Vector position went the other way from this section's "FIFTH named
   element": the four structural names group together, so [`defaultProposalProbs`](../../R/model.R) reads
   `birth_death 0.6, swap 0, change 0.4, perturb 0, birth 0.5` - perturb fourth, ahead of `birth`, not fifth
   after it. Three further points the design left to the implementation: the all-unnamed branch of
   [`dbartsModel`](../../R/model.R)'s `initialize` method needed a `perturb == 0` guard, since without it
   `c(perturb = 0.16)` alone would fall through to the default instead of tripping the "name at least one
   of" refusal; the residual is grouped `1 - (perturb + sum(named))`, one subtraction, rather than
   `1 - perturb - sum(named)`, two; and test-sum-to-one-tolerance.R gained a four-name `makeFullModel` twin
   beside the original `makeModel` rather than an edit to it, so the three-name regression the file guards
   keeps running unchanged. `fillZeroDefaultProposalProbs` (R/model.R) fills the omitted name for
   `resolveSamplerSpec`'s two `all.equal` comparisons, per the trap above.

   Gates (independent run): tests/cpp 279 including `testPerturbMove` (interval invariance, window count, a
   400-step perturb-dominant walk, the all-categorical no-op), ASAN/UBSAN clean; tinytest 7605/0; equivalence
   trio bitwise 50/12/11 against the fbff1989 baselines, which stand; all 22 exact gates quick and
   hazard-exact full (0.0008 / 0.0005); a perturb-dominant one-tree probe (1947 equal-node-count pairs, zero
   variable-set changes) and an explicit-zero-versus-default `identical()` probe; lint 0; R CMD check
   --as-cran OK, zero notes; NEWS 296 entries; API hash unchanged. Not landed: `perturb-balance.R` (slice 2)
   and the Stage 2 benefit run (slice 3).
2. **`perturb-balance.R`.** The prior-only arm on the full factorial, the exact-posterior confirmation arm, both poisons. Roughly
   400 to 500 lines; not startable before slice 1.

   **Landed** (30472110, 2026-09-07). Four design-versus-code points. (1) Section 4's statistic-3 drop rule spoke in 0-based
   split indices; in the 1-based cut numbering [`getTrees`](../../R/dbarts.R) exposes, the retained six states sit at root cuts
   2 and 3 and the dropped seven at 4 and 5, masses matching exactly (0.011281 / 0.005641 / 0.018802 / 0.009401 retained,
   0.003760 / 0.002820 dropped). (2) The confirmation arm scores the within-variable cut law rather than the root-VARIABLE
   marginal this section named, because under a perturb-dominant mixture that marginal mixes on `change`'s timescale -
   [`firstUnder`](../../benchmarks/R/perturb-balance.R)'s ladder puts its first lag under 0.1 at 113 to 138 kept draws, two
   states never getting there within 400 - and scoring it would be refused by the gate's own burn-in rule. (3) Both poisons
   mutate [`poisonRootTarget`](../../benchmarks/R/perturb-balance.R), statistic 1's target, alone. (4) This section's 20 to 40
   minute estimate for a full run was wrong by two orders: the masked sweep costs the move machinery alone on 24 rows,
   microseconds a sweep.

   Gates (independent run): quick PASS (worst `|z|` 1.84 prior-only, 2.66 confirmation); full PASS in 23 seconds wall (worst
   `|z|` 2.04 prior-only at the pair state `x2c2.1`, 1.64 confirmation), 0 of 21 Holm rejections either arm. Poison 1
   (`logProposalCorrection` dropped) FAILS at worst `|z|` 192.8, 8 of 21 rejected; poison 2 (the one-sided window absorbed at
   the boundary) FAILS at 1143.9. The matching engine mutations, m24 and m25 in `benchmarks/R/mutation-battery.R`, are both
   KILLED (worst `|z|` 129.5 and 703.4). lint 0, air clean, yaml parses. Not landed: the Stage 2 benefit run (slice 3).
3. **The Stage 2 harness and run.** Two arms, the three named cells and the two controls at `w = 1, d = 0.16`, twenty matched pairs,
   reusing the landed four-chain C1 arm. Roughly 400 lines; compute is a day for cell 1 and the sham arm plus the 3 unit fits that
   re-measure cell 3's arm A at `useQuantiles = TRUE` (section 5 cell 3, alternative i). Not startable before prerequisites 2 and 4.
   The verdict is recorded here.

   **Run, KILLED** (d73fb4e0, 2026-09-07). `C1-he-hahn.R` gained the sham arm and a `seedBlock` option; cell 1's move-set arms
   (already landed, ba1fb081) supplied arm B's pilot, and this run supplied the sham arm and the mandatory fresh-seed re-run.
   Sham (control (b), same data seeds, sampler seeds offset 1000, 40 fits): summed minimum ESS 14.82 against 12.51,
   -2.3 +/- 9.7 (8/20, t -1.07), inside the +8 bar. P1 house rung re-run (60 fits): default 0.725 (0.682-0.760) held-out,
   matching [10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)
   digit for digit, gate in force. Fresh block (seeds 21 to 40, 80 fits), Trig+poly primary: summed minimum ESS 16.46 against 14.60, -1.9 +/- 6.5 (9/20,
   t -1.27), wall ratio 1.041, held-out RMSE ratio 1.032 - past 6.4's 1.02 margin, one-sided bound excluding the null;
   reproduces the pilot's null and sharpens it negative. Single index, ungated: +9.5 summed (19/20, t 5.35), RMSE ratio
   1.027 (t 6.72), a gain the pilot did not show. Full verdict at
   [5.3 What arm B must produce, and the kill](#53-what-arm-b-must-produce-and-the-kill).
4. **The default share - CLOSED by the slice 3 kill (2026-09-07).** Slice 3 did not pass, so there is no benefit result to carry
   a nonzero default share on; what follows records the shape the question would have taken, for the record only. It would have
   carried 6.4's second kill clause, which needs plateau prediction error in the noise-heavy or large-n stratum; no cell of
   slice 3 measured that, and the battery that does is the grow-from-root harm battery
   ([5. Verdict and consequences](grow-from-root-default.md#5-verdict-and-consequences)), which benchmarks/ does not contain. The
   kernel stays in the tree at default weight 0, an opt-in whose measured benefit on the pre-registered cell is nil; whether it
   is removed before release is the maintainer's call and is not decided here.
