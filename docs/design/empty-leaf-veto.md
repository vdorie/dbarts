# The empty-leaf veto: keep and document (investigation, 2026-07-07)

The no-empty-leaf invariant is enforced one way: the branch
log-likelihood returns a penalty for any branch that contains an empty
leaf, so trees with empty leaves are never accepted into the chain
state. The question this note settles is whether the ordinal proposals
should be made occupancy-aware - so the invariant is enforced by
construction, as categorical rules mostly are - or whether the veto
stays and is documented as deliberate. The conclusion is keep and
document; the reasoning follows.

Correction (2026-07-15): the penalty was a large *finite* constant
(-1e7). That is wrong. A valid branch's log-likelihood is unbounded
below - it carries a -0.5 * centeredSumOfSquares / residualVariance term
that grows with the node's observation count and with a small residual
variance - so at scale (a large fit, or a small sigma during sampling) a
legitimate current branch scores below -1e7, the empty-leaf proposal at
-1e7 wins the finite-vs-vetoed comparison, and the empty leaf enters the
chain state. It then fails the occupancy check on export/restore
(stan4bart's createStoredBARTSampler, observed at n=50000). The penalty
is now -HUGE_VAL, which vetoes unconditionally; the analysis below of why
that stays NaN-free is unchanged.

## Where the constant is read

There is exactly one read of the veto value in the tree: the guard in
logLikelihoodForBranch (moves.hpp), `if numObservations() == 0 return
-10000000.0`. Every conjugate move consumes it through that function:
birth/death score the affected branch before and after and take
exp(newLogLikelihood - oldLogLikelihood); change and swap do the same
with exp(yLogL - xLogL). There is no other literal -1e7 in the engine.

## Which move paths can create an empty leaf

- Ordinal birth draws a cut uniformly over the ancestor-constrained cut
  interval (Tree::splitInterval), which is a function of the cut grid
  and the ancestor splits only - not of occupancy. A cut whose low or
  high side holds no observation reaching the node empties a child.
- The ordinal change move draws a new cut over findGoodOrdinalRules'
  interval (logical descendant satisfiability), again occupancy-blind,
  and can empty any leaf below the changed node once observations
  re-route.
- Swap re-routes observations through the swapped subtree; the validity
  walk (ruleIsValid / ordinalRuleIsValid, categoricalSubtreeIsValid)
  checks logical consistency and the categorical gauge, not occupancy.
- Categorical draws are usually but not always occupancy-safe. The
  canonical gauge keeps at least one *reachable* category on each side
  (drawCategoryPattern rejects the two all-same patterns). "Reachable"
  is the ancestor-filtered category mask (Tree::reachableCategories),
  which is not intersected with the categories actually occupied at the
  node. When a split on some *other* variable has thinned a category's
  support to zero at the node, a reachable-but-unoccupied category placed
  alone on one side empties it. The DART port bug recorded in
  core-generalization.md (empty leaves carried until the veto was
  restored) is direct evidence the moves do emit such proposals.

So the veto is the single, uniform backstop for ordinal and categorical
proposals across birth, change, and swap. Death cannot create an empty
leaf (it collapses two children into their non-empty parent).

## Is vetoed-vs-vetoed reachable?

No. Because the chain state maintains the no-empty-leaf invariant, the
current tree is always non-empty, so its branch score is always finite;
a move can only propose *into* a vetoed state, giving a finite-vs-vetoed
comparison, never vetoed-vs-vetoed. The one path that ever produced two
vetoed scores was the original GP-leaf over-cap guardrail, where the
root itself was over-cap and a birth compared one veto against two. That
path was superseded at the GP stage-1 landing: over-cap leaves now
delegate their marginal and their draw to ConstantGaussianLeaf
(model.hpp, "Leaves larger than maxLeafSize score and draw as constant
leaves"), and there is no -1e7 anywhere in model.hpp. Vetoed-vs-vetoed
is therefore unreachable in the shipped engine.

A corollary: since only finite-vs-vetoed comparisons occur, exp of the
difference is well defined in both orientations, so -HUGE_VAL does not
produce a NaN in the move code (a vetoed proposal against a finite
current branch gives exp(-inf) = 0, a clean rejection). -HUGE_VAL is the
correct penalty precisely because a finite one cannot dominate a branch
score that is itself unbounded below; only an infinite penalty vetoes at
every scale. Should a future non-conjugate path reintroduce a
vetoed-vs-vetoed comparison, that path must supply its own finite
resolution rather than lean on this constant.

## Why not make the proposals occupancy-aware

Full removal of the veto means every ordinal proposal is drawn only from
cuts that leave both sides non-empty, categorical draws use occupied
rather than reachable categories, and the MH ratios carry the matching
correction terms. The cost was assessed and exceeds the item's budget:

- Birth: the occupancy cut range is one O(n_leaf) min/max scan, but
  restricting the proposal to it breaks the current cancellation where
  the rule (and split-variable) proposal density equals the prior
  density. The prior (ruleForVariableLogProbability, growthProbability)
  is grid-based and must stay grid-based to preserve the exact target
  posterior; the proposal becomes occupancy-based. That introduces a
  rule-count correction (occupied vs logical), a split-variable
  correction (occupancy-available vs grid-available variables), and an
  occupancy-based redefinition of node birthability that has to be
  applied consistently in both the birth and the reverse death ratio
  (drawBirthableNode, probabilityOfSelectingNodeForBirth,
  birthableNodeExists, probabilityOfBirthStep). Each term is a place a
  subtle posterior error can hide, catchable only by the exact-posterior
  gates after debugging.
- Change and swap re-route through a whole subtree, so occupancy of a
  deep leaf is not a simple interval. The ordinal change can be made
  occupancy-aware as a rejection sampler whose good set depends only on
  the node's fixed segment and its fixed descendants (invariant to the
  node's own rule, so forward and reverse cancel, mirroring the existing
  categorical flow) - but it must re-route and scan per attempt, and the
  categorical change and both swap validity walks must switch from
  reachable to occupied categories.

Taken together this is a 250-400 line, posterior-changing rewrite of the
move kernels touching moves.hpp, model.hpp, and tree.hpp, plus
regeneration of every RNG-locked snapshot - well past the ~200-line
budget and its 1.5x stop threshold, for the sole benefit of replacing
one finite, well-understood, faithfully-ported guard line. The
risk/reward does not favor it.

## Decision

Keep the veto. It is documented here and inline at its single site as a
deliberate, unconditional (-HUGE_VAL) penalty. If a future consumer
needs occupancy-aware ordinal proposals for a reason beyond the
invariant (e.g. mixing), that work should be scheduled on its own, with
the exact-posterior gates as the arbiter.

## The invariant elsewhere: the transactional predictor surface

The empty-leaf invariant this note keeps (no live tree may hold an unoccupied
bottom) is not confined to the move kernels. `Chain::stateIsValid`'s mean
branch has always re-derived it structurally - build a scratch tree per stored
tree, repartition against the sampler's current data, and refuse unless every
bottom is occupied - as the criterion `$setState` and a warm start
(`installForests`) both gate on. `docs/plans/multiforest-predictor-mutation.md`
made the TRANSACTIONAL predictor surface (`$setPredictor` - whole matrix,
column subset, or per-observation - and the cross-sampler per-observation
session) enforce that same criterion rather than a weaker one: a row installs
only if it empties no leaf in any tree of any forest of any chain, exactly
what `stateIsValid` already required of a
restored state. That arc did not invent a new invariant; it closed a gap
between two paths that were supposed to agree and did not.

One asymmetry survived until that arc's S3: the variance forest's branch of
`stateIsValid` checked well-formedness and strict leaf positivity but not
occupancy, so a heteroscedastic sampler's `$setState` could install a variance
state the mutation veto would have refused. S3 (2026-08-12) closed it, adding
the same scratch-build-and-repartition occupancy check to the variance branch;
see docs/design/heteroscedastic.md section 14.

## What counts as empty: the weight law (2026-08-12)

The veto counted leaf MEMBERS. It now counts POSITIVE-WEIGHT members. A zero
weight is ABSENCE, not reweighting - the shipped contract
(`dbartsSampler-class.Rd`, docs/plans/zero-weight-exactness.md,
docs/plans/sigma-df-zero-weights.md: the leaf suffstats multiply by `w` and the
sigma posterior's df counts positive weights only) - so a leaf all of whose rows
carry weight zero enters no likelihood term of the forest that holds it. Under
the count law such a leaf was legal: it scored exactly `0.0`
(`ConstantGaussianLeaf::logIntegratedLikelihood` returns 0 at `sumWeights == 0`)
and drew its parameter from the prior at posterior precision 0, a state no fit
on the positive-weight subset could produce. It is now vetoed, by the same
`-HUGE_VAL` mechanism at the same site. The mechanism, the penalty value and the
finite-vs-vetoed argument above are unchanged; only the predicate moved.

Two sites carry the predicate, because the branch marginal has two owners:

- `logLikelihoodForBranch` (moves.hpp), the conjugate path every ordinal and
  categorical birth/death/change/swap consumes.
- `MonotoneConstantGaussianLeaf::logLikelihoodForBranchWithParams` (model.hpp),
  which owns the constrained joint outright and therefore returns from
  `logLikelihoodForBranch` BEFORE the loop above it runs. It had its own copy of
  the count test; it now takes the same weight law. Monotone directions compose
  with weights on any family (facade.hpp dispatches on the direction vector
  alone), so leaving it behind would have kept one reachable configuration on
  the old law.

Both call `Tree::leafHasNoWeight(i, weights)`: with `weights == nullptr` it IS
`numObservations() == 0`, so the unweighted path - the overwhelmingly common one
- keeps its decision AND its arithmetic bit for bit; with a weight vector it
scans the leaf's members and stops at the first positive weight, so an ordinary
leaf costs one gather and only a leaf about to be vetoed walks its members.

The obvious cheaper candidate, `Node::sumWeights == 0.0` (exact, since a sum of
nonnegatives is exactly zero iff every addend is), was REJECTED on freshness,
not on arithmetic: `Chain::run` refreshes node statistics only
`if constexpr (leafTracksNodeAverages)`, i.e. `!L::hasVectorParams`, so a
linear-leaf chain never calls `setNodeAverages` and a root-only tree there
carries the field at its `0.0` default. Reading it at the veto would have
vetoed every root branch on that path.

### Which weights the predicate sees

`MoveContext::weights` is the weight vector the forest is actually being scored
against, not the user's: the mean forest under a variance forest sees
`w_i / s^2(x_i)`, a BCF forest sees `composeForestWeights`' product of the
observation weight and the per-forest weight, and a latent family sees its
composed working weights. So a zero per-forest weight (`setForestWeights`) now
also vetoes a leaf of only such rows in THAT forest - stated in
`Chain::setForestWeights`' contract - while the veto for the variance forest
reads the user weights it is handed. Weights ship on gaussian and Student-t
only, and the latent families' own working weights are strictly positive
(a zero Polya-Gamma weight is unreachable, and a zero count is refused at
creation), so no shipped latent configuration changes behavior.

### The sites that still count members, and why that is correct

The fix is deliberately confined to the DRAW LAW. Every other occupancy site
keeps the member count, and each is right to (docs/plans/latent-subset-mask.md,
"Semantics of inactive" rule 2, which depends on this and is written against
it):

- The birth scan's `count` (scan.hpp) sentinels a bin's membership; its
  `sumWeights` is a separate field and the marginal it feeds already returns
  `0.0` at `sumWeights == 0`. Occupancy and weight are two facts there, kept
  apart on purpose.
- `Tree::collapseEmptyNodesBelow` merges on `numObservations() == 0` (its
  weighted merge WEIGHT is a weight sum, but the trigger is the count). It runs
  on structure that must be legal after a data or cut-grid change, where a
  member-empty leaf is unrepresentable and a zero-weight one is merely
  uninformative.
- `Tree::bottomNodesAreOccupied` and `Chain::stateIsValid`'s scratch rebuild -
  the transactional predictor surface and the state-restore criterion. These
  answer "is this partition representable against this data", a question about
  membership; a weight-based criterion there would refuse a state the sampler
  itself could have drawn under a different weight vector, since weights do not
  ride the state block.
- `Tree::numObservations` itself, and the chi-k leaf-count gates that read it.

The fix therefore changes which branches are VETOED, not which are CREATED, and
it does NOT align a masked or zero-weighted sampler's occupancy with a compacted
one's; that claim is struck in the subset-mask plan and is not made here.

### Measured effect

- `benchmarks/R/equivalence.R` against `equivalence-a825263.rds`: 34 of 35
  scenarios reproduce BITWISE ("identical draws (same RNG stream)"), including
  every weighted one whose weights are strictly positive (`weighted`,
  `wtoffset`, `wtgp`, `wtlogistic`, `grouped`, `student`, `logistic`). Only
  `zeroweights` moves - 37 summaries, max |z| = 2.85, so the posterior is
  unmoved and the draw law is not. `bcf-equivalence-a825263.rds` and
  `multinomial-equivalence-1027be5.rds` are bitwise on every channel of every
  scenario (no baseline scenario installs a zero weight on those paths).
- `tests/cpp` non-vacuity, measured: driving 4000 moves on a fixture whose
  lowest cut of x0 isolates a zero-weight block, the count law settles on a leaf
  of only zero-weight rows 546 times and the weight law never does.
- tinytest non-vacuity (`test-empty-leaf-veto-weights.R`): with the zero-weight
  half-space `x1 > 0.5`, a 50-tree sampler under the count law leaves live
  leaves that no positive-weight row reaches; under the weight law every leaf is
  reached by one, for gaussian and for Student-t, and the zero-weight rows still
  receive fits.
- Cost: the added work is one gather and compare per leaf per branch score, and
  only when a weight vector is installed; the unweighted path compiles to the
  same count test it ran before.
