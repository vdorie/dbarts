# The empty-leaf veto: keep and document (investigation, 2026-07-07)

The no-empty-leaf invariant is currently enforced one way: the branch
log-likelihood returns a large finite penalty (-1e7) for any branch that
contains an empty leaf, so trees with empty leaves are never accepted
into the chain state. The question this note settles is whether the
ordinal proposals should be made occupancy-aware - so the invariant is
enforced by construction, as categorical rules mostly are - or whether
the veto stays and is documented as deliberate. The conclusion is keep
and document; the reasoning follows.

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
difference is well defined in both orientations, so -inf would no longer
produce a NaN in the move code today. The finite value is nonetheless
kept deliberately - it is self-documenting as a defensive penalty rather
than a true impossibility, and it makes the invariant robust to any
future path (a non-conjugate strategy, a new guardrail) that could
reintroduce a doubly-penalized comparison without having to re-audit
this reasoning. -1e7 is a load-bearing choice, not an accident.

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
deliberate, finite, load-bearing penalty. If a future consumer needs
occupancy-aware ordinal proposals for a reason beyond the invariant
(e.g. mixing), that work should be scheduled on its own, with the
exact-posterior gates as the arbiter.
