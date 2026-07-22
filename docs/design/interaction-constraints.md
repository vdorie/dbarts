# Interaction constraints for BART

Status: LANDED f455d7c (2026-07-21). Deliverable B shipped in two verified
phases; see the plan's Landing note. Design synthesized from a blind
three-lens research panel (user-value, engine-feasibility, prior-art) and
hardened by an adversarial code-verified critique (SOUND-WITH-CAVEATS); the
two must-fix findings - the cross-variable swap walk and the state-install
containment gate - were folded in and both shipped. Plan:
docs/plans/interaction-constraints.md. P4 follow-on (variant A / blocks(),
the columnMask warm-start bug fix) LANDED aadbbc8/103dbe2/073d3db
(2026-07-21); see "P4 landing" below and plan
docs/plans/interaction-constraints-p4.md.

## The capability

Let a user restrict which predictors may jointly shape the fit: cap the
interaction ORDER (order 1 = additive / GAM, order 2 = bounded pairwise),
and/or forbid named co-occurrences (a DENY form) or confine interactions to
declared groups (an ALLOW form). Per-forest, so BCF / bartCause can make the
treatment (tau) forest additive-or-low-order while the prognostic (mu)
forest stays free - the calibrated-additivity causal use that motivates
this.

## What NOT to copy, and what to skip

bartMachine is the only BART package with interaction constraints: a
per-path ALLOW-list of variable groups mirroring XGBoost / LightGBM,
rJava-backed, with NO max-order knob. Copying its allow-list as the primary
surface misses the two things users actually reach for - CAP complexity
(max order) and FORBID a named effect-modification (deny) - both of which an
allow-list expresses only by painful enumeration. Heredity (Chipman 1996;
hierNet, RAMP) does not fit: a BART "main effect" is an emergent property of
the whole forest, not a local split decision a per-node rule can gate; it
needs an ensemble-level lattice prior. Name it, defer it.

## The design space

Three axes. SCOPE: per-PATH (limit the distinct variables on one
root-to-leaf path), per-TREE (confine a whole tree to one group, giving the
exact block-additive f = sum_G f_G, a clean functional-ANOVA reading),
per-FOREST (already available via the multi-forest Sampler). HARDNESS: hard
(availability ban) vs soft (prior down-weight). VOCABULARY: max-order,
allow-list, deny-list, grouping.

The unifying primitive: a per-node admissibility predicate p(j | A) over the
ancestor split-variable set A. max-order is |unique(A u {j})| <= K;
co-occurrence is "no forbidden (j, a)"; grouped is membership indirection.
All collapse to ONE predicate. Heredity is the one thing outside it.

## Why this is cheap in THIS engine (the feasibility unlock)

Path-dependent availability is ALREADY the exact norm: splitInterval
(tree.hpp) narrows a variable's cut interval by its ancestors,
variableAvailable (tree.hpp) drops a variable once its interval collapses,
and birth / death / change / swap already stay reversible under that
ancestor-dependent narrowing. A hard interaction limit is the SAME
narrowing, keyed by the ancestor variable SET instead of cut exhaustion.
A per-forest static allow-mask already ships (columnMask_ / columnAllowed in
tree.hpp, installed per tree via setColumnMask for BCF) and funnels through
the one availability predicate. So the primitive above EXTENDS an existing
exact mechanism; it is not new machinery.

## Two distinct deliverables (the critique split these)

A. GROUPED block-additive (per-TREE, static). Confine each tree to one
   declared group; the ensemble is exactly f = sum_G f_G. A static per-tree
   group mask is ancestor-INDEPENDENT, so it reuses the already-shipped
   per-tree columnMask / setColumnMask with NO validity walk, NO swap
   detailed-balance subtlety, and is exactly reproducible. Serves
   grouped-GAMI and "these columns are one module" (a factor's dummies, a
   genomics pathway). It does NOT give a global max-order: order-1 over ALL
   variables needs each tree to ADAPTIVELY pick its single variable, which is
   the per-path K = 1 mechanism, not a fixed mask. Nearly free.

B. MAX-ORDER + co-occurrence (per-PATH, the new mechanism). The unifying
   predicate above, applied along the path. This is what the causal
   tau-additive cap, general complexity control, and the deny form need.
   Extend collectAvailableVariables / variableAvailable / hasAnyAvailable-
   Variable to also read the ancestor variable-set: for max-order, once
   |distinct used| >= K keep only the already-used variables; for
   co-occurrence, drop the forbidden partners of used variables. O(p) added
   to the existing single ancestor walk. Birth / death are
   correct-by-construction: growthProbability, drawBirthableNode, and
   drawSplitVariable all route through this ONE predicate, exactly as ordinal
   cut-exhaustion does; a fully-constrained node is a forced leaf, like a
   depth-vetoed one.

## The structure-move exactness (deliverable B), verified

CHANGE is exact and balanced (confirmed against the change move in
moves.hpp): a redraw at node N recomputes N's suffstats under old and new
variable, treeLogProbability re-derives descendant availability against the
LIVE tree so a max-order subtree change IS captured, availability at N
depends only on N's ancestors (unchanged by the move) so forward
p(v | avail N) and reverse p(old-v | avail N) share a normalizer, and the
original subtree is always a live reverse target.

SWAP is the danger case and my first draft got it wrong. Swap exchanges a
parent's and child's rules, altering the ancestor set on BOTH sides, and it
draws no variable, so there is NO selection-density cancellation - its
correctness rests on proposal symmetry plus treeLogProbability(parent)
carrying the full constrained prior, which needs its own argument (not the
change-move one). Worse, the shipped subtree checks
(categoricalSubtreeIsValid, findGoodOrdinalRules) are PER-VARIABLE - complete
for cut-nesting (swapping x1 / x2 cannot affect x3's cuts) but WRONG for
interaction, which couples DIFFERENT variables. Concrete break: forbid the
pair (x2, x3); root splits x1, left child x2, right child x3, both paths
feasible; swapping root and left makes the root x2, so the x3 SIBLING now has
ancestor x2 and the forbidden pair co-occurs - yet neither swapped variable
is x3, so a per-swapped-variable test misses it. FIX: a whole-subtree
feasibility walk (interactionSubtreeIsValid) carrying a running
ancestor-variable bitset that tests EVERY node's variable, unkeyed to the
swapped pair; a violation scores the -1.0 no-op (pi(T') = 0), the auto-reject
categorical and swap moves already use. No leaf-model, marginal, or reduction
change, so within-host thread / SIMD bitwise invariance holds (bitset control
flow, not a numeric reduction).

## Containment: every split consumer must carry the predicate

"Moves-only" is not enough. Trees also enter through buildFromFlat (setState,
the installForest warm start, saved-tree replay), which validates gauge and
cut-correspondence but NOT interaction feasibility - so a donor grown
UNCONSTRAINED (a natural warm start) would install a constraint-violating
tree, which treeLogProbability then mis-scores (splitVariableLogProbability
returns -log(numAvailable) without checking the node's own variable is
available, so it cannot self-detect the violation). FIX: an
interaction-feasibility gate in buildFromFlat / stateIsValid, or refuse a
warm start whose constraint differs from the target. The grow-from-root
producer (growTreeFromRoot in grow.hpp) is a FIFTH split generator; it routes
availability through collectAvailableVariables, so it is correct-by-
construction ONCE the predicate carries the constraint, but it must be named
and tested. So "treeLogProbability needs no change" is true only because it
DELEGATES to the availability primitives that DO change, and only on trees
the constrained sampler itself produced; the state-install path breaks that
assumption and must be gated.

## Exact-posterior gate (cheaper than monotone)

Small p, one or two cuts each, one forbidden pair or a low K; enumerate the
legal structures to small depth; each structure's weight is the constrained
tree prior (path-availability normalizers) times the ORDINARY conjugate
Gaussian marginal - the leaf model is untouched, so there is NO truncated
integral, unlike monotone (docs/design/monotone.md "Exactness gate"). Match
sampler structure-visit frequencies and fitted means to Monte Carlo error.
Component tests: a brute-force path-availability oracle; a feasibility
invariant checked after every accepted move; a change / swap that strands a
descendant scores 0. DE-RISK FIRST: put the swap sibling-stranding
construction above (forbid (x2, x3); root x1, children x2 / x3; swap root and
left) into the gate and drive it through swap BEFORE building the rest - that
one test catches the likeliest implementation error (the per-variable-walk
bug). Also gate monotone-plus-interaction coexistence: they are orthogonal
seams (leaf-model vs split-selection) and should compose, but that is
untested and must be asserted.

## R surface

An interactions() prior object mirroring dart(): interactions(max.order = 2)
and / or interactions(groups = list(c("age", "bmi")), forbid = list(c("z",
"m"))). Groups by name or index, validated at fit time, layered through
parsePriors. Per-forest via bcf(mu.interactions = , tau.interactions = ).
Safe-over-fast R validation (fail at fit time on unknown names, empty groups,
K < 1).

## Cost, phasing, deferred

- P0 (near-zero engine): expose the existing per-variable splitProbabilities
  as SOFT grouped down-weighting. Honest caveat: this is variable SELECTION,
  not an interaction guarantee, and belongs to the grouped-DART story, not
  this item. Separable; listed only for completeness.
- P1 (deliverable A, cheap): per-TREE grouped block-additive via the shipped
  static per-tree columnMask / setColumnMask. No validity walk, no swap DB,
  exactly reproducible. ~150-250 engine LOC (mostly R surface + per-tree
  group assignment + install), Risk 1. Ship this first IF a consumer wants
  block-additive / grouped modules.
- P2 (deliverable B, the higher-value mechanism): per-PATH max-order +
  co-occurrence, with the whole-subtree swap walk, the swap DB argument, and
  the buildFromFlat feasibility gate. Serves the causal tau cap, complexity
  control, and the deny form. Engine-cheaper than monotone (no leaf
  integral) but ~500-700 LOC realistic once the constraint representation,
  the per-forest install at 3+ chain.hpp sites, the state-install gate, and
  the brute-force enumerable gate (test code, not engine) are counted; Risk
  2. The critique's 300-450 was engine-core-only and optimistic.
- P3: grouped-hard co-occurrence as a thin generalization of P2 (membership
  indirection; +~30 LOC).
- P4 / deferred: SOFT path-dependent penalties (no pi = 0 cancellation,
  entangles the DART update, Risk 3); formal heredity (ensemble lattice
  prior).

## The open design call for VD

Two genuinely different targets at different cost / value:
- Deliverable A (per-TREE grouped block-additive): nearly free, clean
  functional-ANOVA, exactly reproducible - but the user must DECLARE groups,
  and it cannot cap global order or express the causal main-effects-only tau.
- Deliverable B (per-PATH max-order + co-occurrence): the flexible,
  higher-value capability (max-order scalar, deny form, causal tau
  additivity, per-forest) - but the real build with the swap walk, the swap
  DB argument, and the state-install gate.
They coincide only at a declared-singleton set (order 1 within it). If the
motivating consumer wants block-additive modules, ship A first (cheap). If
the headline is max-order / causal additivity / deny, that is B and the real
build. Recommendation: B is the differentiated, most-requested capability, so
target B - but if a concrete grouped-module consumer exists, A is a
low-risk down-payment that shares the R surface.

## P4 landing (2026-07-21)

VD reviewed a design memo and an adversarial critique (both grounded in this
worktree; not checked in - referenced by the plan) and directed building all
three P4 deliverables, overriding the critique's recommendation to defer the
consumerless engine. Plan: docs/plans/interaction-constraints-p4.md.

Variant A landed as blocks(groups, trees.per.group) (commit aadbbc8),
confining each WHOLE tree to one declared group via a static per-tree column
mask (Forest::blockMasks / blockOfTree, installed in chain.hpp's
installBlockMasks) - a fixed, deterministic, zero-RNG tree->group
assignment, distinct from the ADAPTIVE per-path allow-list
interactions(groups=) already gives (see "Two distinct deliverables" A
above, and the confirmed idiom below). The partition contract is
REQUIRE-TOTAL: unlike interactions(groups=), which leaves an un-named
predictor unconstrained, a static mask would silently kill an un-named
predictor (masked out of every tree), so resolveBlocks errors on an
un-named predictor, an overlapping name, or a trees.per.group that is the
wrong length, non-positive, or does not sum to n.trees. Applied per forest
(mu.blocks/tau.blocks on the BCF sampler). F3 (the critique's must-fix):
columnMask_ is one pointer per tree and a BCF tau forest's moderator
restriction already owns it, so a tau.blocks row is INTERSECTED with the
existing moderator mask at install rather than overwriting it - tau.blocks
therefore partitions exactly the moderator set, not the full predictor set
(naming a non-moderator errors).

F1 (columnMaskStateFeasible, commit 103dbe2, follow-up 073d3db): the
EXISTING columnMask (BCF moderators, a column-restricted variance forest)
had no warm-start/setState containment gate - buildFromFlatBelow validates
cut/mask gauge, not columnAllowed - so a donor tree splitting outside the
mask installed silently, and splitVariableLogProbability then mis-scored it
against an availability menu that excluded the very column it used (for
BCF, a corrupted causal decomposition). installForests and stateIsValid now
refuse such a donor before touching any live state (a columnMaskMismatch
warm-start result); the same scan covers blocks()'s blockMasks for free,
since a blocks() forest's effective per-tree mask (its block row,
intersected with any base moderator/variance mask) routes through the
identical check.

Idiom confirmed and documented (man/interactions.Rd): interactions(groups =
<a disjoint TOTAL partition of all predictors>) already yields ADAPTIVE
block-additivity - the root's group propagates down every path via
interactionAllows, so f = sum_G f_G holds at every MCMC state, with each
tree's group chosen per sweep rather than fixed. blocks() adds only (a) a
FIXED per-group tree budget and (b) a static, ancestor-independent mask (a
micro-perf save and a stronger "each tree lives in one block forever"
guarantee) - real but narrow value over the free adaptive idiom, built
because VD wanted the fixed-capacity guarantee, not because the adaptive
path was insufficient.

Doors 2 and 3 ("Cost, phasing, deferred" above) are recorded as SHARPENED
DEFER, per the critique:

- Soft path-dependent penalties DEFER: the hard ban is exact because the
  forbidden variable is simply ABSENT from collectAvailableVariables, so the
  forward draw and the reverse density share a normalizer and cancel
  (pi(T') = 0, the same -1.0 no-op cut-exhaustion already uses). A soft
  penalty keeps the variable in the proposal support, breaking that
  cancellation - change/swap would have to compute and carry a real
  penalized-prior ratio, with swap's no-selection-density problem in play
  again. Worse, it entangles DART: forest.treePrior.splitProbabilities IS
  forest.dart.probabilities, and dart.update redraws that Dirichlet from raw
  post-move split counts, so a penalty would bias the very counts the
  update targets - the two mechanisms would fight over one selection
  density. A path-INDEPENDENT penalty is already expressible as a DART/fixed
  weight with no new door; it is specifically the path-DEPENDENT case that
  cannot fold into DART's ancestor-independent count statistic, so no clean
  factorization rescues it.
- Formal heredity DEFER: a BART "main effect" is an emergent marginal of the
  whole ensemble sum (integrating out the other predictors), not a node, a
  path, or a per-tree property; interactionAllows sees only the current
  path's ancestor set, and trees are swept independently (one at a time,
  chain.hpp), so the per-node predicate has no access to the cross-tree
  marginal heredity conditions on - a category mismatch, not an
  implementation gap. A real fix needs a global lattice prior over the
  variable-subset inclusion lattice, coupling every tree and breaking the
  conditional independence per-tree backfitting relies on - effectively a
  different sampler (new global state, prior, moves, exactness gate), with
  severe mixing risk. The two-stage protocol - fit max.order = 1, then
  groups/max.order = 2 (or a BCF mu-additive/tau-restricted split) on the
  screened set - already delivers the practical content of heredity with
  ZERO new engine, and is the recommended substitute.

Gates re-run at landing: default-off equivalence trio BITWISE (27/27 + bcf +
multinomial, no re-record - blocks() only perturbs the stream when active);
tests/cpp incl. block-additive confinement, ASAN clean; full tinytest incl.
the expanded inst/tinytest/test-blocks.R (block-additivity invariant with
the C3 empty-tree special case, fixed per-block tree counts, same-seed
reproducibility, monotone coexistence, and a warm-start refusal that also
checks the target sampler stays usable afterward).
