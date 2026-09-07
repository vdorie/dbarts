# The empty-leaf veto: keep and document (investigation, 2026-07-07)

The no-empty-leaf invariant is enforced one way: the branch
log-likelihood returns a penalty for any branch that contains an empty
leaf, so trees with empty leaves are never accepted into the chain
state. The question this note settles is whether the ordinal proposals
should be made occupancy-aware - so the invariant is enforced by
construction, as categorical rules mostly are - or whether the veto
stays and is documented as deliberate. The conclusion is keep and
document; the reasoning follows.

The penalty is -HUGE_VAL, not a finite constant (0.9-34's `likelihood.cpp`
returns -1e7). A finite penalty is unsound: a valid branch's log-likelihood is
unbounded below - it carries a -0.5 * centeredSumOfSquares / residualVariance
term that grows with the node's observation count and with a small residual
variance - so at scale (a large fit, or a small sigma during sampling) a
legitimate current branch scores below -1e7, the empty-leaf proposal at -1e7
wins the finite-vs-vetoed comparison, and the empty leaf enters the chain
state, where it fails the occupancy check on export/restore (stan4bart's
createStoredBARTSampler, at n = 50000). -HUGE_VAL vetoes unconditionally, and
the analysis below shows it stays NaN-free. (2026-07-15)

## Where the constant is read

Correction (2026-09-03): this section described the pre-RANK mechanism
(a finite -1e7 literal read off a member-count predicate). It now
describes the shipped mechanism instead; see "Is vetoed-vs-vetoed
reachable? Yes; the veto is a RANK (2026-08-18)" for the full argument.

There is no literal -1e7, or any other finite veto constant, anywhere in
`src/bartcore`. The predicate a leaf is scored under is
`Tree::leafHasNoWeight` ([`Tree::leafHasNoWeight`](../../src/bartcore/tree.hpp)): with no
weight vector installed it is the member count
(`numObservations() == 0`); with one installed, it scans the leaf's
members and returns true only if none carries positive weight.
`Tree::leafVetoRank` ([`Tree::leafVetoRank`](../../src/bartcore/tree.hpp)) turns that
predicate into the 2/1/0 rank - 2 for a leaf with no member at all
(checked directly, ahead of and independent of the weight scan), 1 for
`leafHasNoWeight`, 0 otherwise.

The rank is resolved per branch, not per leaf: `logLikelihoodForBranch`
([`logLikelihoodForBranch`](../../src/bartcore/moves.hpp)) walks the branch's bottom nodes and
takes the worst `leafVetoRank` among them. The rank half is leaf-model
independent; the likelihood half is not. Off a `ParamScoringLeafModel` -
the conjugate leaves - it is the per-leaf marginal summed over the rank-0
leaves alone. On one (the monotone constant leaf) the leaf owns the branch
marginal outright and the value returned is
`logLikelihoodForBranchWithParams`
([`logLikelihoodForBranchWithParams`](../../src/bartcore/moves.hpp)) over the whole branch,
with no per-leaf sum running at all. Either way rank and likelihood return
together as a `BranchScore` ([`BranchScore`](../../src/bartcore/moves.hpp)). `resolveVetoRank`
([`resolveVetoRank`](../../src/bartcore/moves.hpp)) then compares a branch's current and
proposed `BranchScore`s lexicographically: the worse-ranked side is
assigned `-HUGE_VAL` (not a finite literal), and ranks equal falls back
to the finite log-likelihoods. Every conjugate move consumes this pair
through `logLikelihoodForBranch` and `resolveVetoRank`: birth/death
score the affected branch before and after and take
exp(newLogLikelihood - oldLogLikelihood); change does the same
with exp(yLogL - xLogL).

## Which move paths can create an empty leaf

- Ordinal birth draws a cut uniformly over the ancestor-constrained cut
  interval (Tree::splitInterval), which is a function of the cut grid
  and the ancestor splits only - not of occupancy. A cut whose low or
  high side holds no observation reaching the node empties a child.
- The ordinal change move draws a new cut over findGoodOrdinalRules'
  interval (logical descendant satisfiability), again occupancy-blind,
  and can empty any leaf below the changed node once observations
  re-route.
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
proposals across birth and change. Death cannot create an empty leaf (it
collapses two children into their non-empty parent).

## Is vetoed-vs-vetoed reachable? Yes; the veto is a RANK (2026-08-18)

The no-empty-leaf invariant holds of the chain STATE, not of the SCORE: the
score's predicate reads WEIGHTS, and weights do not ride the tree. Every
install that can zero the vector a grown forest is scored against reaches a
vetoed CURRENT state - `Chain::setWeights`, `setActiveRows`,
`setForestWeights`, `setForestBasis` (the multiplier is a veto weight),
`setState`/`installForests` and the predictor transaction's revalidation (which
enforce the COUNT law, not this one), `setData` and the donor rebuild through
`collapseEmptyNodes` (count law again), and with no install at all the BCF
per-sweep zero-multiplier snap and `formMeanWeights`' `w_i / s^2(x_i)`
underflow. No install-time gate can be total, and the documented, tested
all-zero mask says so outright.

Left as two infinite penalties, that state is not merely unpriced: the
comparison is `exp(-inf - (-inf)) = NaN`, every comparison against it is false,
and the affected branch is frozen for the rest of the run (measured: a grown
forest under an all-zero mask, structurally unchanged after 300 sweeps, 25 of
25 trees; 766 of 2000 jumps reporting a NaN acceptance probability).

The veto is therefore a lexicographic RANK on branches, compared ahead of the
likelihood (`Tree::leafVetoRank`, `moves.hpp`'s `BranchScore` and
`resolveVetoRank`):

    rank(branch) = max over its leaves of
                   2  the leaf holds no member,
                   1  it holds members but no positive-weight member,
                   0  a likelihood term reaches it.

Worse rank loses outright (ratio exactly 0.0, today's double); better rank wins
outright (+inf inside the same product, today's double); EQUAL rank takes
today's arithmetic on the finite parts, where the finite part SKIPS the
vetoed leaves rather than summing marginals for them - the conjugate leaves
return exactly 0.0 there, but a linear or GP leaf does not, and a leaf with
nothing to estimate should contribute nothing. Equal rank at level 1 is what
keeps a vetoed forest moving: the tree mixes under prior x transition at a
constant likelihood, and any move clearing the veto is a rank decrease,
accepted outright. Level 2 stays separate so the MEMBERSHIP law - what
`bottomNodesAreOccupied`, state export/restore and the predictor surface all
require - is never violated from a vetoed state.

Stationarity. The target is the CGM prior x marginal restricted to the
admissible set S and renormalized, which IS the veto's definition. On S x S the
rank decides nothing and the acceptance is the ordinary one, so the
exact-posterior gates apply; from S no proposal outside S is ever accepted, so
S is absorbing and the target is invariant. Outside S the kernel is off-support
and free, and needs only to return: a death never raises a branch's rank (the
parent's members are the children's union), and every equal-rank move has
strictly positive acceptance, so the tree collapses back into S with positive
probability per sweep. Under an ALL-zero vector S is empty and the kernel
degenerates to the standard CGM structure kernel at constant likelihood - which
is what "the forest sits at its prior" means. A host that rewrites the mask
every sweep is not passing through a burn-in but living off support; the
alternative there is a frozen forest, so the ranked kernel still dominates.

The initializer is asymmetric with the moves, and it is the only place that is.
A forest whose composed vector is entirely zero has no conditional law to draw
from, so `Chain::sampleTreesFromPrior` takes the BARE ROOT there - the unique
structure no later weight restore can strand a member-empty leaf in, since
every row sits in its one leaf - while the same state under the MOVES is a
prior draw over structures. Init and mutate disagree about one law,
deliberately: aligning them is a separate decision (no baseline reaches the
branch).

-HUGE_VAL is the correct penalty because a finite one cannot dominate a branch
score that is itself unbounded below. The one -HUGE_VAL that is NOT the veto -
a constrained leaf model's FEASIBILITY sentinel, an empty monotone cone - can
still meet itself, and `resolveVetoRank` rejects that pair explicitly rather
than reporting NaN.

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
- Change re-routes through a whole subtree, so occupancy of a deep leaf is
  not a simple interval. The ordinal change can be made occupancy-aware as
  a rejection sampler whose good set depends only on
  the node's fixed segment and its fixed descendants (invariant to the
  node's own rule, so forward and reverse cancel, mirroring the existing
  categorical flow) - but it must re-route and scan per attempt, and the
  categorical change's validity walk must switch from reachable to
  occupied categories.

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
(`installForests`) both gate on. `docs/plans/archive/multiforest-predictor-mutation.md`
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

The veto counts POSITIVE-WEIGHT members, not merely members (0.9-34 counts
members: `likelihood.cpp`). A zero weight is ABSENCE, not reweighting - the
shipped contract (`dbartsSampler-class.Rd`,
docs/plans/archive/zero-weight-exactness.md,
docs/plans/archive/sigma-df-zero-weights.md: the leaf suffstats multiply by `w`
and the sigma posterior's df counts positive weights only) - so a leaf all of
whose rows carry weight zero enters no likelihood term of the forest that holds
it. Under a count law such a leaf is legal: it scores exactly `0.0`
(`ConstantGaussianLeaf::logIntegratedLikelihood` returns 0 at `sumWeights == 0`)
and draws its parameter from the prior at posterior precision 0, a state no fit
on the positive-weight subset could produce. The weight law vetoes it, by the
same `-HUGE_VAL` mechanism at the same site.

One site carries the predicate: `logLikelihoodForBranch` (moves.hpp) takes the
rank over the branch's leaves for EVERY leaf model, including the branch-owning
constrained ones, whose marginal is then taken over the whole branch rather
than summed per leaf.
`MonotoneConstantGaussianLeaf::logLikelihoodForBranchWithParams` (model.hpp)
therefore keeps only its own feasibility sentinel and no copy of the weight
law. Monotone directions compose with weights on any family (facade.hpp
dispatches on the direction vector alone), so the shared rank is what keeps
that configuration on one law.

The predicate is `Tree::leafHasNoWeight(i, weights)`: with `weights == nullptr`
it IS `numObservations() == 0`, so the unweighted path - the overwhelmingly
common one - keeps its decision AND its arithmetic bit for bit; with a weight
vector it
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
composed working weights. So a zero per-forest weight (`setForestWeights`)
also vetoes a leaf of only such rows in THAT forest - stated in
`Chain::setForestWeights`' contract - while the veto for the variance forest
reads the user weights it is handed. Weights ship on gaussian and Student-t
only, and the latent families' own working weights are strictly positive
(a zero Polya-Gamma weight is unreachable, and a zero count is refused at
creation), so no USER WEIGHT reaches the law on a latent family. The
active-row mask does: it IS a latent family's working weight vector, so an
inactive-only leaf is weight-empty there exactly as a zero-weighted one is on
gaussian (`Chain::setActiveRows`'s contract).

### The sites that still count members, and why that is correct

The weight law is deliberately confined to the DRAW LAW. Every other occupancy
site keeps the member count, and each is right to
(docs/plans/latent-subset-mask.md, "Semantics of inactive" rule 2, which
depends on this and is written against it):

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

The weight law therefore changes which branches are VETOED, not which are
CREATED, and it does NOT align a masked or zero-weighted sampler's occupancy
with a compacted one's.

### Grow-from-root joins the law (2026-08-18)

The birth scan's occupancy sentinel (scan.hpp) used to belong on the list
above: it read the bin `count`, leaving a zero-weight-only side of a cut
scored rather than vetoed. That was wrong for the one caller it has.
`growTreeFromRoot` builds a LIVE forest the exact sweeps then own, so its
children must be legal exactly where the moves' are, and the scan is the only
place that is decided - the sentinel now reads `sumWeights` on each side. The
vector it reads is the caller's composed one, per forest under a coupling
(`growForestFromRoot` composes before the tree loop), so the two initializers
condition on the same thing. `count` stays as the histogram's own census;
off an installed weight vector it and `sumWeights` are the same numbers, so
the unweighted candidate set - every fit that installs no weight - is
unchanged, and with strictly positive weights the two laws agree row for row.
Nothing else moved: the sentinel has no other caller, and no MH move consumes
the scan.

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

### Measured effect of the rank (2026-08-18)

- `benchmarks/R/equivalence.R` against `equivalence-4a42620a.rds`
  (`--strict-coverage`): 40 of 42 scenarios BITWISE. The two movers are
  `maskprobit` and `maskordinal` - the two that install a PARTIAL mask
  mid-chain on a grown forest - at max |z| 0.48 and 0.65 over 37 and 35
  summaries. `bcf-equivalence-6e3b9fb8` (12 scenarios, its `masked` scenario
  included) and `multinomial-equivalence-4d9a3337` (11) are bitwise on every
  channel. Adjudicated with a counter build: a scenario deviates exactly when
  it REALIZES an equal-rank comparison on its RNG path - 2 in `maskprobit`, 1
  in `maskordinal`, 0 in every bitwise scenario measured (`zeroweights`,
  `wtoffset`, `weighted`, and the whole BCF and multinomial harnesses, which
  score 1603 and 3378 vetoed branches between them without one).
- Oracle: `benchmarks/R/bd-balance.R veto` - the enumerable birth/death gate
  run from OUTSIDE S. Two adjacent cells are zeroed on a grown tree that holds
  them as SIBLING leaves, so the collapse repairing it is an equal-rank
  comparison and is the tree's only death. The chain absorbs in 1 sweep, never
  revisits a vetoed partition in 300000 draws, and matches the exact
  restricted posterior at max |z| = 1.8. Against the pre-rank build the same
  arm never re-enters S at all (20000 sweeps).
- tests/cpp: a tree grown under positive weights and then handed an all-zero
  vector takes no NaN acceptance probability and changes structure ~1300 times
  in 2000 jumps (frozen before); stranding 5 of its 10 leaves, it is absorbed
  back into S by sweep 98. At the sampler surface a 25-tree forest masked
  whole keeps moving, exports and restores its state, and a partial mask's 13
  stranded leaves clear within 200 sweeps.

### Measured occupancy rejection rate (2026-09-06)

A scaffold build put namespace-scope counters at the four
[`resolveVetoRank`](../../src/bartcore/moves.hpp) call sites - birth and death
in [`birthOrDeathMove`](../../src/bartcore/moves.hpp), plus
[`changeMove`](../../src/bartcore/moves.hpp) and the swap move, since deleted
(retired: [`swapMove`](../../src/bartcore/moves.hpp)) - classifying every
scored proposal by the rank pair its two
[`BranchScore`](../../src/bartcore/moves.hpp)s carry
([`Tree::leafVetoRank`](../../src/bartcore/tree.hpp) taken over the branch) and
then by the move's outcome: rejected by the RANK (the proposal's rank
strictly worse, so the likelihood ratio is exactly 0.0), rejected by the
ordinary MH ratio at equal rank, or accepted. Proposals that never reach
a score - pi(T') = 0, an unsatisfiable rule draw, or a tree with no
eligible node - are counted separately as no-ops, so the four shares sum
to the proposals made. Every configuration ran the default prior, 200
burn plus 500 sampled sweeps, one chain, one thread, a fixed seed, and
200 trees unless its row says otherwise. The instrumentation was NOT
landed: it was reverted after the run, and with its environment switch
unset it cost nothing measurable (three repetitions each at n = 2000 and
n = 10000, inside run-to-run noise).

Percent of all proposals made, every move type pooled. "occ 2" is a
member-empty proposed leaf, "occ 1" a proposed leaf of only zero-weight
rows.

    configuration          proposals  occ 2  occ 1  MH rej  accept  no-op
    (a) n = 100               140000   3.46   0.00   52.71   35.12   8.70
    (a) n = 500               140000   0.99   0.00   65.33   24.44   9.24
    (a) n = 2000              140000   0.30   0.00   70.76   18.61  10.32
    (a) n = 10000             140000   0.07   0.00   78.55   10.92  10.46
    (b) 5 factors x 8 lv      140000   0.09   0.00   73.94   15.30  10.67
    (c) n.cuts = 10           140000   0.02   0.00   71.97   17.90  10.10
    (d) probit                140000   0.29   0.00   47.51   43.74   8.46
    (e) 20 pct zero weight    140000   0.27   0.10   69.80   19.85   9.99
    (f) n.trees = 50           35000   0.41   0.00   81.84    9.63   8.13

Configuration (a) at n = 2000, by move type (percent of that move's
proposals):

    move     proposals  no-op   scored  occ 2  MH rej  accept
    birth        37364      0    37364   0.76   76.68   22.56
    death        32502      0    32502   0.00   74.05   25.95
    change       56046   4169    51877   0.23   77.53   14.79
    swap         14088  10284     3804   0.09   20.55    6.36
    all         140000  14453   125547   0.30   70.76   18.61

- Death never fires the veto, in any configuration: it collapses two
  children into their non-empty parent.
- The current branch was never itself vetoed in any of these runs -
  nothing installs a weight or a mask mid-chain - so the rank comparison
  only ever REJECTED, and the rank-improving outright accept was not
  exercised.
- (e) is the only configuration that reaches rank 1 at all: 135 of its
  507 occupancy rejections, 93 of them at birth.
- Occupancy is a small share of the REJECTION budget: 6.2 percent of
  rejections at (a) n = 100, 1.5 at n = 500, 0.43 at n = 2000, 0.09 at
  n = 10000, and between 0.03 and 0.62 for (b) through (f).
- Not measured: what the change rate would be had the proposal not been
  ancestor-conditioned. Reading that off needs a second cut draw, a
  re-route and a second score per proposal, which is neither cheap nor
  RNG-neutral. The cheap descriptor instead - over the ordinal change
  proposals the descendant-valid good set averaged 96.9 to 98.4 percent
  of the variable's whole cut grid, so the conditioning removes only a
  thin tail of cuts at these depths.

The measured rate is one to two orders of magnitude below the roughly 17
percent Pratola (2016) reports for ancestor-conditioned change proposals
in his example, the inefficiency Lakshminarayanan, Roy and Teh (2015)
attribute to the CGM sampler. Only the smallest fit clears 1 percent
(3.46 at n = 100), and the rate falls monotonically with the sample -
0.99, 0.30, 0.07 - because the cut grid is fixed at 100 cuts, so a
larger sample fills every bin and an empty side of a cut becomes rare.
The change move, the one Pratola measures, runs at 0.01 to 2.45 percent
of change proposals, below birth in every configuration; the budget
concentrates in birth, 8.85 percent of births at n = 100 down to 0.19 at
n = 10000. Coarsening the grid to 10 cuts (c) drops the pooled rate to
0.02; unordered factor subset splits (b) and a probit response (d) do
not raise it; the weight law (e) adds a rank-1 stream about a third the
size of its rank-2 one. The alternative priced under "Why not make the
proposals occupancy-aware" - 250 to 400 lines across moves.hpp,
model.hpp and tree.hpp, plus regeneration of every RNG-locked snapshot -
would be recovering between 0.02 and 3.46 percent of proposals, against
an ordinary MH rejection rate of 47 to 82 percent in the same runs.
