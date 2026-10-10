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

There is no literal -1e7, or any other finite veto constant, anywhere in
`src/bartcore`. The predicate a leaf is scored under is
[`Tree::leafIsEmpty`](../../src/bartcore/tree.hpp): the leaf holds no row. It
reads the member count and nothing else - no weight vector, no mask - so it
is the predicate [`Tree::bottomNodesAreOccupied`](../../src/bartcore/tree.hpp)
answers for a state install and the predictor surface, and
"What counts as empty: membership" below is why.

The veto is resolved per branch, not per leaf: `logLikelihoodForBranch`
([`logLikelihoodForBranch`](../../src/bartcore/moves.hpp)) walks the branch's bottom nodes and
records whether any of them is empty. That half is leaf-model independent;
the likelihood half is not. Off a `ParamScoringLeafModel` - the conjugate
leaves - it is the per-leaf marginal summed over the leaves that hold a row.
On one (the monotone constant leaf) the leaf owns the branch
marginal outright and the value returned is
`logLikelihoodForBranchWithParams`
([`logLikelihoodForBranchWithParams`](../../src/bartcore/moves.hpp)) over the whole branch,
with no per-leaf sum running at all. Either way the two return
together as a `BranchScore` ([`BranchScore`](../../src/bartcore/moves.hpp)).
`resolveEmptyLeafVeto`
([`resolveEmptyLeafVeto`](../../src/bartcore/moves.hpp)) then compares a branch's current and
proposed `BranchScore`s: a branch holding an empty leaf, against one holding
none, is assigned `-HUGE_VAL` (not a finite literal), and otherwise the two
finite log-likelihoods stand. Every conjugate move consumes this pair
through `logLikelihoodForBranch` and `resolveEmptyLeafVeto`: birth/death
score the affected branch before and after and take
exp(newLogLikelihood - oldLogLikelihood); change and swap do the same
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

## Is vetoed-vs-vetoed reachable? No

The no-empty-leaf invariant holds of the chain STATE, and the score's
predicate is the state's own: membership rides the tree. Every install that
changes which rows reach a leaf enforces the invariant as it lands -
`setState`/`installForests` and the predictor transaction's revalidation merge
or refuse, `setData`, `setCutPoints` and the donor rebuild merge through
`collapseEmptyNodes` - and the installs that change only what a leaf's rows
WEIGH cannot empty one: `Chain::setWeights`, `setActiveRows`,
`setForestWeights`, `setForestBasis`, the BCF per-sweep zero-multiplier snap
and `formMeanWeights`' `w_i / s^2(x_i)` underflow. So the current branch of a
comparison never carries the penalty, a proposal that does is rejected at a
ratio of exactly 0.0, and two infinite penalties never meet.

From 2026-08-12 to 2026-10-05 that was not so. The veto then counted
positive-weight members, weights do not ride the tree, and every install in
the second list could leave the CURRENT state vetoed. Left as two infinite
penalties, the comparison was `exp(-inf - (-inf)) = NaN`, every comparison
against it false, and the affected branch frozen for the rest of the run
(measured then: a grown forest under an all-zero mask, structurally unchanged
after 300 sweeps, 25 of 25 trees; 766 of 2000 jumps reporting a NaN acceptance
probability). The veto was therefore made a three-level lexicographic rank - a
leaf with no member, a leaf with members but none of positive weight, a leaf a
likelihood term reaches - compared ahead of the likelihood, so that a vetoed
forest kept moving and was absorbed back into the admissible set. With the
middle level gone the rank has nothing to order, and it is gone with it: a
leaf of only zero-weight rows is an ordinary leaf, and the forest an all-zero
vector leaves "at its prior" moves there by the ordinary acceptance.

Stationarity. The target is the CGM prior x marginal restricted to the set S
of trees with no empty leaf and renormalized, which IS the veto's definition.
S is a function of the predictors alone. On S x S the veto decides nothing and
the acceptance is the ordinary one, so the exact-posterior gates apply; from S
no proposal outside S is ever accepted, and no chain starts outside it. Under
an ALL-zero vector the marginal is constant and the kernel is the standard CGM
structure kernel on S - which is what "the forest sits at its prior" means.

The initializer and the moves hold one law.
[`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) draws each tree
from the prior conditioned on S by whole-tree rejection, for every forest and
under every weight vector, an all-zero one included. (While the veto counted
weight an all-zero forest had no admissible tree and took the bare root.)

-HUGE_VAL is the correct penalty because a finite one cannot dominate a branch
score that is itself unbounded below. The one -HUGE_VAL that is NOT the veto -
a constrained leaf model's FEASIBILITY sentinel, an empty monotone cone - can
still meet itself, and `resolveEmptyLeafVeto` rejects that pair explicitly
rather than reporting NaN.

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
bottom) is not confined to the move kernels. The TRANSACTIONAL predictor
surface (`$setPredictor` - whole matrix, column subset, or per-observation -
and the cross-sampler per-observation session) enforces it by refusal: a row
installs only if it empties no leaf in any tree of any forest of any chain. The
forced predictor swap, `setCutPoints` and `setData` enforce it by merging an
emptied bottom into its parent.

Every install of stored trees enforces it by the same merge. A forced
`$setState`, `copy()`, a reload and a warm start build each stored tree against
the sampler's current data and merge a bottom no row reaches into its parent;
an unforced `$setState` declines such a state instead
([install-surface.md](install-surface.md)). Mean
and variance forests alike
([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp),
[`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp)); a tree whose
bottoms are all occupied installs exactly. `Chain::stateIsValid` judges a
stored tree's form, not its occupancy. See docs/design/heteroscedastic.md
section 14 for the variance forest.

## What counts as empty: membership

The veto counts MEMBERS, as 0.9-34 does (`likelihood.cpp`): a leaf is empty
only if no row at all reaches it (dec-B238, 2026-10-05). A row of zero case
weight, a row the active-row mask switches off, a row a per-forest weight or
a zero multiplier leaves weightless in one forest, and a multinomial row with
no trials are all rows of the design. Each occupies its leaf and none enters a
likelihood term. A leaf holding only such rows is legal, contributes nothing
to a branch comparison, and has its value drawn from the prior.

The predicate is [`Tree::leafIsEmpty`](../../src/bartcore/tree.hpp),
`numObservations() == 0`, the test the veto ran before weights existed, so a
single-forest fit that installs no zero weight, no mask and no zero-trial row
keeps its decision AND its arithmetic bit for bit. A two-forest fit does not:
its amplitudes zero rows in one forest with nothing installed, and "Measured
effect" below has what that does to its draws.

### Why: the prior over trees must not read the mask

The mask exists for a larger sampler that redraws which rows are in at every
sweep - principal stratification, a mixture, a compliance class. Such a
sampler alternates a BART sweep given the mask with a draw of the mask given
the fit, from the likelihood alone. That second step is the right conditional
only if the prior over trees does not depend on the mask. A veto that counts
only switched-in rows makes the set of legal trees a function of the mask:
installing a mask can leave the current tree outside the set, and the mask
step's conditional is missing the factor that says so. The combined sampler
then targets neither the model the outer code assumes nor the one with the
restriction written into the prior.

Measured on 2026-10-05 against the exact posterior of such a model, by
enumerating every tree and every mask. Ten rows, one predictor with five
values, two rows each; row i is in the BART component N(f(x_i), 0.5^2) with
probability one half and otherwise in a fixed background N(0.5, 1); one tree
(51 trees x 1024 masks) and two trees (2601 pairs); the mask redrawn every
sweep; 16 million sweeps for one tree and 8 million for two, standard errors
by batch means. Largest difference from the exact joint, in standard errors
and in posterior standard deviations:

                               counting weight     counting members
    membership, one tree       0.066 (319; 0.13)   0.0003 (1.0; 0.001)
    fitted mean, one tree      0.191 (341; 0.27)   0.0006 (0.9; 0.001)
    membership, two trees      0.061 (244; 0.12)   0.0004 (2.3; 0.001)
    fitted mean, two trees     0.172 (258; 0.24)   0.0009 (1.9; 0.001)

- Counting weight, the trees are too small - a single leaf 6.3 percent of the
  time against 1.1 exact, three or more leaves 26 percent against 44 - and the
  fitted means are pulled toward the middle: a region with no switched-in row
  cannot keep its own leaf, so its fit is borrowed from a neighbour. Against
  the joint that writes the restriction into the prior it is further off still
  (membership by up to 0.154).
- An unordered factor is the worst case for it: four levels, subset splits,
  207 trees, where 52 percent of masks leave some level no switched-in row.
  Membership off by up to 0.21 and the fit by 0.86 of a posterior standard
  deviation, against 0.0004 and 0.0011 counting members.
- Zero case weights redrawn every sweep in place of the mask behave as the
  mask does under either rule, and a fixed mask or weight vector is sampled
  exactly under either rule's own target.
- A probit response needs one thing more, the redraw of a reactivated row's
  latent ([The mechanism](active-rows-mask.md#the-mechanism)): counting
  members but keeping stale latents is off by 0.004 in membership and 0.055
  in the fit, with it by 0.0001 and 0.0005.
- At 300 rows, three predictors and 50 trees, the mask redrawn with
  probability one half, the two rules differ by up to 0.044 in a membership
  probability (43 standard errors) and 0.076 in a fitted mean - at the four
  rows with the largest first predictor, all usually switched off - and by
  under 0.005 in membership at 290 of the 300 rows. With the switch-in
  probability rising steeply in the first predictor, so that a whole region
  is mostly off, the fitted means differ by up to 0.111 (0.13 of a posterior
  standard deviation) and counting weight leaves the posterior 9 percent too
  narrow.

`benchmarks/R/mask-redraw-exact.R` is the gate: the one-tree enumeration on a
gaussian, an unordered-factor and a probit arm, in the exact-gates workflow.

### What it costs: a zero weight is not a deleted row

Counting weight, a fixed zero weight was the same as deleting the row.
Counting members it is not: the trees may isolate a region of switched-off
rows, whose fit then comes from the prior and whose split spends a share of
the tree prior. A zero-weight or masked row is IN THE DESIGN AND NOT IN THE
LIKELIHOOD. Measured at 300 rows and 50 trees under one fixed mask, with the
split points, the response scale and the residual prior made equal to those
of a fit on the active rows alone - masked minus deleted, the standard error
of a difference about 0.001:

    mask (active rows)       at active rows            at masked rows
                             largest (in sd)   rms     largest   rms
    random half (150)        0.007 (0.02)      0.001   0.006     0.002
    a region off (169)       0.019 (0.05)      0.002   0.024     0.006

- The residual scale does not move (0.4919 against 0.4922, standard error
  0.0001): its degrees of freedom still count positive-weight rows.
- The posterior standard deviation of the fit rises by 0.6 and 1.9 percent at
  masked rows and not at all at active ones.
- Of about 121 leaves, 0.9 to 1.4 per sweep hold no active row under a fixed
  mask, 1.05 under a mask redrawn at one half and 6.35 under one that
  switches a region off.
- Mixing is not worse: effective sample size per sweep, members over weight,
  has median 0.98 with the mask redrawn and 1.01 and 1.03 with it fixed.

What a zero weight still does is unchanged: the row leaves every sufficient
statistic and every family-level parameter update, the residual variance's
degrees of freedom count positive-weight rows, and its pointwise
log-likelihood is NaN.

### One predicate on every path

Every path that decides whether a branch is legal reads membership:

- The moves, through [`logLikelihoodForBranch`](../../src/bartcore/moves.hpp)
  and [`resolveEmptyLeafVeto`](../../src/bartcore/moves.hpp), for every leaf
  model, the branch-owning constrained ones included. `MoveContext::weights`
  is still the vector the forest is scored against - `w_i / s^2(x_i)` under a
  variance forest, the coupling's product with a per-forest weight, a latent
  family's composed working weights - and it moves a branch's SCORE, never
  whether the branch is one.
- The ordinal cut scan, which reads each side's member count
  ([`scanOrdinalCuts`](../../src/bartcore/scan.hpp)). Grow-from-root asks it
  for the non-missing members of both sides, which is what subsumes the
  ancestor split interval; the rule draw at a nog node holds that interval
  itself and asks for each child's actual members, routed missing rows
  included, which is the moves' own reading
  ([`enumerateNogRuleNeighbourhood`](../../src/bartcore/moves.hpp)).
- The categorical scan, whose enumeration puts a present category on each
  side and so emits no empty child at all
  ([`scanCategoricalPartitions`](../../src/bartcore/scan.hpp)).
- [`growTreeFromRoot`](../../src/bartcore/grow.hpp), through those scans: at a
  node whose rows all carry zero weight every candidate scores 0 and the draw
  is the prior's.
- The prior initializers, mean and variance forests alike
  ([`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp),
  [`Chain::sampleVarianceForestFromPrior`](../../src/bartcore/chain.hpp)): the
  rejection conditions on [`Tree::bottomNodesAreOccupied`](../../src/bartcore/tree.hpp)
  and does not read a weight, so a forest draws the same trees from the same
  seed whatever weights, per-forest weights, basis or mask are installed.
- A multi-forest sampler's per-forest composition. A zero per-forest weight
  ([`Chain::setForestWeights`](../../src/bartcore/chain.hpp)) and a multiplier
  the reparameterization snaps to zero
  ([`AmplitudeForestCombiner::formForestResponse`](../../src/bartcore/combiner.hpp))
  zero a row's precision in that forest and leave it a member there. A
  two-forest construction starts at amplitudes (1, 0, 1), so until the first
  amplitude draw - and for good under a fixed (b0, b1) = (0, 1) - every
  control row is weightless in the treatment forest, which may hold leaves of
  control rows only. "Measured effect" has what that does to an ordinary
  two-forest fit's draws.
- The variance forest's moves, whose leaf counts positive-weight rows in its
  statistic and nowhere else.
- A multinomial zero-trial row, composed into the coupling's mask
  ([`MultinomialForestCombiner::composeEffectiveRows`](../../src/bartcore/combiner.hpp)):
  inactive in every category, a member of its leaf in each.
- `Tree::collapseEmptyNodesBelow`, `Tree::bottomNodesAreOccupied` and the
  chi-k leaf-count gates, which always counted members.

### What a leaf of only switched-off rows does, by leaf model

- Constant: [`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp)
  returns exactly 0 at no weight, and the draw has posterior precision 0, so
  it is the prior's N(0, (scale / k)^2).
- Linear: the formula reaches 0 at no weight only to rounding, so
  [`LinearGaussianLeaf::logIntegratedLikelihoodForNode`](../../src/bartcore/model.hpp)
  returns exactly 0 when the leading entry of U'WU, the members' total weight,
  is not positive. The draw is the ridge alone, the prior on intercept and
  slopes.
- Gaussian process: a leaf with a zero-weight member takes the
  positive-subset paths, which score 0 with no positive member and draw the
  prior process at every member; over the leaf-size cap it is the constant
  leaf.
- Monotone: the leaf owns its branch marginal
  ([`MonotoneConstantGaussianLeaf::logLikelihoodForBranchWithParams`](../../src/bartcore/model.hpp)).
  A weightless leaf's term is the prior mass of its cone given its frozen
  neighbours - finite, and under the "leaf" prior the factor the tree's own
  normalizer divides back out - and its draw is the prior truncated to that
  cone.
- Variance: [`ConstantVarianceLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp)
  returns exactly 0 with no positive-weight row, and the draw is the
  scaled-inverse-chi-squared prior's.

### The rule it replaced (2026-08-12 to 2026-10-05)

For those eight weeks the veto counted positive-weight members: a zero weight
was read as absence from the forest as well as from the likelihood, so that a
weighted fit equalled the fit on the positive-weight rows. The predicate
scanned a leaf's members for a positive weight, the scans' sentinels read the
weight on each side, the prior initializers conditioned each forest on its own
composed vector, and from 2026-08-18 the two failures were ranked apart
("Is vetoed-vs-vetoed reachable? No" above). All of it is gone. What it
measured when it landed still describes the size of the difference between
the two rules on a fixed weight vector: of the equivalence scenarios then
recorded only `zeroweights` moved (37 summaries, max |z| 2.85), and the rank
moved only `maskprobit` and `maskordinal` (0.48 and 0.65).

### Measured effect

- `benchmarks/R/equivalence.R` against `equivalence-36d4ff88.rds`
  (`--strict-coverage`): 52 of 55 scenarios BITWISE. The movers are the three
  that install a zero weight or a mask - `zeroweights` (37 summaries, max |z|
  3.24, one above 3), `maskprobit` (37, 2.43) and `maskordinal` (35, 2.79).
- `bcf-equivalence-d49e2103`: 13 of 15 bitwise on every channel. The two
  that leave their recorded streams are `masked` and `glue_toggle`, which
  holds (b0, b1) at (0, 1) and so every control row weightless in the
  treatment forest for the whole run. That harness's statistic is a
  single-chain one over 40 draws, and reads 12.23 and 6.83 here against 11 to
  21 between two seeds of one build, so it sizes nothing. Over 20 seeds of
  300 draws each the rule moves no summary of either scenario beyond its
  standard error - residual scale 0.194 against 0.195 and 0.197 against
  0.195, leaves per tree within 0.01 in both forests, the treatment effect's
  error and spread at treated and control rows alike. What moves is what the
  rule is about: leaves that no weighted row reaches, none before, now 0.05
  per sweep in `glue_toggle`'s treatment forest (of 31) and 0.25 in
  `masked`'s prognostic forest (of 119).
  `multinomial-equivalence-80b1c8d4`: 11 of 11 bitwise.
- An ordinary two-forest fit, amplitudes drawn, is DRAW-SHIFTING and
  posterior-neutral. Its amplitudes start at (1, 0, 1), and the draws differ
  from the build before exactly when the prior-drawn treatment forest holds a
  leaf with no treated row, which that build never drew, or a first-sweep
  move makes one before b0 is first drawn. Over 60 seeds, `bart()` on
  two-forest data differs in 29 at 50 rows, 7 at 200 and 2 at 1000; a sampler
  started from bare roots, as the BCF harness starts its own, in 3, 0 and 0.
  Bit for bit unchanged are the single-forest fits with no zero weight, mask
  or zero-trial row: 24 seeded fits on each of 19 paths (gaussian, two
  chains, weighted, probit, logistic, Student-t, ordinal, negative binomial,
  DART, linear and GP leaves, a variance forest, monotone, a factor with
  missing values, grow-from-root, a warm start, the five-move mixture), and
  the four seeded snapshot files on the reference build.
- Oracle: `benchmarks/R/mask-redraw-exact.R`. On the build before the change
  its three arms miss the exact joint by 83 to 163 standard errors in quick
  mode (membership by 0.065, 0.208 and 0.058; the fit by 0.19, 0.67 and 0.45);
  after it every one of 44 statistics is within 2.3. With the change but
  without the latent redraw the gaussian and factor arms pass and the probit
  arm fails, membership off by 0.004 at 12 standard errors and the fit by
  0.057 at 28.
- Oracle for a fixed vector: `benchmarks/R/bd-balance.R zeroweight` - the
  enumerable birth/death gate with two adjacent cells zeroed on a grown tree.
  The exact posterior puts 0.155 on partitions holding a leaf of only
  zero-weight rows; the chain is within 0.003 of it in total variation (max
  |z| 1.3 over 8 partitions) and 0.155 from the target restricted to weighted
  leaves, where the build before the change sat.
- tests/cpp: a birth that isolates the zero-weight rows, and one between
  wholly weightless branches, is accepted at prior x transition exactly, on
  the constant, linear, GP and variance leaves
  ([`testZeroWeightLeafContributesNothing`](../../tests/cpp/test_moves.cpp));
  a chain driven under a zero-weight block holds a leaf of only such rows
  after 184 of 4000 moves and never a leaf no row reaches
  ([`testEmptyLeafVetoCountsMembers`](../../tests/cpp/test_moves.cpp)). The
  paths one gaussian forest does not reach have their own: a monotone birth
  into such a leaf and its truncated-prior draw in the first of those tests,
  and the variance forest's prior draw and a coupled sweep under a held zero
  multiplier and under a per-forest weight in
  [`testMembershipAcrossForests`](../../tests/cpp/test_sampler.cpp).
- Cost: none on any path. The weighted path loses the per-leaf scan for a
  positive weight; the unweighted one compiles to the count test it always
  ran.

### Measured occupancy rejection rate (2026-09-06)

This table was taken while the veto counted positive-weight members and
ranked its two failures apart, so the classifier it names is gone
(retired: [`resolveVetoRank`](../../src/bartcore/moves.hpp), whose four call sites are
[`resolveEmptyLeafVeto`](../../src/bartcore/moves.hpp)'s, and
retired: [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp)). Its "occ 2"
column is the veto as it stands; its "occ 1" column is a rejection that no
longer happens, nonzero only in configuration (e).

A scaffold build put namespace-scope counters at the four call sites of the
veto's resolution - birth and death
in [`birthOrDeathMove`](../../src/bartcore/moves.hpp), plus
[`changeMove`](../../src/bartcore/moves.hpp) and
[`swapMove`](../../src/bartcore/moves.hpp) - classifying every scored proposal
by the rank pair its two [`BranchScore`](../../src/bartcore/moves.hpp)s carried
and then by the move's outcome: rejected by the RANK (the proposal's rank
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
- (e) is the only configuration that reached the weight level at all: 135
  of its 507 occupancy rejections, 93 of them at birth.
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
not raise it; counting weight (e) added a second stream of rejections about
a third the size of the member-empty one. The alternative priced under "Why not make the
proposals occupancy-aware" - 250 to 400 lines across moves.hpp,
model.hpp and tree.hpp, plus regeneration of every RNG-locked snapshot -
would be recovering between 0.02 and 3.46 percent of proposals, against
an ordinary MH rejection rate of 47 to 82 percent in the same runs.
