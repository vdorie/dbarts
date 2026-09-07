# rule_gibbs: an exact draw of the split rule at a nog node

Status: PROPOSED, 2026-09-07.

A fifth tree kernel that replaces the Metropolis change proposal at a nog node - an interior node whose two children are both
leaves - with a draw from the rule's own full conditional. The neighbourhood is closed, the acceptance is identically one, and there
is no reverse count.
[15.3 Cross-lens ranking](tree-mixing-proposals.md#153-cross-lens-ranking) ranked it first of ten mechanisms on a validity argument
and left two questions open; both have since been measured, by the nog probe of
[6.1 Stage 0 - the move census (pilot; no kill criterion)](tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)'s
THIRD 2026-09-07 addendum, which is where every census number below comes from. Nog nodes are 48.5 to 98.3 percent of interior nodes
by cell and take 62.7 to 99.1 percent of change's proposals, so the move is not a corner case; and the incumbent rule is far from the
conditional's mode on the cell that matters, holding probability 0.0017 jointly at `c1` (median rank 26.5 of up to 3000) and 0.034 on
the cut axis alone, against 0.737 and 0.746 at `lownoise`, the opposite pole.

**Premise, not reopened here.** The mixture this counts against is the one perturb landed at:
`birth_death 0.6, swap 0, change 0.4, perturb 0, birth 0.5`
([3. The mixture, and the surface](perturb-move.md#3-the-mixture-and-the-surface)), four structural names and a fifth for the
birth/death split. Every count in section 4 is a count of NEW sites beside those, not a replacement of any.

## 1. What change does at a nog node, and what the probe already priced

[`changeMove`](../../src/bartcore/moves.hpp) picks uniformly among interior nodes, redraws the variable from the prior
([`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp)), draws a cut uniformly over the descendant-valid interval, and
scores one candidate. It accepts 3.77 percent of its proposals at the default cell and 1.66 at low noise, and its rejections are not
close calls - the median log-likelihood difference among rejected change proposals is -62.34 and -143.45. At a nog node the whole
conditional is available for the price of a scan: [`scanOrdinalCuts`](../../src/bartcore/scan.hpp) returns the collapsed marginal for
every cut of one variable in a single pass over the node's members, and at a nog node that marginal is EXACT, the two children being a
two-way partition of those members with no skeleton below to reroute through
([15.2 The two facts the lenses agreed on](tree-mixing-proposals.md#152-the-two-facts-the-lenses-agreed-on)).

[`census::nogProbe`](../../src/bartcore/moves.hpp) already assembles exactly the weights this kernel must draw from, and a temporary
assertion cross-checking its incumbent entry against the kernel's own cached branch score agreed to 3.7e-13 over the census run. The
kernel's job is to reproduce that assembly and draw from it.

## 2. The move

### 2.1 The neighbourhood, and why it is closed

At a nog node `u`, the candidate set is every (available ordinal variable `v`, admissible cut `c`), the cut ranging over
[`Tree::splitInterval`](../../src/bartcore/tree.hpp)'s ancestor-constrained `SI_v(u)` and the variable over
[`Tree::collectAvailableVariables`](../../src/bartcore/tree.hpp)'s availability set at `u`. Both read ANCESTORS only, and at a nog node
[`findGoodOrdinalRules`](../../src/bartcore/moves.hpp) coincides with `splitInterval` because there are no descendants to keep
satisfiable. So the neighbourhood does not depend on the incumbent rule, and neither does the node's member set - a rule change
repartitions those members, it does not change which rows are there. Three consequences, and they are the whole correctness argument:
the candidate set is identical from every state in it; its normalizer is the same before and after the draw; and no proposal count
survives into the acceptance.

The weight of a candidate is the tree posterior restricted to `u`'s rule, every factor that does not read it having cancelled:

    log w(v, c) = scanScore(v, c)                                  the two children's collapsed marginals
                + log P_splitvar(v)                                VARIES only under DART
                - log |SI_v(u)|  (- log 2 if the node routes NAs)  VARIES across variables
                + log(1 - growth(left)) + log(1 - growth(right))   VARIES across candidates

What cancels: `growth(u)` itself, every prior factor at or above `u`, the marginals of every other leaf, and the residual sum of
squares the scan deliberately omits, which is additive over any partition of a fixed member set. What varies: the marginal; the rule
prior `1/|SI_v|`, constant within a variable and different between them, which is the factor
[The gate](change-move-balance.md#the-gate) exists to protect; and the two
[`CGMTreePrior::growthProbability`](../../src/bartcore/model.hpp) terms, which take one of two values per child - `0` when the
candidate leaves that child no available variable at all, `log(1 - base/(1 + depth)^power)` otherwise. The split-variable factor
[`CGMTreePrior::splitVariableLogProbability`](../../src/bartcore/model.hpp) is `-log(numAvailable)` and cancels under the default
prior; under DART it is `log(p_v / total)` and does NOT, so the kernel carries a per-variable `log splitProbabilities[v]` term the
probe omits. The census identity of section 5 is therefore stated at `dart = FALSE`, which is what the census ran.

Draw `(v, c)` from the normalized weights. Acceptance is one; no [`resolveVetoRank`](../../src/bartcore/moves.hpp) comparison, no
snapshot-and-restore, no `logProposalCorrection`. Node selection is uniform over the eligible nog nodes and its reciprocal cancels
because that set is invariant under the move: shape is preserved, so `fillNoGrand`'s set does not move
([`Tree::fillNoGrand`](../../src/bartcore/tree.hpp)), and NO NOG NODE IS AN ANCESTOR OF ANOTHER - a node with an interior descendant is
not nog - so changing `u`'s rule cannot move any other nog node's availability set either.

Three of `changeMove`'s guards drop outright. No mask pool, an ordinal rule allocating no words. No stranding walk. And no interaction
walk: `tree.interactionSubtreeIsValid` exists because a redrawn variable can strand a descendant SPLIT, and a nog node has none, so
`collectAvailableVariables`'s own interaction test at `u` is the whole constraint.

### 2.2 The veto, and the missing direction

`scanOrdinalCuts` scores a cut with a zero-weight side as [`cutScanEmptySentinel`](../../src/bartcore/scan.hpp), `-inf`, so such a
candidate carries weight exactly zero and is excluded. That is the right restriction and not merely a convenient one: the shipped
kernel's veto is LEXICOGRAPHIC on [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp), so from a rank-0 state every rank-worsening
proposal is rejected outright and the chain is reversible with respect to the posterior TRUNCATED to occupancy-admissible trees
([Which move paths can create an empty leaf](empty-leaf-veto.md#which-move-paths-can-create-an-empty-leaf)). The Gibbs draw over the
surviving candidates is the exact full conditional of that truncated target. The incumbent is always among them, its own two children
being rank-0 by the invariant every other site enforces
([`Tree::bottomNodesAreOccupied`](../../src/bartcore/tree.hpp)), so the neighbourhood is never empty.

**Where the sentinel is not enough, and it decides section 5.** The veto is rank-RELATIVE. A weight or mask install can strand a whole
branch at rank 1 - members but no positive weight - and from there the shipped moves compare finite parts and mix under prior x
transition. The scan's sentinel tests WEIGHT, so under an all-zero weight vector every candidate is excluded, the neighbourhood is
empty, and the kernel is inert exactly in the configuration the house's cheap balance gates run in. **A, weight occupancy only**: the
kernel no-ops from a rank-1 state; free, correct, and it costs the prior-only gate arm. **B, a rank-aware scan**: give
`scanOrdinalCuts` an optional per-candidate rank out-parameter, defaulted null so [`growTreeFromRoot`](../../src/bartcore/grow.hpp)'s
call is untouched and bitwise unchanged, and let the kernel keep the candidates whose rank is no worse than the current branch's.
About fifteen lines. **RECOMMEND B**: it makes the kernel's veto semantics identical to `changeMove`'s rather than merely equal on the
common case, and it is the only thing that lets section 5's prior-only arm run at all.

**The missing direction is part of the candidate exactly when the node routes a missing row.** `scanOrdinalCuts` returns `2 * numCuts`
entries then, scoring the direction rather than leaving it to a coin, and the rule prior widens by the same factor two the candidate
count does - which is what the probe's `doubled` branch does. Where the column declares missing values but the NODE holds none, both
directions score identically, so enumerating the cut once at prior `1/|SI_v|` is the exact collapsed weight and the direction is drawn
from its own conditional, a fair coin, after the cut. Reproducing this convention is what puts the kernel on the probe's numbers.

### 2.3 A categorical variable at the node

**A, skip the node** - require every available variable ordinal, else no-op or fall back to a change proposal at that node. Simple,
and it is what the probe does; cost, on any mixed design a nog node is usually ineligible and the move is close to inert.
**B, enumerate the categorical partitions too.** Refused, and not on cost:
[`scanCategoricalPartitions`](../../src/bartcore/scan.hpp) is marked INIT-ONLY and says why in the file - above the cap it is
truncated to a sorted-prefix family and REWEIGHTED so the family carries the whole variable's prior mass, and below it the enumeration
is over the partitions of the categories PRESENT at the node with the absent reachable positions filled by post-draw coins. Neither is
the conditional over the `2^R - 2` masks [`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp) normalizes, so
drawing from it targets the wrong law. **C, a mixed rule: enumerate the ORDINAL available variables and let the move act only when the
incumbent rule is itself ordinal.** The restricted candidate set is still ancestor-determined, so the draw is a Gibbs step on the rule
conditional on the rule being ordinal, and a node whose rule is categorical is a fixed point of this component - a valid restriction,
not an approximation. **RECOMMEND C.** It keeps the move alive on the mixed designs A gives up on, at the cost of never proposing a
categorical rule; `changeMove` still does, and
[`drawCategoricalRuleFromPrior`](../../src/bartcore/moves.hpp) remains the only mechanism that reaches a categorical rule at all. On
an all-categorical design the move is inert, exactly as perturb is.

**At every non-nog interior node the ordinary change move stays**, unchanged and at its own share. Nothing here proposes below the nog
frontier: a moved cut there reroutes members through a fixed skeleton and no prefix scan exists for it.

### 2.4 Cost, and the two restricted variants

One scan per available ordinal variable over the nog node's members, against `changeMove`'s single pass over the same members. In
cut-scan units - one `scanOrdinalCuts` pass over a node's members for one variable, so a full pass over `n` is about `L` units at `L`
leaves per tree - a nog node holds two leaves' worth of members, about `2n/L`, so a proposal costs `2 p_avail` units against change's
2. Per sweep, at the full change share (`m = 75`, one proposal per tree per sweep, 0.4 of them change) and the census's own
leaves-per-tree and target-nog shares:

    cell      p    L     target-nog%   units/sweep at d = 0.4   a sweep's own traffic
    default   10   2.83     71.9              431                       637
    lownoise  10   3.79     62.7              376                       853
    wide      50   2.53     76.5             2295                       569
    c1        30   2.52     77.3             1391                       567

The last column is `3 m L`, three full passes over `n` per tree for the residual, the leaf statistics and the fit rebuild - an
estimate, not a measurement. **At C1 the move costs about 1390 cut-scan units a sweep at the full change share and 557 at
`d = 0.16`, so it roughly doubles a sweep at the dosage section 6 runs.** `bcf` is absent because its treatment forest refuses a
non-default `proposal.probs` outright.

**Two restricted variants, both priced by the same probe.** *Cut-only Gibbs* holds the incumbent variable and enumerates its cuts: ONE
scan, the pass change already makes, acceptance still one - the variable is a deterministic function of the state and the restricted
set is again ancestor-determined. *Randomised variable* draws `v` from the split-variable prior as `changeMove` does and then draws
the cut from its exact conditional. **It is NOT a Gibbs step and the brief's "one scan" is wrong.** Writing `R_v` for the sum of
`w(v, .)` over `v`'s cuts with the variable prior divided out, the forward density is `P(v') w(v',c')/Z_{v'}` and the reverse
`P(v) w(v,c)/Z_v`, so the variable prior cancels and `alpha = min(1, R_{v'} / R_v)` - the incumbent variable's scan sum is required,
which is a SECOND scan. Two scans and a Hastings term, against thirty scans and none.

**The entropy columns say the variable axis is nearly free of content in four cells of five.** Joint against cut-restricted median
weight entropy is 0.911 / 0.894 at `default`, 0.334 / 0.328 at `lownoise`, 1.066 / 1.059 at `wide` - at `p = 50` the fifty-fold cost
buys 0.007 nats - and 1.659 / 1.579 at `bcf f1`. Only `c1` separates: 6.425 against 3.426, `P(incumbent)` 0.0017 against 0.034.
**RECOMMEND building the FULL rule draw first**, because `c1` IS the benefit cell of section 6 and is the one place the extra
`p_avail`-fold price has anything to buy, because the cut-only variant is the same code with the variable loop deleted and can be
taken on a private `-D` build the way perturb's width arms are, and because the census identity of section 5 exists only for the joint
neighbourhood. The randomised-variable middle is recorded as a door and not built: it needs its own correction and its own balance
gate for a saving the cut-only variant already has without either.

## 3. Leaf models and families

The scan is templated on [`ScalarLeafModel`](../../src/bartcore/model.hpp) and the probe gates on
[`ScannableLeafModel`](../../src/bartcore/moves.hpp), a scalar leaf carrying the four-argument `(k, sigma^2, sum w, sum wz)` marginal.
That admits [`ConstantGaussianLeaf`](../../src/bartcore/model.hpp) and nothing else the engine ships, so the move is available for
every response family whose mean forest carries a constant leaf: gaussian directly, the latent families (probit, logistic, ordinal,
nbinom, multinomial, aft, hazard) through the working response and weights the scan reads exactly as `computeLeafStats` does, and
heteroscedastic through the composed weights.

Three fallbacks, all of which `changeMove` serves today and must keep serving - it is templated on `MoveScorableLeafModel` and runs
for all of them. [`LinearGaussianLeaf`](../../src/bartcore/model.hpp) is a vector leaf and
[`GPGaussianLeaf`](../../src/bartcore/model.hpp) a function leaf; neither carries the scan's scalar marginal.
[`ConstantVarianceLeaf`](../../src/bartcore/model.hpp) is a `ScaleLeafModel` and NOT a `ScalarLeafModel`, so the variance forest cannot
run this move - a difference from perturb, which reaches [`sweepVarianceForest`](../../src/bartcore/chain.hpp) free. The variance
forest reads the mean forest's mixture, so at a nonzero share it would spend that share on a no-op; at the shipped default of zero
nothing is spent, and folding the share into change for the variance forest alone is a door, not this design's business.

**The monotone leaf needs an explicit guard the probe does not have.**
[`MonotoneConstantGaussianLeaf`](../../src/bartcore/model.hpp) satisfies `ScalarLeafModel` - it forwards the constant leaf's marginal
- and therefore satisfies `ScannableLeafModel`, but it is also a [`ParamScoringLeafModel`](../../src/bartcore/model.hpp), so the
branch score is the constrained joint over the touched leaves given frozen neighbours and NOT the sum of two unconstrained scan
entries. The kernel's predicate must be `ScannableLeafModel<L> && !ParamScoringLeafModel<L>`. In practice
[`resolveSamplerSpec`](../../R/spec.R) rewrites a monotone fit to birth/death-only before this can fire, which is why the probe never
saw it; the guard is what makes that a convenience rather than the correctness argument.

## 4. The mixture, and the surface

**The name.** `gibbs` alone is refused: a user reads it as the whole sampler. `rule` under-specifies against `change`, which also
redraws the rule. `informed_change` is the literature's word for the locally-balanced APPROXIMATIONS this move does not need, and
naming an exact draw "informed" misleads. **RECOMMEND `rule_gibbs`**, the qualifier ruling out the whole-sampler reading and matching
the survey's own name for the mechanism; `informed_change` is the runner-up.

**A, a new structural name at default zero**, additive, bitwise neutral, entering the dispatch FOURTH at threshold
`birthOrDeath + swap + perturb + rule_gibbs` with `changeMove` still the `else`, so the added test IS the perturb test in IEEE at a
share of exactly zero (section 7). **B, an option flag making the existing change move take the Gibbs draw at nog nodes.** Cheaper in
surface, but it changes the DEFAULT kernel's stream and takes a bundled re-record of every baseline and every hardcoded snapshot, and
it fuses two kernels under one name so section 6's arms cannot be a mixture contrast at all. **C, both.** Pays B's re-record for A's
flexibility. **RECOMMEND A.** The share comes from change, as perturb's does, the proposal count per tree per sweep staying one.
**Default share: 0**, which satisfies
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
absolute null-control gate by construction.

### 4.1 Every site, enumerated

The pattern is [2. What is removed](swap-removal.md#2-what-is-removed) and
[3. The surface](swap-removal.md#3-the-surface) read as an addition, and
[3.1 Every site, enumerated](perturb-move.md#31-every-site-enumerated) counted the same tree one name ago.

**R, five files, six default vectors.** [`defaultProposalProbs`](../../R/model.R); the [`dbarts`](../../R/dbarts.R), `bart2`
(R/bart.R) and [`dbartsSpec`](../../R/spec.R) formals; and the two literals in the monotone branch of
[`resolveSamplerSpec`](../../R/spec.R), the comparison default and the birth/death-only rewrite. Both `all.equal` comparisons - the
monotone one and the treatment-forest one - route through
[`fillZeroDefaultProposalProbs`](../../R/model.R), which gains the new name or every caller passing the documented default is refused
spuriously. [`dbartsModel`](../../R/A_class.R) gains a `p.rule_gibbs` slot, a prototype of 0 and a fifth term in validity's sum.
**The fill rule takes a SECOND zero-default name, and the rule perturb landed generalizes to it verbatim**: `rule_gibbs` resolves
ahead of the three-name fill and never enters it, the residual becoming `1 - (perturb + rule_gibbs + sum(named))` as one subtraction;
the all-unnamed guard becomes `perturb == 0 && rule_gibbs == 0`; the frozen test gains the same term; and the refusal
(["name at least one of"](../../inst/tinytest/test-proposal-probs.R)) now names three zero-default moves rather than two.

**C++, five files.** In src/bartcore/moves.hpp: the kernel beside [`perturbMove`](../../src/bartcore/moves.hpp), a fifth probability
on [`MoveContext`](../../src/bartcore/moves.hpp), a sixth enumerator on [`StepType`](../../src/bartcore/moves.hpp), the dispatch
branch, a fifth argument to [`structureIsFrozen`](../../src/bartcore/moves.hpp), and the census hooks and legend the other kernels
carry. src/bartcore/scan.hpp takes section 2.2's rank out-parameter. Then TWELVE in
[`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp) - three struct fields, four copies into a
forest, the variance forest's copy, the two `MoveContext` initializers and the two `structureIsFrozen` call sites - and THREE in
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp), without which the move is unreachable
from BCF and multinomial fits and silently zero. EIGHT in
[`parseModel`, `refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp): the parsed field, the slot read, the sum
check's fifth term, the creation printout, the `SamplerOptions` copy, the two-forest refusal's hard-coded mixture, the forest spec
copy and `setModel`'s copy. Twenty-three probability sites outside the kernel file.

**tests/cpp, two files: SEVENTEEN positional `MoveContext` initializers**, fourteen in tests/cpp/test_moves.cpp and three in
tests/cpp/test_interaction.cpp. The probability block sits ahead of `const double* weights`, so a short initializer binds that pointer
to a `double` - a hard compile error, which is what makes the count safe. **Four Rd files** (man/dbarts.Rd, man/bart2.Rd,
man/dbartsSpec.Rd usage and argument text; man/bart.Rd argument text alone). **SIX tinytest files**: test-proposal-probs.R,
test-argument-surface.R's pinned default, test-bcf-creation.R's refusal literal, test-spec.R and test-monotone.R's slot reads, and
test-sum-to-one-tolerance.R, whose three-name regression must keep failing where it fails today. **Plus inst/NEWS.Rd.**
**Twenty-three files.** Roughly 200 lines of kernel and scan, 60 across the surface sites, 150 of tests.

**The flat C header does not move** and no stored state does:
[`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) takes the model as a `SEXP`, so no proposal
probability crosses the ABI and no `LinkingTo` consumer recompiles.

## 5. Correctness: `rule-gibbs-balance.R`

Two arms, and unlike perturb's the second is mandatory: the prior-only arm exercises the node selection, the enumeration and the prior
factors but NOT one scan entry, and the scan-weighted draw is the whole new mechanism.

**The prior-only arm and its target.** Under an all-zero weight mask every leaf holding rows is rank 1, the rank-0 set is empty, the
likelihood difference is exactly 0 and - given section 2.2's rank-aware neighbourhood - every candidate's weight collapses to its
prior factors alone. The kernel is then reversible with respect to the CGM prior TRUNCATED to member-occupied trees and renormalized,
the same target [4. Correctness: perturb-balance.R](perturb-move.md#4-correctness-perturb-balancer) names. **The design makes the
truncation vacuous, and is chosen to make both poisons resolvable**: two ordinal columns as a FULL FACTORIAL, x1 on 6 distinct values
and x2 on 2, at least one row per cell of the 12, `useQuantiles = TRUE` so the induced grids are 5 cuts and 1, and
`tree.prior = cgm(0.95, 0.5)`. The lopsided cut counts are what give poison (ii) its size; the low `power` is what gives poison (i)
its size, deepening trees so a candidate can exhaust a child.

**Three statistics.** (1) The root's (variable, cut) marginal plus the stump, SEVEN states, closed form `P(grow) x P(v) x 1/|SI_v|`:
0.095 for each of x1's five cuts, 0.475 for x2's one, 0.05 for the stump. (2) The leaf count, from a dynamic program whose state is
(remaining x1 cuts, remaining x2 cuts, depth), `growthProbability` being depth-dependent; support 1 to 12. (3) The left child's cut
CONDITIONAL on the root splitting x2 and that child being itself a nog node, FIVE states - the child then holds only x1, the
conditioning makes the closed form exact, and it is the one place the `1 - growth` terms vary. All three draw the split variable
uniformly among the variables that STILL have a cut, which is what `drawSplitVariable` does with no `splitProbabilities` set.

**Two poisons, both sized in advance.** (i) **Drop the two `log(1 - growth(child))` terms.** At the depth-1 left child, splitting x1 at
either end leaves a grandchild with no available variable and a factor of 1, while a middle cut leaves both splittable at
`1 - 0.95/3^0.5 = 0.4515` each: the true law over statistic 3's five states is (0.2981, 0.1346, 0.1346, 0.1346, 0.2981) and the
poisoned one is uniform at 0.2 - the ends fall 33 percent, the middles rise 49. Statistic 1 is BLIND to this poison by construction,
the root's children always retaining the other variable, which is why statistic 3 exists. (ii) **Drop the `1/|SI_v|` rule factor.**
The root's six candidates then carry equal weight: x1's cuts go 0.095 to 0.158333 (+67 percent) and x2's single cut 0.475 to 0.158333
(-67 percent), the same low-cardinality bias in mirror image that
[The gate](change-move-balance.md#the-gate) repaired. Both effects are diluted by the arm's other shares; the arm is
`birth_death 0.10, change 0.10, rule_gibbs 0.80`, change retained because statistic 1 needs a mechanism at a non-nog root and
birth/death because statistic 2 moves through nothing else. Batch-means z per retained state, Holm at family alpha 0.05, and the run
length taken from the same ladder perturb's gate uses ([`firstUnder`](../../benchmarks/R/perturb-balance.R)); the masked sweep costs
the move machinery alone on 12 rows, so length is nearly free and statistic 3's conditioning, not the run, is what has to be paid for.

**The confirmation arm** runs the same grid with positive weights and no mask, scored against `change-balance.R`'s region dynamic
program ([The gate](change-move-balance.md#the-gate)), and it is what tests the scan entries themselves.

**The census identity, as a tests/cpp assertion rather than a script.** Assemble the neighbourhood twice at a fixed tree: once by the
kernel, once by a reference that installs each candidate rule and calls
[`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp) and `logLikelihoodForBranch` directly. Require agreement to 1e-12
on every candidate. That is stronger than agreeing with `census::nogProbe`, which the census build's own 3.7e-13 already records, and
it does not need the census build to run.

## 6. Benefit, pre-registered

**Prerequisites, both outside this design and both met.** benchmarks/R/surfaces/P1-friedman.R's 90 percent coverage must return near
0.71 in the control arm before any verdict is valid; it reads 0.725
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)).
And the chain configuration is settled: the shipped four chains at 500 + 500 on C1's independent design at 75 trees, which is what
turns the primary from coverage into minimum ESS
([5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic)).

**Two arms, matched seeds, twenty pairs.** **A**, the shipped mixture, which reads 15(8-31) summed minimum ESS on Trig+poly. **B**,
`rule_gibbs d` taken from change. **Dosage: `d` in `{0.16, 0.32}`.** `d = 0.40` is excluded, and the reason is sharper here than it
was for perturb: change lands on a nog node 77.3 percent of the time at `c1`, so a full share does not merely make change the Gibbs
draw - it leaves the 34.8 percent of `c1`'s interior nodes that are NOT nog with no mechanism that moves their rule at all. **The primary is the summed minimum ESS over C1's 25 fixed points at the +8 bar**, the bar
[5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic) derives from the cell's own spread and which
the eight move-set arms already recorded on that cell calibrate at a paired standard error of 1.8 to 2.5. Four secondaries at
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s own margins, all must-not-degrade: coverage at
-0.010 absolute against arm A's 0.961, held-out RMSE at a ratio above 1.02, interval length reported, and cost as below.

**Controls, as [5.3 What arm B must produce, and the kill](perturb-move.md#53-what-arm-b-must-produce-and-the-kill) has them.** The
bitwise null at `d = 0`, which slice 1 runs and which is 6.4's own first absolute gate. A sham arm, A against A at twenty FRESH
sampler seeds, whose paired difference must sit inside the +8 bar. And P1's rung re-run as the absolute control at the time slice 3
runs. A flagged cell takes a mandatory fresh-seed re-run before the flag counts
([6.1 The rule, stated operationally](benchmark-surfaces.md#61-the-rule-stated-operationally)).

**Cost, and why it is not wall time.** 6.4's "wall time per sweep at a ratio above 1.05" cannot be met and it would be dishonest to
pretend otherwise: section 2.4 prices the move at roughly a doubled sweep at `d = 0.16` on `c1`, and 10.4's host carried a load of 9
to 67 throughout its runs, so ESS per second measured there would be noise
([10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)). **RECOMMEND measuring cost as the
CUT-SCAN COUNT**, a deterministic count the census build already produces - nog proposals times available ordinal variables - reported
per sweep beside the primary; and, as the cost-adjusted secondary, ONE equal-cost arm on a quiet machine in which arm A is given the
sweep count that ratio buys it. The draws-based primary answers whether the kernel mixes better; the equal-cost arm answers whether it
should ever be a default, which is slice 4's question anyway.

**The rooting-lock probe.** P2's duplicate-column cell at `m = 1`, must-not-degrade, and it is where the move's reach can be stated
exactly. The shipped mixture parks 5 of 40 chains on an x3 root, change being vetoed once a child splits on x1
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)). At `m = 1` the ROOT IS NOG ONLY
WHEN THE TREE HAS ONE SPLIT. So the Gibbs draw resolves the rooting exactly and in one step whenever the tree is down to two leaves -
it draws the root's variable from its full conditional rather than proposing one - and it can do NOTHING for a parked deeper tree,
where the root is not nog and the move cannot reach it. It does not rotate a rule up the tree; swap remains the only move that does.
The threshold is arm B's mean switches per chain inside arm A's seed range, chains parked no more than 10 of 40, and pooled
`p(root on x1)` within Monte Carlo error of 0.5.

**Kill criterion. KILL if, at `d = 0.16` on the shipped four-chain configuration, arm B does not improve C1's summed minimum ESS over
arm A by more than +8, over at least 20 matched pairs, with a mandatory fresh-seed re-run of any flagged cell before a flag counts.**
Departures from [6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered), each stated: the
statistic is minimum ESS rather than coverage, which has no headroom at 0.961; the cell is C1 rather than the low-noise cell, three
shipped mixtures being indistinguishable on P1; the wall-time conjunct is replaced by the cut-scan count above, which is a LOOSER
gate on this move than 6.4's and is why the equal-cost arm is mandatory rather than optional; and 6.4's plateau-error clause needs a
stratum no cell here measures and belongs to slice 4.

**What the design does not predict.** Perturb's pilot moved this exact statistic by +0.1 +/- 8.1, and the frozen-structure run says
the worst C1 point's minimum ESS rises only from 1.6 to 4-21 with structure held fixed entirely - so most of that deficit survives any
structural kernel. Nothing connects a wider rule neighbourhood to a move in the minimum ESS over 25 points, and the entropy the probe
measured is a one-step statistic on a chain this move never ran. Slice 3 can fail honestly.

## 7. RNG and baselines

At a share of exactly zero the added dispatch test is `u < (bd + swap + perturb) + 0.0`, which is the perturb test exactly in IEEE for
any finite sum, so it fails wherever that one failed and control reaches `changeMove` at the same stream position - the identity the
swap restore and the perturb landing both used
([9. Reversal: the move returns at default zero](swap-removal.md#9-reversal-the-move-returns-at-default-zero)). The move-type
selection consumes one uniform per tree per sweep whichever branch it takes. Section 2.2's scan out-parameter is defaulted null, so
`growTreeFromRoot` passes nothing and the grow path is byte-identical. **No baseline is re-recorded**; benchmarks/R/equivalence.R
stays green by construction, which is what makes the correctness gate runnable before anything changes for users. Section 6's arms
consume the stream differently, which is fine off the default, and a nonzero default is slice 4's stream shift to pay for.

## 8. Slices

**Prerequisites and order.** (1) The perturb kernel at weight zero, LANDED (ab49f83a), together with the structural fill rule
(6934e487): section 4 counts against that tree and the fill generalizes rather than being rewritten. (2) The nog probe, LANDED and
MEASURED (6.1's third 2026-09-07 addendum): it supplies the weights the kernel must reproduce, the nog and target-nog shares section
2.4 prices off, and the entropy columns section 2.4's variant recommendation rests on. (3) P1's absolute gate, MET. (4) The C1
four-chain configuration, SETTLED. Nothing waits on perturb's slice 3.

1. **The kernel, at default weight 0.** The move, section 2.2's rank-aware scan parameter, the dispatch branch and enumerator,
   `MoveContext`'s fifth probability, `structureIsFrozen`'s fifth argument; twelve chain.hpp and three combiner.hpp sites; eight
   bridge sites; six R default vectors, the slot, prototype and validity term, `rule_gibbs` resolved ahead of the three-name fill,
   both `all.equal` comparisons; seventeen positional initializers in tests/cpp; four Rd, six tinytest and NEWS - twenty-three files.
   Tests: bitwise neutrality; the census identity of section 5 against a reference assembly; an all-categorical no-op and a
   mixed-design test that the move fires on the ordinal-rule nodes and is a fixed point on the categorical ones; the monotone and
   variance-forest guards; and **a Gibbs-dominant walk on a two-column design small enough for the conditional to be enumerated in
   the test, where the realized draw frequencies must match it** - itself a correctness gate at the kernel level, and the one that
   would catch a mis-assembled weight before any script runs.
2. **`rule-gibbs-balance.R`.** The prior-only arm on the 6-by-2 factorial at `power = 0.5`, the exact-posterior confirmation arm, both
   poisons. Roughly 450 to 550 lines; not startable before slice 1, and it is slice 1's rank-aware scan that makes the prior-only arm
   exist at all.
3. **The Stage 2 harness and run.** Two arms at `d` in `{0.16, 0.32}`, C1's four-chain cell, the sham arm, P1's rung, P2's
   duplicate-column cell and the equal-cost arm. Roughly 400 lines; compute on the order of a day plus the equal-cost arm's quiet
   machine. The verdict is recorded here.
4. **The default share, if slice 3 passes** - POST-RELEASE, and its gate does not exist: 6.4's second kill clause needs plateau
   prediction error in the noise-heavy or large-n stratum, which the grow-from-root harm battery measured and benchmarks/ does not
   contain ([5. Verdict and consequences](grow-from-root-default.md#5-verdict-and-consequences)). The equal-cost arm of section 6 is
   the other thing slice 4 needs and slice 3 supplies. A nonzero default is also a stream shift and pays for its own re-record.
