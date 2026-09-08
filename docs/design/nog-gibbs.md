# rule_gibbs: an exact draw of the split rule at a nog node

Status: PROPOSED, 2026-09-07; AMENDED 2026-09-07 (the veto's real law and the neighbourhood as a rank stratum, the cost table at 1 - stump%, the cost instrument, the balance gate sized, the surface at twenty-four files); SLICE 1 LANDED 2026-09-07 (the kernel at weight zero, 7fb166ca); SLICE 2 LANDED 2026-09-07 (rule-gibbs-balance.R, d888c9f3); SLICE 3 RUN 2026-09-07: NOT KILLED at d = 0.16, the coverage secondary fails (50032833); CUT-ONLY PILOT 2026-09-08: the private cut-only variant keeps about half the Trig+poly gain and all of the Single index one at 1.04 sweep-equivalents against 2.21 (a6f44e12); coverage flag dissolved by the reference arm 2026-09-08 (e002c10d); DOSE RESPONSE 2026-09-08 (17505c50).

A fifth tree kernel that replaces the Metropolis change proposal at a nog node - an interior node whose two children are both
leaves - with a draw from the rule's own full conditional. The neighbourhood is closed, the acceptance is identically one, and there
is no reverse count - on the branch-rank stratum the empty-leaf veto's own lexicographic law makes current, which section 2.2 states
and which is the whole of the correctness argument under a weight mask or a routed missing row.
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
repartitions those members, it does not change which rows are there. Three consequences: the candidate set is identical from every
state in it; its normalizer is the same before and after the draw; and no proposal count survives into the acceptance. Section 2.2
restricts that set to one branch-rank stratum, which is ancestor-determined in the same way and carries the same three consequences;
everything below reads on the stratum.

The weight of a candidate is the tree posterior restricted to `u`'s rule, every factor that does not read it having cancelled:

    log w(v, c) = S(v, c)                                          the children's collapsed marginals, 2.2's rank-0 sum
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
probe omits. The census identity of section 5 is therefore stated at `dart = FALSE`, which is what the census ran. `S` is the
scan's entry wherever both children carry positive weight, which is every candidate with no mask installed and no missing member at
the node; section 2.2 gives it in general.

Draw `(v, c)` from the normalized weights. Acceptance is one; no pairwise [`resolveVetoRank`](../../src/bartcore/moves.hpp)
comparison - 2.2 does that law's work in the enumeration instead - no snapshot-and-restore, no `logProposalCorrection`. Node
selection is uniform over the eligible nog nodes and its reciprocal cancels because that set is invariant under the move: shape is
preserved, so `fillNoGrand`'s set does not move ([`Tree::fillNoGrand`](../../src/bartcore/tree.hpp)), and NO NOG NODE IS AN ANCESTOR
OF ANOTHER - a node with an interior descendant is not nog - so changing `u`'s rule cannot move any other nog node's availability
set either.

Three of `changeMove`'s guards drop outright. No mask pool, an ordinal rule allocating no words. No stranding walk. And no interaction
walk: `tree.interactionSubtreeIsValid` exists because a redrawn variable can strand a descendant SPLIT, and a nog node has none, so
`collectAvailableVariables`'s own interaction test at `u` is the whole constraint.

### 2.2 The veto's law, and what the scan must emit

**The scan's occupancy test is not the veto's, and the gap is what defines the neighbourhood.**
[`scanOrdinalCuts`](../../src/bartcore/scan.hpp) writes [`cutScanEmptySentinel`](../../src/bartcore/scan.hpp), `-inf`, when either
side's weight over the NON-MISSING bins is non-positive, and on that branch it writes the sentinel to both missing directions and
computes neither side's marginal. The veto reads something else. [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp) is 2 when a leaf
holds no member, 1 when it holds members but no positive weight and 0 otherwise, all three off the leaf's ACTUAL index span with any
routed missing rows in it; [`logLikelihoodForBranch`](../../src/bartcore/moves.hpp) takes a branch's rank as the maximum over its
leaves and its log-likelihood as the marginal summed over the RANK-0 leaves alone; and
[`resolveVetoRank`](../../src/bartcore/moves.hpp) applies that pair LEXICOGRAPHICALLY - a rank-improving proposal takes `-HUGE_VAL`
on the current side and is accepted outright, a rank-worsening one takes it on the proposal side and is refused, and at equal ranks
the finite parts are compared as they always were.

Two states where the two laws disagree, and both are reachable.

- **A routed missing row.** At a node holding missing members, a cut whose non-missing members all fall one way with the missing rows
  routed the other leaves both children positive weight and branch rank 0, while the scan sentinels both of its directions.
  [`changeMove`](../../src/bartcore/moves.hpp) installs exactly that rule - at a nog node
  [`findGoodOrdinalRules`](../../src/bartcore/moves.hpp) collapses to [`Tree::splitInterval`](../../src/bartcore/tree.hpp), ancestors
  only and no member test, and the missing direction is a fair coin - so a scan-defined neighbourhood need not contain the incumbent,
  and a draw that leaves a positive-target state with probability one and cannot return is reversible with respect to nothing.
- **An installed weight mask.** Weights do not ride the tree, so a mask install strands whole branches at rank 1 - members but no
  positive weight - and from there the shipped moves compare finite parts and mix under prior x transition. The sentinel tests
  WEIGHT, so under an all-zero mask every candidate is sentinelled and the neighbourhood is empty, which is exactly the configuration
  the house's cheap balance gates run in.

**So the scan has to emit per SIDE, not per branch.** A single branch-rank flag cannot reconstruct the score: a MIXED candidate, one
child rank 0 and the other rank 1, has branch rank 1 and a log-likelihood equal to the rank-0 child's marginal alone, and the
sentinel path never computes it. **A, weight occupancy only**: read the scan as it stands and no-op from a rank-1 branch. Free, and
still wrong under routing unless the kernel also detects an incumbent it did not enumerate and no-ops there too; it costs section 5's
prior-only arm outright. **B, a rank-aware scan**: an optional out-parameter, defaulted null so
[`growTreeFromRoot`](../../src/bartcore/grow.hpp)'s call is untouched and bitwise unchanged, carrying per candidate the two SIDES'
ranks under the ROUTED occupancy and the marginal summed over the rank-0 sides. About forty lines rather than fifteen: the routed
count and weight are carried per side, and the sentinel branch computes the surviving side's marginal instead of skipping both.
**RECOMMEND B.** It is what puts the kernel on the veto's own law rather than on a law that merely agrees with it in the common case,
and it is the only thing that lets section 5's prior-only arm run at all.

**The neighbourhood is a rank STRATUM, and that is what makes the draw exact.** For a candidate `(v, c, s)` - `s` the missing
direction where the node routes missing rows - write `r_L` and `r_R` for its two sides' ranks and `S(v, c, s)` for the marginal
summed over its rank-0 sides. Then `r = max(r_L, r_R)` is `logLikelihoodForBranch`'s rank for that candidate and `S` is its
log-likelihood, candidate for candidate. Candidates of rank 2 are dropped absolutely: no move may install a member-empty leaf even
from a vetoed state, the membership law [`Tree::bottomNodesAreOccupied`](../../src/bartcore/tree.hpp) every site outside the move
kernels enforces. Let `r*` be the smallest rank the survivors carry. **The kernel draws over the `r*` stratum ALONE**, weighting each
of its members by section 2.1's `log w` with `S` as the marginal term. The law has two cases and they are not the same law:

- **`r*` equals the incumbent's rank.** The ordinary case, and the only one that occurs with no mask installed and no missing member
  at the node. The stratum contains the incumbent and is ancestor-determined, so it is identical from every state in it, and the draw
  IS the exact full conditional of the shipped chain's own target restricted to it: at `r* = 0` the posterior truncated to
  occupancy-admissible trees
  ([Which move paths can create an empty leaf](empty-leaf-veto.md#which-move-paths-can-create-an-empty-leaf)), at `r* = 1` the prior
  times `exp S`, which is what the shipped moves compare when ranks are equal. Acceptance is one.
- **`r*` is strictly better than the incumbent's.** The incumbent is outside the stratum and this is NOT a Gibbs step. It is a
  Metropolis-within-Gibbs step with acceptance one, which is what `resolveVetoRank` already gives every rank-improving proposal, and
  it is valid for the reason the shipped moves are: the better stratum is absorbing, no rank-worsening proposal ever being accepted,
  so the chain enters it in one step and is stationary there. No stationarity is claimed for the stratum it leaves.

`r*` is never WORSE than the incumbent's rank, so the stratum is never empty: `collectAvailableVariables` and `splitInterval` both
ignore `u`'s own rule, so the incumbent is always a candidate, and its own rank is at most 1 because every site outside the move
kernels enforces `bottomNodesAreOccupied`. With no weight vector installed the incumbent is rank 0 and the stratum is the whole
occupied candidate set - the rank-0 law is [`Tree::bottomNodesHaveWeight`](../../src/bartcore/tree.hpp), and tree.hpp's own note is
that it and the membership law agree exactly there.

**The missing direction is part of the candidate exactly when the node routes a missing row.** `scanOrdinalCuts` returns `2 * numCuts`
entries then, scoring the direction rather than leaving it to a coin, and the rule prior widens by the same factor two the candidate
count does - which is what the probe's `doubled` branch does. Where the column declares missing values but the NODE holds none, both
directions score identically, so enumerating the cut once at prior `1/|SI_v|` is the exact collapsed weight and the direction is
drawn from its own conditional, a fair coin, after the cut. That convention is the kernel's and the probe's;
[`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp) takes its `- log 2` from the COLUMN instead, so the two
differ by `log 2` per variable at such a node - which is what section 5's identity has to be stated around rather than a defect in
either.

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
2. The kernel picks its node uniformly among the ELIGIBLE NOG NODES, so it scans on every proposal a non-stump tree gets - every
binary tree carrying a split has a nog node - and the multiplier is `1 - stump%`, NOT the census's `target-nog%`, which is the share
of change proposals landing on a nog node and so a statistic of `changeMove`'s uniform-over-interior-nodes selection rather than of
this kernel's. Per sweep, at `m = 75` and one proposal per tree per sweep,

    units/sweep = m x d x (1 - stump%) x 2 x p_avail

against the census's own leaves-per-tree and stump shares:

    cell      p    L     stump%   units/sweep at d = 0.4   a sweep's own traffic
    default   10   2.83    2.4            586                      637
    lownoise  10   3.79    0.5            597                      853
    wide      50   2.53    7.5           2775                      569
    c1        30   2.52    6.1           1690                      567

The last column is `3 m L`, three full passes over `n` per tree for the residual, the leaf statistics and the fit rebuild - an
estimate, not a measurement. A cut-scan unit is a pass over `n/L` rows, so `L` cancels out of the move's own column and the units are
not commensurable ACROSS cells; only the ratio to a cell's own traffic is. **At C1 the move costs about 1690 cut-scan units a sweep at
the full change share and 676 at `d = 0.16`, so a sweep at the dosage section 6 runs costs 2.19 times what it costs today.** `bcf` is
absent because its treatment forest refuses a non-default `proposal.probs` outright.

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
gate for a saving the cut-only variant already has without either. **The cut-only variant was built and measured after the full
draw was**, on the private `-D` build this section named: one scan a proposal, 1.04 sweep-equivalents against the full
draw's 2.21, and about half of Trig+poly's gain with all of Single index's - the cut-only pilot of
[6. Benefit, pre-registered](#6-benefit-pre-registered).

## 3. Leaf models and families

The scan is templated on [`ScalarLeafModel`](../../src/bartcore/model.hpp) and the probe gates on
[`ScannableLeafModel`](../../src/bartcore/moves.hpp), a scalar leaf carrying the four-argument `(k, sigma^2, sum w, sum wz)` marginal.
That admits [`ConstantGaussianLeaf`](../../src/bartcore/model.hpp) and, because it forwards the same four-argument marginal,
[`MonotoneConstantGaussianLeaf`](../../src/bartcore/model.hpp), which the guard below excludes on a different conjunct and not on the
scan predicate. So the move is available for every response family whose mean forest carries a constant leaf: gaussian directly,
the latent families (probit, logistic, ordinal, nbinom, multinomial, aft, hazard) through the working response and weights the scan
reads exactly as `computeLeafStats` does, and heteroscedastic through the composed weights.

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
carry. src/bartcore/scan.hpp takes section 2.2's per-side rank and
rank-admitted marginal out-parameter. Then TWELVE in
[`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp) - three struct fields, four copies into a
forest, the variance forest's copy, the two `MoveContext` initializers and the two `structureIsFrozen` call sites - and THREE in
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp), without which the move is unreachable
from BCF and multinomial fits and silently zero. EIGHT in
[`ParsedModel`, `parseModel`, `printInitialSummary`, `optionsFromParsed`, `refuseUnsupportedAmplitudeComposition`, `buildMultinomialSampler`, `bartcore_setModel`](../../src/R_interface_bartcore.cpp),
spread over seven declarations and not the two a shorter list implies: `ParsedModel`'s field; `parseModel`'s slot read and its sum
check's fifth term; `printInitialSummary`'s creation printout; `optionsFromParsed`'s `SamplerOptions` copy;
`refuseUnsupportedAmplitudeComposition`'s hard-coded two-forest mixture; `buildMultinomialSampler`'s forest-spec copy; and
`bartcore_setModel`'s parameter copy. **Twenty-three probability SITES outside the kernel file** - twelve chain.hpp, three
combiner.hpp, eight bridge. That is a count of SITES; the count of FILES below is a different quantity, and no longer the same
number.

**tests/cpp, THREE files.** SEVENTEEN positional `MoveContext` initializers, fourteen in tests/cpp/test_moves.cpp and three in
tests/cpp/test_interaction.cpp: the probability block sits ahead of `const double* weights`, so a short initializer binds that pointer
to a `double` - a hard compile error, which is what makes the count safe. And FOUR positional `structureIsFrozen` calls in
[`testFrozenForest`](../../tests/cpp/test_sampler.cpp), which the same mechanism breaks - `structureIsFrozen` takes four required
doubles and no defaulted argument - so they are four more edits, not four to dodge: giving the fifth argument a default would give up
the compile-error safety this paragraph rests on everywhere else. **Four Rd files** (man/dbarts.Rd, man/bart2.Rd,
man/dbartsSpec.Rd usage and argument text; man/bart.Rd argument text alone). **SIX tinytest files**: test-proposal-probs.R,
test-argument-surface.R's pinned default, test-bcf-creation.R's refusal literal, test-spec.R and test-monotone.R's slot reads, and
test-sum-to-one-tolerance.R, whose three-name regression must keep failing where it fails today. **Plus inst/NEWS.Rd.**
**Twenty-four files.** Roughly 230 lines of kernel and scan, 60 across the surface sites, 150 of tests.

**The flat C header does not move** and no stored state does:
[`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) takes the model as a `SEXP`, so no proposal
probability crosses the ABI and no `LinkingTo` consumer recompiles.

## 5. Correctness: `rule-gibbs-balance.R`

Two arms, and unlike perturb's the second is mandatory: the prior-only arm exercises the node selection, the enumeration and the prior
factors but NOT one scan entry, and the scan-weighted draw is the whole new mechanism.

**The prior-only arm and its target.** Under an all-zero weight mask every leaf holding rows is rank 1 and every leaf holding none is
rank 2, so `r*` is 1, section 2.2's stratum is exactly the member-occupied candidate set, `S` is 0 across all of it and every
candidate's weight collapses to its prior factors alone. The kernel is then reversible with respect to the CGM prior TRUNCATED to
member-occupied trees and renormalized, the same target
[4. Correctness: perturb-balance.R](perturb-move.md#4-correctness-perturb-balancer) names. **The design makes the truncation
vacuous, and is chosen to make both poisons resolvable**: two ordinal columns as a FULL FACTORIAL, x1 on 6 distinct values
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
birth/death because statistic 2 moves through nothing else.

**The run, sized against [4. Correctness: perturb-balance.R](perturb-move.md#4-correctness-perturb-balancer)'s template, clause for
clause.** States of prior mass below 0.004 are dropped by the same pre-stated rule, the leaf count's tail pooled into a `>= j` bin
the dynamic program fixes ahead of the run: statistic 1 keeps all seven states and statistic 3 all five, statistic 2 at most twelve,
so the family is at most **`m = 24`**. Batch-means z per retained state, Holm at family alpha 0.05, thresholds running from
`|z| = 1.96` (least strict) to `|z| = 3.08` (strictest, at `m = 24`). **Run length 4 chains x 250,000 kept draws at `n.thin = 20`**,
batch means over 500 batches of 500 consecutive kept draws per chain, the four independent chains pooled as perturb's gate pools them
([`batchMeanSE`](../../benchmarks/R/change-balance.R) is the same estimator at one chain). **Burn-in 20,000 sweeps per chain**,
discarded before the batching and sized off the leaf count, the slowest statistic here as there: birth/death holds 0.10 of the one
proposal a tree gets per sweep, so that is roughly 2,000 dimension proposals against a support of 1..12. The script does not assume
it - it reads the first lag at which each statistic's kept-draw autocorrelation falls under 0.1
([`firstUnder`](../../benchmarks/R/perturb-balance.R), the BURN-IN adequacy ladder and not a run-length rule) and refuses to score a
run whose burn-in is under fifty of those lags, doubling and re-running instead. **The undiluted detection floors**, at the family's
strictest threshold: about 210 effective draws for poison (i) on statistic 3's end states, about 205 for poison (ii) on statistic 1's
x1 cuts and 24 on its x2 cut. Statistic 3's are CONDITIONAL draws, so the same refuse-and-double rule covers them: the script takes
the conditioning event's mass from the same dynamic program and refuses to score statistic 3 unless the realized conditional
effective count clears ten times its floor. The masked sweep costs the move machinery alone on 12 rows, so length is nearly free;
statistic 3's conditioning, not the run, is what has to be paid for.

**The confirmation arm** runs the same grid with positive weights and no mask, scored against `change-balance.R`'s region dynamic
program ([The gate](change-move-balance.md#the-gate)), and it is what tests the scan entries themselves.

**The census identity, as a tests/cpp assertion rather than a script.** Assemble the neighbourhood twice at a fixed tree: once by the
kernel, once by a reference that installs each candidate rule, calls [`Tree::refreshSubtree`](../../src/bartcore/tree.hpp) - without
it [`ConstantGaussianLeaf`](../../src/bartcore/model.hpp)'s per-node marginal reads the STALE cached `sumWeights` and
`sumWeightedResponse`, which is why `changeMove` refreshes before it scores - and only then calls
[`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp) and `logLikelihoodForBranch`. Two things cannot be asserted per
candidate as they stand. The `- log 2` for a missing direction is the column's in `ruleForVariableLogProbability` and the node's in
the kernel, so on a column declaring missing values at a node holding none the two differ by `log 2` per variable; and a rank-1
candidate's reference score has to be `logLikelihoodForBranch`'s partial sum, not a full branch marginal. **The assertion's fixture is
therefore a design with no missing values and no mask installed**, where neither qualification bites and the identity holds candidate
by candidate at 1e-12; the missing-data convention is asserted separately, PER CUT summed over its two directions, which is the form
in which it is exact. That is stronger than agreeing with `census::nogProbe`, which the census build's own 3.7e-13 already records,
and it does not need the census build to run.

## 6. Benefit, pre-registered

**Prerequisites, both outside this design and both met.** benchmarks/R/surfaces/P1-friedman.R's 90 percent coverage must return near
0.71 in the control arm before any verdict is valid; it reads 0.725
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)).
And the chain configuration is settled: the shipped four chains at 500 + 500 on C1's independent design at 75 trees, which is what
turns the primary from coverage into minimum ESS
([5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic)).

**Two arms, matched seeds, twenty pairs.** **A**, the shipped mixture, which reads 15(8-31) summed minimum ESS on Trig+poly. **B**,
`rule_gibbs d` taken from change. **Dosage: `d` in `{0.16, 0.32}`.** `d = 0.40` is excluded, and the reason is sharper here than it
was for perturb - stated in its two DIFFERENT denominators, which the census keeps apart. `target-nog%` is 77.3 at `c1`, the share of
change PROPOSALS landing on a nog node, so a full share strands the other 22.7 percent of that traffic; `nog%` is 65.2, the share of
`c1`'s INTERIOR NODES that are nog, so the 34.8 percent that are not would be left with no mechanism that moves their rule at all.
It is the second that decides the exclusion. **The primary is the summed minimum ESS over C1's 25 fixed points at the +8 bar**, the bar
[5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic) derives from the cell's own spread and which
the eight CELLS already recorded there - four move-set arms crossed with two mean functions, not eight arms - calibrate at a paired
standard error of 1.8 to 2.5. TWO secondaries at
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s own margins, both must-not-degrade: coverage at
-0.010 absolute against arm A's 0.961, and held-out RMSE at a ratio above 1.02. Interval length is REPORTED and carries no margin,
6.4 setting none for it; cost is not a secondary here at all but a conjunct of the kill, below.

**Controls, as [5.3 What arm B must produce, and the kill](perturb-move.md#53-what-arm-b-must-produce-and-the-kill) has them.** The
bitwise null at `d = 0`, which slice 1 runs and which is 6.4's own first absolute gate. A sham arm, A against A at twenty FRESH
sampler seeds, whose paired difference must sit inside the +8 bar. And P1's rung re-run as the absolute control at the time slice 3
runs. A flagged cell takes a mandatory fresh-seed re-run before the flag counts
([6.1 The rule, stated operationally](benchmark-surfaces.md#61-the-rule-stated-operationally)).

**Cost, its instrument, and why it is not wall time.** 6.4 carries TWO cost-bearing metrics, "wall time per sweep" at a ratio above
1.05 and "minimum ESS over 25 fixed points, per second" at a ratio below 0.90, and neither can be read here. Section 2.4 prices the
move at 2.19 sweeps at `d = 0.16` on `c1`, so both would fire by construction; and 10.4's host carried a load of 9 to 67 throughout
its runs, so a per-second statistic measured there would be noise in any case
([10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)). **Cost is measured as the CUT-SCAN
COUNT instead, and the census build does not produce it today.** The `g` record's legend carries `scanned` and `jointCandidates`,
which counts finite-weight (variable, cut) PAIRS and is not `p_avail`; there is no available-ordinal-variable column; and the hook
that writes `g` sits inside `changeMove`, which at a nonzero `rule_gibbs` share is the branch the dispatch does not take. The
instrument is therefore a NEW census hook in the new kernel and a new column beside
[`nogNames`](../../benchmarks/R/move-census.R) - proposals reaching an eligible nog node, ordinal variables scanned, and their
product per sweep - and it is a slice 1 deliverable, not a free read. The census build draws nothing and restores every byte it
touches, so a census run at slice 3's own seeds and mixture reproduces the arm's chain exactly and its scan count IS the arm's.
Beside it, as the cost-adjusted secondary, ONE equal-cost arm on a quiet machine in which arm A is given the sweep count that ratio
buys it; that arm is where 6.4's per-second question is actually answered. The draws-based primary answers whether the kernel mixes
better; the equal-cost arm answers whether it should ever be a default, which is slice 4's question anyway.

**The rooting-lock probe, as a control and not as a place the move can act.** P2's duplicate-column cell at `m = 1`,
must-not-degrade. The shipped mixture parks 5 of 40 chains on an x3 root, change being vetoed once a child splits on x1
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)). At `m = 1` the ROOT IS NOG ONLY
WHEN THE TREE HAS ONE SPLIT, and 10.1 records every one of those 5 chains at an x3 root with 3 to 4 interior nodes and a child
already split on x1 - four to five leaves, root not nog. **The move reaches none of the recorded parked chains**, so this probe
cannot pass by improvement and is not evidence for the move. The only cell in which a rooting effect could show is a chain that
passes through a one-split tree, where the draw resolves the rooting exactly and in one step rather than proposing at it; none of the
5 does. The move does not rotate a rule up the tree either; swap remains the only move that does. What the probe is for is the
mechanism the move DISPLACES: arm B's mean switches per chain inside arm A's seed range, chains parked no more than 10 of 40, and
pooled `p(root on x1)` within Monte Carlo error of 0.5, all as must-not-degrade controls on the change share `d` takes away.

**Kill criterion. KILL if, at `d = 0.16` on the shipped four-chain configuration, arm B does not improve C1's summed minimum ESS over
arm A by more than +8 over at least 20 matched pairs; OR if arm B's cut-scan cost, read off the instrument above, exceeds 2.4
sweep-equivalents against section 2.4's own predicted 2.19. A mandatory fresh-seed re-run of any flagged cell before a flag counts.**
FIVE departures from [6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered), each stated.
(a) The statistic is minimum ESS rather than coverage, which has no headroom at 0.961. (b) The cell is C1 rather than the low-noise
cell, three shipped mixtures being indistinguishable on P1. (c) 6.4's second conjunct, arm C against arm D, is dropped with arm C, so
the kill fires on B's failure alone rather than on B's and C's together - STRICTER than 6.4, and the same departure
[5.3 What arm B must produce, and the kill](perturb-move.md#53-what-arm-b-must-produce-and-the-kill) names. (d) BOTH of 6.4's
cost-bearing metrics go: wall time per sweep and minimum ESS per second would each fire by construction on a move priced at 2.19
sweeps, and neither is readable on 10.4's host, so the cut-scan ratio above stands in their place as the kill's cost conjunct and the
equal-cost arm carries the per-second question. That is a LOOSER cost gate than 6.4's, which is why the equal-cost arm is mandatory
rather than optional. (e) 6.4's plateau-error clause needs a stratum no cell here measures and belongs to slice 4.

**What the design does not predict.** Perturb's pilot moved this exact statistic by +0.1 +/- 8.1, and the frozen-structure run says
the worst C1 point's minimum ESS rises only from 1.6 to 4-21 with structure held fixed entirely - so most of that deficit survives any
structural kernel. Nothing connects a wider rule neighbourhood to a move in the minimum ESS over 25 points, and the entropy the probe
measured is a one-step statistic on a chain this move never ran. Slice 3 can fail honestly.

**Verdict (2026-09-07).** The Stage 2 harness ran both dosages on the shipped four-chain cell over twenty matched pairs, both mean
functions, with the control re-run in the same session; it reproduces
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s own `independent75pool4` rows digit
for digit on both. Absolute readings, mean over seeds (min-max):

    mean fn      arm                             seeds  95% coverage        length  RMSE  min ESS (sum)  per chain  between
    trigpoly     independent75pool4              1-20   0.961(0.945-0.977)  4.61    1.12  15(8-31)       2(1-2)     0.78
    trigpoly     independent75pool4ruleGibbsB    1-20   0.939(0.911-0.958)  3.98    1.09  36(19-53)      3(2-4)     0.58
    trigpoly     independent75pool4ruleGibbs32   1-20   0.939(0.917-0.957)  3.96    1.09  52(23-97)      3(2-6)     0.55
    trigpoly     independent75pool4              21-40  0.964(0.950-0.977)  4.58    1.10  16(9-28)       2(1-2)     0.80
    trigpoly     independent75pool4ruleGibbsB    21-40  0.938(0.917-0.951)  3.96    1.08  39(15-55)      2(2-3)     0.59
    singleindex  independent75pool4              1-20   0.895(0.878-0.915)  6.45    1.92  21(9-32)       2(1-2)     0.68
    singleindex  independent75pool4ruleGibbsB    1-20   0.900(0.872-0.921)  6.57    1.97  35(16-57)      2(2-3)     0.40
    singleindex  independent75pool4ruleGibbs32   1-20   0.904(0.872-0.923)  6.61    1.96  36(18-74)      3(2-5)     0.37
    singleindex  independent75pool4              21-40  0.891(0.868-0.911)  6.43    1.93  15(9-27)       2(1-2)     0.67
    singleindex  independent75pool4ruleGibbsB    21-40  0.898(0.879-0.914)  6.55    1.97  39(20-60)      3(2-6)     0.39

Paired differences against each row's own control on its own seeds, mean +/- sd (seeds positive of 20):

    mean fn      arm                             seeds  d min ESS (sum)               d per-chain min ESS    d 95% coverage                     d RMSE, ratio
    trigpoly     independent75pool4ruleGibbsB    1-20   +21.5 +/- 12.8 (19/20) t 7.51  +1.14 +/- 0.59 t 8.66  -0.022 +/- 0.008 (0/20) t -12.46   -0.028 +/- 0.071, ratio 0.975
    trigpoly     independent75pool4ruleGibbs32   1-20   +37.0 +/- 19.6 (20/20) t 8.45  +1.37 +/- 0.96 t 6.37  -0.023 +/- 0.008 (0/20) t -13.20   -0.026 +/- 0.070, ratio 0.977
    trigpoly     independent75pool4ruleGibbsB    21-40  +22.1 +/- 12.2 (19/20) t 8.10  +0.86 +/- 0.53 t 7.21  -0.026 +/- 0.012 (0/20) t -9.93    -0.016 +/- 0.045, ratio 0.986
    singleindex  independent75pool4ruleGibbsB    1-20   +14.4 +/- 13.2 (17/20) t 4.87  +0.73 +/- 0.47 t 6.94  +0.005 +/- 0.010 (13/20) t 2.38    +0.046 +/- 0.035, ratio 1.024
    singleindex  independent75pool4ruleGibbs32   1-20   +15.1 +/- 12.8 (18/20) t 5.29  +1.09 +/- 0.81 t 5.99  +0.009 +/- 0.010 (17/20) t 3.93    +0.042 +/- 0.036, ratio 1.022
    singleindex  independent75pool4ruleGibbsB    21-40  +23.9 +/- 11.0 (20/20) t 9.74  +1.14 +/- 0.86 t 5.91  +0.006 +/- 0.008 (18/20) t 3.51    +0.043 +/- 0.027, ratio 1.022

Held-out RMSE ratio, in the order of the paired table: 0.979, 0.980, 0.987, 1.028, 1.027, 1.023. Wall ratio: 3.12, 6.20, 3.05,
3.28, 6.86, 4.15. Paired standard error of the summed minimum ESS runs 2.5 to 4.4 across the six cells.

**The kill does not fire.** Both of its clauses are met at the confirmatory dosage. Arm B improves the summed minimum ESS by +21.5
at the first block and +22.1 at the fresh one, each over twenty matched pairs and each about eight paired standard errors, against
a +8 bar; and the cut-scan cost is 2.21 sweep-equivalents against the 2.4 conjunct, below. The per-chain minimum moves with it -
1.51 to 2.65 and 1.55 to 2.42 - which no four-chain arm recorded on this mean function did before; the summed gain is four chains
mixing better, not one chain carrying three.

**What fails is a secondary, and it fails twice.** 95 percent coverage on the primary cell moves -0.022 (t -12.46) and -0.026
(t -9.93), past
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
-0.010 margin with the one-sided bound excluding the null in both blocks, so the mandatory fresh-seed re-run
([6.1 The rule, stated operationally](benchmark-surfaces.md#61-the-rule-stated-operationally)) confirms the flag rather than
dissolving it. The other secondary passes: held-out RMSE reads 0.979 and 0.987, better than the control. The mechanism is in the
same table. Arm B's pooled interval is 14 percent shorter at slightly better point accuracy, and its between-chain ratio falls from
0.78 to 0.58; 10.4 reads this cell as four chains each sitting in its own place with pooling supplying the width, and the move makes
them agree, so part of that width goes. What the width was carrying was OVER-coverage: the control sits at 0.961 against a nominal
0.95 and arm B at 0.939, the same 0.011 miss with the sign reversed. Against nominal the two arms are equally far out. Against arm
A, which is what the margin is stated against and what this design pre-registered, arm B regresses, and the pre-registration stands
as written.

**Single index, reported and not gated**, gains as much and pays for it in RMSE: +14.4 and +23.9 summed minimum ESS
(t 4.87 and 9.74), coverage +0.005 and +0.006, held-out RMSE ratio 1.028 and 1.023, past the 1.02 margin at both blocks. The two
mean functions move coverage in OPPOSITE directions, which is what a narrowing interval does when the control over-covers on one
and under-covers on the other.

**`d = 0.32` buys more ESS and nothing else.** +37.0 against +21.5 on Trig+poly and +15.1 against +14.4 on Single index, at the
same coverage move and twice the wall. It is reported; it is not the kill's cell and nothing here recommends it.

**The controls.** P1's rung re-ran as the absolute gate, 60 fits: default 0.725 (0.682-0.760) held-out, birthdeath 0.709, swap
0.728, identical to
[10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07),
so the gate is in force. P2's duplicate-column cell took a fourth arm at arm B's mixture: 146.8 (109.6-159.0) root switches per
chain against arm A's 70.8 (57.9-81.2), outside arm A's seed range and ABOVE it, which is not a degradation; 6 of 40 chains parked
against arm A's 5, inside the 10 allowed; pooled p(root on x1) 0.437 (0.373-0.524) against arm A's 0.457, and 0.514 against 0.522
among on-pair draws, so the null is intact. On the confounded design, which the move cannot act on and which gates nothing, the
share on x1 given the pair reads 0.354 (0.267-0.500) against arm A's 0.523; at five seeds and half a switch per chain that is forty
locked chains, and [10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)'s own
birth/death arm reads 0.504 over a 0.167-0.833 range.

**Cost, on the instrument this section named.** A census run of the `c1` cell reports 337.32 cut scans per sweep, 2.47 leaves per
tree and a 5.25 percent stump share, the kernel reaching an eligible nog node on 95.05 percent of its proposals and scanning all 30
variables there. One variable scanned at a nog node is TWO of
[2.4 Cost, and the two restricted variants](#24-cost-and-the-two-restricted-variants)'s units, the node holding two leaves' worth of
members, so 674.6 units a sweep against that section's predicted 676; against `3 m L` = 556 at the measured leaf count, or 567 at
2.4's own 2.52, the move costs **2.21 sweep-equivalents** (2.19 at 2.4's `L`), under the 2.4 conjunct and within one percent of
what the design predicted. Two departures from what this section asked for, both stated. The census `c1` cell runs at its own data
and sampler seed and one chain of 200 + 500 sweeps, not the arm's twenty seeds at four chains, so this is not literally the arm's
own chain. And benchmarks/R/move-census.R's rule_gibbs mixture takes its 0.16 from change AND birth/death rather than from change
alone, so the rule_gibbs share - the quantity the scan count is proportional to - is arm B's exactly while the tree population is
not; no scaling is applied and that is the assumption. Wall does not follow the unit count: arm B measures 3.05 to 3.12 on
Trig+poly against a cut-scan 2.21, a scan pass doing more per row than the residual pass the unit is calibrated on. Wall gates
nothing here, and the host carried a load of 22 to 85 throughout.

**Residue, and what is owed.** The equal-cost arm - arm A given the sweep count the 2.21 ratio buys it - was NOT run here; it is
where 6.4's per-second question is answered and it is what slice 4 needs, and the dose-response paragraph below runs it. `rule-gibbs-balance.R`, slice 2,
landed at d888c9f3 while this run was in flight, so the correctness argument rests on that prior-only detailed-balance gate as well as on slice 1's own tests.
Slice 4 is NOT closed by a kill here: what blocks it is the coverage secondary, which a nonzero default would carry onto every core
cell, and its own missing gate. What slice 3 settles is narrower, and is worth stating plainly - this is the first kernel change
measured on this cell that moves its gated mean function's mixing statistic, by a wide margin and at a price the design predicted
to within one percent, and it cannot be defaulted until the interval it narrows is understood against nominal rather than against
arm A.

**Cut-only pilot (2026-09-08).**
[2.4 Cost, and the two restricted variants](#24-cost-and-the-two-restricted-variants)'s first restricted variant, built and measured
on a private `-DBARTCORE_RULE_GIBBS_CUT_ONLY` build the way the census scaffolding is built, and off in every shipped build. With the
macro defined [`enumerateNogRuleNeighbourhood`](../../src/bartcore/moves.hpp) holds the node's incumbent variable and enumerates that
variable's cuts alone - ONE scan a proposal - and every other term of the law stands: the branch-rank stratum, the rank-admitted
marginal, the missing coin, the categorical fixed point, and the prior factors, of which the split-variable prior and `1/|SI_v|` now
cancel outright over one variable's cuts and are KEPT rather than dropped, so a candidate's weight stays on the joint kernel's own
scale. It is still an exact Gibbs step: the restricted set is a deterministic function of state the move cannot change, so it is the
same set from every state in it and no proposal count survives. What it gives up is the variable axis, which then moves only at
change's own share.

**The gates hold on the restricted kernel.** The macro-OFF build is bitwise the shipped engine and the macro-ON build is bitwise the
same at the shipped defaults, the rule_gibbs share being zero there: 50 / 12 / 11 identical against the `fbff1989` baselines on both.
[`rule-gibbs-balance.R`](../../benchmarks/R/rule-gibbs-balance.R) PASSES on the cut-only kernel in both modes at its own 100,000-sweep
burn-in - the autocorrelation ladder lengthens from 52 kept-draw lags to 69 against a need of 100, so the refuse-and-double rule does
not fire - the prior-only arm's worst `|z|` 2.93 at `root x1c3` (0 of 23 Holm rejections, strictest threshold 3.07) against the joint
kernel's 1.74, and the confirmation arm's 1.09 at `x2c1` (0 of 7). Both poisons still FAIL, at `|z|` 69.02 (5 of 23) and 103.38
(6 of 23).

**Cost, on the same instrument.** A census run of the `c1` cell at the same mixture reports 10.97 cut scans per sweep against the full
draw's 337.32 - one scan at the 93.7 percent of proposals reaching an eligible nog node, where the full draw scanned all 30 variables -
which is 21.9 cut-scan units a sweep against `3 m L` = 549 at the run's own 2.44 leaves per tree: **1.04 sweep-equivalents against the
full draw's 2.21**. The variable axis is 97 percent of the move's own price.

**What it keeps.** The same harness, the same four-chain cell, the same twenty matched pairs and the same arm-B mixture, the arm name
selecting the mixture and the library the kernel. The control reproduces
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s own rows digit for digit at both seed
blocks and on both mean functions, which is the bitwise null read in situ.

    mean fn      arm                             seeds  95% coverage        length  RMSE  min ESS (sum)  per chain  between
    trigpoly     independent75pool4              1-20   0.961(0.945-0.977)  4.61    1.12  15(8-31)       2(1-2)     0.78
    trigpoly     ruleGibbsB, cut-only build      1-20   0.955(0.931-0.969)  4.33    1.11  26(8-53)       2(1-3)     0.61
    trigpoly     independent75pool4              21-40  0.964(0.950-0.977)  4.58    1.10  16(9-28)       2(1-2)     0.80
    trigpoly     ruleGibbsB, cut-only build      21-40  0.952(0.926-0.968)  4.36    1.11  26(9-61)       2(1-3)     0.64
    singleindex  independent75pool4              1-20   0.895(0.878-0.915)  6.45    1.92  21(9-32)       2(1-2)     0.68
    singleindex  ruleGibbsB, cut-only build      1-20   0.897(0.869-0.916)  6.51    1.97  37(14-68)      3(2-5)     0.41

Paired differences against each row's own control on its own seeds, mean +/- sd (seeds positive of 20), with the full draw's own
reading beside each:

    mean fn      seeds  d min ESS (sum)               full draw  d 95% coverage                   full draw  d RMSE, ratio
    trigpoly     1-20   +11.3 +/- 16.2 (16/20) t 3.13  +21.5     -0.007 +/- 0.009 (5/20) t -3.41  -0.022     -0.006 +/- 0.075, ratio 0.995
    trigpoly     21-40  +9.1 +/- 12.9 (15/20) t 3.14   +22.1     -0.012 +/- 0.012 (2/20) t -4.69  -0.026     +0.016 +/- 0.059, ratio 1.014
    singleindex  1-20   +16.2 +/- 11.3 (20/20) t 6.43  +14.4     +0.002 +/- 0.008 (12/20) t 1.36  +0.005     +0.049 +/- 0.024, ratio 1.025

Per-chain minimum ESS moves +0.36, +0.40 and +0.88 (t 4.78, 4.40, 5.82). Held-out RMSE ratio: 0.992, 1.012, 1.033. Wall carries no
claim - the host ran a load of 12 to 88 and the same two arms read a wall ratio of 1.84 at one block and 0.62 at the other.

**The variable axis earns about half its keep on Trig+poly and none on Single index.** Trig+poly holds +11.3 and +9.1 of the full
draw's +21.5 and +22.1, above the +8 bar at both blocks and at a third of the paired t; Single index holds +16.2 against +14.4, which
is the whole of it. The secondaries shrink with the gain on the gated cell - coverage moves -0.007 and -0.012 against the full draw's
-0.022 and -0.026, straddling
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
-0.010 margin rather than clearing it by two to three times, and the interval narrows 6 and 5 percent against the full draw's 14 - so
the mechanism is the same one at a smaller amplitude, not a different one. Single index pays what it paid before and slightly more:
held-out RMSE 1.033 against 1.028, past the 1.02 margin. **So the thirty-fold price of the variable axis buys, on the one cell whose
entropy columns said it carried content at all, about half of one mean function's ESS gain and none of the other's, and it buys most
of the coverage regression that flagged the full draw.** The cut-only kernel is not proposed as a default here and nothing above
changes slice 4, which the coverage secondary and its missing gate still block; what the pilot settles is the price of the axis,
which section 2.4 could only bound.

**Reference arm (2026-09-08).** The coverage secondary above is stated against arm A, and arm A's own 0.961 is now measured to be
an artefact of its chain configuration. A well-mixed reference on this cell - `rule_gibbs` at `d` = 0.32, four chains of 1000 +
2500, between-chain ratio 0.48 against arm A's 0.78 and the same coverage at half its kept length, 0.939 against 0.941 - reads
0.941 on Trig+poly against a nominal 0.95, within 0.002 of arm B's 0.939 and 0.020 below arm A's 0.961; and an eight-chain arm at
one kernel and one length moves coverage 0.939 to 0.956 where five times the length moves it 0.939 to 0.941, so what carries arm
A's extra coverage is chain COUNT and not the posterior
([10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)). Residual disagreement widens a
pooled interval, so 0.941 is an upper bound on a fully mixed reading and the direction is not in doubt.
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)
now reads coverage against that reference rather than against the shipped-configuration control, margin unchanged at -0.010, and
arm B's -0.002 against it is inside the margin: **the coverage flag is dissolved**. Single index goes the same way - the reference
reads a held-out RMSE ratio of 1.024 against arm A, so arm B's 1.028 is 1.004 against the reference - and that mean function is
reported, not gated, either way. What slice 4 still lacked at that point was 6.4's own missing plateau-error gate and the
equal-cost arm, the second of which the paragraph below runs.

**Dose response and the equal-cost arm (2026-09-08).** The cut-only variant of
[2.4 Cost, and the two restricted variants](#24-cost-and-the-two-restricted-variants) run at a third dosage below the pilot's and at
the pilot's own two, and the equal-cost arm the verdict above owed. Same private `-DBARTCORE_RULE_GIBBS_CUT_ONLY` library, same
four-chain cell, same twenty matched pairs. The harness gains two arms (17505c50):
[`independent75pool4ruleGibbs08`](../../benchmarks/R/surfaces/C1-he-hahn.R), the half dose taken from change alone, and
[`independent75pool4equalCost`](../../benchmarks/R/surfaces/C1-he-hahn.R), the shipped mixture at four chains of 1105 burn-in and
1105 kept. Every Trig+poly dose cleared the +8 bar at the first block, so every dose took the fresh-seed block on seeds 21 to 40,
control and all. The control reproduces
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s own rows digit for digit at both
blocks and on both mean functions and the `d` = 0.16 arm reproduces the cut-only pilot's, which is the bitwise null read in situ
twice over.

    mean fn      arm                             seeds  95% coverage        length  RMSE  min ESS (sum)  per chain  between
    trigpoly     independent75pool4              1-20   0.961(0.945-0.977)  4.61    1.12  15(8-31)       2(1-2)     0.78
    trigpoly     independent75pool4ruleGibbs08   1-20   0.955(0.927-0.975)  4.38    1.11  23(9-52)       2(1-3)     0.63
    trigpoly     independent75pool4ruleGibbsB    1-20   0.955(0.931-0.969)  4.33    1.11  26(8-53)       2(1-3)     0.61
    trigpoly     independent75pool4ruleGibbs32   1-20   0.956(0.934-0.974)  4.58    1.17  28(14-68)      2(1-3)     0.61
    trigpoly     independent75pool4              21-40  0.964(0.950-0.977)  4.58    1.10  16(9-28)       2(1-2)     0.80
    trigpoly     independent75pool4ruleGibbs08   21-40  0.954(0.937-0.971)  4.31    1.10  18(8-32)       2(1-2)     0.68
    trigpoly     independent75pool4ruleGibbsB    21-40  0.952(0.926-0.968)  4.36    1.11  26(9-61)       2(1-3)     0.64
    trigpoly     independent75pool4ruleGibbs32   21-40  0.955(0.919-0.981)  4.62    1.16  26(9-66)       2(1-4)     0.63
    singleindex  independent75pool4              1-20   0.895(0.878-0.915)  6.45    1.92  21(9-32)       2(1-2)     0.68
    singleindex  independent75pool4ruleGibbs08   1-20   0.894(0.864-0.910)  6.50    1.97  29(12-54)      2(2-3)     0.44
    singleindex  independent75pool4ruleGibbsB    1-20   0.897(0.869-0.916)  6.51    1.97  37(14-68)      3(2-5)     0.41
    singleindex  independent75pool4ruleGibbs32   1-20   0.898(0.859-0.920)  6.54    1.97  48(9-88)       3(2-5)     0.40

Paired differences against each row's own control on its own seeds, mean +/- sd (seeds positive of 20), with the full draw's own
summed-ESS reading and the dose's cut-scan cost beside each:

    mean fn      d     seeds  d min ESS (sum)                full   d per-chain            d 95% coverage                    d RMSE, ratio            held-out  cost
    trigpoly     0.08  1-20   +8.1 +/- 10.3 (15/20) t 3.54   -      +0.29 +/- 0.29 t 4.43  -0.006 +/- 0.009 (6/20) t -3.02   -0.005 +/- 0.069, 0.995  0.999     1.02
    trigpoly     0.16  1-20   +11.3 +/- 16.2 (16/20) t 3.13  +21.5  +0.36 +/- 0.34 t 4.78  -0.007 +/- 0.009 (5/20) t -3.41   -0.006 +/- 0.075, 0.995  0.992     1.04
    trigpoly     0.32  1-20   +13.4 +/- 13.7 (20/20) t 4.37  +37.0  +0.42 +/- 0.44 t 4.23  -0.005 +/- 0.011 (6/20) t -2.16   +0.057 +/- 0.118, 1.052  1.057     1.08
    trigpoly     0.08  21-40  +1.8 +/- 6.1 (11/20) t 1.29    -      +0.19 +/- 0.28 t 3.04  -0.010 +/- 0.010 (2/20) t -4.38   +0.002 +/- 0.049, 1.002  0.994     1.02
    trigpoly     0.16  21-40  +9.1 +/- 12.9 (15/20) t 3.14   +22.1  +0.40 +/- 0.41 t 4.40  -0.012 +/- 0.012 (2/20) t -4.69   +0.016 +/- 0.059, 1.014  1.012     1.04
    trigpoly     0.32  21-40  +9.3 +/- 16.8 (15/20) t 2.48   -      +0.34 +/- 0.60 t 2.58  -0.009 +/- 0.014 (5/20) t -2.83   +0.062 +/- 0.100, 1.057  1.048     1.08
    singleindex  0.08  1-20   +8.3 +/- 12.4 (16/20) t 2.98   -      +0.59 +/- 0.47 t 5.60  -0.000 +/- 0.009 (10/20) t -0.19  +0.051 +/- 0.033, 1.027  1.033     1.02
    singleindex  0.16  1-20   +16.2 +/- 11.3 (20/20) t 6.43  +14.4  +0.88 +/- 0.67 t 5.82  +0.002 +/- 0.008 (12/20) t 1.36   +0.049 +/- 0.024, 1.025  1.033     1.04
    singleindex  0.32  1-20   +27.5 +/- 21.8 (19/20) t 5.63  +15.1  +1.12 +/- 0.76 t 6.62  +0.003 +/- 0.010 (13/20) t 1.26   +0.051 +/- 0.034, 1.027  1.032     1.08

Paired standard error of the summed minimum ESS runs 1.35 to 4.88 across the nine cells. Cost is in sweep-equivalents off the census
instrument: a run of the `c1` cell on the cut-only census build reports 10.97 cut scans a sweep at its own 0.16 share, reach 93.7
percent and one variable scanned, against 2.44 leaves per tree - the pilot's numbers exactly - so 21.9 units against `3 m L` = 549
and 1.04. The 0.08 and 0.32 rows scale that measured scan count LINEARLY in the share, which holds the tree population fixed and is
an assumption rather than a measurement; the full draw's own 337.32 scans a sweep scale the same way, to 1.61 and 3.43.

**The dose response is not monotone on the gated mean function, and the restriction is why.** Trig+poly rises +8.1, +11.3, +13.4 at
the first block and +1.8, +9.1, +9.3 at the second, so `d` = 0.08 does not survive its own fresh block (t 1.29 against a +8 bar) and
`d` = 0.32 buys two ESS points over `d` = 0.16 while regressing a gated secondary at BOTH blocks: held-out RMSE 1.057 and 1.048
against
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
1.02. The mechanism is the restriction itself. The cut-only kernel never moves the split variable, so its share comes out of change,
the only move that does; at `d` = 0.32 change is left at 0.08, the pooled interval goes back UP to 4.58 and 4.62 from `d` = 0.16's
4.33 and 4.36, and in-sample RMSE with it, 1.17 and 1.16 against 1.11 - chains that agree better on a worse fit. The full draw at the
same dose does the opposite, +37.0 at a held-out ratio of 0.980, because it carries the variable axis itself. Single index shows none
of it: one rotated ridge gives the variable axis little to buy, and the cut-only kernel reads +8.3, +16.2, +27.5 monotonically at a
flat 1.03 held-out ratio throughout. Coverage flags nowhere. Read against
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s well-mixed reference of 0.941 rather
than against the control, every cut-only Trig+poly dose sits 0.011 to 0.015 ABOVE the reference and the -0.010 margin is nowhere
approached.

**Cut-only at `d` = 0.32 does NOT match the full draw at `d` = 0.16, so the equal-cost question is not closed by dominance.** +13.4
and +9.3 against +21.5 and +22.1 on the gated mean function, at a held-out RMSE the full draw did not pay, against a cost of 1.08
sweep-equivalents to its 2.21. On Single index it more than matches, +27.5 against +14.4, but that mean function is reported and not
gated. The equal-cost arm is therefore the read, and it was run.

**The equal-cost arm.** Arm A given the sweep count the 2.21 ratio buys it, four chains of 1105 + 1105 against the shipped 500 + 500,
on the shipped library, with the control and arm B in the same session; both of those reproduce their recorded rows digit for digit.
Trig+poly, seeds 1 to 20.

    arm                           chains  95% coverage        length  RMSE  held-out  min ESS (sum)  per chain  between
    independent75pool4            4x500   0.961(0.945-0.977)  4.61    1.12  1.15      15(8-31)       2(1-2)     0.78
    independent75pool4ruleGibbsB  4x500   0.939(0.911-0.958)  3.98    1.09  1.13      36(19-53)      3(2-4)     0.58
    independent75pool4equalCost   4x1105  0.952(0.925-0.970)  4.13    1.05  1.08      23(8-67)       2(1-2)     0.73

**2.21 times the sweeps do not buy what the kernel buys.** At one cut-scan budget - 2210 units either way, the full draw's 1000
sweeps at 2.21 against the shipped kernel's 2210 at 1.00 - the summed minimum ESS reads 36.3 against 22.9, so the sweeps buy 63
percent of what the kernel buys. Per 500 kept draws the equal-cost arm reads 10.3 against the control's own 14.8, this statistic
growing sub-linearly in chain length; per second, which carries NO timing claim on this host, 1.31 against the control's 1.75 and the
full draw's 1.19, arm B's wall ratio being 3.6 against its cut-scan 2.21 as the verdict above already recorded. Two readings beside
the ratio. The arm's between-chain ratio is 0.73 against the control's 0.78 where the kernel's is 0.58, so length barely moves what
the kernel moves; and its coverage falls to 0.952 from 0.961 with no kernel change at all, a third of the way to the reference's
0.941, which is the reference arm's finding confirmed from the length side. **This is 6.4's per-second question answered on the
instrument the design chose: at equal cut-scan cost the kernel wins.**

**The default-share candidate is the cut-only kernel at `d` = 0.16, recommended and not decided.** It is the only cut-only dose that
clears the +8 bar at both seed blocks, +11.3 and +9.1, with every secondary clean at 6.4's margins read at the reference: Trig+poly
coverage 0.955 and 0.952 against the reference's 0.941, held-out RMSE 0.992 and 1.012 against 1.02, and Single index +16.2 at a
held-out 1.033 that is 1.009 against that mean function's own reference ratio of 1.024. `d` = 0.08 fails its fresh block and `d` =
0.32 fails a gated secondary at both. And it is the cheapest thing on the table per unit of what it buys: 26.1 summed minimum ESS for
1040 cut-scan units against the control's 14.8 for 1000, the full draw's 36.3 for 2210 and the equal-cost arm's 22.9 for 2210, which
is 25.1 ESS per thousand units against 14.8, 16.4 and 10.4. Taking the equal-cost arm's own measured length scaling, 2.21 times the
sweeps buying 1.55 times the ESS, the cut-only kernel run out to the full draw's 2210-unit budget would read about 40 against its
36.3 - an extrapolation off a scaling measured on the shipped kernel and not on this one, and the reason the full draw's larger raw
gain does not settle the question the other way. Against the full draw the trade is the whole of section 2.4's: about half of
Trig+poly's gain and more than all of Single index's, for 4 percent added cost instead of 121.

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
MEASURED (6.1's third 2026-09-07 addendum): it supplies the weights the kernel must reproduce, the stump share section 2.4 prices
off, the nog and target-nog shares section 6's dosage argument reads, and the entropy columns section 2.4's variant recommendation
rests on. (3) P1's absolute gate, MET. (4) The C1
four-chain configuration, SETTLED. Nothing waits on perturb's slice 3.

1. **The kernel, at default weight 0.** The move, section 2.2's per-side rank and rank-admitted marginal on the scan's optional
   out-parameter, the dispatch branch and enumerator, `MoveContext`'s fifth probability, `structureIsFrozen`'s fifth argument -
   required, not defaulted; twelve chain.hpp and three combiner.hpp sites; eight bridge sites across seven declarations; six R
   default vectors, the slot, prototype and validity term, `rule_gibbs` resolved ahead of the three-name fill, both `all.equal`
   comparisons; seventeen positional `MoveContext` initializers and four positional `structureIsFrozen` calls across three tests/cpp
   files; four Rd, six tinytest and NEWS - **twenty-four files**. Section 6's cut-scan census hook and its column land here too:
   slice 3 cannot measure cost without them, and the ordinary build carries neither.
   Tests: bitwise neutrality; the census identity of section 5 against a reference assembly; an all-categorical no-op and a
   mixed-design test that the move fires on the ordinal-rule nodes and is a fixed point on the categorical ones; the monotone and
   variance-forest guards; and **a Gibbs-dominant walk on a two-column design small enough for the conditional to be enumerated in
   the test, where the realized draw frequencies must match it** - itself a correctness gate at the kernel level, and the one that
   would catch a mis-assembled weight before any script runs.

   **Landed** (7fb166ca, 2026-09-07). Three design-versus-code points. (1) [`scanOrdinalCuts`](../../src/bartcore/scan.hpp)
   emits the branch rank alone - the max of the two sides - with the marginal summed over the rank-0 sides, not the two
   sides' own ranks as section 2.2 specified; that pair is all section 2.2's law consumes. (2) The cut-only variant of
   section 2.4 was not built - not free, only cheaper than the full draw, and slice 1 took 2.4's own recommendation to
   build the full rule draw first. (3) [`census::nogProbe`](../../src/bartcore/moves.hpp) now calls the shared enumerator
   for its weights, but its own eligibility gate is left at section 2.3's option A (every available variable ordinal)
   rather than moved to the kernel's option C, so the census numbers already recorded stay comparable while the kernel
   itself ships option C.

   Gates (independent run): tests/cpp 281, including a rule_gibbs test covering the census identity at 1.25e-16 against
   an independent reference assembly, candidate-by-candidate rank agreement on plain and weight-stranded fixtures, an
   all-categorical no-op, and a Gibbs-dominant walk on an enumerable two-column design (chi-square 8.41 on 7 df against
   24.32, the same counts scoring 88 against a uniform draw), ASAN/UBSAN and census builds clean; tinytest 7947/0 (a
   rule_gibbs-dominant one-tree run that changes rules only at nog nodes and never a categorical rule, and an
   explicit-zero-versus-default `identical()` probe); equivalence trio bitwise 50/12/11 against the fbff1989 baselines,
   which stand; 20 quick exact gates green; the census package build reproduces the default cell's own 7.98 percent and
   -62.34; lint 0; R CMD check --as-cran OK, zero notes; NEWS 298 entries; API hash unchanged. Not landed:
   `rule-gibbs-balance.R` (slice 2) and the Stage 2 benefit run (slice 3).
2. **`rule-gibbs-balance.R`.** The prior-only arm on the 6-by-2 factorial at `power = 0.5`, the exact-posterior confirmation arm, both
   poisons.

   **Landed** (d888c9f3, 2026-09-07). Six design-versus-code points. (1) Section 5's `cgm(0.95, 0.5)` reads as (base, power)
   while [`cgm`](../../R/model.R)'s own signature is `cgm(power, base)`; the script calls `cgm(power = power, base = base)`
   at `power = 0.5, base = 0.95`, the pairing section 5's own arithmetic assumes. (2) Family size is `m = 23`, not the
   section's `m = 24`, strictest Holm threshold `|z| = 3.07`: the twelve-leaf state carries prior mass 0.00147, under the
   pre-stated 0.004 floor, so statistic 2 keeps ten singletons plus a pooled `leaves >= 11` bin. (3) Burn-in is 100,000
   sweeps per chain, not the section's 20,000, and is sized off the autocorrelation ladder rather than the leaf count: at
   `power = 0.5` the slowest series is the root's x2 cut indicator, at about 52 kept-draw lags to `acf < 0.1`. (4) Quick
   mode batches 125 x 400 rather than full mode's 500 x 500, batch LENGTH rather than count being what has to clear the
   roughly 30-draw integrated autocorrelation. (5) Section 5's constants otherwise reproduce exactly against the engine:
   0.095 / 0.475 / 0.05, the `(0.2981, 0.1346, 0.1346, 0.1346, 0.2981)` child law, the leaf-count dynamic program, and the
   undiluted detection floors. (6) Length 1012 lines against the section's 450 to 550 estimate,
   [`perturb-balance.R`](../../benchmarks/R/perturb-balance.R)'s shape. The gate is a nominal alpha-0.05 family-wise test:
   an off-seed probe produced one Holm rejection in one of four alternate seeds (`|z|` 3.31), no bias when the four are
   pooled, and the shipped seed passes both modes with margin.

   Gates (independent run): quick mode 16 seconds, prior-only arm worst `|z|` 1.17 (`x1c4|nog`), 0 of 23 Holm rejections;
   confirmation arm worst `|z|` 1.11 (`x2c1`), 0 of 7. Full mode 32 seconds, burn-in 100,000 sweeps (5,000 kept), slowest
   `acf < 0.1` lag 52 (root:x2c1) against the need of `<= 100`; prior-only worst `|z|` 1.74 (root `x1c3`), 0 of 23;
   confirmation worst `|z|` 1.64 (`x1c4`), 0 of 7, power against the prior marginal `|z|` 65.9; statistic 3's conditioning
   mass 0.0967 (realized 0.0964). Poison 1 (the two `1 - growth` factors dropped) FAILS as required, worst `|z|` 66.69, 5
   of 23; poison 2 (`1/|SI_v|` dropped) FAILS, worst `|z|` 111.66, 6 of 23. The matching engine mutations, m26 and m27 in
   `benchmarks/R/mutation-battery.R`, are both KILLED with this gate as sole killer (m26's confirmation arm alone would
   have passed at `|z|` 2.25, which is why the prior-only arm exists). lint 0, air clean, yaml parses, verify-anchors
   clean. Not landed: the Stage 2 benefit run (slice 3).
3. **The Stage 2 harness and run.** Two arms at `d` in `{0.16, 0.32}`, C1's four-chain cell, the sham arm, P1's rung, P2's
   duplicate-column cell as a control, the equal-cost arm, and the cut-scan count read off a census build at the arms' own seeds.
   Roughly 400 lines; compute on the order of a day plus the equal-cost arm's quiet machine. The verdict is recorded here.

   **Run** (50032833, 2026-09-07). Two arms on benchmarks/R/surfaces/C1-he-hahn.R's shipped four-chain configuration and one on
   benchmarks/R/surfaces/P2-confounded-step.R, 64 added lines across the two scripts and no engine change: 200 C1 fits over
   both mean functions and both seed blocks, 60 P1 fits, 40 P2 fits and one census cell. Not run: the sham arm, which perturb's
   own slice 3 already measured on this cell and control (-2.3 +/- 9.7 inside the +8 bar, at
   [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)), and the equal-cost arm, which
   needs a quiet machine and is owed. The cut-scan count came off the census build's own `c1` cell rather than off a replay of the
   arms' seeds; section 6's verdict paragraph states both departures.
4. **The default share.** Slice 3's kill did NOT fire; section 6's coverage secondary no longer blocks this item, the reference
   arm of 2026-09-08 having dissolved that flag with 6.4 now reading coverage against the cell's own well-mixed reference rather
   than against the shipped four short chains, which over-cover it; and the equal-cost arm this item was waiting on is run, at one
   cut-scan budget the kernel buying 36.3 summed minimum ESS where 2.21 times the sweeps buy 22.9. **The candidate is the CUT-ONLY
   kernel at `d` = 0.16**, not the full draw: the dose response of 2026-09-08 clears the +8 bar at both seed blocks on the gated
   mean function with every secondary clean at the reference-read margins, at 1.04 sweep-equivalents against the full draw's 2.21,
   and it is the best per unit of cost of anything measured on this cell. THREE things stand between that candidate and a default,
   none of them another measurement of benefit. (a) The decision itself, which is the maintainer's; nothing above is more than a
   recommendation, and it carries the surface question with it, the restricted kernel being a private `-D` build today rather than
   a mode anything can select. (b) 6.4's second kill clause, whose gate does not exist: it needs plateau prediction error in the
   noise-heavy or large-n stratum, which the grow-from-root harm battery measured and benchmarks/ does not contain
   ([5. Verdict and consequences](grow-from-root-default.md#5-verdict-and-consequences)). (c) The re-record. A nonzero default is a
   stream shift and pays for every RNG-locked baseline in
   [7. RNG and baselines](#7-rng-and-baselines)'s sense. POST-RELEASE.
