# level: an exact Gibbs draw on the level fibre

Status: PROPOSED, 2026-09-07.

A leaf-value step, not a tree kernel. Add a constant `c_t` to every occupied leaf of tree `t`, with the constants summing to zero
across the forest's trees: the fitted function is unchanged exactly, so the conditional of the shift vector on that subspace is the
leaf prior alone and the draw is closed-form.
[16.3 Ranking](tree-mixing-proposals.md#163-ranking) ranks it first of eight mechanisms in the second brainstorm round, with the
derivation confirmed and its one quantitative claim (a "1.7 sweeps" timescale) refuted as a scale mismatch. Its falsifier has since
run: [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s frozen-structure paragraph
freezes the structures and reads the leaf Gibbs alone, and at the MEDIAN point the deficit is overwhelmingly structural (ESS 14.9
unfrozen against about 670 of 2500 frozen), while at the WORST point it is not (minimum ESS 1.6 against 4 to 21 frozen, three
orders below the 2500 kept, the worst frozen point carrying two to three times the median posterior spread). So the leaf Gibbs
itself is slow at the worst coordinate, which is the coordinate the benefit stage's primary statistic is taken over.

**Premise, not reopened here.** The mixture is the one perturb landed at, `birth_death 0.6, swap 0, change 0.4, perturb 0,
birth 0.5`, and perturb's own benefit run was KILLED on this cell
([5.3 What arm B must produce, and the kill](perturb-move.md#53-what-arm-b-must-produce-and-the-kill)). Nothing here takes a share
from any structural move: the step is an addition to the sweep, not a redistribution inside it, and it changes no rule.

**The house already ships two moves of this shape, and they are the template.**
[`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) draws the common level the softmax cannot identify
from the leaf prior restricted to that direction, absorbs it uniformly over each forest's trees, skips empty leaves, and carries no
flag at all; [`rescaleAmplitudeRidge`](../../src/bartcore/combiner.hpp) travels the multiplicative amplitude ridge with a GIG draw
on the same principle. This is the ADDITIVE within-forest sibling of the first, over the `m - 1` directions a single forest carries,
and section 1's arithmetic reduces to the multinomial one under its own direction - the cheapest check that it is right.

## 1. The draw

**The subspace.** Let `Z` be the `n x sum_t L_t` leaf-indicator design over the forest's OCCUPIED leaves, so `f = Z mu`. Condition
on everything else: structures, occupancy, `sigma`, `k`, weights, latents, and `f` itself. On the affine set `{mu : Z mu = f}` the
likelihood is constant, so the conditional is the leaf prior restricted to `ker(Z)`. Each tree's indicator columns sum to `1_n`, so
every difference of two trees' block-indicators lies in `ker(Z)`: the level fibre is the `m - 1` dimensional set
`{(c_1..c_m) : sum_t c_t = 0}`, acting as `mu_{t,l} -> mu_{t,l} + c_t`.

**The conditional.** With independent leaf priors `mu_{t,l} ~ N(0, tau_{t,l}^2)`, write `P_t = sum_l 1/tau_{t,l}^2` and
`Q_t = sum_l mu_{t,l}/tau_{t,l}^2` over tree `t`'s occupied leaves. The log prior in `c` is
`-0.5 sum_t (P_t c_t^2 + 2 Q_t c_t) + const`, which factorizes, so the UNCONSTRAINED conditional is independent per tree,

    u_t ~ N(m_t, v_t),   m_t = -Q_t / P_t,   v_t = 1 / P_t

and the constrained draw is the gaussian conditioning of `u` on `1'u = 0`,

    c = u - v (1'u) / (1'v)

componentwise in `v`, giving `c ~ N(m - v (1'm)/(1'v), diag(v) - v v' / (1'v))`, rank `m - 1`, supported on `1'c = 0`. At the
homogeneous constant leaf every `tau_{t,l}` is `tau = node.scale / (k sqrt(m))`
([`ConstantGaussianLeaf`](../../src/bartcore/model.hpp): `scale` already carries the `1/sqrt(m)`), so `P_t = L_t/tau^2`,
`Q_t = S_t/tau^2` for `S_t` the tree's occupied leaf sum, and `m_t = -S_t/L_t`, `v_t = tau^2/L_t` - which is
[16.3 Ranking](tree-mixing-proposals.md#163-ranking)'s row 1 verbatim, re-derived. The general `(P_t, Q_t)` form is what section 3's
monotone leaf needs, its per-leaf prior sd differing within a tree.

**What it changes, and what it does not.** `f` is unchanged EXACTLY in the algebra and to rounding in the arithmetic
(`fl(mu + c_t)` rounds, and `sum_t c_t` is zero only to the accuracy of the projection). Every occupied leaf value changes. The tree
prior is untouched: [`CGMTreePrior`](../../src/bartcore/model.hpp) reads rule indices and grid shape and no leaf value at all
([16.2 The two structural facts this round rests on](tree-mixing-proposals.md#162-the-two-structural-facts-this-round-rests-on)).
Sigma, the latents and the response families see only `f`, and see it unchanged.

**Where the next sweep feels it, precisely.**
[`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp) reduces to
`0.5 log(P/(P+Q)) + 0.5 b^2 / (s^4 (P+Q))` with `P` the prior precision, `Q = W/s^2` the node's posterior precision, `W` its weight
sum and `b` its weighted residual sum. Tree `t`'s residual next sweep is its old one plus `c_t` (section 2), so every node of that
tree has `b -> b + c_t W`. `P` and `Q` read `W` only, so for a birth splitting `v` into `(l, r)` the likelihood ratio picks up

    linear:     (2 c_t / s^4) [ b_l W_l/(P+Q_l) + b_r W_r/(P+Q_r) - b_v W_v/(P+Q_v) ]
    quadratic:  (c_t^2 / s^4) [ W_l^2/(P+Q_l) + W_r^2/(P+Q_r) - W_v^2/(P+Q_v) ]

with `b_v = b_l + b_r`, `W_v = W_l + W_r`. `W^2/(P + W/s^2)` is convex in `W` and vanishes at zero, hence superadditive, so the
quadratic bracket is at most zero: a nonzero `c_t` PENALIZES births at tree `t` by `c_t^2` times a structure-dependent constant, and
the linear term reweights which cut wins. Nothing cancels. The magnitude is not a rounding perturbation - `c_t` is of the order of
one leaf value (sd `tau/sqrt(L_t)`, and `tau = 0.028868` on the internal scale at C1's shape), and `c_t W_v` is therefore comparable
to `b_v` itself - but no timescale is claimed here, 16.3 having refuted the one that was.

**Is the intercept all of `ker(Z)`?** No, generically yes but not always. `rank Z = dim(sum_t col(Z_t)) <= sum_t L_t - (m - 1)`, with
equality exactly when the only dependencies among the leaf indicators are the `m` copies of `1_n`. It fails whenever two trees induce
the SAME row partition (each matched column pair contributes one more direction) or one tree's partition coarsens another's (a coarse
cell's indicator is the sum of the fine ones). At 75 trees of 2.5 leaves over 30 columns and 100 cuts that is plausible and
UNMEASURED; nothing logs it. **A, the intercept directions only**: always valid, since a Gibbs step on a subspace of `ker(Z)` chosen
by a function of the CONDITIONED variables leaves the target invariant whether or not the subspace is all of it; cost as section 2.
**B, the full kernel**: needs `rank Z` per sweep, an `n x sum_t L_t` factorization, about `n (sum L_t)^2 = 3.5e8` flops at C1's shape
against a sweep's own 2.3e6 row touches - 150x a sweep. **C, the duplicate-partition extension**: hash each tree's `leafOf` row to
detect identically-partitioned pairs (one pass over `n` per tree, about a third of a sweep) and add the per-cell antisymmetric
directions, whose conditional at equal prior variances is `d_l ~ N(-(mu_{j,l} - mu_{k,l})/2, tau^2/2)` per cell, closed form.
**RECOMMEND A**, with C recorded as a door and its own generator (count duplicate partitions per sweep off `getTrees`) as the cheap
evidence that would justify it. B is refused on cost.

## 2. Placement and cost

**Once per sweep, not once per tree**: the constraint couples all `m` trees, so there is no per-tree form of the step. Per FOREST,
each forest's own zero-sum shift leaving its own fits invariant whatever multiplier a combiner applies.

**Where in the sweep.** **A, at the END of the sweep, after the leaf draws.** Costs three fixups. `forest.kSumSquaredParams` is
accumulated during the leaf draws and consumed by the `k` hyperprior at the tail of
[`Chain::run`](../../src/bartcore/chain.hpp), so a shift landing before it stales the statistic; the running residual `treeY` is
then stale by `-c_last`, which breaks the second identity [`runEnsembleTests`](../../tests/cpp/test_ensemble.cpp) asserts after every
sweep; and [`storeSavedTreeRecord`](../../src/bartcore/chain.hpp) flattens each tree INSIDE the tree loop, so a `keepTrees` record
would have to be patched or re-flattened. **B, at the TOP of the sweep, before the forest loop.** Every one of those disappears:
`kSumSquaredParams` is zeroed at the top of the forest loop and re-accumulated against the shifted leaves;
[`rollTreeResidual`](../../src/bartcore/chain.hpp)'s `t == 0` branch rebuilds the residual as `y - totalFits + mu_0[leafOf]` rather
than continuing the stale one, so `treeY` is dead across the sweep boundary and needs no pass; and every recorded channel of the
sweep is written after the shift. **RECOMMEND B.** It also composes cleanly with the two existing orbit moves, both of which fire in
the previous sweep's tail.

**No data pass, and no fits update.** `totalFits` is the cached `sum_t f_t` and is stale by `sum_t c_t`, which is zero, so it is
already correct for the shifted state. The roll then reads it and reconstructs tree 0's residual as `y - sum_{t>0} f_t^{new}` =
`y - sum_{t>0} f_t^{old} + c_0` - exactly right, and the reason the shift reaches the next sweep's structural scores with no
bookkeeping. Test fits are untouched, `f` being unchanged. So the step writes `forest.muByTree[t]` and nothing else, the same
wholesale-leaf-table write [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) makes.

**Cost.** One pass over the leaf tables to accumulate `(P_t, Q_t)`, `m` standard normals, one pass to apply: `2 sum_t L_t` leaf
touches plus `m` draws plus `O(m)` arithmetic. At C1's shape (`m = 75`, 2.52 leaves per tree by
[16.2 The two structural facts this round rests on](tree-mixing-proposals.md#162-the-two-structural-facts-this-round-rests-on)'s
third fact) that is about 450 operations against a sweep's own `3 m n = 2.3e6` row touches: **1/5000 of a sweep**. Against ONE cut
scan at a nog node (about `2n/L = 7900` members) it is 1/18, not 16.3's "under 1/100 of a cut scan"; the sweep-level conclusion is
unchanged and the row's cost column is optimistic by about 5x. Under multiple chains it is per chain, inside each chain's own loop
and on its own generator, so no barrier and no thread-count dependence.

**Degenerate cases.** `m = 1` gives a zero-dimensional fibre and the step is a no-op; a tree with no occupied leaf drops out and the
projection runs over the rest; fewer than two eligible trees is a no-op, mirroring `rescaleAmplitudeRidge`'s own `numLeaves < 2`
guard. The step fires in the sampling loop only, not in the grow-from-root warm start, which is an initializer and not MH-exact.

## 3. Leaf models and response families

**Occupied leaves only, and that is not a detail.** An empty leaf is forced to `0` by
[`sampleParametersAndSetFits`](../../src/bartcore/chain.hpp) rather than drawn, and a TEST row routed to it reads that zero. Shifting
it would move predictions at new `x` while leaving the training fit alone. So `L_t`, `S_t`, `P_t` and `Q_t` are over occupied leaves
and the empty ones keep their convention - which is also what
[`rescaleAmplitudeRidge`](../../src/bartcore/combiner.hpp) and
[`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) already do, in the same words.

- **Constant gaussian leaf: YES**, unchanged, the derivation of section 1.
- **Latent families (probit, logistic, multinomial, ordinal, nbinom, hurdle, t, AFT, hazard): YES**, unchanged. The leaves fit the
  working response and the shift leaves `f` fixed, so nothing outside the forest can observe it: an offset, an intercept held
  outside the forest, the ordinal thresholds, the nbinom dispersion, the Student-t degrees of freedom and every latent refresh read
  `f` and only `f`. Enumerated the other way, the only channels in the engine that read LEAF VALUES rather than fits are the `k`
  hyperprior (section 2's placement), the amplitude ridge, the monotone truncation, and the empty-leaf convention above. Multinomial
  gains a SECOND level step beside the one it ships: its cross-forest direction shifts every forest by a common `c` (softmax-
  invariant, `f` moves), this one is within-forest (`f` fixed). They are orthogonal and compose.
- **Heteroscedastic variance forest: SKIPPED, and not because its level is identified.**
  [`ConstantVarianceLeaf`](../../src/bartcore/model.hpp) is not a log-scale gaussian leaf: it is a MULTIPLICATIVE scaled-inverse-
  chi-squared factor, `s^2(x) = prod_j h_j(x)` through
  [`applyLeafFactor`](../../src/bartcore/chain.hpp). The fibre exists - multiply tree `j`'s leaves by `g_j` with `prod_j g_j = 1`
  and `s^2(x)` is exactly fixed - but the restricted conditional is a product of inverse-gamma densities on the log-sum-zero
  manifold, with no closed-form projection and no gaussian to condition. **A, skip.** **B, a PAIRWISE variant**: pick two variance
  trees, multiply one by `g` and the other by `1/g`; the conditional of `g` is then
  `g^{(L_k - L_j) nu'/2} exp(-b_j/g - b_k g)`, a generalized inverse gaussian, drawable with the
  `ext_rng_simulateGeneralizedInverseGaussian` the amplitude ridge already calls. **RECOMMEND A for this design**, B recorded as a
  door with its draw written out: the variance forest is a second forest whose mixing nothing has measured, and adding an untested
  kernel there would be priced against no deficit.
- **Monotone leaves: YES, with the per-leaf precision, and 15.2's exclusion does not bite.**
  [15.2 The two facts the lenses agreed on](tree-mixing-proposals.md#152-the-two-facts-the-lenses-agreed-on) scopes monotone forests
  out of fibre moves because the truncation reads leaf boxes and frozen neighbour values, both of which a REPARTITIONING move
  changes. This move changes neither. [`monotoneNeighborBounds`](../../src/bartcore/model.hpp) returns bounds that are max and min of
  NEIGHBOUR values in the same tree and never an absolute constant, so the cone is a set of pure difference constraints, invariant
  under a common shift of the whole tree, and its normalizer is a function of the structure alone and cancels. So the conditional of
  `c` is section 1's exactly - but the per-leaf prior sd is `cInflation * scale / k` for a leaf with a constrained neighbour and
  `scale / k` otherwise, so the homogeneous `(L_t, S_t)` form is WRONG here and the `(P_t, Q_t)` form is required. One eligibility
  test: a tree carrying an EMPTY leaf is skipped, because that leaf sits at 0 as a hard bound on its occupied neighbours and an
  occupied-only shift can leave the cone. Truncating the shift instead is refused: a truncated `c` is not a draw from the restricted
  prior and breaks the Gibbs law, which is the same reason the constrained prior draw uses rejection rather than a sweep of
  truncated conditionals.
- **Linear leaves: YES, restricted to the intercept.** [`LinearGaussianLeaf`](../../src/bartcore/model.hpp)'s
  `fitForObservation` is `params[0] + sum_j beta_j u_ij`, so `params[0]` adds to every row of the leaf and is the level; the slope
  block is not. Every coordinate shares the `scale / k` prior sd and the prior is independent across coordinates, so conditioning on
  the slopes leaves the intercepts' restricted conditional exactly section 1's, with `L_t` the occupied leaf count.
- **GP leaves: SKIPPED, and not on cost.** The parameters ARE the fitted values at the leaf's rows, and
  [`GPGaussianLeaf`](../../src/bartcore/model.hpp)'s `fitForTestObservationForNode` predicts a test row as
  `sum_r exp(-d^2/2) alpha_r` with `alpha = K^{-1}` times the drawn values. Kriging weights do not sum to one, so adding `c` to every
  training value does NOT add `c` to the prediction: the shift is not `f`-invariant out of sample and there is no level fibre to
  draw on. The in-sample conditional would additionally need `1'K^{-1}1` per leaf per sweep.

## 4. The surface

`proposal.probs` is the wrong home: it is a probability vector over STRUCTURAL proposals that must sum to one, and this is a logical
that adds a step rather than redistributing one. [`dbartsControl`](../../R/A_class.R) is the home - it already carries
creation-time sampler settings that are not run mechanics (`useQuantiles`, `n.cuts`, `storage`) - and the slot is `levelGibbs`,
default `FALSE`.

**A, always on now.** It is an exact Gibbs step that cannot hurt the target and has no tunable, and the two moves of its shape the
engine ships carry no flag; the cost is the whole re-record of section 7 taken before any benefit is measured, and a second one if
the benefit turns out nil. **B, an opt-in control flag at default off, bitwise neutral, flipped in slice 4 with the re-record then.**
**B', a private compile-time switch and no R surface at all**, the shape perturb's width arm was to be run in: slices 2 and 3 then
need a private library per arm, and nothing has to be un-shipped if the step is killed - which is exactly the residue perturb left,
a shipped knob whose measured benefit is nil and whose removal is now a maintainer's call.
**C, on by default above some chain count**: refused, a hidden rule that makes the kernel a function of a run setting.
**RECOMMEND B.** B' is cheaper if the step dies, but a compile-time switch is unreachable by the embedding consumers whose cells are
not C1 - a driver loop that freezes structure between response swaps is exactly where a level step should be reachable - and one
logical at `FALSE` is a far smaller residue than a fifth element in a probability vector with a fill rule. A is right only after the
benefit is in hand, which is what slice 4 is.

**Every site.** **R, three files.** [`dbartsControl`](../../R/A_class.R)'s slot, prototype and a validity length check;
[`dbartsControl`](../../R/dbarts.R)'s formal and its `newValidated` argument; [`bart2`](../../R/bart.R)'s formal and the control it
builds. `bart()` takes the default and only its Rd moves; `dbartsSpec` resolves the (control, model, data) triple and passes the
control through unchanged, so nothing in R/spec.R moves. **C++, two files.** A `SamplerOptions` field and the step itself in
[`Chain::run`](../../src/bartcore/chain.hpp); a `ParsedControl` field, the [`parseControl`](../../src/R_interface_bartcore.cpp) read,
the `SamplerOptions` copy and the verbose creation printout in the bridge. **It is a CREATION-TIME setting**, like `useQuantiles`
and `n.cuts`: [`bartcore_setControl`](../../src/R_interface_bartcore.cpp) pushes four settings through `SamplerBase` setters and
ignores the rest, so the flag needs NO new facade virtual and no `--preclean` hazard. **Rd, three files**: man/dbartsControl.Rd,
man/bart2.Rd (usage line and argument), man/bart.Rd's mention. **tinytest, three**: a new gate file, the control round-trip in
test-control-valuesAreUsed.R, and test-argument-surface.R's pinned formals. **tests/cpp, one**: test_ensemble.cpp. **Plus
inst/NEWS.Rd.** Thirteen files, roughly 60 lines of engine and 40 of surface.

**Neither the ABI nor the state moves.** The flat C header takes the control as a `SEXP`, so `DBARTS_C_API_HASH` is unchanged and no
`LinkingTo` consumer recompiles. `storeState` writes forests, sigma, scale, latents, DART, RNG, glue and digests and no sampler
option among them, so `stateFormatVersion` does not move; a shifted leaf value rides in [`FlatNode`](../../src/bartcore/tree.hpp)'s
`value` field, the same shape it had, and the round trip is bit-identical in structure and equal in value to the shifted state.

## 5. Correctness

The step is a Gibbs draw from a conditional in closed form: acceptance is identically one, there is no proposal, no reverse count and
no involution. **No detailed-balance gate is proposed, and that is a decision, not an omission**: `bd-balance.R` and
`change-balance.R` exist to test a Metropolis ratio against a target, and there is no ratio here. What replaces them is a direct
test of the conditional's law, which is strictly more informative about a Gibbs step than a marginal is.

**(a) Invariance and law, tests/cpp.** [`runEnsembleTests`](../../tests/cpp/test_ensemble.cpp) already asserts
`totalFits[i] == sum_t fits_t[i]` after every sweep at 200 trees, at a tolerance of `1e-11` against a measured worst deviation of
`1.2e-15`; re-run with the step on it IS the invariance test, the shift adding `m` roundings of about `eps` times a leaf value
(some `7e-16` at 200 trees), so it must hold at the same tolerance and within an order of magnitude of the recorded worst. A second
test freezes one forest, calls the step some `2e5` times against a seeded generator, and
scores the empirical mean and covariance of `c` against section 1's closed form, per-coordinate `z` at a pre-stated threshold, plus
`|1'c| < 1e-12` on every draw. **(b) tinytest**: with the flag on, the per-draw sum over trees of `getTrees` leaf values must
reproduce the recorded training fit (the deterministic precondition
["the gaussian backfit at ensemble scale"](../../benchmarks/R/backfit-exact.R) already runs), `predict` on held-out rows must agree
between the flag's two settings within Monte Carlo error, and a `storeState`/`setState` round trip must return the same leaf values.

**(c) Poisons, both sized.** (i) Drop the zero-sum projection and apply `u` directly: `sum_t u_t` is then
`N(sum_t m_t, sum_t v_t)` with `sum_t v_t = m tau^2/L`, sd about `tau sqrt(m/L)` = 0.16 at C1's shape, so `f` moves by that much at
every row and the invariance identity fails by 14 orders. (ii) Drop the prior mean term and centre `u` at zero: the shift's law keeps its covariance and loses its mean, so the
empirical mean fails its `z` while the covariance passes - which is the poison that distinguishes a correct conditional from a
correctly-shaped one.

**(d) The exact-posterior gates, and one that cannot pass.** The 22 exact gates run at the default and are untouched at the flag off.
At the flag ON, most are single-tree and the fibre is empty there, so they are inert rather than passing.
["the gaussian backfit at ensemble scale"](../../benchmarks/R/backfit-exact.R) is the exception and **it cannot be exact with the
step on, for a bookkeeping reason and not a correctness one**: its reference reconstructs tree `j`'s residual from the RECORDED
trees, taking trees `j+1..m` from the PREVIOUS sweep's record, which is pre-shift, while the engine drew tree `j` against those
values plus this sweep's shift. Trees `1..j-1` match, being recorded after it. So the reference is short by
`sum_{t>j} c_t` and could only be repaired by recording the shift as an `m`-vector channel of its own - a door, not this design's.
The gate therefore names the flag off explicitly and stays exactly as it is; the distributional test in (a) covers what it would
have covered, and slice 4 must state this rather than discover it.

## 6. Benefit, pre-registered

**The pilot runs FIRST, and it is the direct test.**
["How much of C1's autocorrelation lives in the leaf values"](../../benchmarks/R/surfaces/C1-frozen-ess.R) already reads the leaf
Gibbs alone at five seeds and two freeze points, minimum ESS 1.6 unfrozen against 4 to 21 frozen. Re-run it with the step on, same
seeds and same freeze points, paired. **Pilot kill: the paired ratio of frozen minimum ESS (on over off) must exceed 5 at both freeze
points on at least 4 of the 5 seeds, with the median-point frozen ESS not falling below its recorded ~600.** Below that, the level
fibre is not what pins the worst coordinate, the reading 16.3's row 1 rests on is refuted, and slice 3 does not run. It costs 10
short fits against slice 3's day.

**Slice 3, the confirmatory stage.** Arm A is the shipped kernel at the flag off; arm B is arm A plus the step, matched seeds,
paired. **Primary: the minimum ESS over C1's 25 fixed points, summed over four chains, on the independent design at 75 trees,
Trig+poly, at the shipped `n.chains = 4` at 500 + 500**, arm A reading 15(8-31), **at the +8 bar**
([5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic)).
Secondaries, each must-not-degrade at
[6.4 What "no regression on the core" means numerically](benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)'s
own margin: 95 percent coverage at -0.010 absolute (arm A 0.961), held-out RMSE at a ratio above 1.02, per-chain minimum ESS
reported beside the primary, and wall time per sweep at a ratio above 1.05 - which section 2's cost makes a formality, but it is the
clause that stops an ESS win bought with compute. Controls as perturb has them: **the sham arm**, arm A against itself at sampler
seeds offset 1000 on the same twenty data seeds, whose paired difference must sit inside the bar (it read -2.3 +/- 9.7, paired SE
2.17, so the bar is four times 2.17 = 8.7 and +8 is marginally optimistic rather than wrong); and **the P1 rung as the absolute
gate**, whose control arm must return 90 percent coverage near 0.71 and 0.725 held-out
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)),
or no verdict from any cell is valid.

**Kill criterion.** **KILL if, at the shipped four-chain configuration, arm B does not improve C1's summed minimum ESS over arm A by
more than +8, over at least 20 matched pairs, with wall per sweep inside the 1.05 ratio, with the sham control inside the bar, and
with a mandatory fresh-seed re-run of any flagged cell before a flag counts**
([6.1 The rule, stated operationally](benchmark-surfaces.md#61-the-rule-stated-operationally)). A must-not-degrade secondary past its
own margin kills independently. The design does not predict the bar will be met: the pilot measures the leaf half in ISOLATION, and
nothing connects a frozen-chain gain to the same gain under a live structural kernel, where 10.4 reads the median point's deficit as
overwhelmingly structural.

**How to read it against the He-Hahn pooling finding.** 10.4's coverage of 0.961 comes from four chains sitting in DIFFERENT places
(between-chain ratio 0.78) and pooling widening the interval. A level step mixes within a chain, so each chain should cover more of
the posterior and the between-chain ratio should FALL - but only together with a rise in summed minimum ESS. **Between falling with
ESS flat is the perturb signature** (0.67 and 0.50 against 0.78, with paired ESS differences of +3.1 and +0.1) and reads as lost
chain diversity, not gained exploration: it is a fail. Read the three columns together, and read coverage against its -0.010 margin
in every arm.

**P2 is untouched.** [10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)'s rooting
lock is a connectivity failure in TREE space - two rootings separated by rule changes - and the level fibre lies entirely in leaf
space, so no leaf step can cross it. P2 stays a must-not-degrade cell (switches per chain within arm A's own seed range, parked
chains no worse) and is not a cell the step can win on.

## 7. RNG and baselines

**At the flag off, bitwise neutral by construction**, not by coincidence: the step is guarded before any generator call, so no
stream advances and benchmarks/R/equivalence.R stays green with the three `fbff1989` baselines standing. That is what lets slices 1
and 2 land and the correctness gate run before anything changes for users.

**At the flag on, every draw moves**, from the first sweep of the first chain, and slice 4 pays
[4. RNG, and the one re-record](swap-removal.md#4-rng-and-the-one-re-record)'s list in one bundle: the three equivalence baselines
`equivalence-fbff1989.rds`, `bcf-equivalence-fbff1989.rds` and `multinomial-equivalence-fbff1989.rds`, each a new file with a
MANIFEST row marked current and the workflow pins moved; feature-matrix.md's Evidence paragraph, whose counts
[`makeScenarios`](../../benchmarks/R/equivalence.R) and its two siblings recompute; the gate ledger's baseline footnote; and the four
seeded-drift tripwires test-reproducibility-continuousResponse-singleThreaded.R, -multithreaded.R, -binaryResponse.R and -xbart.R,
regenerated by ["Regenerates the reference values"](../../tools/regenerate-snapshots.R). bench-sampler is a COMPARE, not a re-record:
section 2's cost cannot move a 1.05 ratio. Consumer snapshots on the lockstep branches re-record with it. The bitwise oracle is the
usual one: record the baselines at a throwaway build with the flag defaulted on, then `compare` from the landed tip and require
identical draws.

## 8. Slices

**Prerequisites and order.**

1. **The perturb kernel and its kill, LANDED** (`ab49f83a`, `d73fb4e0`): section 4's surface counts and section 6's arm A are both
   against that tree, and the C1 four-chain arm slice 3 reuses is already in the harness.
2. **The P1 absolute gate, MET** (0.725 held-out in the control arm). Slice 3 re-runs the rung as its own control.
3. **The C1 chain configuration, SETTLED**: four chains at 500 + 500, summed minimum ESS, the +8 bar and its sham calibration.

Then, in order:

1. **The step behind the flag, at default off.** The draw and its projection at the top of the sweep, the `SamplerOptions` field,
   the control slot and its four R sites, the bridge parse, three Rd files, NEWS; tests/cpp's ensemble invariance arm and the
   frozen-forest distributional test with both poisons; tinytest's fit-sums-leaves precondition, predict agreement, state round trip
   and control round trip. Roughly 60 lines of engine, 40 of surface, 150 of tests, thirteen files. Gates: tests/cpp green,
   tinytest 0 fail, the equivalence trio bitwise against the standing baselines, the 22 exact gates unchanged, ASAN/UBSAN clean.
2. **The frozen-structure pilot.** `C1-frozen-ess.R` re-run with the step on, five seeds, two freeze points, paired against the
   recorded table. Its kill is section 6's; nothing below runs if it fires. Ten fits.
3. **The benefit run.** Two arms on C1's four-chain cell at twenty matched pairs, the sham arm, the P1 control rung, the P2
   must-not-degrade cell. About a day of compute; the verdict is recorded here.
4. **The default flip, and the one re-record.** Only on a slice 3 pass. It carries section 7's bundle, the
   backfit-exact.R clause of section 5 (d), and
   [6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered)'s independent default clause, which
   needs plateau prediction error in the noise-heavy or large-n stratum that no cell of slice 3 measures. On a slice 3 fail the step
   stays at default off as an opt-in whose measured benefit on the pre-registered cell is nil, and whether it is removed before
   release is the maintainer's call.
