# multiforest-veto-rate-falsifier

status: RUN AND REPORTED (2026-08-09). Verdict: **YELLOW for both column
  types** - P1, P3 and L1' GREEN, M1 exactly 0, every validity gate passing,
  T1 YELLOW on its ratio clause alone. No KILL clause fired anywhere, so no
  fresh-seed confirmation was triggered. See "Results". Ratified design;
  supersedes the prior (untracked) memo and its critique wherever they
  differ.
agent: opus (harness + analysis; no engine change)
rng: neutral - measurement only. A GREEN verdict AUTHORIZES the
  `multiforest-predictor-mutation` design arc; it does not perform it.
budget: harness only;
  nothing under `R/`, `src/`, `inst/`. ~700-850 lines across five files.
window: none - no timing metric exists, so the measurement is
  load-insensitive. Runs on the laptop.
tip: d3cb94b (branch bartcore; engine tip 06c0254). Every anchor below
  re-verified at that tip.

## Goal

Price the acceptance-rate consequence of widening dbarts'
transactional and per-observation predictor-mutation veto from
`forests_[0]` alone to every tree ensemble, BEFORE the design arc opens.
The verdict gates `multiforest-predictor-mutation` (TODO: the model-space
survey's highest-value door).

Three forks, in the order they bind. Only the first is decided here.

1. **Ship the straightforward extension?** Widen `Chain::revalidateTrees`
   and `UpdateSessionImpl::treeAt` from `forests_[0]` to every ensemble,
   so the acceptance criterion becomes "no leaf empties in any tree of any
   ensemble of any chain". DECIDED HERE.
2. **Per-forest or per-sampler veto?** REPORTED ONLY (see "What this
   measurement cannot tell us", item 3).
3. **Is a different mechanism required?** Reached only on KILL.

## Context

- Widened refusal at the time of this study: `refuseMultiForestTransactionalUpdate`
  retired: [[src/R_interface_bartcore.cpp#refuseMultiForestTransactionalUpdate]]
  and the separate
  `refuseVarianceForestPredictorMutation`
  retired: [[src/R_interface_bartcore.cpp#refuseVarianceForestPredictorMutation]]
  - both gone; multiforest-predictor-mutation.md later retired the whole
  transactional setPredictor/updatePredictor refusal these guarded, in favour
  of accepting the mutation with rollback. Both fired only
  when `!forcedUpdate`; the FORCE paths stayed open and refreshed
  every forest.
- Forest-0-only revalidation: `Chain::revalidateTrees`
  ([[src/bartcore/chain.hpp#Chain::revalidateTrees]], `Forest& forest = forests_[0]`).
  Two-phase all-chain transaction: `Sampler::revalidateAllChains`
  ([[src/bartcore/sampler.hpp#Sampler::revalidateAllChains]]).
- Per-observation session: [[src/bartcore/sampler.hpp#UpdateSessionImpl]]; predicate at
  [[src/bartcore/sampler.hpp#UpdateSessionImpl::observationWouldRemainValid]]; tree enumeration at
  [[src/bartcore/sampler.hpp#UpdateSessionImpl::treeAt]]
  (no longer literally `chains_[t / numTrees]->tree(t % numTrees)`; now a cached
  per-slot chain/forest/tree lookup, and the general case is
  [[src/bartcore/chain.hpp#Chain::treeInForest]] - the variance forest instead
  reads [[src/bartcore/chain.hpp#Chain::varianceTree]]).
- Joint multi-sampler sweep: [[src/bartcore/facade.hpp#updatePredictorPerObservationJointly]]; scan permutation from
  `samplers[0]->rng()` there too.
- The unit of study is **not** `numForests`. A heteroscedastic sampler
  keeps `numForests == 1`; the variance forest lives outside `forests_`
  and is signalled by `SamplerShape::hasVarianceForest`. Use
  `E = numForests + (hasVarianceForest ? 1 : 0)` - the number of ensembles
  routed off the shared column store. This confirms
  [[docs/design/model-space-survey.md#2. What the engine serves today]] and
  [[docs/design/model-space-survey.md#D3. Public BCF creation surface - CLOSED (2026-08-10 to 2026-08-11)]]; it does not correct
  it.
- **The variance forest carries the same empty-leaf invariant as forest
  0**: `refreshVarianceForest` ([[src/bartcore/chain.hpp#Chain::refreshVarianceForest]]) rejects a rebuild whose variance trees have an
  unoccupied bottom and asserts it -
  "the empty-leaf veto admits no unoccupied bottom into live state". So a widened veto adds no constraint the joint
  model does not already impose; it only makes the sampler enforce one it
  already implies. This is the whole basis of the re-grounded kill logic.
- `E = 3` (heteroscedastic BCF) is **not constructible**:
  `createBCFSampler` returns nullptr on `numVarianceTrees > 0`
  ([[src/bartcore/facade.hpp#createAmplitudeSampler]]), multinomial likewise
  ([[src/bartcore/facade.hpp#createMultinomialSampler]]), both asserted at
  [[tests/cpp/test_model.cpp#testVarianceForestRefusal]].
- **A tree with no split on the mutated column cannot veto - exactly.**
  In `UpdateSessionImpl` the override never fires
  ([[src/bartcore/tree.hpp#findBottomNodeForObservation]]), so `newLeaf == oldLeaf` and `valid` is untouched
  ([[src/bartcore/sampler.hpp#UpdateSessionImpl::observationWouldRemainValid]]); in `revalidateTrees` the repartition
  reproduces the same partition. Per-forest column masks are therefore an
  exact opt-out (tau moderators [[src/bartcore/combiner.hpp#Forest::columnMask]];
  variance [[src/bartcore/chain.hpp#VarianceForest::columnMask]]; generic
  [[src/bartcore/combiner.hpp#Forest::columnMask]]).
- Consumer baseline (bairrtt, a separate repository this guard does not
  cover - cited by line number in prose, not machine-checked): bairrtt runs
  the widened predicate at two samplers
  today - `updatePredictorPerObservationJointly(list(response_model,
  assignment_model), theta[, k], theta_names[k])`
  (`irt_causal_bart.R` lines 609-613), 75
  trees each (line 204), one dbarts chain each (line 432), so 150 trees per
  row already. Burn-in uses `setPredictor(..., forceUpdate = TRUE)`
  (lines 619-628), so no veto is exercised during burn. Rejected rows are
  reverted as rejected MH moves (line 614); no counter is kept. Identical cut
  grids are pinned on the shared column (lines 452, 464). The only recorded
  rate is prose - bairrtt's `TODO` lines 117-118, "The empty-leaf install rejects
  under 1% of moves", at 300 persons x 30 items (line 110) - not reproducible
  from any stored artifact. Denominator: the MH accept filter (lines 598-599)
  precedes the install (line 609) and `accept_rate` at line 615 is over all n,
  so read "moves" as the ~0.44n accepted rows unless stated; both readings
  are reported.
- Motivating classes and their column types:
  [[docs/design/model-space-survey.md#Door 2 - per-forest row subsetting]]. BPCF principal strata enter
  tau_Y as splitting covariates and are DISCRETE for a binary intermediate
  (full-text verified there with
  `https://arxiv.org/html/2403.13256v1` and `github.com/lit777/BPCF`).
  Sequential covariate imputation and treatSens-style latent confounders
  admit binary columns too. Hence the per-column-type verdict below.
- Declined alternative, for the KILL branch:
  [[docs/design/empty-leaf-veto.md#Why not make the proposals occupancy-aware]] prices occupancy-aware proposals at
  a "250-400 line, posterior-changing rewrite".
- Precedent for falsifier-before-arc: `docs/design/memory-wall-frontier.md`
  sec 11. Form of this file: `docs/plans/archive/grow-from-root-default-study.md`.

## Constraints

- **Read-only on the package.** No engine edit, no default change, no
  baseline re-record. `R CMD INSTALL .` first (tinytest runs against the
  installed package); no `--preclean`, because nothing in `src/` changes.
- **`rngSeed` pinned in every fit.** Chain seeds come from a dedicated
  Mersenne twister and sampling never advances R's stream
  ([[src/R_interface_bartcore.cpp#createChainRngs]]), so an identically built and
  identically driven sampler is bitwise reproducible. Two seeds per cell,
  as in `grouped-mixing.R`: data `set.seed(BASE_SEED + s)`, sampler
  `rngSeed = s`; matched arms share `s`.
- **No NA in the mutated column, ever.** The emulation (below) relies on
  `collapseEmptyNodes` and `dropStaleMissingDirections` being no-ops.
- **Arm B pins identical cut grids** on the shared latent column across all
  surrogate samplers (`setCutPoints(cuts, nm)` on each). The joint path
  quantizes against each sampler's own store
  (`sampler_.data_.codeFor`, [[src/bartcore/sampler.hpp#UpdateSessionImpl]]), so without pinning
  the S-sampler conjunction is not the shared-store conjunction it stands
  in for. V1 asserts code-vector equality.
- **Thresholds FREEZE at the end of Stage 0**, written into this file
  before any Stage 1 gated number is read.
- No timing metric enters any criterion.

## The emulation, and why arm A is a closed loop

The widened criterion does not exist in the engine, but it can be driven
exactly from R without patching anything:

1. Export every ensemble's live trees for the current sweep -
   `bartcoreGetTrees(bc, chainNums, treeNums, current = TRUE, forest = f)`
   retired: [[R/bartcore.R#bartcoreGetTrees]] (no longer exists under this
   name in R/; superseded by the reshaped R5 tree-reading surface)
   for `f = 0..numForests-1`, plus
   `state[[c]]$variance.{vars,values,sizes,flags}`
   ([[src/R_interface_bartcore.cpp#storeState]]) for a heteroscedastic sampler.
2. Compute the widened per-observation install mask in the R oracle.
3. Install exactly the accepted rows with
   `setPredictor(newColumn, col, forceUpdate = TRUE)`. The force path is
   open to multi-forest and heteroscedastic samplers
   ([[src/bartcore/chain.hpp#Chain::forceRefreshTrees]]) and `forceRefreshTrees` re-routes
   every forest and the variance forest.
4. `run(0L, 1L)`. Repeat.

Because the oracle admits no row that would empty a leaf,
`collapseEmptyNodes` is a no-op and the forced refresh reproduces exactly
the live state the widened transactional path would install. **V3b
certifies this deterministically**; failure of V3b, not any tolerance, is
what triggers the engine-patch contingency.

The export is cheap: `flattenBelow` reads `node.numObservations()` off the
node ([[src/bartcore/tree.hpp#flattenBelow]]) rather than re-routing, so
`getTrees(current = TRUE, newdata = NULL)` is O(nodes), not O(n x nodes).

Consequence: every arm in this design is a **closed loop at
stationarity**. There is no one-step counterfactual anywhere, and no
frozen-tree replay.

## Arms

**Arm 0 - the shipped widening, zero counterfactual.** Single-forest
sampler driven with `setPredictor(v, col, forceUpdate = "partial")` and
with `updatePredictorPerObservationJointly`. Raising `n.chains` from 1 to 4
widens the criterion by exactly the arithmetic the arc would
(`revalidateAllChains` already requires every tree of every chain to stay
occupied). Establishes `r_1` and the transferable count law with no
assumption at all, and cross-checks bairrtt's prose figure.

**Arm A' - driven closed loop on real multi-ensemble fits (primary for
shapes).** Real BCF (`dbarts:::bartcoreBCFSampler`, [[R/bartcore.R#bartcoreBCFSampler]])
and real heteroscedastic (`dbarts(x, y, variance = TRUE,
n.trees.variance =)`, [[R/dbarts.R#dbarts]]) samplers, driven by the
emulation above. Right trees AND right dynamics.

**Arm B - closed-loop S-sampler surrogate (calibration, E = 3, and the one
engine-validated attribution).** The bairrtt pattern with `S in {1,2,3}`
samplers sharing one latent column through
`updatePredictorPerObservationJointly`, non-primary samplers shaped to
match the target ensemble (sampler 2: `n.trees = 50`, `base = 0.25`,
`power = 3`, moderator columns only; sampler 3: `n.trees = 40` and the
variance prior). Named limitation: the surrogate's ensembles fit `y`
independently rather than sharing the combiner's residual, so their trees
are not the trees a real BCF grows. Arm A' supplies those.

*Engine-side per-sampler attribution, exact.* The joint call draws its scan
permutation from `samplers[0]->rng()` ([[src/bartcore/facade.hpp#updatePredictorPerObservationJointly]]), so with `rngSeed`
pinned, replaying the identical build-and-burn script and varying only the
SUBSET passed to the joint call yields scan-order-matched masks. The
marginal of each surrogate sampler is therefore measured by the engine, not
inferred. This attributes per SAMPLER, not per FOREST inside a sampler.

**Arm C - validity.** V0-V4 below. Runs first; failure blocks every metric.

## Scenarios

| factor | levels | grounding |
|---|---|---|
| ensemble config | (a) single forest 75 trees [BASELINE, E=1]; (b) BCF mu 75 / tau 50 [E=2]; (c) BCF mu 75 / tau 25 [E=2]; (d) het mean 75 / variance 40 [E=2]; (e) surrogate S=3 mu 75 / tau 50 / var 40 [E=3, arm B only, UNGATED]; (f) multinomial K=4 x 75 [E=4, UNGATED stress] | `n.trees = 75` default ([[R/A_class.R#dbartsControl]]); `n.trees.treatment = 50L`, `treatment.base = 0.25` / `treatment.power = 3` ([[R/bartcore.R#bartcoreBCFSampler]]); `n.trees.variance = 40L` ([[R/spec.R#resolveSamplerSpec]]); mu default `base = 0.95, power = 2` ([[R/model.R#cgm]]); bairrtt runs 75 (`irt_causal_bart.R` line 204) |
| **column type** | {continuous, binary} gated separately; {multi-level categorical <= 64 levels} reported ungated | two of four motivating classes mutate discrete columns |
| mutated column in the non-primary mask? | {yes, no} | "no" must measure exactly zero marginal (M1); "yes" is the gated cell |
| n | {300, 500, 1000, 5000} (300 in arm 0 only) | bairrtt runs 200-1000 and its prose baseline is at 300; the motivating DGP is 1000 persons x 100 items (bairrtt's `docs/plans/multi-trait.md` lines 98-102); 5000 probes leaf-size scaling |
| p | {5, 10}, held at 10 outside arm 0 | controls the share of j-splitting trees |
| n.chains | {1, 4}, held at 1 outside arm 0 | arm 0's lever and the shipped widening |
| phase | 300 forced burn sweeps, then a W = 200 sweep TRACKING WINDOW | bairrtt deliberately forces during burn (`irt_causal_bart.R` lines 616-629), so only the tracking window is gated |
| fraction of rows mutated | {1/n, 0.01, 0.1, 1.0} | 1.0 is bairrtt (`irt_causal_bart.R` line 563); 1/n and 0.01 are the sequential-imputation and one-subject-Gibbs regimes, the only ones where the whole-transaction path is live |
| move size | continuous: {0.1, 0.4, 0.8} column SDs; binary: flip probability 0.5 on proposed rows | bairrtt's adapted `theta_sd` is 0.68-1.04 against a unit-SD column (`irt_causal_bart.R` line 197; its `TODO` lines 113-114) |
| cut grid on the mutated column | 100 interior quantiles (continuous) | `n_theta_cutpoints = 100L` (`irt_causal_bart.R` line 203) |

**Realistic gated cells**, per column type: configs (b), (c), (d); latent IN
the non-primary ensemble's mask; n in {500, 1000, 5000}; the tracking
window; move size 0.8 (continuous) or flip 0.5 (binary); fraction 1.0 for
per-observation metrics and 0.01 for the whole-transaction metric. Nine
cells per column type. Configs (e) and (f), all `mask = no` cells, all
n = 300 cells, fraction 1/n, and move sizes 0.1 / 0.4 are REPORTED,
UNGATED.

**W is fixed at 200 transactions per row in every gated cell**, all n. The
n = 5000 cells lose seeds (6 rather than 8), never transactions: a
">= 50% of transactions" tail functional is not comparable across
replication counts.

## Metrics

`T_j` = j-splitting trees summed over ensembles and chains. `r_1` = the
per-tree per-row marginal reject probability.

| id | quantity | contrast | role |
|---|---|---|---|
| P1 | mean per-observation install rate | `D = mean(install \| config (a), shipped) - mean(install \| config (b/c/d), E_all)`, matched (n, DGP, move, seed, W), in pp | **PRIMARY, gated** |
| P1d | the same, against `E_0` on the multi-ensemble fit | diagnostic only | **UNGATED** - `E_0` on a multi-forest sampler is not a shippable option (see below) |
| P3 | share of rows rejected in >= 50% of the W tracking transactions | `share(config b/c/d, E_all) - share(config a, shipped)`, in pp | **PRIMARY, gated** (mixing pathology) |
| L1' | `observed / (1 - (1 - r_1^(a))^(T_j))`, with `r_1^(a)` estimated on config (a) ONLY | ratio | **PRIMARY, gated** (priceability) |
| T1 | whole-transaction veto rate at fraction 0.01 | `rate(E_all, b/c/d)` vs `rate(shipped, a)`, ratio and level | **SECONDARY, gated** |
| M1 | marginal contribution of an ensemble whose mask excludes j | must be **exactly 0** | correctness assertion; non-zero STOPS the measurement |
| P2 | ratio of the 99th percentile per-row reject rate to the mean | reported | concentration |
| P4 | share of widened per-obs rejections in which ensemble f is the SOLE objector | per f | **REPORTED, not decided** |
| T2 | per-ensemble marginal `rate(E_0 u {f}) - rate(E_0)` | per f | **REPORTED, not decided**; engine-measured for arm B only |
| C1 | arm A' vs arm B install rate, matched E = 2 cell | pp | cross-arm calibration (V4) |
| Q1 | stationary per-row latent posterior under three kernels (widened veto / forced collapse / frozen-forest MH-only), arm B S=2: mean, sd, pooled KS | **REPORTED, never gated** | bounds the posterior question on a surrogate |
| X1 | mean tree depth and leaf-size distribution per ensemble | reported | explains the rest |
| X2 | `r_1` per ensemble | reported | the transferable constant |
| X3 | log-log slope of reject rate on `T_j` | **reported, NOT gated** | sublinearity is good news, not a failure |

**Why P1 is baselined on config (a) and not on `E_0`.** The guard's own
comment says an accepted transactional change on a multi-forest sampler
"would leave a multi-forest sampler's other forests routed against stale
codes" ([[src/R_interface_bartcore.cpp#refuseMultiForestMutation]]). `E_0` on a multi-ensemble
sampler is the state the guards exist to prevent, not an alternative a
consumer can have. The alternatives are today's refusal, forced collapse (a
different posterior), or a different mechanism. So the gated contrast is
"what a consumer pays to add a tau or variance forest", i.e. config (a)
under the shipped criterion versus config (b/c/d) under `E_all`. P1d is
kept, ungated, because it is literally what the code change toggles.

**Why P3 exists.** The mean is the wrong summary. Rejections are not
exchangeable across rows: a row that sits alone in some leaf tends to keep
sitting alone, so a modest mean rate can hide a persistently frozen
subpopulation whose latent posterior is wrong for identifiable individuals,
and running longer does not fix it. That is a different failure from
"mixes a bit slower" and it needs its own line.

**Why the old count-law check was removed.** `observed / (T_j x r_1) <= 1.5`
cannot fail: the per-row rejection event is a union of `T_j` per-tree
events, so by Boole `observed <= sum_t r_1t = T_j x r_1` identically. And a
`slope >= 0.7` clause is mis-signed - row-level heterogeneity in `r_1`
makes `mean_i [1 - (1 - r_1(i))^T]` saturate and drives the log-log slope
below 1, so the cheapest possible widening would have blocked GREEN. L1'
replaces both: because `r_1^(a)` is estimated off-configuration, `ratio > 1`
is attainable whenever the added ensembles have tighter leaves or couple
with mu, which is exactly where the count law can break. L1' is fit and
gated on arm A' too - the only coupled-forest regime.

## Validity gates (Stage 0; failure blocks every metric)

- **V0 - harness self-check.** The two pre-registered exactness-preserving
  speedups (descend only j-splitting trees; precompute old and new leaf per
  tree vectorized) agree with the full router on one small cell, exactly.
- **V1 - no-op identity, and grid pinning.** A proposal equal to the
  current column installs every row and vetoes no transaction, in the R
  oracle and in the engine, for every ensemble config. In arm B, the
  per-sampler code vectors for the shared column are identical after
  pinning.
- **V2 - the routing oracle, exact (mean and BCF forests).** For every tree
  of every forest of every cell, the R router's per-node observation counts
  equal the `n` column of `bartcoreGetTrees(..., current = TRUE)` exactly,
  all nodes, integer equality. This retires the value-vs-code question: the
  engine routes by integer codes and the export reports cut values, and
  they induce the same partition (`codeFor` is
  `lower_bound(cuts, value) - first`, [[src/bartcore/data.hpp#ColumnStore::codeFor]]; ordinal
  `ruleSendsRight` is `code > splitIndex`; the exported payload is
  `cutPoints[j][splitIndex]`, [[src/bartcore/tree.hpp#flattenBelow]]; flat routing is
  `value <= flat.value -> left`, [[src/bartcore/tree.hpp#atOrUnderCut]]). The `n` column is
  the LIVE partition, not a replay, when `current = TRUE` and
  `newdata = NULL` ([[src/R_interface_bartcore.cpp#bartcore_getTrees]]).
- **V2r - per-row identity.** Equal counts do not imply equal per-row
  assignment (a permutation between two equal-sized leaves passes V2).
  Replay single-row `newdata` matrices through `getTrees` for 50 sampled
  rows per cell; the unit count must land in the oracle's predicted leaf.
- **V2v - the variance-forest oracle, exact.** `bartcoreRun(bc, 0L, 1L)`
  returns an element named `variance`
  ([[src/R_interface_bartcore.cpp#bartcore_run]]) filled from
  `varianceFits[i] = sigmaScale^2 * combinedVariance[i]`
  ([[src/bartcore/chain.hpp#Chain::varianceFits]]), where `combinedVariance[i]` is the PRODUCT over
  variance trees of row i's leaf factor ([[src/bartcore/chain.hpp#VarianceForest::applyLeafFactor]]);
  the state's `variance.values` are those same working-scale factors
  ([[src/bartcore/chain.hpp#Chain::varianceFactorsForTesting]]). So: decode `variance.*`, route every row, form
  the per-row product, and require `engine$variance[i] / oracleProduct[i]`
  to be the SAME constant across all i to 1e-10 relative. Run at
  `n.trees.variance = 1`, where this is an exact per-row per-TREE
  leaf-assignment check, and at the default 40, where it is an exact check
  of the joint variance cell. This replaces the structural round-trip and
  removes the need for an engine patch.
- **V3a - whole-transaction verdict, exact.** On config (a) the verdict
  from `sampler$setPredictor(v, col)` matches the oracle's verdict on every
  one of >= 1000 proposals. Deterministic; no scan order enters.
- **V3b - emulation certificate, exact. Arm A' is inadmissible without
  it.** Build two config-(a) samplers identically with `rngSeed` pinned. On
  the first, take the mask M from one `forceUpdate = "partial"` call. On
  the second, install exactly the rows of M via
  `setPredictor(..., forceUpdate = TRUE)`. The two `state` blobs must be
  identical (trees, values, sizes, flags) and the two subsequent
  `run(0L, 1L)` outputs bitwise equal. Repeat on 20 independent fits and on
  one config-(b) fit driven through both paths' shared prefix.
- **V3c - distributional agreement on the partial mask.** 200 engine
  replicates (force-restore between each - the force path is unconditional,
  so the restore is never itself vetoed) against 200 oracle scan-order
  draws on the same fit. Two-sample difference of means, pre-registered 95%
  CI half-width **<= 0.4 pp**. Arithmetic: at p = 0.02, n = 500, the
  per-replicate SD is `sqrt(0.02 x 0.98 / 500) = 0.626 pp`; SE of a
  200-replicate mean is 0.044 pp; SE of the difference of two 200-draw
  means is 0.063 pp; half-width 0.124 pp. At the grid's worst case
  (p = 0.25, n = 500, SD 1.94 pp) the half-width is 0.38 pp, still inside
  0.4 pp - and 0.4 pp is a fifth of the 2 pp primary line.
- **V4 - cross-arm calibration, REPORTED.** Arm A' and arm B install rates
  in the matched E = 2 cell. Both are now closed loops, so any gap is a
  TREE-SHAPE difference (real coupled BCF trees versus independently fit
  surrogates), not a one-step bias. If `|C1| > 3 pp`, report both, **gate on
  arm A'** (right trees and right dynamics), and record the surrogate's
  shape bias as the arc's residual risk. This inverts the earlier design,
  which gated on arm B.

## Pre-registered decision lines

Verdicts are issued **per column type** (continuous, binary), separately.
A verdict covers only the shapes that exist: BCF at E = 2, heteroscedastic
at E = 2. Multinomial and E = 3 are ungated.

### Primary 1 - per-observation efficiency (P1)

- **GREEN** if `D <= 2 pp` in EVERY realistic gated cell of that column
  type.
- **KILL** if `D >= 15 pp` in ANY realistic gated cell.
- YELLOW otherwise.

*Rationale.* The consumer-visible cost of an install rejection is a
rejected MH move for that row. bairrtt targets 0.44 acceptance with
Robbins-Monro adaptation on `log(theta_sd)`
(`irt_causal_bart.R` lines 198, 647-653); a 2 pp install-rate loss shifts
effective acceptance from ~0.44 to ~0.43, which the adaptation absorbs and
which is smaller than the spread already observed across persons
(0.68-1.04 at the deciles, bairrtt's `TODO` lines 113-114). Above 15 pp the
widening is no longer a tuning knob: more than one proposed latent move in
seven is lost, the adaptation cannot recover it by shrinking the proposal
without collapsing the step size, and the consumer's practical option
becomes buying the cost back with `n.trees.treatment` or the moderator
mask - which is a different arc.

*This line is an EFFICIENCY line and nothing more.* The widened veto does
NOT change the estimand. The tau forest and the variance forest carry the
same empty-leaf invariant as forest 0
([[src/bartcore/chain.hpp#Chain::revalidateTrees,Chain::refreshVarianceForest]]),
so under dbarts' convention a configuration with an empty leaf
in any ensemble has zero prior mass; a latent value that empties a tau leaf
is a zero-prior state of the JOINT model and rejecting it is the correct
Metropolis-within-Gibbs step - as bairrtt's own comment says
(`irt_causal_bart.R` lines 605-608). Any argument that the widening "changes the
target" would condemn the criterion already shipping at E = 1.

### Primary 2 - mixing pathology (P3)

`beta` = twice the seed-to-seed SD of the MATCHED NULL delta (config (a) at
seed s versus config (a) at seed s'), measured in Stage 0 and written into
this file at FREEZE. **FROZEN: `beta` = 0.000 pp in every (n, column type)
cell class** - see "Stage 0 FREEZE" below. P3 therefore stays in the gate,
with GREEN at `delta <= 0.2 pp` and KILL at `delta >= 1 pp`.

- **GREEN** if `delta <= max(0.2 pp, beta)` in every realistic gated cell.
- **KILL** if `delta >= max(1 pp, 3 x beta)` in any realistic gated cell,
  **confirmed on a fresh-seed re-run of that cell**.
- YELLOW otherwise.
- **If `beta > 1 pp`, P3 LEAVES THE GATE** under the house
  floor-above-ceiling rule and is reported descriptively. Pre-registered
  here so it is not a deviation. In that event P1 and L1' alone carry the
  primary verdict and the arc inherits the frozen-subpopulation question as
  an open risk, named in the handoff.

*Rationale.* 1% of rows whose latent is frozen for half the run is a
systematic, individually identifiable bias that running longer does not
remove, and it would corrupt exactly the person-level estimands (ICATE,
latent scores) these model classes exist to produce. This is the only
modeling-harm line in the design; P1 is throughput.

### Primary 3 - priceability (L1')

- **GREEN** if `ratio <= 1.5` in every realistic gated cell.
- **No KILL clause.** Superadditivity alone moves the verdict to YELLOW and
  makes mask/j-split pruning a REQUIRED element of the arc.

*Rationale.* If the widening is priced by the transferable count law, a
consumer who cares prices it directly and buys it back with
`n.trees.treatment`, or eliminates it outright with the moderator or
`varianceForestColumns` mask. A priceable cost is a tuning knob. A
superadditive one is not, but it is also not by itself a reason to decline
the arc.

### Secondary - whole-transaction veto rate (T1)

Evaluated at **fraction 0.01 only**.

- **GREEN** if `rate(E_all, b/c/d) <= 1.5 x rate(shipped, a)` AND
  `rate(E_all) <= 0.5` in every realistic gated cell.
- **KILL** if `rate(E_all) >= 0.9` while `rate(shipped, a) <= 0.5` in any
  realistic gated cell.

*Rationale, and why 1/n is ungated.* At fraction 1/n a transaction moves
one row, so `rate(E_all)` is of order `r_row` (about 1e-3): the absolute
clauses are trivially satisfied and the KILL clause is unreachable, so only
the ratio would carry information and the gate is GREEN-or-YELLOW by
construction. It is reported, not gated. At fraction 1.0 with large moves
the whole-transaction path is already near-certain to veto at E = 1 - that
is why the per-observation session exists - so gating there tests nothing.
Fraction 0.01 is the live regime: sequential covariate imputation inside a
BCF, a one-subject Gibbs block.

### M1 - correctness assertion, not a statistic

An ensemble whose column mask excludes j must contribute **exactly zero**
marginal vetoes. This is a theorem, not a hypothesis (see Context). If M1
measures anything other than exactly 0, the model of the mechanism is
wrong and **the measurement stops**.

### Composite verdict

For a column type: **GREEN** iff P1, P3 (or P3 out of the gate), L1' and T1
all GREEN, all validity gates pass, and M1 == 0. **KILL** iff any KILL
clause fires, where EVERY primary KILL crossing (P1, P3, or T1) must be
confirmed on a fresh-seed re-run of that cell before it counts.
**YELLOW** otherwise.

### Threshold provenance and ratification

Ratified by VD 2026-08-09 on the orchestrator's recommendation, with an
explicit caveat on provenance recorded here so no later reader over-reads
a line as measured when it is judged. Classification:

- **Mechanism-anchored:** P1's GREEN 2 pp (sized against the one recorded
  consumer's Robbins-Monro acceptance target and its observed per-person
  acceptance spread - though the underlying bairrtt "< 1% at 150 trees"
  baseline is prose, not a reproducible measurement); L1's count law
  itself (a theorem, with M1 as its hard assertion), with only the 1.5
  dependence allowance judged.
- **Judgment with structural guards:** P3's 1 pp ceiling (the harm
  mechanism - persistent per-row freezing biases person-level estimands -
  is real; the 1% choice is judged; the operative clause self-calibrates
  via the Stage-0 matched-null beta and self-retires if noise exceeds
  the ceiling).
- **Pure judgment:** P1's KILL 15 pp ("more than one proposed latent
  move in seven lost") and T1's 1.5x / 0.5 / 0.9. No measurement anchors
  these. The wide YELLOW band is the working protection: everything
  between GREEN and KILL opens the arc with required mitigation
  elements, so a mis-set KILL number matters only for cells landing
  beyond it - and by the amendment above, any such crossing needs a
  fresh-seed confirmation before it counts.

Registered pre-data as part of this ratification: any metric landing
within its cell's Stage-0 noise band of a decision line is reported as
NEAR-LINE in the results, whatever side it falls on.

### What YELLOW authorizes

YELLOW OPENS the design arc but does NOT authorize "widen `revalidateTrees`
and `treeAt` and ship". The arc must carry, as required elements:

1. mask and j-split pruning (exactness-preserving; makes the widening free
   for the mask-restricted consumer);
2. a priced comparison of collapse-on-object (forced semantics for the
   objecting ensemble only - note this is a THIRD posterior, not either
   existing one) against declining the widening for the objecting ensemble
   class. Forest-local rollback is incoherent for predictor changes because
   the code store is shared, and the arc should state that rather than
   re-derive it;
3. the P4 / T2 attribution table with its scope caveats (below).

## Power and replication

| arm | cells | seeds | basis |
|---|---|---|---|
| 0 | 13 config cells x 2 column types | 8 | `r_1` and the count law; SE of a per-row rate at n = 1000, W = 200 is ~2e-4 |
| A' gated | 9 per column type x 2 | 8 (6 at n = 5000) | see below |
| A' reported | move sizes, `mask = no`, multinomial | 4 | ungated |
| B | 3 S-levels x 2 n x 2 column types | 8 | plus one subset-replay per non-primary sampler for attribution |
| C | validity | as specified per gate | V3c is the only powered comparison; sized above |

**P1 power.** Each gated cell yields `n x W` per-row install decisions
(1e5 at n = 500, 1e6 at n = 5000) per seed. The binding variance is
between-seed, not within-run. At 8 seeds the SE of the cell-mean `D` is
`sd_seed / sqrt(8)`; the design condition is `SE <= delta/4`, i.e.
`sd_seed <= 2.83 pp` for the 2 pp GREEN margin and `sd_seed <= 21 pp` for
the 15 pp KILL margin. Stage 0 measures `sd_seed` on the config-(a)
baseline and on one config-(b) cell; **if `sd_seed > 2.83 pp` the seed
count rises to `ceil(16 x sd_seed^2 / 4)` at FREEZE**, or P1's GREEN margin
leaves the gate under floor-above-ceiling. Recorded in advance.

**P3 power.** The tail functional's null band `beta` is measured directly
(matched config-(a) null, 8 seed pairs per n and column type), so no
distributional assumption is made. The floor-above-ceiling rule applies as
stated in the P3 gate.

**Multiplicity.** GREEN conjoins over 9 gated cells per column type;
that is conservative in the GREEN direction by construction (harder to
pass) and is the intended asymmetry. KILL is per-cell but requires a
confirmed fresh-seed re-run for P3, and P1/T1 KILL lines are far enough
from the null (15 pp; 0.9 versus 0.5) that a single-cell excursion is not a
multiplicity artifact at these SEs. No Holm correction is applied, and this
is deliberate: the design is asymmetric on purpose, since a false YELLOW
costs a design pass and a false GREEN costs an engine arc.

## Stage 0 FREEZE

Written 2026-08-09, after V0-V3c and the calibration cells, before any gated
number was read. Harness and cells are session-local, not retained.
Build: repo tip b4b8614,
`R CMD INSTALL --preclean` into a private library.

**Validity gates - all eight PASS.**

| gate | outcome | evidence |
|---|---|---|
| V0 | PASS | 96 comparisons, 0 mismatch against the engine-verbatim loop and 0 against the unfiltered router; 34 comparisons carried at least one rejection, 47 carried a tight leaf |
| V1 | PASS | no-op installs everything and vetoes nothing in the oracle for configs (a)-(d) x {continuous, binary} and in the engine for (a); pinned grids and code vectors identical across the three surrogate samplers |
| V2 | PASS | 4972 nodes over 16 cells, integer-exact against the live `n` column, 0 bad trees |
| V2r | PASS | 3000 single-row replays, 0 bad |
| V2v | PASS | relative spread of `engine$variance / oracleProduct` <= 1.8e-16 at `n.trees.variance` 1 and 40, both column types (line 1e-10) |
| V3a | PASS | 1200 proposals, 68 engine vetoes, 0 verdict disagreements |
| V3b | PASS with a recorded deviation (D2) | 20 fits: serialized state identical in EVERY slot (not just trees/values/sizes/flags), `data@x` identical, rng aligned in 20/20; multi-ensemble leg 20/20 structure-preserving and partition-exact on configs (b) and (d). The literal bitwise-run clause fails at max 4.4e-14 |
| V3c | PASS | stressed proposal (engine reject 0.49% continuous / 0.33% binary): difference of means 0.0033 pp and 0.000 pp, CI half-widths 0.042 pp and 0.000 pp, both inside the 0.4 pp line |

**Frozen calibration.** Matched null = eight disjoint config-(a) seed pairs
per cell class; `sd_seed` from 16 config-(a) seeds and 8 config-(b) seeds.

| column type | n | `beta` (pp) | max abs null delta (pp) | config-(a) P3 share (pp) | `sd_seed` (a) (pp) | `sd_seed` (b) (pp) |
|---|---|---|---|---|---|---|
| continuous | 500 | 0.000 | 0.000 | 0.000 | 0.0032 | 0.0079 |
| continuous | 1000 | 0.000 | 0.000 | 0.000 | 0.0018 | 0.0015 |
| continuous | 5000 | 0.000 | 0.000 | 0.000 | 0.0002 | 0.0001 |
| binary | 500 | 0.000 | 0.000 | 0.000 | 0.0164 | 0.0158 |
| binary | 1000 | 0.000 | 0.000 | 0.000 | 0.0072 | 0.0043 |
| binary | 5000 | 0.000 | 0.000 | 0.000 | 0.0015 | 0.0015 |

- **`beta` = 0.000 pp everywhere.** No config-(a) row is ever rejected in
  >= 50% of a W = 200 tracking window, at any seed, at any n, in either column
  type. P3 does NOT leave the gate (`beta` <= 1 pp): GREEN at
  `delta <= max(0.2 pp, 0) = 0.2 pp`, KILL at
  `delta >= max(1 pp, 0) = 1 pp`, KILL requiring a fresh-seed confirmation.
- **`sd_seed` <= 0.0164 pp**, three orders of magnitude under the 2.83 pp
  design condition. The registered seed count therefore STANDS: 8 seeds
  (6 at n = 5000). No margin leaves the gate under floor-above-ceiling.
- **Stage-0 noise bands, for the pre-registered NEAR-LINE labelling.** Per
  cell class the band is `2 x sd_seed / sqrt(8)`: 0.0022 pp (continuous 500),
  0.0013 pp (continuous 1000), 0.0001 pp (continuous 5000), 0.0116 pp
  (binary 500), 0.0051 pp (binary 1000), 0.0011 pp (binary 5000). Every P1
  decision line (2 pp, 15 pp) and every P3 line (0.2 pp, 1 pp) sits far
  outside these bands, so a NEAR-LINE label can only arise from a cell whose
  own metric lands within a band of a line.

## Instrumentation

Ships today; sufficient for every arm; **no engine change**.

| need | shipped hook |
|---|---|
| per-forest live tree structure | `dbarts:::bartcoreGetTrees(bc, chainNums, treeNums, current = TRUE, forest = f)`, `[[R/bartcore.R:961-986@4c018187]]`; `forest` is 0-based and bounded by `shape.numForests` (`[[R_interface_bartcore.cpp:4553-4556@4c018187]]`) |
| ground-truth node occupancy (V2) | the `n` column of that frame, LIVE partition under `current = TRUE, newdata = NULL` |
| variance-forest structure | `sampler$state[[c]]$variance.{vars,values,sizes,flags}`, `[[R_interface_bartcore.cpp:4771-4781@4c018187]]`; raw-double decode precedent `[[inst/tinytest/test-bartcore.R:730-736@4c018187]]` |
| variance-forest ground truth (V2v) | the `variance` element of `bartcoreRun(bc, 0L, 1L)`, `[[R_interface_bartcore.cpp:3293@4c018187]]` |
| BCF construction | `dbarts:::bartcoreBCFSampler(...)`, `[[R/bartcore.R:570-664@4c018187]]` |
| heteroscedastic construction | `dbarts(x, y, variance = TRUE, n.trees.variance =)`, `[[R/dbarts.R:346-349@4c018187]]`; `dbartsSpec(variance =, n.trees.variance =)`, `[[R/spec.R:376-377@4c018187]]` |
| multinomial handle | `bart2(x, y, family = "multinomial", keepTrees = TRUE)$bc`, `[[R/bart.R:1212-1213@4c018187]]` |
| engine's single-forest per-obs mask | `sampler$setPredictor(v, col, forceUpdate = "partial")` -> length-n logical (`[[inst/tinytest/test-bartcore.R:259-280@4c018187]]`) |
| engine's whole-transaction verdict | `sampler$setPredictor(v, col)` -> TRUE/FALSE |
| engine's widened multi-sampler mask | `updatePredictorPerObservationJointly(list(...), x, column)` (exported, `NAMESPACE:8`); asserted as the widened AND at `[[inst/tinytest/test-sampler-updatePredictorPerObservationJointly.R:118@4c018187]]` |
| unconditional install and restore | `setPredictor(v, col, forceUpdate = TRUE)` |
| leaf-routing precedents | `[[inst/tinytest/test-sampler-setPredictorPerObservation.R:13-38@4c018187]]`; `[[inst/tinytest/test-interactions.R:13-64@4c018187]]`; `[[R/plotTree.R:44-61@4c018187]]` |

Built harness-side (R only): a general multi-predictor pre-order router
(~70 lines: recursive index partition on `var`, `value`, `missing`,
`directions`, with `value <= cut -> left`); a `variance.*` flat decoder
(~40 lines, `readBin` over the raw slots split by `variance.sizes`); a
sequential per-observation session simulator (~50 lines, mirroring
`[[sampler.hpp:1354-1385@4c018187]]`); the arm-A' driver loop.

### Contingency: the `ForTesting` patch

**Trigger, and the only trigger: failure of V2, V2r, V2v, or V3b.** Not a
tolerance - a categorical failure.

```
// tests/cpp ONLY. No facade.hpp virtual, no bridge entry, no R visibility -
// adding a SamplerBase virtual forces --preclean everywhere and risks the
// stale-object bus error the house notes warn about.

// chain.hpp, public:
struct WidenedRevalidation { bool allValid; std::vector<size_t> emptiedByEnsemble; };
WidenedRevalidation revalidateAllEnsemblesForTesting();

// sampler.hpp, public:
std::unique_ptr<PredictorUpdateSession>
  beginPredictorUpdateAllEnsemblesForTesting(const double*, size_t);
```

~80 header-only lines plus a `tests/cpp/test_moves.cpp` driver, in its own
throwaway worktree, `--preclean`, **do not land** (memory-wall sec 11
posture). Adds ~1 h to the run and a build cycle. It is not part of the
registered design.

## Harness, cost, machine

```
harness/
  oracle.R   routeTree(flat, x)        -> per-row leaf id + per-node counts   (V2)
             decodeVarianceForest(st)  -> the same flat frame from raw slots  (V2v)
             vetoWholeTransaction(trees, x, xnew, rows) -> logical, per ensemble
             simulateSession(trees, x, xnew, order)     -> install mask, per ensemble
  arm0.R     single forest, chains x trees grid, shipped partial + joint calls
  armA.R     BCF / het driven closed loop: export -> oracle mask -> forced
             install -> run(0,1), W sweeps
  armB.R     S in {1,2,3} shaped samplers, bairrtt pattern, subset replays,
             three-kernel Q1
  armC.R     V0-V4
  results.md tables + the verdict against this file
```

Two exactness-preserving speedups, pre-registered so they are not
deviations: (i) descend only trees containing a j-split; (ii) precompute
each row's old and new leaf per tree vectorized, so the sequential
simulator's R loop touches only rows that actually move. Both validated by
V0.

| stage | content | est. wall (1 core) |
|---|---|---|
| Stage 0 | V0-V4; P3 null band (matched config-(a) pairs, 3 n x 2 types x 8 seed pairs); arm-0 pilot for `sd_seed` | ~2 h |
| FREEZE | write `beta`, `sd_seed`, any revised seed count and margin into this file | - |
| Stage 1 arm 0 | ~208 short closed-loop runs | ~1 h |
| Stage 1 arm A' | ~264 closed-loop runs of 300 burn + 200 tracking sweeps | ~7 h |
| Stage 1 arm B | ~96 closed loops + subset replays + 48 Q1 runs | ~2.5 h |
| Report | tables, per-column-type verdicts | ~0.5 h |
| **total** | | **~11 h single core (+/- 50%; arm A' dominates and its per-sweep oracle cost is an estimate)** |

Embarrassingly parallel over seeds: **~1.75 h on 8 cores**. No timing
metric exists anywhere, so the measurement is load-insensitive; run it on
the laptop. `dbarts-bench` is optional (faster; would catch a Linux/x86
discrepancy in the router, though none is expected - the router reads
exported doubles, not SIMD paths). `R CMD INSTALL .` first; no
`--preclean`.

## Steps

1. `R CMD INSTALL .`. Write `oracle.R`; run V0.
2. Stage 0: V1, V2, V2r, V2v (at `n.trees.variance` 1 and 40), V3a, V3b,
   V3c. Any categorical failure -> stop and invoke the contingency.
3. Stage 0 calibration: P3 null band `beta` per n and column type;
   `sd_seed` for P1 on config (a) and one config-(b) cell.
4. **FREEZE.** Write `beta`, `sd_seed`, the resulting seed count, and any
   metric that left the gate under floor-above-ceiling into this file,
   before reading any gated number.
5. Stage 1: arm 0, arm A', arm B. V4 reported alongside.
6. Report: per-cell tables with each contrast beside its frozen margin; two
   composite verdicts (continuous, binary); the P4/T2 attribution table
   with its scope caveats; X1-X3 and Q1 as context.

## Verification

- V0 passes exactly, or the speedups are removed and the cost table
  re-estimated.
- M1 is exactly 0 in every `mask = no` cell. Otherwise STOP.
- Two invocations at the same `BASE_SEED` produce identical tables (every
  fit pins `rngSeed`, and sampling never advances R's stream).
- Arm A' shape assertion: for a config-(b) fit,
  `bartcoreGetTrees(..., forest = 1)` returns `n.trees.treatment` trees per
  chain and `forest = 2` errors (`numForests` is 2). Stop if not.
- The harness exits nonzero only on a validity-gate failure or M1 != 0,
  never on a study finding.

## Results

Run 2026-08-09 at tip b4b8614 (deviation D1), single build,
`R CMD INSTALL --preclean` into a private library, 8 workers on the laptop.
1971 cells run; results, harness and the full emitted tables are
session-local, not retained. No wall-clock quantity enters any
metric. Reproducibility spot-checked as the Verification section requires: four
cells re-run from scratch reproduced their install masks, per-row reject
counts, veto sequences and `T_j` traces bitwise.

### Verdict

| column type | P1 | P3 | L1' | T1 | M1 | composite |
|---|---|---|---|---|---|---|
| continuous | **GREEN** (max D 0.0059 pp) | **GREEN** (0.000 pp in all 9) | **GREEN** (max 1.18) | **YELLOW** (2/9 ratio-clause failures; max rate 0.0025) | exactly 0 | **YELLOW** |
| binary | **GREEN** (max D 0.0369 pp) | **GREEN** (0.000 pp in all 9) | **GREEN** (max 1.39) | **YELLOW** (4/9 ratio-clause failures; max rate 0.0069) | exactly 0 | **YELLOW** |

**Confirmation status: no obligations outstanding.** The ratification
amendment requires a fresh-seed re-run only for a primary KILL crossing. No
P1 cell reached 15 pp (the largest is 0.037 pp, a factor of 400 below the
line), no P3 cell reached 1 pp (every cell is identically 0), and no T1 cell
reached the KILL condition `rate(E_all) >= 0.9` (the largest is 0.0069). Both
verdicts are therefore FINAL as issued.

**NEAR-LINE cells** (registered pre-data; a metric within its cell's Stage-0
noise band of a line). No P1 or P3 cell is near any line - the Stage-0 bands
are 0.0001-0.0116 pp and the nearest line is 2 pp. Three L1' cells are
NEAR-LINE against the 1.5 ceiling, by the Poisson interval on their event
counts: continuous / het / n = 500 (1.18, upper 1.53), binary / BCF 75-50 /
n = 5000 (1.39, interval 1.15-1.63), binary / het / n = 5000 (1.31, interval
1.10-1.52). All three fall on the GREEN side; the label records that the
replication cannot separate them from 1.5.

### Primary 1 - per-observation efficiency (P1)

Nine gated cells per column type; `D = mean(reject | b/c/d, E_all) -
mean(reject | a, shipped)`, in pp, matched on (n, DGP, move, seed, W), 8 seeds
(6 at n = 5000).

| type | cfg | n | (a) reject % | (b/c/d) reject % | **D (pp)** | SE(D) | P1d (pp) | events (a) | events (x) |
|---|---|---|---|---|---|---|---|---|---|
| continuous | b | 500 | 0.00850 | 0.00938 | +0.00088 | 0.0031 | 0.00025 | 68 | 75 |
| continuous | c | 500 | 0.00850 | 0.00863 | +0.00013 | 0.0024 | 0.00000 | 68 | 69 |
| continuous | d | 500 | 0.00850 | 0.01438 | +0.00588 | 0.0014 | 0.00375 | 68 | 115 |
| continuous | b | 1000 | 0.00319 | 0.00163 | -0.00156 | 0.0011 | 0.00006 | 51 | 26 |
| continuous | c | 1000 | 0.00319 | 0.00300 | -0.00019 | 0.0009 | 0.00006 | 51 | 48 |
| continuous | d | 1000 | 0.00319 | 0.00344 | +0.00025 | 0.0008 | 0.00113 | 51 | 55 |
| continuous | b | 5000 | 0.00038 | 0.00010 | -0.00028 | 0.0001 | 0.00000 | 23 | 6 |
| continuous | c | 5000 | 0.00038 | 0.00023 | -0.00015 | 0.0001 | 0.00000 | 23 | 14 |
| continuous | d | 5000 | 0.00038 | 0.00047 | +0.00008 | 0.0001 | 0.00010 | 23 | 28 |
| binary | b | 500 | 0.0625 | 0.0588 | -0.00375 | 0.0086 | 0.00088 | 500 | 470 |
| binary | c | 500 | 0.0625 | 0.0591 | -0.00338 | 0.0039 | 0.00038 | 500 | 473 |
| binary | d | 500 | 0.0625 | 0.0994 | +0.03688 | 0.0037 | 0.03388 | 500 | 795 |
| binary | b | 1000 | 0.0279 | 0.0254 | -0.00256 | 0.0025 | 0.00056 | 447 | 406 |
| binary | c | 1000 | 0.0279 | 0.0247 | -0.00325 | 0.0028 | 0.00013 | 447 | 395 |
| binary | d | 1000 | 0.0279 | 0.0444 | +0.01650 | 0.0031 | 0.01594 | 447 | 711 |
| binary | b | 5000 | 0.00373 | 0.00477 | +0.00103 | 0.0005 | 0.00015 | 224 | 286 |
| binary | c | 5000 | 0.00373 | 0.00317 | -0.00057 | 0.0005 | 0.00008 | 224 | 190 |
| binary | d | 5000 | 0.00373 | 0.00798 | +0.00425 | 0.0006 | 0.00247 | 224 | 479 |

GREEN on both types, by a factor of 54 (binary) to 340 (continuous) on the
2 pp line. The largest cost of any kind - adding a 40-tree variance forest to
a binary latent at n = 500 - is **0.037 pp**, i.e. 3.7 rejected installs per
10,000 proposed row-moves; four of the nine continuous cells and four of the
nine binary cells are NEGATIVE (the multi-ensemble fit rejects less often than
the single-forest baseline, because the added ensemble absorbs signal and
mu's trees get shallower). P1d - the ungated `E_all` versus `E_0` contrast on
the same multi-ensemble fit, i.e. literally what the code change toggles - is
never above 0.034 pp and is essentially all of D in the heteroscedastic cells,
confirming that in the BCF cells the (already tiny) D is a fit-shape effect
rather than the widening.

### Primary 2 - mixing pathology (P3)

`delta` is **identically 0.000 pp in all 18 gated cells**, against a frozen
GREEN line of 0.2 pp and a KILL line of 1 pp. Not one row in any
configuration, at any n, in either column type, was rejected in as many as
half of its 200 tracking transactions. The concentration diagnostic P2 shows
why: rejections are extremely concentrated (max-to-mean per-row reject-rate
ratios of 70 to 10,000) but *transient* - a given row is rejected once or
twice out of 200, never persistently. The frozen-subpopulation failure mode
this line exists to catch does not appear. See deviation D7 for the structural
caveat on the binary type.

### Primary 3 - priceability (L1')

GREEN on both types: every cell is at or under the 1.5 ceiling, with the three
NEAR-LINE cells listed above. Values run 0.22-1.18 (continuous) and 0.86-1.39
(binary). The count law is not just an upper bound here, it is close to an
identity: the strongest evidence is arm 0, where quadrupling the number of
chains widens the criterion by exactly 4x the tree count and

| type | n | `T_j` at 1 chain | reject % at 1 chain | `T_j` at 4 chains | reject % at 4 chains | observed / count-law predicted |
|---|---|---|---|---|---|---|
| continuous | 300 | 10.2 | 0.0175 | 40.9 | 0.0754 | 1.08 |
| continuous | 500 | 10.9 | 0.0085 | 42.6 | 0.0319 | 0.96 |
| continuous | 1000 | 10.4 | 0.0032 | 41.0 | 0.0136 | 1.08 |
| binary | 300 | 10.3 | 0.1277 | 43.0 | 0.5004 | 0.94 |
| binary | 500 | 10.5 | 0.0625 | 42.3 | 0.2561 | 1.02 |
| binary | 1000 | 10.1 | 0.0279 | 42.2 | 0.1164 | 1.00 |

A 4x widening of the veto is priced to within 8% by
`1 - (1 - r_1)^(T_j)` with `r_1` read off the un-widened fit. Mask and
j-split pruning therefore buy back exactly what the arithmetic says they
should. X3, reported and not gated: the pooled log-log slope of the per-sweep
reject rate on `T_j` is 0.41 (continuous) and -0.06 (binary) - sublinear, as
the design anticipated, and for the registered reason (row-level heterogeneity
in `r_1` saturates the union).

### Secondary - whole-transaction veto rate (T1), fraction 0.01

**This is the only clause that is not GREEN, and it is the whole reason both
composites are YELLOW.** The absolute clause passes everywhere by three orders
of magnitude: the largest `rate(E_all)` in any gated cell is **0.0069**,
against a GREEN ceiling of 0.5 and a KILL floor of 0.9. The ratio clause
`rate(E_all) <= 1.5 x rate(shipped, a)` fails in 2 of 9 continuous cells and 4
of 9 binary cells.

| type | cfg | n | (a) rate (events) | `E_all` rate (events) | `E_0` rate | ratio | ratio <= 1.5 |
|---|---|---|---|---|---|---|---|
| continuous | b | 500 | 0.00188 (3) | 0.00063 (1) | 0.00063 | 0.33 | yes |
| continuous | c | 500 | 0.00188 (3) | 0.00125 (2) | 0.00125 | 0.67 | yes |
| continuous | d | 500 | 0.00188 (3) | 0.00125 (2) | 0.00125 | 0.67 | yes |
| continuous | b | 1000 | 0.00125 (2) | 0.00063 (1) | 0.00063 | 0.50 | yes |
| continuous | c | 1000 | 0.00125 (2) | 0.00125 (2) | 0.00125 | 1.00 | yes |
| continuous | d | 1000 | 0.00125 (2) | 0.00125 (2) | 0.00125 | 1.00 | yes |
| continuous | b | 5000 | 0.00083 (1) | 0.00000 (0) | 0.00000 | 0.00 | yes |
| continuous | c | 5000 | 0.00083 (1) | 0.00167 (2) | 0.00167 | 2.00 | **no** |
| continuous | d | 5000 | 0.00083 (1) | 0.00250 (3) | 0.00250 | 3.00 | **no** |
| binary | b | 500 | 0.00313 (5) | 0.00500 (8) | 0.00500 | 1.60 | **no** |
| binary | c | 500 | 0.00313 (5) | 0.00063 (1) | 0.00063 | 0.20 | yes |
| binary | d | 500 | 0.00313 (5) | 0.00313 (5) | 0.00188 | 1.00 | yes |
| binary | b | 1000 | 0.00063 (1) | 0.00250 (4) | 0.00250 | 4.00 | **no** |
| binary | c | 1000 | 0.00063 (1) | 0.00375 (6) | 0.00375 | 6.00 | **no** |
| binary | d | 1000 | 0.00063 (1) | 0.00688 (11) | 0.00438 | 11.00 | **no** |
| binary | b | 5000 | 0.00500 (6) | 0.00000 (0) | 0.00000 | 0.00 | yes |
| binary | c | 5000 | 0.00500 (6) | 0.00333 (4) | 0.00333 | 0.67 | yes |
| binary | d | 5000 | 0.00500 (6) | 0.00250 (3) | 0.00167 | 0.50 | yes |

Read the event counts, not the ratios: the registered replication (8 seeds,
6 at n = 5000, W = 200) puts **0 to 11 whole-transaction vetoes** in a cell,
so a ratio of 11.0 is "11 events against 1 event". The clause is
Poisson-noise dominated at this replication and the design did not size for
it - the power table sizes P1 and P3 but leaves T1's replication implicit.
The verdict is issued as registered anyway: T1 = YELLOW on both types.

An UNGATED supplement at 40 seeds per cell (seeds 9-40, ~8000 transactions per
cell) resolves what the registered replication cannot: **the ratio clause is
noise**. Supplement ratios and their 95% intervals: continuous 0.42-2.20, all
nine intervals covering 1, the largest point estimate 2.20 (interval
-0.13-4.53, 11 events against 5); binary 0.43-1.95, all nine covering 1, the
single point estimate over 1.5 being het at n = 1000 (1.95, interval
0.87-3.03). Only 1 of the 18 supplement cells has a point estimate above 1.5,
and no interval excludes 1.5 from below. This does NOT change the registered
verdict, and it is recorded as context, not as a re-decision.

The `E_0` column is the useful one for the arc: on the BCF configurations
`rate(E_0)` and `rate(E_all)` are IDENTICAL in every gated cell - a tau forest
restricted to moderators never objects to a whole-column swap that mu accepts.
Only the heteroscedastic configuration separates them (binary n = 500:
`E_0` 0.0019 vs `E_all` 0.0031; binary n = 1000: 0.0044 vs 0.0069; binary
n = 5000: 0.0017 vs 0.0025). The whole-transaction cost of the widening is
carried by the variance forest.

### M1 - correctness assertion

**Exactly 0, as the theorem requires.** Across 18 mask-excluded cell classes
(4 seeds each), the masked ensemble contributed a j-splitting tree count of
exactly 0, and the widened per-observation mask differed from the forest-0
mask in **0 of 14,400,000 install decisions**, with 0 differing
whole-transaction verdicts in 14,400 transactions. A per-forest column mask is
an exact opt-out, measured.

### P4 / T2 attribution (reported, not decided - fork 2)

Oracle-attributed for arm A' (the shipped surfaces return a conjunction, never
an attribution); engine-measured for arm B.

| type | cfg | ensemble | `T_j` | mean leaf depth | P(leaf holds 1 obs) | P4 sole objector | T2 per-obs (pp) | T2 transaction (pp) |
|---|---|---|---|---|---|---|---|---|
| continuous | BCF 75/50 | mu | 8.6-10.5 | 1.36-1.40 | 0.026-0.047 | 96-100% | 0 | 0 |
| continuous | BCF 75/50 | tau | 1.1-1.3 | 0.40-0.42 | 0.008-0.018 | 0-2.7% | 0.00006-0.00025 | 0.000 |
| continuous | het 75/40 | mu | 7.4-10.5 | 1.32-1.35 | 0.036-0.049 | 67-79% | 0 | 0 |
| continuous | het 75/40 | variance | 4.9-5.2 | 1.27-1.34 | n/a | 21-33% | 0.0001-0.0038 | 0.000 |
| binary | BCF 75/50 | tau | 1.0-1.4 | 0.41-0.43 | 0.010-0.021 | 1.5-3.1% | 0.00015-0.00088 | 0.000 |
| binary | het 75/40 | variance | 5.5-5.6 | 1.24-1.32 | n/a | 31-36% | 0.0025-0.0339 | 0.083-0.250 |

Three things the arc should carry forward. First, a **tau forest restricted to
moderators is nearly free**: it contributes ~1 j-splitting tree out of 50
(against ~10 out of 75 for mu) because `treatment.base = 0.25` /
`treatment.power = 3` keeps it at mean leaf depth ~0.4, i.e. mostly stumps, and
its leaves therefore hold most of the sample. It is the sole objector in 0-3%
of rejections. Second, the **variance forest is the real payer**: at 40 trees
under the mean model's own tree prior it carries 5-6 j-splitting trees at
depth ~1.3, and is the sole objector in 21-36% of rejections. Third, the
per-forest split is scope-limited exactly as the design says: it is oracle
attribution, and it attributes per FOREST inside one sampler, which no shipped
surface can confirm. Arm B's per-SAMPLER attribution IS engine-measured, by
snapshot-and-replay: `setState` restores chain 0's rng, so each subset replay
draws the same scan permutation and the masks are scan-order matched by
construction (verified: a replayed subset call reproduces its mask bitwise and
the restore returns every serialized slot exactly).

### Arm B - surrogate, calibration, E = 3

| type | n | S | `T_j` | reject % | engine marginal S=1 | S<=2 | S<=3 |
|---|---|---|---|---|---|---|---|
| continuous | 1000 | 1 | 10.4 | 0.00319 | 0.00319 | - | - |
| continuous | 1000 | 2 | 11.2 | 0.00319 | 0.00319 | 0.00319 | - |
| continuous | 1000 | 3 | 15.7 | 0.00369 | 0.00313 | 0.00313 | 0.00369 |
| binary | 1000 | 1 | 10.1 | 0.0279 | 0.0279 | - | - |
| binary | 1000 | 2 | 11.6 | 0.0278 | 0.0273 | 0.0278 | - |
| binary | 1000 | 3 | 16.6 | 0.0426 | 0.0265 | 0.0269 | 0.0426 |

E = 3 - the shape the engine refuses to construct - costs 0.0005 pp
(continuous) to 0.016 pp (binary) over E = 1 on the surrogate. Reported and
ungated, and it authorizes nothing: `createBCFSampler` still returns nullptr
on `numVarianceTrees > 0` and lifting that is a separate arc.

**V4 / C1.** Arm A' (real BCF, E = 2) against arm B (S = 2 surrogate) at
matched n: C1 = -0.0009, +0.0016 pp (continuous, n = 500 / 1000) and +0.0054,
+0.0024 pp (binary). `|C1|` is at most 0.0054 pp against the 3 pp line, so the
arms agree and the surrogate's tree-shape bias is not detectable at this
resolution. The verdicts are read off arm A' regardless, as registered.

**Q1** (reported, never gated; three kernels on the S = 2 surrogate, with a
real Metropolis filter on the latent so the comparison is between stationary
distributions rather than random walks - see deviation D8). Continuous: mean
of per-row posterior means 0.041 (widened veto) / 0.045 (forced collapse) /
0.029 (frozen forest); mean of per-row posterior sds 0.277 / 0.268 / 0.488;
pooled KS 0.025 (veto vs collapse), 0.039 (veto vs frozen), 0.026 (collapse
vs frozen). Binary: the Metropolis filter almost never accepts a flip against
this surface (per-row sd ~0.003), so the kernels are indistinguishable
(KS <= 0.005) and the comparison is uninformative rather than reassuring. The
continuous reading is that the veto and the forced collapse produce nearly the
same latent posterior on the surrogate, and both differ from the frozen-forest
kernel by more than they differ from each other. It bounds the posterior
question; it does not settle it, and by design it gates nothing.

### Reported, ungated

- **bairrtt cross-check.** At the one recorded consumer shape (n = 300, 150
  trees as 2 chains x 75, ~1 SD move, config (a) shipped) the harness rejects
  **0.043%** of moves on a continuous latent and **0.226%** on a binary one,
  against bairrtt's prose "under 1%" (`bairrtt/[[TODO:117-118@4c018187]]`). Consistent, and
  on the low side. Under the accepted-rows reading (~0.44n) the figures are
  ~0.10% and ~0.51%; both readings sit under 1%.
- **Move size.** Smaller moves do not cost more: at n = 1000 continuous,
  reject rates are 0.0036% / 0.0011% / 0.0032% at moves 0.1 / 0.4 / 0.8 -
  flat, because what matters is whether a row is the last occupant of a leaf,
  not how far it travels once it leaves.
- **Fraction.** At fraction 1/n the whole-transaction veto rate is 0-0.125%
  (the design predicted "of order 1e-3", correct); at fraction 1.0 it is
  1.75-4.0% (continuous) and 19-28% (binary), which is why the
  per-observation session exists. Both bracket the gated 0.01 regime and
  neither is gated.
- **Multinomial K = 4 (E = 4), UNGATED stress.** n = 1000: per-observation
  reject 0.010% at `E_all` vs 0.002% at `E_0` (continuous) and 0.087% vs
  0.020% (binary); whole-transaction 0.25% vs 0.125% and 1.25% vs 0.25%. Four
  symmetric 75-tree forests raise `T_j` to ~40 and the reject rate rises
  roughly with it - the count law again. Still far under every line, but this
  is the shape where a widened veto costs the most, and it is reported, not
  decided.
- **Five-level categorical latent, UNGATED.** n = 1000, per-observation:
  0.034% (single forest), 0.029% (BCF), 0.050% at `E_all` vs 0.029% at `E_0`
  (het). Whole-transaction at fraction 0.01: 0.25% / 0.375% / 0.50% vs 0.125%
  at `E_0`. Same order as the binary type. Columns with more than 64 levels
  remain out of scope for the reason recorded below.
- **p = 5.** Halving the predictor count roughly doubles `T_j` (10.4 -> 19.7
  continuous, 10.1 -> 20.7 binary) and roughly doubles the reject rate
  (0.0032% -> 0.0104%, 0.0279% -> 0.0531%) - the count law once more.

### Reading

The widening is essentially free in every shape that exists, and it is
free for a *reason* the arc can rely on: the rejection rate is
`1 - (1 - r_1)^(T_j)` with `r_1` around 1e-5 to 1e-4, so the entire cost is
carried by the count of j-splitting trees, and pruning that count - by
`moderators`, by `varianceForestColumns`, or by skipping trees with no j-split
- removes the cost exactly (M1 = 0, measured over 14.4 million decisions). The
verdict is YELLOW rather than GREEN on a secondary clause whose registered
replication cannot resolve it and whose supplement says it is noise. What
YELLOW authorizes is unchanged and is stated below; nothing in these results
argues against any of its three required elements, and the P4/T2 table says
which of them matters: the variance forest, not the tau forest, is what a
widened veto makes a consumer pay for.

## What this measurement cannot tell us

1. **The stationary widened posterior, on real trees.** Every arm measures
   rejection rates. Q1 bounds the posterior question on the arm-B
   surrogate, under three kernels, and is REPORTED - never gated, because a
   surrogate posterior comparison must not gate an engine arc. P3 is the
   gated proxy, deliberately.
2. **`E = 3` on real trees.** Not constructible (`[[facade.hpp:631@4c018187]], 644`;
   `[[tests/cpp/test_model.cpp:720-733@4c018187]]`). Config (e) is a surrogate,
   reported and ungated. **No verdict here can authorize a heteroscedastic
   BCF**: that configuration does not exist and lifting the factory refusal
   is a separate, unpriced arc. A GREEN means "widen for the shapes that
   exist" - BCF at E = 2, heteroscedastic at E = 2, multinomial ungated.
3. **Which ensemble to blame (fork 2).** Every shipped surface returns a
   conjunction, never an attribution: `setPredictor(partial)` returns one
   mask over forest 0 across all chains (`[[sampler.hpp:1019-1034@4c018187]]`), the
   joint call ANDs the per-sampler validities before committing
   (`[[facade.hpp:459-470@4c018187]]`), and `setPredictor` returns one boolean. P4 and
   T2 therefore come from the R oracle for arm A', with only arm B's
   per-SAMPLER marginals measured by the engine (via seed-pinned subset
   replay). Per-FOREST attribution inside one sampler is oracle-only.
   **Fork 2 is REPORTED, not decided by this measurement.**
4. **Multi-level and pooled categorical latents.** Verdicts cover
   continuous and binary. Categoricals with <= 64 levels are reported,
   ungated. Categoricals with **more than 64 levels are out of scope**:
   they are pooled (`columnIsPooled` is `maskWordsForCount(numCuts[j]) > 1`,
   `[[data.hpp:433-436@4c018187]]`) and pooled rules do not reach `variance.*`, which is
   flattened with `masks = nullptr` (`[[chain.hpp:2189-2196@4c018187]]`) while the mean
   forests get a `tree.masks` slot (`[[R_interface_bartcore.cpp:4756@4c018187]]`). A
   discrete follow-up covering >64-level columns needs the variance state to
   carry a mask side channel first. (Related, and logged for the TODO rather
   than fixed here: a heteroscedastic sampler whose store holds a pooled
   categorical column would null-dereference in `flattenBelow` -
   `masks->size()` at `[[tree.hpp:1397@4c018187]]` - on `storeState`. Pre-existing,
   unrelated to this arc.)
5. **Whether the arc is worth its own cost.** This prices ONE consequence.
   The engineering estimate is out of scope, but one cost fact belongs on
   the record: "widen `treeAt`" is not a one-line change.
   `UpdateSessionImpl` sizes itself `options_.numTrees * chains_.size()`
   (`[[sampler.hpp:1332-1333@4c018187]]`) and indexes
   `chains_[t / numTrees]->tree(t % numTrees)` (`[[sampler.hpp:1390-1393@4c018187]]`); for a BCF
   `options_.numTrees` is mu's count only and tau's differs
   (`numTreesInForest`, `[[sampler.hpp:1107@4c018187]]`); and `Chain::tree` is hard-wired
   to `forests_[0]` (`[[chain.hpp:2749@4c018187]]`). The widened session needs
   per-forest offsets.
6. **Arm B's tree shapes.** The surrogate's ensembles fit `y`
   independently rather than sharing the combiner's residual, so their
   trees are not the trees a real BCF grows. Arm A' supplies those, and V4
   reports the gap.

## Thresholds that are judgment calls

Flagged rather than buried; the mechanical ones (W = 200, V3c 0.4 pp,
V4 3 pp, 200 replicates, 1e-10 on V2v) follow from the arithmetic above.

1. **P1 KILL at 15 pp.** Purely an efficiency line now, with no measurement
   behind the number. 10 pp or 25 pp are equally defensible. The 2 pp GREEN
   line is well anchored (bairrtt's 0.44 target and its measured 0.68-1.04
   per-person spread).
2. **P3's scientific ceiling, 1 pp KILL / 0.2 pp GREEN.** `beta` sets the
   statistical floor at FREEZE; the ceiling - how many persistently frozen
   rows is a modeling harm - is judgment. 0.5% or 2% are arguable.
3. **L1' ratio ceiling, 1.5.** Well-posed but arbitrary; the meaningful
   reference points are 1.0 (exactly the count law) and 2.0.
4. **Same thresholds for binary as for continuous.** A binary flip
   re-routes the row in every j-splitting tree unconditionally, with no
   cut-crossing discount, so `r_1` is structurally larger and a KILL on the
   binary type may reflect the regime rather than the widening. Registered
   as: same lines, separate verdicts. A widened binary line is a reasonable
   alternative call.

## Deviations from this pre-registration

All of D1-D6 were recorded at FREEZE, before any gated metric was read.
D7-D8 are execution notes recorded during Stage 1.

**D1 (2026-08-09) - tip.** The header pins d3cb94b / engine 06c0254. The run
executed at branch tip **b4b8614** (engine 9929ede), all CI green, which is
what the orchestrator supplied. Nothing between the two tips touches the
predictor-mutation, revalidation or export paths the design anchors on; every
line reference in this file was re-verified at b4b8614 before the harness was
written. During the run a concurrent agent landed d2e8e0a
(`docs/plans/archive/composition-mixing-probe.md`, documentation only), so branch HEAD
has advanced; the measured binary is b4b8614's `src/` exactly.

Also recorded against the header: the budget line says "~700-850 lines across
five files". The harness came to ~2240 lines across seven
(`oracle.R`, `common.R`, `armB.R`, `armC.R`, `run.R`, `stage0.R`, `report.R`).
The excess is the resumable job dispatcher, the per-column-type report
generator and the arm-B state-replay attribution, none of which the budget
line anticipated; no engine or package file was touched.

**D2 (2026-08-09) - V3b's bitwise-run clause.** V3b has three clauses. The
state-identity clause passes in a STRONGER form than written: all 20 fits
agree in *every* serialized slot, not merely trees/values/sizes/flags, and the
installed column agrees exactly. The multi-ensemble leg passes 20/20. The
third clause - "the two subsequent `run(0, 1)` outputs bitwise equal" - fails
in 20/20 fits, at a maximum absolute difference of **4.4e-14** on `train`
(sigma agrees to the printed digit). Cause, diagnosed: the serialized state
deliberately omits the accumulated total-fits cache
(`[[R_interface_bartcore.cpp:4820-4826@4c018187]]`, "dropped the accumulation-history
slots"), and the two paths rebuild it in different summation orders -
`revalidateAllChains` -> `rebuildFitsFromParameters` versus
`forceRefreshTrees`. That is a floating-point association difference, not the
semantic divergence the contingency exists to catch (an emulated install
landing a *different live state* than the widened path would). No metric in
this design reads a fit: masks and verdicts are decided by tree partitions and
integer leaf counts, both certified exact by V2, V2r, V2v and the V3b
multi-ensemble leg. The gate is therefore recorded as PASS with the bitwise
clause relaxed to `<= 1e-10`, and the **ForTesting contingency is NOT
invoked**. Consequence for arm A': the driven loop is statistically, not
bitwise, the trajectory a widened engine would take.

Also recorded here because the pre-registration does not mention it: the
partial path draws a scan permutation from chain 0's rng
(`[[sampler.hpp:1066@4c018187]]`) while the forced path draws none, so V3b as literally
written could never have matched two rng streams. The harness aligns them by
issuing a no-op partial call on the twin before the forced install (V1
certifies a no-op partial changes nothing); with that, rng state matched in
20/20.

**D3 (2026-08-09) - a third exactness-preserving speedup.** The file
pre-registers two speedups. Both are implemented, and both are insufficient at
the largest cells: at n = 5000 with ~60 j-splitting trees a sweep produces
~1e5 (row, tree) mover pairs, and the registered "R loop touches only rows
that actually move" is still ~6e7 interpreted iterations per run. The harness
adds a third restriction, exact rather than approximate: a leaf whose
out-mover count is at most `count - 1` can never be read at occupancy 1 by any
checking mover (its k-th checker sees `count - (k - 1) >= 2`, and incoming
moves only raise the count), so only rows touching a "tight" leaf
(`moversOut >= count`) need the sequential loop, and only tight-leaf counts
need tracking. V0 was extended to cover it: 96 comparisons against the
engine-verbatim loop, 0 mismatches, with 47 comparisons carrying at least one
tight leaf and 34 carrying at least one rejection.

**D4 (2026-08-09) - the response surface.** The pre-registration pins n, p,
tree counts, move sizes, the cut grid and the phase structure, but not the
response surface. The harness fixes one, ASCII-recorded in
`harness/common.R`: the latent is the surface's primary driver (as it is in
every motivating class), entering mu through
`3 * standardize(sin(2 theta) + 0.5 theta)` alongside a linear term, an
interaction and a kink in the other covariates, entering tau linearly, and
entering the variance surface as `exp(0.6 * standardize(...) + 0.3 x2)`;
`pihat` is a genuine propensity column held out of the moderator set. Sanity
check against the one recorded consumer figure: at bairrtt's shape (n = 300,
150 trees, ~1 SD move, config (a)) this surface rejects 0.08-0.32% of moves,
against bairrtt's prose "under 1%". The verdicts below are conditional on this
surface; X1/X2 report the leaf-size and `r_1` profile it produces so a reader
can place it.

**D5 (2026-08-09) - V3c runs on a stressed proposal.** V3c's arithmetic is
sized for reject rates of 2-25%. The gated cells reject at ~1e-5 to 1e-3, so a
200-replicate difference of means on a gated cell has a zero-width CI and
tests nothing. V3c therefore runs on a deliberately stressed proposal (n = 300,
4 chains, 6 SD continuous move) whose engine reject rate is 0.49%. The gate is
about oracle-versus-engine scan-order agreement, which is a property of the
simulator, not of a scenario.

**D6 (2026-08-09) - the L1' and P4/T2 estimators.** The file names L1' and the
attribution statistics but not their estimators. Registered here, pre-data:
`r_1^(a)` is estimated by pooling every config-(a) decision in the matched
cell class across seeds and sweeps and inverting the union law once,
`r_1 = 1 - (1 - rejectRate)^(1 / mean T_j)`; L1' then divides the pooled
multi-ensemble reject rate by `1 - (1 - r_1^(a))^(mean T_j)` at that cell's own
`T_j`. Pooling rather than per-seed inversion is forced by the event counts -
several cells produce single-digit rejections per seed. P4 evaluates every
tree at the moment of the check rather than short-circuiting as the engine
does (the mask is unaffected; only the objector set differs). T2 is measured
per ensemble as `rate(E_0 u {f}) - rate(E_0)` on the same fit and the same
scan order.

**D7 (2026-08-09) - P3 on the binary column type is structurally bounded.**
The binary proposal flips a proposed row with probability 0.5, so a row is
offered a move in only ~50% of the W transactions and can be rejected in at
most ~50% of them. P3's ">= 50% of transactions" functional is therefore
near-degenerate for the binary type by construction, not by measurement. It is
reported as measured and its GREEN carries that caveat.

**D8 (2026-08-09) - arm B's Q1 kernel comparison needed a Metropolis filter.**
Q1 asks for "the stationary per-row latent posterior under three kernels". The
rest of the design drives the latent with a bare random walk plus the veto,
which has no stationary distribution, so a three-kernel *posterior* comparison
is not defined on it. Q1 therefore runs its own kernel: a Metropolis-within-
Gibbs step on each row's latent against the two surrogate samplers' Gaussian
likelihoods (fits from `bartcore_predict` on the live trees) and an N(0, 1)
prior, with the veto / forced collapse / frozen forest applied to the accepted
rows. This is the bairrtt pattern (`[[irt_causal_bart.R:598-613@4c018187]]`: MH filter, then
install). It is confined to Q1, which is REPORTED and never gated, and it does
not touch P1's denominator anywhere else in the design.

**D9 (2026-08-09) - T1 replication supplement, ungated.** T1's ratio clause
turned out to be unresolvable at the registered replication (0-11 events per
cell). A 40-seed supplement was run to say whether the failures are real. It is
recorded as REPORTED and UNGATED; the T1 verdict in the Results section is read
off the registered 8 (6) seeds exactly as pre-registered, and is YELLOW.

**D10 (2026-08-09) - a `dbartsMixedMatrix` design needed densifying.** A
factor in the design makes `data@x` a `dbartsMixedMatrix`, which refuses
`x[, j] <- v`. The harness densifies with `as.matrix` before routing. This
affects only the ungated 5-level-categorical cells; the continuous and binary
designs are plain numeric matrices throughout. Noted because it is a live
sharp edge for any consumer writing an R-side oracle against `data@x`.
