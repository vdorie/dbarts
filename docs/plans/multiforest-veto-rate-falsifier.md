# multiforest-veto-rate-falsifier

status: PRE-REGISTERED, NOT RUN. Ratified design; supersedes
  `.claude/d1-veto-rate-falsifier/memo.md` and its critique wherever they
  differ. Adjudication record: `.claude/d1-veto-rate-falsifier/synthesis.md`.
agent: opus (harness + analysis; no engine change)
rng: neutral - measurement only. A GREEN verdict AUTHORIZES the
  `multiforest-predictor-mutation` design arc; it does not perform it.
budget: harness only, under `.claude/d1-veto-rate-falsifier/harness/`;
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

- Widened refusal today: `refuseMultiForestTransactionalUpdate`
  (`src/R_interface_bartcore.cpp:2087-2095`) and the separate
  `refuseVarianceForestPredictorMutation` (`:2069-2075`). Both fire only
  when `!forcedUpdate` (`:2090`); the FORCE paths stay open and refresh
  every forest.
- Forest-0-only revalidation: `Chain::revalidateTrees`
  (`src/bartcore/chain.hpp:1533-1543`, `Forest& forest = forests_[0]`).
  Two-phase all-chain transaction: `Sampler::revalidateAllChains`
  (`src/bartcore/sampler.hpp:1138-1151`).
- Per-observation session: `sampler.hpp:1326-1403`; predicate at
  `:1354-1373`; tree enumeration at `:1390-1393`
  (`chains_[t / numTrees]->tree(t % numTrees)`, and `Chain::tree` is
  `forests_[0].trees[t]`, `chain.hpp:2749`).
- Joint multi-sampler sweep: `facade.hpp:446-476`; scan permutation from
  `samplers[0]->rng()` at `:458`.
- The unit of study is **not** `numForests`. A heteroscedastic sampler
  keeps `numForests == 1`; the variance forest lives outside `forests_`
  and is signalled by `SamplerShape::hasVarianceForest`. Use
  `E = numForests + (hasVarianceForest ? 1 : 0)` - the number of ensembles
  routed off the shared column store. This confirms
  `docs/design/model-space-survey.md:73-76, 441-450`; it does not correct
  it.
- **The variance forest carries the same empty-leaf invariant as forest
  0**: `chain.hpp:2666-2668` rejects a rebuild whose variance trees have an
  unoccupied bottom, and `refreshVarianceForest` asserts it -
  "the empty-leaf veto admits no unoccupied bottom into live state"
  (`chain.hpp:3286-3291`). So a widened veto adds no constraint the joint
  model does not already impose; it only makes the sampler enforce one it
  already implies. This is the whole basis of the re-grounded kill logic.
- `E = 3` (heteroscedastic BCF) is **not constructible**:
  `createBCFSampler` returns nullptr on `numVarianceTrees > 0`
  (`facade.hpp:631`), multinomial at `:644`, both asserted at
  `tests/cpp/test_model.cpp:720-733`.
- **A tree with no split on the mutated column cannot veto - exactly.**
  In `UpdateSessionImpl` the override never fires
  (`tree.hpp:991-1004`), so `newLeaf == oldLeaf` and `valid` is untouched
  (`sampler.hpp:1354-1373`); in `revalidateTrees` the repartition
  reproduces the same partition. Per-forest column masks are therefore an
  exact opt-out (tau moderators `combiner.hpp:219,237` + `chain.hpp:649,3952`;
  variance `chain.hpp:3050-3060`; generic `chain.hpp:593-600`).
- Consumer baseline: bairrtt runs the widened predicate at two samplers
  today - `updatePredictorPerObservationJointly(list(response_model,
  assignment_model), theta[, k], theta_names[k])`
  (`/Users/vdorie/Repositories/bairrtt/R/irt_causal_bart.R:609-613`), 75
  trees each (`:204`), one dbarts chain each (`:432`), so 150 trees per
  row already. Burn-in uses `setPredictor(..., forceUpdate = TRUE)`
  (`:619-628`), so no veto is exercised during burn. Rejected rows are
  reverted as rejected MH moves (`:614`); no counter is kept. Identical cut
  grids are pinned on the shared column (`:452, 464`). The only recorded
  rate is prose - `bairrtt/TODO:117-118`, "The empty-leaf install rejects
  under 1% of moves", at 300 persons x 30 items (`:110`) - not reproducible
  from any stored artifact. Denominator: the MH accept filter (`:598-599`)
  precedes the install (`:609`) and `accept_rate` at `:615` is over all n,
  so read "moves" as the ~0.44n accepted rows unless stated; both readings
  are reported.
- Motivating classes and their column types:
  `docs/design/model-space-survey.md:570-635`. BPCF principal strata enter
  tau_Y as splitting covariates and are DISCRETE for a binary intermediate
  (`:184-193`, full-text verified there with
  `https://arxiv.org/html/2403.13256v1` and `github.com/lit777/BPCF`).
  Sequential covariate imputation and treatSens-style latent confounders
  admit binary columns too. Hence the per-column-type verdict below.
- Declined alternative, for the KILL branch:
  `docs/design/empty-leaf-veto.md:115` prices occupancy-aware proposals at
  a "250-400 line, posterior-changing rewrite".
- Precedent for falsifier-before-arc: `docs/design/memory-wall-frontier.md`
  sec 11. Form of this file: `docs/plans/grow-from-root-default-study.md`.

## Constraints

- **Read-only on the package.** No engine edit, no default change, no
  baseline re-record. `R CMD INSTALL .` first (tinytest runs against the
  installed package); no `--preclean`, because nothing in `src/` changes.
- **`rngSeed` pinned in every fit.** Chain seeds come from a dedicated
  Mersenne twister and sampling never advances R's stream
  (`R_interface_bartcore.cpp:1811-1817`), so an identically built and
  identically driven sampler is bitwise reproducible. Two seeds per cell,
  as in `grouped-mixing.R`: data `set.seed(BASE_SEED + s)`, sampler
  `rngSeed = s`; matched arms share `s`.
- **No NA in the mutated column, ever.** The emulation (below) relies on
  `collapseEmptyNodes` and `dropStaleMissingDirections` being no-ops.
- **Arm B pins identical cut grids** on the shared latent column across all
  surrogate samplers (`setCutPoints(cuts, nm)` on each). The joint path
  quantizes against each sampler's own store
  (`sampler_.data_.codeFor`, `sampler.hpp:1334-1336`), so without pinning
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
   for `f = 0..numForests-1` (`R/bartcore.R:961-986`), plus
   `state[[c]]$variance.{vars,values,sizes,flags}`
   (`R_interface_bartcore.cpp:4771-4781`) for a heteroscedastic sampler.
2. Compute the widened per-observation install mask in the R oracle.
3. Install exactly the accepted rows with
   `setPredictor(newColumn, col, forceUpdate = TRUE)`. The force path is
   open to multi-forest and heteroscedastic samplers
   (`R_interface_bartcore.cpp:2088-2095`) and `forceRefreshTrees` re-routes
   every forest and the variance forest.
4. `run(0L, 1L)`. Repeat.

Because the oracle admits no row that would empty a leaf,
`collapseEmptyNodes` is a no-op and the forced refresh reproduces exactly
the live state the widened transactional path would install. **V3b
certifies this deterministically**; failure of V3b, not any tolerance, is
what triggers the engine-patch contingency.

The export is cheap: `flattenBelow` reads `node.numObservations()` off the
node (`tree.hpp:1372-1374`) rather than re-routing, so
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
shapes).** Real BCF (`dbarts:::bartcoreBCFSampler`, `R/bartcore.R:570-664`)
and real heteroscedastic (`dbarts(x, y, variance = TRUE,
n.trees.variance =)`, `R/dbarts.R:346-349`) samplers, driven by the
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
permutation from `samplers[0]->rng()` (`facade.hpp:458`), so with `rngSeed`
pinned, replaying the identical build-and-burn script and varying only the
SUBSET passed to the joint call yields scan-order-matched masks. The
marginal of each surrogate sampler is therefore measured by the engine, not
inferred. This attributes per SAMPLER, not per FOREST inside a sampler.

**Arm C - validity.** V0-V4 below. Runs first; failure blocks every metric.

## Scenarios

| factor | levels | grounding |
|---|---|---|
| ensemble config | (a) single forest 75 trees [BASELINE, E=1]; (b) BCF mu 75 / tau 50 [E=2]; (c) BCF mu 75 / tau 25 [E=2]; (d) het mean 75 / variance 40 [E=2]; (e) surrogate S=3 mu 75 / tau 50 / var 40 [E=3, arm B only, UNGATED]; (f) multinomial K=4 x 75 [E=4, UNGATED stress] | `n.trees = 75` default (`R/A_class.R:243`); `n.trees.treatment = 50L` (`R/bartcore.R:573`), `treatment.base = 0.25` / `treatment.power = 3` (`:574-575`); `n.trees.variance = 40L` (`R/spec.R:377`); mu default `base = 0.95, power = 2` (`R/model.R:871`, `:1002-1003`); bairrtt runs 75 (`irt_causal_bart.R:204`) |
| **column type** | {continuous, binary} gated separately; {multi-level categorical <= 64 levels} reported ungated | two of four motivating classes mutate discrete columns |
| mutated column in the non-primary mask? | {yes, no} | "no" must measure exactly zero marginal (M1); "yes" is the gated cell |
| n | {300, 500, 1000, 5000} (300 in arm 0 only) | bairrtt runs 200-1000 and its prose baseline is at 300; the motivating DGP is 1000 persons x 100 items (`bairrtt/docs/plans/multi-trait.md:98-102`); 5000 probes leaf-size scaling |
| p | {5, 10}, held at 10 outside arm 0 | controls the share of j-splitting trees |
| n.chains | {1, 4}, held at 1 outside arm 0 | arm 0's lever and the shipped widening |
| phase | 300 forced burn sweeps, then a W = 200 sweep TRACKING WINDOW | bairrtt deliberately forces during burn (`irt_causal_bart.R:616-629`), so only the tracking window is gated |
| fraction of rows mutated | {1/n, 0.01, 0.1, 1.0} | 1.0 is bairrtt (`irt_causal_bart.R:563`); 1/n and 0.01 are the sequential-imputation and one-subject-Gibbs regimes, the only ones where the whole-transaction path is live |
| move size | continuous: {0.1, 0.4, 0.8} column SDs; binary: flip probability 0.5 on proposed rows | bairrtt's adapted `theta_sd` is 0.68-1.04 against a unit-SD column (`irt_causal_bart.R:197`; `bairrtt/TODO:113-114`) |
| cut grid on the mutated column | 100 interior quantiles (continuous) | `n_theta_cutpoints = 100L` (`irt_causal_bart.R:203`) |

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
codes" (`R_interface_bartcore.cpp:2077-2086`). `E_0` on a multi-ensemble
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
  `lower_bound(cuts, value) - first`, `data.hpp:588-597`; ordinal
  `ruleSendsRight` is `code > splitIndex`; the exported payload is
  `cutPoints[j][splitIndex]`, `tree.hpp:1391-1393`; flat routing is
  `value <= flat.value -> left`, `tree.hpp:1695-1697`). The `n` column is
  the LIVE partition, not a replay, when `current = TRUE` and
  `newdata = NULL` (`R_interface_bartcore.cpp:4596-4598, 5616-5634,
  5821-5830`).
- **V2r - per-row identity.** Equal counts do not imply equal per-row
  assignment (a permutation between two equal-sized leaves passes V2).
  Replay single-row `newdata` matrices through `getTrees` for 50 sampled
  rows per cell; the unit count must land in the oracle's predicted leaf.
- **V2v - the variance-forest oracle, exact.** `bartcoreRun(bc, 0L, 1L)`
  returns an element named `variance`
  (`R_interface_bartcore.cpp:3293`) filled from
  `varianceFits[i] = sigmaScale^2 * combinedVariance[i]`
  (`chain.hpp:4022-4029`), where `combinedVariance[i]` is the PRODUCT over
  variance trees of row i's leaf factor (`chain.hpp:722-726, 2789-2792`);
  the state's `variance.values` are those same working-scale factors
  (`chain.hpp:2189-2196`). So: decode `variance.*`, route every row, form
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
(`irt_causal_bart.R:198, 647-653`); a 2 pp install-rate loss shifts
effective acceptance from ~0.44 to ~0.43, which the adaptation absorbs and
which is smaller than the spread already observed across persons
(0.68-1.04 at the deciles, `bairrtt/TODO:113-114`). Above 15 pp the
widening is no longer a tuning knob: more than one proposed latent move in
seven is lost, the adaptation cannot recover it by shrinking the proposal
without collapsing the step size, and the consumer's practical option
becomes buying the cost back with `n.trees.treatment` or the moderator
mask - which is a different arc.

*This line is an EFFICIENCY line and nothing more.* The widened veto does
NOT change the estimand. The tau forest and the variance forest carry the
same empty-leaf invariant as forest 0 (`chain.hpp:1533-1543, 2666-2668,
3286-3291`), so under dbarts' convention a configuration with an empty leaf
in any ensemble has zero prior mass; a latent value that empties a tau leaf
is a zero-prior state of the JOINT model and rejecting it is the correct
Metropolis-within-Gibbs step - as bairrtt's own comment says
(`irt_causal_bart.R:605-608`). Any argument that the widening "changes the
target" would condemn the criterion already shipping at E = 1.

### Primary 2 - mixing pathology (P3)

`beta` = twice the seed-to-seed SD of the MATCHED NULL delta (config (a) at
seed s versus config (a) at seed s'), measured in Stage 0 and written into
this file at FREEZE.

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
clause fires (P3's confirmed on re-run). **YELLOW** otherwise.

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

## Instrumentation

Ships today; sufficient for every arm; **no engine change**.

| need | shipped hook |
|---|---|
| per-forest live tree structure | `dbarts:::bartcoreGetTrees(bc, chainNums, treeNums, current = TRUE, forest = f)`, `R/bartcore.R:961-986`; `forest` is 0-based and bounded by `shape.numForests` (`R_interface_bartcore.cpp:4553-4556`) |
| ground-truth node occupancy (V2) | the `n` column of that frame, LIVE partition under `current = TRUE, newdata = NULL` |
| variance-forest structure | `sampler$state[[c]]$variance.{vars,values,sizes,flags}`, `R_interface_bartcore.cpp:4771-4781`; raw-double decode precedent `inst/tinytest/test-bartcore.R:730-736` |
| variance-forest ground truth (V2v) | the `variance` element of `bartcoreRun(bc, 0L, 1L)`, `R_interface_bartcore.cpp:3293` |
| BCF construction | `dbarts:::bartcoreBCFSampler(...)`, `R/bartcore.R:570-664` |
| heteroscedastic construction | `dbarts(x, y, variance = TRUE, n.trees.variance =)`, `R/dbarts.R:346-349`; `dbartsSpec(variance =, n.trees.variance =)`, `R/spec.R:376-377` |
| multinomial handle | `bart2(x, y, family = "multinomial", keepTrees = TRUE)$bc`, `R/bart.R:1212-1213` |
| engine's single-forest per-obs mask | `sampler$setPredictor(v, col, forceUpdate = "partial")` -> length-n logical (`inst/tinytest/test-bartcore.R:259-280`) |
| engine's whole-transaction verdict | `sampler$setPredictor(v, col)` -> TRUE/FALSE |
| engine's widened multi-sampler mask | `updatePredictorPerObservationJointly(list(...), x, column)` (exported, `NAMESPACE:8`); asserted as the widened AND at `inst/tinytest/test-sampler-updatePredictorPerObservationJointly.R:118` |
| unconditional install and restore | `setPredictor(v, col, forceUpdate = TRUE)` |
| leaf-routing precedents | `inst/tinytest/test-sampler-setPredictorPerObservation.R:13-38`; `inst/tinytest/test-interactions.R:13-64`; `R/plotTree.R:44-61` |

Built harness-side (R only): a general multi-predictor pre-order router
(~70 lines: recursive index partition on `var`, `value`, `missing`,
`directions`, with `value <= cut -> left`); a `variance.*` flat decoder
(~40 lines, `readBin` over the raw slots split by `variance.sizes`); a
sequential per-observation session simulator (~50 lines, mirroring
`sampler.hpp:1354-1385`); the arm-A' driver loop.

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
.claude/d1-veto-rate-falsifier/harness/
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

## What this measurement cannot tell us

1. **The stationary widened posterior, on real trees.** Every arm measures
   rejection rates. Q1 bounds the posterior question on the arm-B
   surrogate, under three kernels, and is REPORTED - never gated, because a
   surrogate posterior comparison must not gate an engine arc. P3 is the
   gated proxy, deliberately.
2. **`E = 3` on real trees.** Not constructible (`facade.hpp:631, 644`;
   `tests/cpp/test_model.cpp:720-733`). Config (e) is a surrogate,
   reported and ungated. **No verdict here can authorize a heteroscedastic
   BCF**: that configuration does not exist and lifting the factory refusal
   is a separate, unpriced arc. A GREEN means "widen for the shapes that
   exist" - BCF at E = 2, heteroscedastic at E = 2, multinomial ungated.
3. **Which ensemble to blame (fork 2).** Every shipped surface returns a
   conjunction, never an attribution: `setPredictor(partial)` returns one
   mask over forest 0 across all chains (`sampler.hpp:1019-1034`), the
   joint call ANDs the per-sampler validities before committing
   (`facade.hpp:459-470`), and `setPredictor` returns one boolean. P4 and
   T2 therefore come from the R oracle for arm A', with only arm B's
   per-SAMPLER marginals measured by the engine (via seed-pinned subset
   replay). Per-FOREST attribution inside one sampler is oracle-only.
   **Fork 2 is REPORTED, not decided by this measurement.**
4. **Multi-level and pooled categorical latents.** Verdicts cover
   continuous and binary. Categoricals with <= 64 levels are reported,
   ungated. Categoricals with **more than 64 levels are out of scope**:
   they are pooled (`columnIsPooled` is `maskWordsForCount(numCuts[j]) > 1`,
   `data.hpp:433-436`) and pooled rules do not reach `variance.*`, which is
   flattened with `masks = nullptr` (`chain.hpp:2189-2196`) while the mean
   forests get a `tree.masks` slot (`R_interface_bartcore.cpp:4756`). A
   discrete follow-up covering >64-level columns needs the variance state to
   carry a mask side channel first. (Related, and logged for the TODO rather
   than fixed here: a heteroscedastic sampler whose store holds a pooled
   categorical column would null-dereference in `flattenBelow` -
   `masks->size()` at `tree.hpp:1397` - on `storeState`. Pre-existing,
   unrelated to this arc.)
5. **Whether the arc is worth its own cost.** This prices ONE consequence.
   The engineering estimate is out of scope, but one cost fact belongs on
   the record: "widen `treeAt`" is not a one-line change.
   `UpdateSessionImpl` sizes itself `options_.numTrees * chains_.size()`
   (`sampler.hpp:1332-1333`) and indexes
   `chains_[t / numTrees]->tree(t % numTrees)` (`:1390-1393`); for a BCF
   `options_.numTrees` is mu's count only and tau's differs
   (`numTreesInForest`, `sampler.hpp:1107`); and `Chain::tree` is hard-wired
   to `forests_[0]` (`chain.hpp:2749`). The widened session needs
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

(none yet - append here, with date and reason, before reading any gated
metric)
