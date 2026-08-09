# variance-forest-mutation-routing

agent: opus, every slice (each one can move a posterior). Serialized: one
  implementer, each slice lands before the next starts.
rng: neutral for every homoscedastic path, EVERY slice - the trio bitwise is
  the gate at each tip. Heteroscedastic draws move deliberately in S3-S5;
  there is no heteroscedastic baseline to re-record
  (benchmarks/R/heteroscedastic-exact.R is a falsifier script, not a MANIFEST
  entry). A trio deviation is a bug in the slice, never a re-record.
window: S1 is a memory-safety stop-loss and does not wait on anything. The
  refused transactional door (see "Doors held open") rides
  multiforest-predictor-mutation and must not be reopened here.
budget: S1 ~90 engine/bridge + ~70 R/test; S2 ~130 test + ~15 chain.hpp hook;
  S3 ~150 engine + ~90 test; S4 ~70 engine + ~60 test; S5 ~100 engine + ~80
  test.

## Goal

A heteroscedastic sampler stops reporting `s^2(x)` from a partition of data it
no longer holds. Three memory-safety holes close (an out-of-bounds heap write
on `setData` with a new n, an out-of-bounds cut-grid read on the DEFAULT
`run()` path after `setCutPoints`, and the same read after `setData`), the
global sigma stops being unpinnable behind the variance forest's back, and the
predictor-mutation surface a heteroscedastic sampler is allowed to use is
routed correctly rather than accepted and silently wrong. No `dbarts.h` change.

## Mechanism (the one fact the whole arc rests on)

`Chain` holds the mean forests in `forests_` and the variance forest in a
separate `varianceForest_` member, so `Chain::numForests()` reports 1 for a
heteroscedastic sampler. Both bridge guards that protect a second forest -
`refuseMultiForestMutation` and `refuseMultiForestTransactionalUpdate` - key on
`numForests >= 2`, so a heteroscedastic sampler passes both; and every
mutation helper (`revalidateTrees`, `rebuildFitsFromParameters`,
`recoverTreeParameters`, `applyNewData` on `forests_[0]`; `repartitionTrees`,
`dropStaleMissingDirections`, `forceRefreshTrees` over `forests_`) misses
`varianceForest_`. The only path that re-anchors a variance tree to the store
is `rebuildVarianceForest`, reached from `setState` alone. Full census, probes
and provenance: `.claude/varianceforest-mutation-design/` (memo, blind
critique, synthesis). Read the synthesis before starting - it carries the
conflict resolutions and two overturns this plan encodes.

## Binding decisions inherited (do not reopen)

1. **The transactional and per-observation entries are REFUSED, not
   repaired.** No consumer combines them with a variance forest today
   (bairrtt's per-observation IRT path is homoscedastic), and the repair is
   not the shared wiring the design memo priced: `UpdateSessionImpl` caches
   occupancy over mean trees only (`treeAt` resolves to `Chain::tree` =
   `forests_[0]`), so extending `revalidateTrees` alone lets a variance leaf
   empty with no cell-guard veto, `finalize()` returns false, and the bridge
   raises a hard error AFTER `commitObservation` has written cells with no
   journal - a reachable unrecoverable state where today there is a recoverable
   wrong answer. See "Doors held open".
2. **Sigma is not a parameter under a variance forest** - one semantic, all
   four doors. See S1 step 2.
3. **The variance forest's collapse merges geometrically**, and
   `recoverVarianceLeafValues` gives an empty bottom 1.0, not 0.0. See S3
   steps 1-2 and the rationale in "Collapse semantics".
4. **`refreshVarianceForest` is called explicitly per site**, never from
   inside `dropStaleMissingDirections` (which has three callers spanning three
   slices).
5. `n.trees.variance` / `base.variance` / `power.variance` stay
   creation-only. `applyVarianceAttributes` has one reader and one call site;
   no mutation entry re-reads them, and none is being added.

## Context (seams, all read at 84d2594)

- `VarianceForest` (chain.hpp) holds plain `Tree`s, not the templated mean
  leaf, so `mapOldCutPointsOntoNew`, `repartitionSubtree`, `collapseEmptyNodes`,
  `dropStaleMissingDirections` and `resetObservations` (tree.hpp) all apply
  verbatim. Both halves of the recover/scatter pair already exist -
  `recoverVarianceLeafValues` is the scale analogue of
  `recoverParametersFromFits`, and `rebuildVarianceForest` is the scatter. The
  note in docs/design/heteroscedastic.md section "Commit decomposition" item C2
  ("templated on the mean L and does not transfer to a differently-typed
  forest") is over-cautious and must not be read as a reason to refuse.
- Seven n-sized allocations are pinned at CREATION n and never resized:
  `Chain::meanWeights_` plus `VarianceForest`'s `indexBuffer`, `factorByTree`,
  `combinedVariance`, `meanResidual`, `divisor`, `treeResidual`
  (`VarianceForest::initialize`, `buildVarianceForest`). `resizeTestStorage`
  DOES cover `combinedVarianceTest`, and `initializeSavedTrees` /
  `setSavedSlotBase` DO cover the saved variance trees - the gap is the
  train-side mutation surface only.
- `Tree::flattenBelow` reads `data.cutPoints[j][rule.splitIndex()]` with no
  bound check. `ruleIsUnrepresentable` exists for exactly this and its only
  caller is `collapseEmptyNodesBelow`, which no variance tree reaches.
  `setCutPointsForColumn` shrinks by `assign`, retaining capacity, so the read
  is stale rather than faulting. `updateState = TRUE` is the `dbartsControl`
  default and `getState` flattens every variance tree, so the read is on the
  plain default `run()` path - `keepTrees` (default FALSE) and
  `predictVariance` are two further doors, not the only ones.
- The forced-predictor path cannot shrink a grid: `bartcore_setPredictor`
  returns `invalidCutPoints` ("number of induced cut points in new predictor
  less than previous"). `setCutPoints` and `setData` can. That asymmetry is why
  those two are the memory-safety entries and the forced swap is not.
- The empty-leaf veto (`logIntegratedLikelihoodForBranch`,
  docs/design/empty-leaf-veto.md) keeps every LIVE variance-tree bottom
  occupied, which is what makes a recover-before-repartition ordering safe.
  `rebuildVarianceForest` has no such guarantee: it repartitions a flat donor
  tree with no occupancy check, and `stateIsValid` checks well-formedness and
  strict leaf positivity but not occupancy.
- `columnMaskStateFeasible` loops `forests_` only. `installForests`'
  containment comment already claims it covers "a column-restricted variance
  forest"; the code does not, and it is vacuously true only because no variance
  tree is ever installed. S5 makes it live.
- E10/E11 (`setResponse` / `setOffset`, both `updateScale` flavors) are
  OK-WITH-CAVEAT, not OK: the scale leaf is calibrated once at creation
  (`ConstantVarianceLeaf::calibrated` from the sigma-prior triple on the
  working scale) and nothing recalibrates it or rescales `combinedVariance`
  when `sigmaScale` moves. The surface self-heals over the next sweeps; the
  PRIOR does not. Recorded, not scheduled - see "Doors held open".
- Design docs: docs/design/heteroscedastic.md (the standing contract),
  docs/design/empty-leaf-veto.md, docs/design/data-store.md (read first -
  data-adjacent work).

## Defects fixed, and where

| defect | lands | failing-first test |
|---|---|---|
| `setData` with a new n: OOB heap write in `formMeanWeights`, segfault on the next `run` | S1 refuses, S4 repairs | E1 probe: n 200 -> 5000 then `run(0, 5)`; refusal form in S1, survive-and-finite form in S4 |
| `setCutPoints`: OOB read in `flatten` on the default `run()` path; 21/21 serialized thresholds illegal; the sampler cannot restore its own state | S1 refuses, S3 repairs | flat-state legality: every serialized variance threshold is a member of the new grid |
| forced `setPredictor` / `updatePredictor`: variance forest routes rows by the OLD predictors | S1 refuses, S3 repairs | distinct-value bound on an all-identical design (correct bound 1) |
| unforced transactional + per-observation: same, through `revalidateAllChains` | S1 refuses (door recorded) | the refusal, both entries |
| `setModel` unpins sigma; `installForest` overwrites the pin; `Chain::setSigma` unguarded at the engine | S1 | pinned-sigma test with a homoscedastic non-vacuity arm |
| warm start silently cold-starts the variance forest and accepts a shape mismatch | S1 refuses the mismatch, S5 repairs the match | state-level: donor and destination variance trees agree immediately after install |
| `samplePriorPredictive(type = "ppd")` returns a homoscedastic prior predictive | S1 | the refusal, with a `type = "ev"` arm proving it still works |
| an installed variance tree with an empty bottom reads factor 0.0: `s.test == 0`, and a flat state that fails its own positivity check | S3 (source guard), S5 (install guard) | tests/cpp: hand-edit a valid state to empty one variance bottom, round-trip it |
| `growFromRoot` scans against unit precision instead of `w/s^2` | S5 | init-quality test vs the same sampler with a following sweep |
| `numVarianceTrees > 0` silently dropped by the BCF and multinomial factories | S1 | tests/cpp: both factories return null |

## Constraints

- **Append-only draw neutrality.** Every added statement sits inside a null
  test on `varianceForest_`, which is null for all 35 trio scenarios, so
  neutrality is compile-time visible and measured anyway. APPEND after an
  existing `for (Forest& forest : forests_)` loop; never move, merge or
  reorder that loop, and never change what it writes.
- **No refusal message names `setState` as the recovery.** After a grid
  change `setState` fails with "state is not consistent with this sampler" -
  `stateIsValid` rejects the flattened garbage before `rebuildVarianceForest`
  is reached. Messages name the variance forest and say to build a new
  sampler.
- Refuse, never coerce. A temporary refusal says so in this plan, not in the
  message text; the message a user sees must not promise a future.
- The mean-forest instantiation of any templated helper must compile to
  today's code. A codegen or bitwise difference on a homoscedastic path is a
  falsifier, not a tolerance to widen.
- **WEAK-GATE discipline (standing, all slices).** When a slice widens what is
  compared, `statesAgree` must widen with it, and every new invariant must be
  shown to FAIL against a deliberately broken build before it counts as a
  gate. A green suite that cannot go red is not evidence.
- Out of scope: the transactional and per-observation repair (door below);
  post-creation `n.trees.variance` / `base.variance` / `power.variance`; a
  prior draw of the variance forest (which is what a correct heteroscedastic
  `type = "ppd"` would need); any `inst/include/dbarts/dbarts.h` edit.

## S1. Stop-loss refusals and the sigma pin

No repair lands here. Every tip after this slice is a coherent surface: a
heteroscedastic sampler is IMMUTABLE on the predictor side, and S3/S4 reopen
it entry by entry.

1. Refuse on `shape.hasVarianceForest`, one message text reused, each naming
   its caller: `bartcore_setData` (TEMPORARY, S4 lifts), `bartcore_setCutPoints`
   (TEMPORARY, S3 lifts), and the four transactional entries by extending
   `refuseMultiForestTransactionalUpdate`'s predicate (`bartcore_setPredictor`,
   `bartcore_updatePredictor`, `bartcore_updatePredictorPerObservation` and its
   joint variant - the forced flavors TEMPORARY, S3 lifts; the transactional
   flavors held, door below). Refusing the forced flavors here too is
   deliberate: it costs one slice of capability and buys a uniform surface
   instead of a graded one, and no tip ships a known silent-wrong answer.
2. The sigma pin, at the ENGINE, once. `Chain::setSigma` returns early when
   `varianceForest_ != nullptr` - that alone closes `installForest`'s
   `setSigma(state.sigma)` and the flat C API. `Chain::setModel` skips its
   WHOLE gaussian sigma clause under a variance forest: both branches move the
   pin, the `sigmaIsFixed` branch by overwriting the working-scale 1.0 with
   `model.sigmaEstimate / sigmaScale()` and the other by re-prioring and
   unpinning. `refusePinnedSigmaChange` stays as the user-facing refusal at
   `bartcore_setSigma`; the engine guard is the backstop for internal callers
   that have no error channel. Add the contract comment; do not add one per
   call site.
3. `Chain` retains the creation residual-prior triple (three doubles).
   `bartcore_setModel` refuses a model whose triple differs from it under a
   variance forest, so `resid.prior` is never SILENTLY inert. S5 replaces the
   refusal with a recalibration.
4. `Sampler::installForests` returns `WarmStartResult::shapeMismatch` when
   donor and destination disagree on `hasVarianceForest`, both directions.
   PERMANENT.
5. `samplePriorPredictive` refuses a heteroscedastic sampler at
   `type = "ppd"` ONLY. `type = "ev"` is a legitimate mean-surface prior draw
   that never reaches the chisq draw (the gate is
   `type == "ppd" && !responseIsBinary`) and must keep working. Re-run the
   probe at `type = "ppd"` before writing the assertion. PERMANENT.
6. `createBCFSampler` and `createMultinomialSampler` return null when
   `options.numVarianceTrees > 0`, mirroring `createSampler`'s own refusal.
   This is an asymmetry closure, not a user fix: `applyVarianceAttributes` has
   one call site (`createHolder`, which reaches only `createSampler`), so the
   option is identically 0 at both factories on every shipped surface. No R
   behavior changes; the tests/cpp case is the whole gate.

Tests: the refusals in `test-sampler-errors.R` beside the existing `setSigma`
case; the pinned-sigma test (record `getSigmas()`, `setModel(sampler$model)`,
run, assert every recorded sigma identical) WITH its homoscedastic non-vacuity
arm showing sigma moving; the E1 five-repetition probe turning five
segfaults/SIGBUS into five clean errors; a `type = "ev"` arm; the two factory
cases in tests/cpp.
Gate: trio bitwise 27/27 + 5x6 + 3x5; tests/cpp from clean plus ASAN/UBSAN
(the sigma guard is engine code on a live path); full tinytest, no snapshot
regenerated; `air format --check .`; rchk next scheduled run (the bridge moves).

## S2. Gate hardening

Test-side, plus one chain.hpp accessor. This is NOT a zero-engine-change slice
and must not be described as one. It comes BEFORE any repair because every
repair slice's gate depends on it.

The residual gap, precisely - tests/cpp is not variance-free (test_model.cpp
carries nine heteroscedastic tests including a state round trip and a
`varianceTrees` corruption negative, and test_shape.cpp asserts
`hasVarianceForest`), and none of them mutate:

1. `statesAgree` omits `varianceTrees` - the ONLY `ChainStateData` field it
   omits - so every state round trip outside that one hand-written test is
   vacuous on the variance side. Add `sameFlatTrees(x.varianceTrees,
   y.varianceTrees)`. Falsify by perturbing one leaf value.
2. No accessor exposes a variance tree to a test: the whole `ForTesting`
   inventory and `Chain::tree` are forest-0, and the public variance surface
   is value-only. Add `Chain::varianceTreeForTesting(j)` in that style, with a
   `Sampler` passthrough.
3. `fuzzInvariantViolation` covers forest 0 only, and its coverage check is
   interval arithmetic that a STALE partition satisfies (`fuzzSubtreeCovers`
   never consults a split rule). Add routing agreement: for every tree of every
   forest AND of the variance forest, every observation must sit in the leaf
   `findBottomNodeForObservation` routes it to. This is the invariant the whole
   defect violates, and it catches the mean-forest analogue for free. Its
   `isfinite(s.sigma(c))` check already exists - do not re-add it.
4. Variance product invariant: `combinedVariance[i] == prod_j
   factorByTree[j*n + i]` to the sweep's own recompute tolerance, every factor
   strictly positive.
5. `ConfigSpec` cannot express a variance forest (its fields are `name`,
   `family`, `cols`, `numChains`, `opMask`). Add `numVarianceTrees` and thread
   it into the sampler build, then add ONE heteroscedastic config with the SAFE
   mask only: `OP_SET_RESPONSE | OP_SET_WEIGHTS | OP_SET_OFFSET | OP_SET_TEST |
   OP_RUN | OP_STATE`. `OP_SET_CUTS` joins at S3, `OP_SET_DATA` at S4.
   `OP_SET_PREDICTOR | OP_UPDATE_COLUMNS | OP_PER_OBS | OP_SESSION_ABANDON`
   stay masked out for the whole arc (refused at the bridge). `OP_SET_SIGMA`
   stays masked out (bridge-refused, and after S1 an engine no-op).

Gate: item 3 must be shown to FAIL when the mean forest's repartition is
commented out, and item 1 when a variance leaf value is perturbed - otherwise
neither is a gate. tests/cpp from clean plus ASAN/UBSAN; trio bitwise (a
formality that catches an accidental hook-related edit); full tinytest.

## S3. `refreshVarianceForest`, the geometric merge, and the forced paths

1. `Tree::collapseEmptyNodes` and `Tree::mapOldCutPointsOntoNew` gain a
   merge-policy template parameter defaulting to today's weight-weighted
   arithmetic merge, so every existing call site moves by zero lines and the
   mean instantiation compiles to today's code. The variance forest passes a
   weight-weighted GEOMETRIC merge, `exp(sum_b w_b log h_b / sum_b w_b)`, with
   the `weightTotal == 0` arm the plain `exp(mean log h_b)`.
2. `recoverVarianceLeafValues` gives an empty bottom 1.0 - the multiplicative
   identity, i.e. that tree abstains where it has no training support - instead
   of 0.0, with an NDEBUG assert that a live tree never needs the fallback (the
   empty-leaf veto guarantees it). This also fixes the pre-existing
   `s.test == 0` and self-rejecting-state pair on the `setState` path.
3. New `Chain::refreshVarianceForest(const std::vector<std::vector<double>>*
   oldCutPoints)`: recover node-indexed factors from the live slab through the
   LIVE partition (this ordering is load-bearing - recovering after any
   repartition reads zeros), drop the variance trees' stale missing directions,
   optionally remap onto the rebuilt grid, optionally `resetObservations` on an
   n change, `repartitionSubtree`, collapse with the geometric merge, scatter
   through the new partition, recompute `combinedVariance` as the product.
4. Call it from `forceRefreshTrees` ONLY, appended after the existing
   `forests_` loop. Never from `dropStaleMissingDirections` - it has three
   callers spanning three slices, and step 3 already drops the variance trees'
   directions internally, so the per-site discipline is structural.
5. Lift S1's `setCutPoints` and forced-predictor refusals.

Covers the forced `setPredictor` / `updatePredictor` entries, `setCutPoints`,
and the derived `predict()$variance` / `s.test` corruption.
Tests: the distinct-value bound on an all-identical design (correct bound 1)
with the fresh-fit control proving the bound is achievable and the mean-channel
arm proving non-vacuity; the `setCutPoints` bound (2 cuts x 3 columns, ceiling
27 = 3^3 quantization cells - a combinatorial ceiling with no tolerance to
tune); flat-state legality (every serialized variance threshold a member of the
new grid) and a `setState` round trip after the grid change; the empty-bottom
state round trip in tests/cpp; a statistical-agreement test (a mutated-then-run
sampler's posterior mean `s(x)` agrees with a from-scratch fit on the
post-mutation data, per docs/plans/post-mutation-assertions.md) - the only test
that catches a repair that re-routes correctly but scatters the wrong factors.
Widen the fuzz mask with `OP_SET_CUTS`.
Gate: trio bitwise - this is the neutrality-critical slice, since
`forceRefreshTrees` is on the `setdata` and `quants` scenario paths and the
merge template touches `tree.hpp`; tests/cpp from clean plus ASAN/UBSAN; full
tinytest. `R CMD INSTALL --preclean` (tree.hpp and chain.hpp are headers) and
delete the `benchmarks/kernels` binaries (no header dependency tracking there).
Falsifier: any trio scenario moving by one ULP means the merge template did
not instantiate to today's code on the mean path - abort, do not re-record.

## S4. `setData`

1. `Sampler::setData`: recover the variance leaf factors alongside
   `recoverTreeParameters`, before the store moves.
2. `Chain::applyNewData`: resize all seven n-sized allocations
   (`meanWeights_`, and the variance forest's `indexBuffer`, `factorByTree`,
   `combinedVariance`, `meanResidual`, `divisor`, `treeResidual`), then call
   `refreshVarianceForest` with the old grid, appended after the existing
   forest-0 body.
3. Lift S1's `setData` refusal; convert the E1 test from its refusal form to
   its repair form.

Tests: E1 at n 200 -> 5000 -> 200, then `run(0, 5)`, asserting the run
completes AND every reported variance is finite AND the distinct-value bound
holds - a test that only asserts `setData` does not error is vacuous, because
`setData` returns cleanly today and the fault is on the following `run`
(4 of 5 repetitions segfault, the fifth survives with non-finite variance);
the fixed-n all-identical bound; statistical agreement against a from-scratch
fit on the replacement data.
Widen the fuzz mask with `OP_SET_DATA`.
Gate: as S3, plus a MANDATORY ASAN/UBSAN leg on the fuzz run and on the E1
probe - a plain run survives this overflow 1 time in 5, so a green non-ASAN run
is not evidence.

## S5. Warm start, grow-from-root, and the setModel recalibration

1. `Sampler::installForests` carries `varianceTrees` into the reassembled
   `ChainStateData`; `Chain::installForest` rebuilds them - same grid through
   `rebuildVarianceForest`, cross grid through `refreshVarianceForest`'s remap
   arm, mirroring `rebuildLiveForestRemapped`.
2. `columnMaskStateFeasible` extends to the variance forest, using
   `varianceForest_->columnMask`. The containment comment in `installForests`
   already claims this coverage and the code does not deliver it; it is
   vacuously true only while no variance tree is installed, and step 1 ends
   that. Without it a donor's variance tree can split on a column the
   destination's `variance = ~ subset` forbids, mis-scored by
   `splitVariableLogProbability`.
3. Reject an install whose rebuilt variance tree leaves a bottom unoccupied,
   for the same reason `stateIsValid` requires strict positivity: S3's 1.0
   fallback keeps it safe, this keeps it honest.
4. `growForestFromRoot` calls `formMeanWeights` and scans against
   `meanWeights_`, mirroring `run()`'s heteroscedastic pre-step. The variance
   forest stays un-swept during grow - grow-from-root is an initializer and the
   following run's first sweep fits it. Document that in the engine comment.
5. `Chain::setModel` recalibrates `varianceForest_->leaf` from the model's
   residual-prior triple, converted to the working scale
   (`(model.sigmaEstimate / sigmaScale())^2 * model.sigmaRawScale`), and
   updates the retained triple. This REPLACES S1 step 3's refusal.

Tests: the warm start at the STATE level - donor and destination variance trees
identical immediately after install - not behaviorally; a behavioral probe of
this shape was run twice during design and could not separate warm from cold
(one sweep of 20 variance trees on a strong signal already reaches
`cor(s, x1) = 0.75` cold), so do not build the gate on it. Plus the shape
mismatch both directions (already landed in S1, re-run here), the column-mask
refusal, the empty-bottom install refusal, a grow-from-root init test, and a
`setModel` recalibration test replacing S1's refusal test.
Gate: as S3.

## Collapse semantics (decided; do not re-derive)

`collapseEmptyNodesBelow` fires on three disjuncts: left child empty, right
child empty, or `ruleIsUnrepresentable`. The third fires on a grid SHRINK with
both children occupied - reachable from `setCutPoints` (S3) and `setData` (S4)
- and merges the whole subtree's bottoms by weight-weighted ARITHMETIC mean.
For a multiplicative scale leaf that is AM over GM-natural quantities. In the
one-occupied-child case GM and AM agree exactly (the empty sibling carries
weight 0 in both), so the change is confined to the third disjunct.

Draw consequence: the merged factor becomes the geometric mean, which is `<=`
the arithmetic mean, so the arithmetic rule would bias the merged `s^2` high.
Both are approximations to a partition the new grid cannot express, and
`sweepVarianceForest` redraws every factor from the residuals on the next
sweep, so the difference is a one-sweep transient on heteroscedastic samplers
only. Positivity becomes structural rather than asserted: `exp` of a finite
mean of logs of positive numbers cannot be zero, so `stateIsValid`'s strict
positivity and `formMeanWeights`' division are both safe by construction.

Two alternatives were evaluated and rejected. A **zero-fit reset** (the
docs/plans/sampletreesfromprior-midchain.md pattern - return every variance
tree to root and discard the surface) is right there because forgetting is the
caller's stated intent; here the caller asked to change data, not to forget the
scale surface, and the enabling idiom is one sweep per outer-Gibbs iteration,
so a per-mutation reset leaves the variance forest permanently flat - strictly
worse than today's stale-but-informed surface. An **explicit positive floor**
on the arithmetic merge is a magic constant to choose and defend, and converts
a semantic question into a clamp.

## Doors held open (recorded, not scheduled)

- **Transactional and per-observation mutation under a variance forest.**
  Refused in S1. The enabled class is heteroscedastic latent-covariate models:
  a latent ability or confounder imputed by an outer sampler that enters BOTH
  the mean surface and the variance surface, so the noise scale depends on the
  imputed quantity (IRT response-time or careless-responding variance keyed on
  latent ability; latent-confounder sensitivity where the residual scale moves
  with the unobserved confounder). The door rides
  `multiforest-predictor-mutation`, which is the SAME work - generalizing
  `UpdateSessionImpl` and `revalidateTrees` past `forests_[0]` - and which
  already owns the unpriced modeling consequence both share: widening the
  acceptance criterion from "no leaf empties in any tree" to "... of any
  container" makes the sampler reject strictly more transactions and drops the
  partial-install rate by an unmeasured amount. Price it once, there, with one
  veto-rate measurement covering both.
- **Scale-leaf calibration staleness under `updateScale = TRUE`.** The leaf is
  calibrated once at creation and nothing rescales it when `sigmaScale` moves,
  so a response or offset swap leaves the PRIOR on the old scale while the
  surface self-heals. Same defect class as S5 step 5. Not scheduled: the
  mechanism is certain, the rescale factor is not verified against
  `GaussianResponse`, and scheduling it would be scheduling unverified algebra.
- **A prior draw of the variance forest**, which is what a correct
  heteroscedastic `samplePriorPredictive(type = "ppd")` needs. S1 refuses
  instead.
- **setState column-mask check for variance trees** (found at S5,
  pre-existing, unlisted): `stateIsValid`'s variance branch checks count,
  well-formedness, and strict positivity only, and `rebuildVarianceForest`
  has no `columnMaskSubtreeIsValid` backstop - the mean path has both. The
  INSTALL path is covered by S5's `installForests` pre-flight; a
  column-restricted variance forest restored through `setState` is not.
  Small, same scratch-pass shape as the S5 install check.

## Verification (every slice)

- `R CMD INSTALL --preclean` for S2-S5 (chain.hpp and tree.hpp are headers);
  plain `R CMD INSTALL .` suffices for S1's bridge-and-R work only if no header
  moved, which it will have - assume `--preclean`. Delete the
  `benchmarks/kernels` binaries after any header edit; `tests/cpp` tracks
  headers via `-MMD`.
- `cd tests/cpp && make clean && make && ./test_bartcore` - all pass. ASAN/UBSAN
  leg for every slice that makes new engine code reachable: S1 (the sigma
  guard), S3, S4 (mandatory - see S4's gate), S5.
- Full `tinytest::test_package("dbarts")` from a preclean install. New tests
  ADD; no snapshot is regenerated. `test-heteroscedastic.R` is 103 lines of
  statistical and structural assertions with no hardcoded draws, calls no
  sampler method, and must keep passing verbatim; `test-sampler-errors.R` gains
  refusal cases; `test-spec.R` builds a spec only. If a snapshot is forced,
  that is a signal the slice changed more than intended - stop and report.
- The trio, EVERY slice, expecting no deviation:
  `benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-c8f661a.rds`
  -> 27/27 "identical draws (same RNG stream)";
  `bcf-equivalence-99205ee.rds` -> 5 scenarios x 6 channels bitwise;
  `multinomial-equivalence-ec2a3d0.rds` -> 3 scenarios x 5 channels bitwise.
  No max-|z| line anywhere. No trio scenario carries a variance forest -
  `equivalence.R`'s only `variance` hit is a prose comment, and for BCF and
  multinomial the create-side option is silently DROPPED rather than refused
  (`applyVarianceAttributes` reaches only `createSampler`), so exposure is nil
  because no scenario sets the option, not because setting it is refused. THE
  TRIO IS NECESSARY, NOT SUFFICIENT: it proves only that the homoscedastic
  paths did not move. The new tests are the real oracle.
- `air format --check .` on any slice touching R/. rchk on the next scheduled
  run for S1 (the bridge moves).
- Speed: no slice touches `run()`'s inner loop, and no mutation entry is inside
  `bench-sampler.R`'s timed arms, so `bench-sampler.R compare
  benchmarks/baselines/bench-sampler-ab1dc52.csv` is a formality. Run it only
  if a slice ends up touching `run()` - none should. S5's
  `growForestFromRoot` change is outside the timed arms; confirm before
  skipping.

Stop conditions per docs/plans/README.md: a step fails twice, the diff exceeds
1.5x budget, or a needed change is out of scope - report and stop.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S1: a heteroscedastic sampler (`variance = TRUE`) now refuses predictor
  mutation rather than silently reporting a variance surface routed by the
  previous predictors; `setData` and `setCutPoints` previously risked memory
  corruption. Refusals on `setPredictor`, `updatePredictor`,
  `updatePredictorPerObservation` and `setData`/`setCutPoints` are lifted in
  later releases as each path is routed correctly.
- S1: `setModel` no longer unpins the residual standard deviation of a
  heteroscedastic sampler, and a warm start no longer overwrites its pinned
  value; `setModel` refuses a changed `resid.prior` rather than ignoring it.
- S1: `installForests` refuses a donor and destination that disagree on
  whether they carry a variance forest.
- S1: `samplePriorPredictive(type = "ppd")` refuses a heteroscedastic sampler
  rather than returning a homoscedastic prior predictive; `type = "ev"` is
  unaffected.
- S3: `setPredictor(forceUpdate = TRUE)`, `updatePredictor` and `setCutPoints`
  re-route the variance forest, so `predict()$variance` and `s.test` follow the
  new predictors.
- S4: `setData` works on a heteroscedastic sampler, including a change in the
  number of observations.
- S5: a warm start carries the donor's variance trees instead of cold-starting
  the variance surface, and `n.grow.sweeps` initializes against the variance-
  weighted precisions.

## Open items

- The transactional/per-observation door is the one decision for VD (see
  "Doors held open"). Recommended: refuse now, ride
  `multiforest-predictor-mutation`. Deciding evidence: a consumer that needs
  per-observation mutation under a variance forest before that arc runs.
- Design artifacts (memo, blind critique, synthesis) are durable at
  `.claude/varianceforest-mutation-design/`. The critique's verdict was STANDS
  WITH AMENDMENTS; A1-A10 are adopted here with two overturns recorded in the
  synthesis (A5's zero-variance mechanism, and A6's reachability), plus three
  findings neither document has (the default-path `getState` flatten, the
  `setModel` fixed branch as a fourth sigma door, and
  `columnMaskStateFeasible`'s missing variance arm).

## Landing notes

All five slices landed 2026-08-08/09, each behind an independent
orchestrator battery (tests/cpp plain + ASAN/UBSAN from clean, full
tinytest, the trio bitwise, air), trio bitwise at every tip.

- S1 51378c0: every stop-loss refusal + the engine sigma pin (all four
  doors); E1's five segfault/SIGBUS repetitions became five clean errors.
  Deviations: `refuseMultiForestTransactionalUpdate` gained a
  `forcedUpdate` parameter (predicate-only could not refuse the forced
  flavors); `readWarmStartState` parses donor `variance.vars`
  structurally so the shape check sees them.
- S2 23050b9: gate hardening, 93 append-only lines; both new invariants
  falsified-then-reverted (the routing-agreement check caught a disabled
  repartition in every predictor-mutating config; the `statesAgree`
  variance arm caught a 1.001 leaf perturbation AND went green with the
  widening removed - the weak-gate proof). Two extra `ForTesting`
  accessors beyond the prescribed one.
- S3 cd6af7d: merge-policy template (`ArithmeticMerge` default compiles
  to the prior code - trio bitwise is the codegen proof), geometric
  variance merge, empty bottom -> 1.0 (healing the setState flat-state
  hole), `refreshVarianceForest` from `forceRefreshTrees` only;
  setCutPoints + forced refusals lifted. The intermediate
  refusals-lifted-no-repair build proved the tests detect the original
  corruption (79 distinct values against bound 1). Honest weakness: the
  statistical-agreement test does not discriminate the repair - the
  variance forest re-fits within a few sweeps; the bounds and the
  routing invariant are the discriminating gates.
- S4 4a4a0b1: setData repaired - factors recovered BEFORE the store
  moves (the slab stride is the old n), the seven pinned allocations
  resized, refresh through the old grid appended after the forest-0
  body; refusal lifted. The mandatory ASAN leg caught a defect S4 made
  reachable: S3's recovery looped the node arena, where a free-listed
  pair reads as a bottom with a stale index range that a shrinking n
  puts out of bounds - fixed by recursing the live tree
  (`recoverVarianceLeafValuesBelow`).
- S5 dbcd255: warm starts carry variance trees (same-grid rebuild /
  cross-grid remap), with `WarmStartResult::varianceMismatch` for the
  unoccupied-bottom and out-of-mask refusals (a reused `shapeMismatch`
  text would have lied); `columnMaskStateFeasible` variance pass;
  `growForestFromRoot` scans against `w/s^2` (het grow-init correlation
  with truth 0.139 -> 0.962, homoscedastic control bit-identical);
  `setModel` recalibrates the scale leaf, replacing the S1 refusal (the
  working-scale formula verified equal to `buildVarianceForest`'s own).
  New door recorded above: the setState variance column-mask gap.

Budget note for future plans: every test-heavy slice overran its test
budget while staying at or under the engine budget; size test budgets
to the oracle actually mandated, not to the engine delta.
