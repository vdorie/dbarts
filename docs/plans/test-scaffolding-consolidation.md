# Test scaffolding consolidation

agent: Sonnet (R and test slices S1-S3, S5), Opus (S4 shifting fold, S6 engine)
rng: neutral, except S4 (shifting: unseeded ordinal and nbinom fits)
window: before the merge (dec-A52, dec-A51)
budget: ~3,300 changed lines over 6 serial slices, about half deletions

## Goal

Tests and benchmarks reach the package through exported functions and
`dbartsSampler` methods. The handle layer in R/bartcore.R keeps only
what production calls (xbart's fold views, the quiet run). The creators
that only tests call are gone, and so is the test-only wrapper file
inst/common/bartcoreHandle.R, except one per-forest tree reader. The
two test hooks on the response base class leave its vtable. The other
test hooks are non-virtual and truly additive, and they stay where they
are.

## Context

- Decision: docs/decisions.md dec-A52 and dec-A51; TODO "test scaffolding
  consolidation". Starting list: VD-H in
  docs/plans/review-2026-08-24/consolidated-report.md. That list is stale:
  four of its handle functions have since moved into
  inst/common/bartcoreHandle.R, and `dataSlotOrNULL`, `Tree::rightChildOf`
  and `Sampler::setCurrentSampleNum` are gone. Everything below was
  re-derived at 7e9e0a71 by parse-token caller counts. Comments were
  excluded. Callers were counted across R/, inst/tinytest, inst/common,
  benchmarks, tools and .github/workflows. An independent critique
  (scafcrit-*.R probes) corrected revision 1, and this revision carries
  its findings.
- The rule, applied strictly. KEEP-ADDITIVE means the helper reaches an
  engine state or a code path that no exported function, method or
  run() channel reaches. "Faster for tests" or "fewer lines" does not
  count. A function with a production caller is PROD, not scaffolding.
  Tests that call a production helper directly (`splitRhat`,
  `channelMeans`, `walkFormulaTerms` and about 40 others) are unit tests
  of production code; no second validation path exists, so they stay.
  So do tests that hand-build a spec and pass it to the exported
  `new("dbartsSampler", control, model, data)`, which is class API.
- A handle is not "the sampler". [`bartcoreSampler`](../../R/bartcore.R)
  rebuilds an engine from the sampler's current (control, model, data)
  plus an optional family. It therefore sees any slot a test mutated
  after `dbarts()` ran, which the sampler's own engine never sees.
  Appendix A5 lists every such site (the sweep: every assignment to a
  `$control`, `$model@`, `$data@` slot or a `bartcore.*` attribute,
  and every `family =`, before a creator call). Each site gets its
  public spelling, all of them verified to exist on the installed
  package:
  - aft: `dbarts(x, cbind(time, status), family = "aft")`;
  - leaf scale: `normal(scale = )`, checked through `$getCalibration`;
  - categorical and ordered columns: `factor` and `ordered` columns
    (`varTypes` 1 and 2);
  - cut counts: a per-column `n.cuts` vector on `dbartsControl`;
  - binary families: `family = "logistic"`;
  - test offset: `$setTestOffset(NULL)`.
- Creation seeding: [`createChainRngs`](../../src/R_interface_bartcore.cpp)
  draws one R uniform per chain at every unseeded creation, so removing
  a second creation moves the chain seeds.
- Public routes already bitwise (critique-confirmed):
  [`bartcoreBCFSampler`](../../R/bartcore.R) equals
  `dbarts(forests = list(forest(), forest(basis = ~ factor(z), ...)))`.
  Every creator argument has a `forest()` field. The multinomial
  creators equal `dbarts(x, factor(labels), family = "multinomial")` and
  `dbarts(dbartsData(x, counts = counts), family = "multinomial")`; the
  second spelling was checked on the installed package.
  [`predict.bartMultinomial`](../../R/generics.R) already predicts
  through `object$fit$predict`.

## Constraints

- dbarts.h and src/C_interface.cpp do not change. The response classes
  are all `final`, and no LinkingTo consumer includes an engine header,
  so S6 is not an ABI event.
- The equivalence baselines stay bitwise for bcf-equivalence.R and
  multinomial-equivalence.R. Only S4 re-records, and only equivalence.R's
  `ordinal` and `nbinom` scenarios.
- A slice migrates a file WHOLESALE: every creator, wrapper and handle in
  it. A name is deleted only in a slice after its last caller migrated.
  The slice's first gate is a grep for each name it deletes, which must
  return nothing outside the definition.
- Mapping traps (critique-found):
  - `forest` counts from 0 in the wrappers and from 1 in the methods.
  - `bartcoreSetPredictor(bc, x)` maps to
    `$setPredictor(x, forceUpdate = FALSE)`: the method's own default
    for a whole-matrix update is `forceUpdate = TRUE`.
  - `bartcoreUpdatePredictor(bc, x, cols)` maps to
    `$setPredictor(x, cols)` (both default FALSE).
- Out of scope: [`bartcore_runWithCallback`](../../src/R_interface_bartcore.cpp),
  which has no caller but is kept by the per-draw-callbacks design; the
  inst/common fixtures; the production `bartcoreSampler*` delegates.

## Steps (serial slices, each leaving a green tree)

A/B protocol (S1-S3). Before editing a site whose creator or spelling
changes, run the old and new routes side by side in a scratch script
with the same `set.seed` placement, and require `identical()` on every
recorded channel.
- For an exact-gate benchmark, which is statistical, the criterion is
  instead the same `$getCalibration` row and the same data coding
  (`data@varTypes`, `data@n.cuts`, cut points), followed by a gate
  pass.
- A test site with no bitwise public spelling is not re-recorded. It
  becomes KEEP-ADDITIVE on a raw `.Call` if it reaches something named
  in A5, or it is dropped.
- A check that becomes vacuous under the fold is rewritten against the
  public spelling, not kept. Example: "uncensored aft equals gaussian"
  is identical by construction once the handle's family is lost
  (scafcrit-aft.R).

S1. BCF-centred files (R only; neutral). No deletion.
  The tinytest files are test-bcf.R, test-bcf-creation.R,
  test-bcf-family.R, test-bcf-mutation-pins.R,
  test-bcf-zero-multiplier.R, test-blocks.R, test-interactions.R,
  test-forest-weights.R (also its multinomial creator and handle),
  test-multi-forest-seam.R (the same) and test-bcf-r5-surface.R. The
  benchmarks are bcf-equivalence.R, bcf-exact.R, bcf-exact-weak.R,
  bcf-exact-restricted.R and bcf-latent-exact.R.
  - Benchmarks create the public sampler right after their existing
    engine `set.seed`.
  - test-bcf-creation.R's public-vs-internal oracle and its
    unseeded-differs block are deleted as vacuous.
  - test-bcf-family.R's creator-family checks (`family = "logistic"`
    on the creator) move to `dbarts(..., family = "logistic")`.
  - Per-forest tree reads (forest > 0) use `forestTrees()` (A2), which
    this slice adds to inst/common/bartcoreHandle.R.
  ~800 lines.
  Gates: full tinytest; bcf-equivalence.R compare --bitwise on the
  reference build (count the identical-draws lines: all scenarios, no
  |z|); exact-gates quick for the four bcf gates.

S2. Multinomial-centred files (R only; neutral). No deletion.
  The tinytest files are test-multinomial-surface.R,
  test-multinomial-category-offset.R and
  test-multinomial-test-offset.R (both also use the BCF creator and a
  handle), test-multinomial-counts-mutation.R (the same),
  test-calibration-creation.R, test-calibration-midchain.R,
  test-composition-sequences.R, test-fits-without-offset.R,
  test-forest-basis-r5.R and test-bcf-reporting.R. The benchmarks are
  multinomial-equivalence.R, multinomial-exact.R (its
  `family = "logistic"` host and varTypes flip, A5),
  composition-matrix.R and sbc.R (its BCF and multinomial arms).
  Handle-level refusals that the methods now raise R-side are pinned by
  their R message.
  test-forest-basis-r5.R edits `attr(control, "bartcore.forests")$params`
  before calling `new()`. Re-spell that through `forest(sd = ,
  amplitude.prior.variance = )` if the A/B holds; otherwise keep it and
  name it in A5.
  ~650 lines.
  Gates: full tinytest; multinomial-equivalence.R compare --bitwise on
  the reference build; exact-gates quick for multinomial-exact.R.

S3. Single-forest files (R only; neutral in tests). No deletion.
  The tinytest files are test-bartcore.R, test-aft.R,
  test-aft-heteroscedastic.R, test-active-rows-pins.R (also its count
  creator), test-bartcore-keepfits.R, test-data-handle.R,
  test-sparse-factor.R, test-linear-leaves.R, test-gp-leaves.R,
  test-prior-init-composed-law.R, test-proposal-probs.R,
  test-sampler-bridge-errors.R and test-engine-constants.R, plus
  leafPriorChecks.R. The benchmarks are aft-exact.R, aft-hetero-pit.R,
  categorical-exact.R, logistic-reference.R, t-exact.R, negbin-exact.R
  and ordinal-exact.R.
  - Every A5 site takes its public spelling under the A/B protocol.
  - test-bartcore.R's bridge family refusals (`"cauchit"`, `"logistic"`
    on a continuous response) move to `dbarts(family = )` and pin the
    R-side message.
  - test-data-handle.R keeps its view handles, which are PROD. Its
    oracle becomes `set.seed(7); ref <- dbarts(...)`. Its view-mutation
    refusal blocks follow the A3 rule, because no production path
    mutates a view.
  - Files where `sampler` and `bartcoreSampler(sampler)` were created
    with no `set.seed` between now consume fewer uniforms. No migrated
    file hardcodes draws, but each is replayed whole.
  ~800 lines.
  Gates: full tinytest; exact-gates quick for the seven benchmarks.

S4. Production folds (R; shifting for unseeded ordinal and nbinom fits).
  - `bart2Ordinal` and `bart2Negbin` run the sampler's own engine, and
    the `bartcoreSampler` plus `$adoptPointer` second creation goes.
    Ordinal calls `sampler$run(n.burn, n.samples, updateState = FALSE)`.
    The nbinom per-sample loop keeps `bartcoreRun` on
    `list(ptr = sampler$getPointer())`, for its summed, once-per-fit GP
    warning.
  - The `bartcoreRun` sites in `bart2Multinomial` and
    `bart2MultinomialCounts` become `sampler$run(...,
    updateState = FALSE)`, and the manual `warnOnGPFallback` goes.
    keepFits is guaranteed TRUE there by the earlier refusal.
  - `predict.bartOrdinal` and `predict.bartNegbin` call
    `object$fit$predict(newdata, offset, n.threads)`. For negbin this
    changes scalar-offset handling: `$predict` recycles a length-1
    offset, and treats a lone `NA` as no offset, where
    `bartcorePredict` did neither. Pin both behaviours in
    test-negbin*.R and state them in `predict.bart`'s Rd, where the
    offset is documented. No NEWS entry: never on main.
  - Delete `bartcoreSampler`, `bartcorePredict`, the `bartcorePredict`
    alias in inst/common/bartcoreHandle.R, and the RC method
    `adoptPointer`, with its sentence in man/dbartsSampler-class.Rd
    (`getPointer`'s paragraph and the `reapply*` "called only from"
    list) and its entry and header comment in test-host-shell-pins.R's
    method census.
  Seeded fits are bitwise unchanged; prove that first on both families.
  Unseeded fits change their draws but not the posterior. ~350 lines.
  Gates (shifting): full tinytest, with RNG-locked snapshots replayed
  per whole file; tests/cpp; tools/check-rc-codoc.R. Re-record
  equivalence.R and diff against the old baseline: only `ordinal` and
  `nbinom` may move, because every fit is seeded at `set.seed(seed)`.
  Run z-mode compare against the old baseline. P17 oracle for the
  MANIFEST row: ordinal-exact.R and negbin-exact.R full-mode reruns,
  with their gaps recorded.

S5. Test-only creators, wrappers and bypass tests (R + bridge; neutral).
  - Delete `bartcoreBCFSampler`, `bartcoreMultinomialSampler`,
    `bartcoreMultinomialCountSampler`, `bartcoreMultinomialDataSampler`
    and `validateCategoryTestOffset`.
  - Delete [`bartcore_createBCF`](../../src/R_interface_bartcore.cpp)
    and its only callee `createBCFHolder` (about 110 lines), with the
    registration and declaration. `applyAmplitudeSpec` stays; the
    public route uses it.
  - Trim inst/common/bartcoreHandle.R to `forestTrees()`.
  - Apply A3 to the direct `.Call` sites. A block moves to the method
    if the method raises the same refusal. It moves to test-capi.R if
    C_interface.cpp reaches the same guard (the `refuseMultiForest*`
    family is shared). It stays raw only for a memory-safety backstop.
    Otherwise it is dropped.
  - Retire or re-point the doc cites (Verification).
  ~550 lines, mostly deletion.
  Gates: `R CMD INSTALL --preclean`; full tinytest; tests/cpp;
  lint_package; check-doc-freshness.

S6. Engine: the two response-base virtuals, plus small deletions
  (engine; neutral). Agent-made, recorded as such: the orchestrator
  adjudicated the engine scope under VD's "leave the internal functions
  when they're truly additive". Only the two virtuals shape a production
  vtable, so only they move. The non-virtual accessors (Appendix A4)
  have no vtable or layout effect and stay where they are.
  - `ResponseModel::varianceSurfaceForTesting` and
    `ResponseModel::sigmaDegreesOfFreedomForTesting` leave the base
    class. The overrides on the `final` classes `GaussianResponse` and
    `AFTResponse` lose `override` and become plain members.
  - `Chain::installedVarianceSurfaceForTesting` and
    `Chain::sigmaDegreesOfFreedomForTesting` `dynamic_cast`
    `response_.get()` to those two classes, returning null or 0 for
    every other family, as the base defaults did.
  - This beat the friend-struct route of revision 1, which moved about
    30 hook bodies into a new test header for no vtable gain; this
    route is about 40 lines.
  Small deletions in the same slice:
  - `Sampler::varianceTreeForTesting`: its test uses the public
    `chain(c).varianceTree(j)`.
  - `TResponse::estimatesResidualDfForTesting` and
    `NBResponse::estimatesDispersionForTesting`: their tests check
    behaviour instead (nu or r fixed across sweeps under a fixed spec,
    moving under a grid).
  - `Chain::setLevelGibbsForTesting`: the orchestrator's note said "no
    callers", but it has one, test_sampler.cpp's
    testLevelGibbsAutomatic. That test builds its samplers with
    `SamplerOptions::levelGibbs` at construction.
  - The `Chain::totalFits()` and `totalFitsInForest` sites that hold a
    sampler move to [`forestTotalFits`](../../src/bartcore/facade.hpp).
    The accessors stay for bare-chain sites (non-virtual).
  Update docs/design/aft-status-setter.md
  (`GaussianResponse::sigmaDegreesOfFreedomForTesting`, which still
  resolves but whose sentence says "virtual") and
  docs/design/negative-binomial.md
  (`TResponse::estimatesResidualDfForTesting`, which is deleted and so
  becomes `retired:`). ~250 lines.
  Gates: `R CMD INSTALL --preclean`; `make clean && make` in tests/cpp;
  tests/cpp plain and under ASan/UBSan; full tinytest; the three
  bitwise equivalence compares on the reference build; the four
  seeded-drift snapshot files.

Ordering: S1, S2 and S3 in any order but serial, since each has to
leave a green tree. Then S4, which needs S3 because `bartcoreSampler`
and `bartcorePredict` callers must be gone. Then S5, which needs S1-S4.
S6 is file-disjoint (src/bartcore, tests/cpp, and the two docs above,
which no R slice touches) and may run in a parallel worktree, stacking
after S5 by rebase with one merged-tree battery.

## Verification

- `R CMD INSTALL --preclean -l <lib> .` for S5 and S6, and a plain
  install otherwise. Check provenance with `dbarts:::buildInfo()$mode`.
- `tinytest::test_package("dbarts")` against the slice's library.
  Warnings are counted with `withCallingHandlers` where `$run` now warns
  on GP fallback and `bartcoreRun` did not (test-gp-leaves.R,
  test-engine-constants.R).
- `Rscript -e 'lintr::lint_package()' && air format --check . &&
  Rscript tools/check-rc-codoc.R . && Rscript tools/check-win-drift.R . &&
  Rscript tools/check-doc-freshness.R .`, each on its own exit status.
  Deleted-symbol cites become `retired:` in the slice that deletes them.
  - S4: docs/architecture.md, docs/design/threaded-predict.md,
    docs/plans/predict-surface.md.
  - S5: docs/design/bcf.md, docs/design/multinomial.md,
    docs/design/multinomial-mutation-arc.md,
    docs/design/public-surface.md,
    docs/design/core-generalization.md,
    docs/design/bart-as-a-component.md,
    docs/plans/bcf-latent-evidence.md.
- Per-slice deletion grep:
  `grep -rnw <name> R inst benchmarks src tests man docs`, which must be
  empty except for `retired:` cites.
- After S5: `grep -rn "dbarts:::bartcore" inst benchmarks` shows only
  the PROD handle functions xbart uses (`bartcoreDataHandle`,
  `bartcoreSamplerFromHandle`, `bartcoreSetModel`, `bartcoreRun`).

## Budget and risk

S1 ~800, S2 ~650, S3 ~800, S4 ~350 plus the re-record, S5 ~550, S6 ~250.
Six implementer runs, each with a second-reader diff review. S4 adds
one full-mode exact-gate pair.

1. A handle site whose engine differs from the sampler's (A5) is
   silently folded to the sampler. The fold then passes a check that has
   become vacuous. The A5 table is the list, and the reviewer checks
   each row against its A/B record.
2. A public route is not bitwise for some scenario. The A/B catches it
   before deletion. The fix is the spelling, never a re-record.
3. Mapping traps: forest indexing, and `forceUpdate` defaults. The
   reviewer diffs every `forest` and `forceUpdate` argument.
4. Method-side mirroring (`$setPredictor` updates data@x; `$setData`
   re-installs the mask) changes what a test sees. A test that relied on
   the engine and the R object diverging tested an unreachable state;
   rewrite it or drop it under A3.
5. GP-fallback warnings fire on migrated `run` calls.
6. S4's baseline push cancels an in-flight exact-gates run. Push the
   records separately.
7. S6 changes virtuals in model.hpp. A stale object bus-errors, so use
   `--preclean` and `make clean`.

## Decisions

D1, settled as an agent-made call (2026-09-26; recorded in
docs/decisions.md section A for the maintainer's mark): `$getTrees`
always reads forest 1, so a BCF treatment forest's trees and a
multinomial fit's categories after the first are unreachable from R.
The test-only `forestTrees(sampler, forest, ...)` stays, as truly
additive under the maintainer's rule, and TODO's
sampler-gettrees-forest door holds a public `forest =` for after the
merge. The alternative, adding it now (about 60 lines plus Rd), widens
the surface before the merge.

## Appendix A1. R/bartcore.R and other package-internal R

R = production occurrences outside the definition; tt, cm, bn =
tinytest, inst/common and benchmark files.

| function | R | tt | cm | bn | verdict |
|---|---|---|---|---|---|
| samplerCarriesAmplitudes, refuseAmplitudeMutation, samplerCarriesCounts, refuseCountsMutation, requireCountsCapability | dbarts.R, bartcore.R | 0-1 | 0 | 0 | PROD |
| bartcoreSamplerSetCounts, -SetCategoryOffset, -SetCategoryTestOffset, -Run, -SetPredictor, -SetResponse, -SetOffset, -SetData, -SetCutPoints, -SetTestPredictor | dbarts.R (method bodies) | 0-1 | 0 | 0 | PROD |
| validateCallback, warnOnGPFallback, resolveColumnIndex, resolveForestIndex, bartcoreNumForests, refuseMultiForestWarmStart | bart.R, xbart.R, dbarts.R | 0 | 0 | 0 | PROD |
| asCountMatrix, validateCategoryOffset | data.R, generics.R | 0 | 1 | 0 | PROD |
| bartcoreSampler | bart.R (ordinal, nbinom double creation) | 13 | 1 | 8 | FOLD -> the sampler's own engine, with A5 sites re-spelled (S1-S3); DELETE in S4 |
| bartcoreBCFSampler | none | 12 | 0 | 6 | FOLD -> `dbarts(forests = list(forest(...), forest(basis = ~ factor(z), ...)))` (S1-S2); DELETE in S5 |
| bartcoreMultinomialSampler | none | 11 | 0 | 4 | FOLD -> `dbarts(x, factor(labels), family = "multinomial")`; DELETE in S5 |
| bartcoreMultinomialCountSampler | none | 5 | 0 | 2 | FOLD -> `dbarts(dbartsData(x, counts = ), family = "multinomial")`; DELETE in S5 |
| bartcoreMultinomialDataSampler, validateCategoryTestOffset | only the two above | 0 | 0 | 0 | DELETE (S5) |
| bartcorePredict | generics.R (predict.bartOrdinal, predict.bartNegbin) | 6 | 1 (alias) | 0 | FOLD -> `$predict` (S4) |
| bartcoreRun | bart.R, xbart.R | 19 | 2 | 15 | PROD for xbart's views and the quiet nbinom loop; bart.R's multinomial and ordinal sites FOLD -> `$run` (S4); test callers -> `$run` |
| bartcoreDataHandle, bartcoreSamplerFromHandle | xbart.R | 2-3 | 1 | 0 | PROD; their view tests are KEEP-ADDITIVE, since a fold view's draws are what xbart never returns |
| bartcoreSetModel | xbart.R | 3 | 1 | 0 | PROD; test callers -> `$setModel` |
| dbartsSampler$adoptPointer | bart.R (the double creations) | 1 (census) | 0 | 0 | DELETE (S4), with the Rd sentence and the census entry |
| buildInfo | none | 5 | 0 | 0 (tools 1, gh 2) | KEEP-ADDITIVE: build mode and SIMD level, read by CI and the reviewer checklist |

## Appendix A2. inst/common/bartcoreHandle.R wrappers

All 29 are test and benchmark only. They go in S5, except `forestTrees()`; the `bartcorePredict` alias goes in S4.

| wrapper | tt / bn files | verdict |
|---|---|---|
| bartcoreSetCounts, -SetCategoryOffset, -SetCategoryTestOffset | 2/1, 1/1, 1/1 | FOLD -> `$setCounts`, `$setCategoryOffset`, `$setCategoryTestOffset` |
| bartcoreSetForestBasis, -SetForestWeights | 4/1, 2/0 | FOLD -> `$setForestBasis`, `$setForestWeights` (forest + 1) |
| bartcoreForestAmplitudes, -ForestFits, -ForestVariableCounts | 8/5, 10/8, 7/2 | FOLD -> `$getForestAmplitudes`, `$getForestFits`, `$getForestVariableCounts` (forest + 1) |
| bartcoreFitsWithoutOffset, -ForestCalibration, -SetForestPriorScale | 1/0, 2/0, 1/0 | FOLD -> `$getFitsWithoutOffset`, `$getCalibration`, `$setCalibration` |
| bartcoreSetOffset, -SetResponse, -SetActiveRows, -SetWeights, -SetTestOffset, -SetData, -SetTestPredictor | 5, 7, 1, 3, 4, 7, 6 | FOLD -> the same-named methods |
| bartcoreSetPredictor, -UpdatePredictor, -UpdatePredictorPerObservation | 6, 3, 5 | FOLD -> `$setPredictor(x, forceUpdate = FALSE)` (the method defaults to TRUE for a whole matrix; scafcrit-txn.R), `$setPredictor(x, column)`, `$setPredictor(x, column, forceUpdate = "partial")` |
| bartcoreUpdatePredictorPerObservationJointly | 2/0 | FOLD -> exported `updatePredictorPerObservationJointly()` |
| bartcoreSetCutPoints, -GetLatents, -StoreState, -SetState | 3, 2/1, 10/3, 6/1 | FOLD -> `$setCutPoints`, `$getLatents`, `$storeState`/`$state`, `$setState` |
| bartcorePredictPerForest | 0 | DELETE (duplicates `$predictForests`; no caller) |
| bartcoreGetTrees, forest = 0 | 8/0 | FOLD -> `$getTrees(..., current = )` |
| bartcoreGetTrees, forest > 0 | 6 files | KEEP-ADDITIVE as `forestTrees(sampler, forest, ...)`, calling `C_dbarts_bartcore_getTrees` on `sampler$getPointer()`. No method reads a non-first forest's trees (D1) |
| aliases bartcoreSetModel, bartcoreRun, bartcorePredict | - | go with the file |

## Appendix A3. Direct .Call sites in tests and benchmarks

| entry | files | verdict |
|---|---|---|
| lastPredictPartition, lastTestFitPartition, columnStorageIsSparse | test-engine-constants.R, test-generics-multithreaded.R | KEEP-ADDITIVE: they record which routing path ran, and the results are bitwise identical by design, so nothing else can tell the paths apart |
| setPredictParallelCutoff | test-generics-multithreaded.R, constant-predict-parallel-cutoff.R | KEEP-ADDITIVE: the constant's own tuning benchmark must vary it |
| setSIMDInstructionSet, getMaxSIMDInstructionSet | test-simd.R, test-build-info.R, check-standard.yaml | KEEP-ADDITIVE: the only way to run each SIMD kernel level in the R-loaded build |
| isValidPointer | test-binary-weight-mask.R, test-capi.R, test-sampler-state-format.R | KEEP-ADDITIVE: pointer death after serialize is not observable through a method, because `$getPointer` re-creates |
| assignInPlace | test-assignInPlace-bounds.R | KEEP: memory-safety bounds guard on a production utility |
| setForestWeights, setWeights, installForests, setPredictor, setCounts, setCategoryOffset, setCategoryTestOffset, setTestPredictorAndOffset, predict, create | test-forest-weights-r5.R, test-sampler-errors.R, test-multiforest-warmstart-refusal.R, test-mutate-sparse-valued.R, test-multinomial-*.R, test-predict-sparse.R, test-predict-code-channel.R, test-sum-to-one-tolerance.R, test-sampler-bridge-errors.R | S5 rule, per block: move to the method; else to test-capi.R where C_interface.cpp shares the guard (the `refuseMultiForest*` family is shared); else keep only a memory-safety backstop (for example negative counts, which allocate about 1.8e19); else drop |
| sampleTreesFromPrior, sampleNodeParametersFromPrior, getForestFits, storeState (sbc.R, test-active-rows-pins.R, test-bcf-mutation-pins.R) | 3 | FOLD -> methods |
| makeModelMatrixFromDataFrame (sparse-indicator-cutoff.R) | 1 | FOLD -> exported `makeModelMatrixFromDataFrame()` |

## Appendix A4. Engine test-only accessors

None has a production caller (grep-proved), and all ship. Verdicts
follow the orchestrator's adjudication (S6): only a virtual moves.

| hook (class) | verdict | where / folds to |
|---|---|---|
| varianceSurfaceForTesting, sigmaDegreesOfFreedomForTesting (ResponseModel, virtual; GaussianResponse, AFTResponse overrides) | KEEP-ADDITIVE (installed-surface pointer identity; the nu_0 + #{w > 0} df), off the vtable | S6: plain members of the two `final` classes; Chain's forwarders `dynamic_cast` |
| Chain::installedVarianceSurfaceForTesting, Chain::sigmaDegreesOfFreedomForTesting | KEEP (forwarders) | S6: body becomes the cast |
| Sampler::varianceTreeForTesting | DELETE | S6: public chain(c).varianceTree(j) |
| TResponse::estimatesResidualDfForTesting, NBResponse::estimatesDispersionForTesting | DELETE | S6: behavioural checks |
| Chain::setLevelGibbsForTesting | DELETE (one caller, testLevelGibbsAutomatic) | S6: SamplerOptions::levelGibbs at construction |
| Chain::totalFits, totalFitsInForest (unsuffixed, "consistency read for tests") | FOLD where a sampler is held | S6: SamplerBase::forestTotalFits; stays for bare-chain sites |
| Chain: leafOf, leafOfStale, residual, workingResponse, varianceFactors, varianceLeaf, accountStrandedLeafKStats, drawLevelShift, muByTree, fusedSuffstatRuns, checkFusedSuffstatAgainstStock (+ FusedSuffstatCheck), treeFits, forestTreeFits | KEEP-ADDITIVE, in place (non-virtual): internal caches, per-tree fits, and the states the distributional tests build | - |
| OrdinalResponse: ordinalThresholdLogAcceptance, computeScales, updateOrdinalThresholds, drawLatents; NBResponse: dispersionKernel, collapsedStatistic; NBDispersionPrior::kernelValue | KEEP-ADDITIVE, in place: kernel-level conditionals | - |
| ColumnStore: testColumnIsSparseForTesting, testSparseColumnForTesting | KEEP-ADDITIVE, in place | - |
| LinearGaussianLeaf::statisticsCacheResidentBytes (test_model.cpp only; unlisted in rev 1) | KEEP-ADDITIVE, in place: the crossproduct cache's resident capacity, which no read exposes | - |
| globals testFitPartition, predictPartition + bridge lastTestFitPartition, lastPredictPartition, setPredictParallelCutoff, columnStorageIsSparse | KEEP-ADDITIVE, in place (A3). Cost: relaxed atomic stores per test-row routing, and an O(slabs) fill per predict | - |
| bridge setSIMDInstructionSet, getMaxSIMDInstructionSet | KEEP-ADDITIVE (A3) | - |
| bridge bartcore_createBCF + createBCFHolder | DELETE | S5 |
| bridge bartcore_runWithCallback | out of scope | - |
| Chain::fusedSuffstatRuns_ member | stays: one non-atomic increment per fused tree-sweep | - |

## Appendix A5. Handle sites whose engine differs from the sampler's

The sweep covered every assignment to a `$control`, `$model@` or
`$data@` slot or a `bartcore.*` attribute, and every `family =`, before
a creator call. Each row is migrated under the A/B protocol in its
slice. Rows marked (gate) are statistical: they need matching
calibration and coding, then a pass.

| file | what the handle saw | public spelling | slice |
|---|---|---|---|
| test-aft.R, test-aft-heteroscedastic.R, test-active-rows-pins.R | gaussian host + `attr(ctrl, "bartcore.survival")` + `family = "aft"` | `dbarts(x, cbind(exp(log.t), status), family = "aft")`; A/B, since exp/log may not round-trip bitwise; "uncensored aft == gaussian" is re-stated against it | S3 |
| aft-exact.R, aft-hetero-pit.R (gate) | the same, plus `model@node.scale <-` (aft-exact) | aft two-column response + `normal(scale = )` matching `$getCalibration` | S3 |
| t-exact.R (gate) | `model@node.scale <-` | `normal(scale = )` | S3 |
| logistic-reference.R (gate) | `model@node.scale <-`, `data@offset.test <- NULL`, `family = family` | `dbarts(family = family, node.prior = normal(scale = ))`, `$setTestOffset(NULL)` | S3 |
| categorical-exact.R (gate) | `data@varTypes[1] <- 1L` or `2L`, `n.cuts[1] <-`, `offset.test <- NULL` | a `factor` or `ordered` column; per-column `dbartsControl(n.cuts = )`; `$setTestOffset(NULL)` | S3 |
| negbin-exact.R, ordinal-exact.R (gate) | `data@varTypes[1] <- 1L` + `family = "nbinom"` or `"ordinal"` | factor column + `dbarts(family = )` | S3 |
| multinomial-exact.R (gate) | `family = "logistic"` host; `data@varTypes[1] <- 1L` | `dbarts(family = "logistic")`; factor column | S2 |
| test-bartcore.R | `data@varTypes[1] <- 1L` on four hosts (category fits, bad codes, wide, over-cap); `family = "cauchit"` and `"logistic"` refusals | factor columns, with the bad-code and over-cap refusals pinned at `dbartsData`; family refusals pinned at `dbarts(family = )` | S3 |
| test-bcf-family.R | `family = "logistic"` on the BCF creator | `dbarts(forests = , family = "logistic")` | S1 |
| test-forest-basis-r5.R | `attr(control, "bartcore.forests")$params` edited before `new()` | `forest(sd = , amplitude.prior.variance = )` if bitwise; else KEEP (it reaches a raw param slot) | S2 |
| test-data-handle.R, leafPriorChecks.R | `data@n.cuts`, `model@node.prior` edited before PROD view creation | unchanged (PROD handle) | - |

## Landing

Slice S1 LANDED (pending hash): all ten BCF-centred tinytest files
(test-bcf.R, test-bcf-creation.R, test-bcf-family.R,
test-bcf-mutation-pins.R, test-bcf-zero-multiplier.R, test-blocks.R,
test-interactions.R, test-forest-weights.R, test-multi-forest-seam.R,
test-bcf-r5-surface.R) and all five BCF benchmarks
(bcf-equivalence.R, bcf-exact.R, bcf-exact-weak.R, bcf-exact-restricted.R,
bcf-latent-exact.R)
migrated off `dbarts:::bartcoreBCFSampler` and the
inst/common/bartcoreHandle.R wrappers, onto
`dbarts(forests = list(forest(), forest(basis = ~ factor(z), ...)))`
and the `dbartsSampler` methods. `forestTrees()` added to
inst/common/bartcoreHandle.R (KEEP-ADDITIVE, Appendix A2, decision D1).
No production R or C++ change.

Dropped as vacuous, per the Steps rule: test-bcf-creation.R's
public-vs-internal oracle (the positive-half comparison against
`internalSampler`), its unseeded-differs arm, and its pinned-amplitude
internal comparison (restated instead as a public-only glue check,
matching test-bcf.R's `bcFixed` pin); test-bcf-r5-surface.R's
low-level-vs-method comparison (both sides call the same `.Call`);
test-bcf-family.R's internal-creator calibration-anchor
section (duplicate of the public probit/logistic checks earlier in the
same file). Restated rather than dropped: test-bcf.R's all-zero-column
`setForestBasis` probe, which `$setForestBasis` now refuses
(`validateForestBases`) and which is pinned as that refusal; and its
`moderators = NULL` neutrality check, vacuous on either route (the
default is NULL), now `vars = colnames(x.mod)` against the omitted
default.

bcf-exact-restricted.R needed a per-column `n.cuts` (`c(2, 1)`) so
the uniform grid put one cut between each pair of adjacent cells, and
the public multi-forest route refuses per-column `n.cuts`. That
refusal is a guard, not an engine limit: the amplitude sampler's
constructor passes per-column caps to `ColumnStore::build` as the
single-forest one does (TODO multiforest-per-column-ncuts). The gate
migrated instead through the quantile grid (`useQuantiles = TRUE`,
scalar `n.cuts = 2`), which cuts at the midpoints of adjacent observed
values: the same partitions and per-column cut counts, so the exact
oracle is unchanged. Against the internal route: identical
`$getCalibration` rows for both forests, identical `varTypes`, and
bitwise-identical draws over 1200 sweeps at matched seed placement
(`data@n.cuts` reads `c(2, 2)` where the internal route read `c(2, 1)`;
the grid built from it is the same). quick and full mode pass.

Traps found beyond the plan's own list: `$setPredictor(x, forceUpdate =
TRUE)` suppresses its return value where the low-level route always
returned the logical; `$setData`/`$setModel`/`$setResponse(updateScale =
TRUE)`/`$setOffset(updateScale = TRUE)` on an amplitude-carrying sampler
refuse R-side (`refuseAmplitudeMutation`) before the bridge's generic
wording, except `updateScale = NA`, uncaught by the R-side `isTRUE()`
check, which still reaches the bridge unchanged; `$setForestWeights`'s
own length and multinomial-refusal wording differ from the bridge's;
`$setTestOffset`'s "test matrix is NULL" precondition fires before the
bridge's "have no off-sample basis" wording on a sampler with no test
predictor at all; `updatePredictorPerObservationJointly()` needs column
names on the shared design even for a single sampler.

Gates: full tinytest (shipped) 9222/9222; bcf-equivalence.R
`--bitwise` on the reference build, 15/15 scenarios identical (no
`max |z|` line), against benchmarks/baselines/bcf-equivalence-d49e2103.rds
- confirmed also in FULL (non-quick) mode, since every scenario places
`set.seed()` immediately before its (now single) creation call, matching
the exact stream position the old two-creation code left it at; exact-gates
quick for bcf-exact.R, bcf-exact-weak.R, bcf-exact-restricted.R
and bcf-latent-exact.R, all OK; `lintr::lint_package()` clean;
`air format --check .` clean; `tools/check-doc-freshness.R` OK (one stale
quoted-fragment cite in docs/plans/forest-cache-drift.md repointed to the
respelled test-bcf.R line).
