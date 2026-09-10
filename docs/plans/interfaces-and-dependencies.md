# interfaces-and-dependencies

agent: sonnet for all four slices; opus only if the R-hat/ESS
  implementation (S1) needs a correctness review.
rng: neutral throughout - no engine change - EXCEPT S3's wide-factor
  auto-sparse switch (step 12), which changes the design-matrix
  REPRESENTATION for a `factors = "indicators"` fit with a wide factor;
  an indicator column's density (1/K) always lands the switch in the
  engine's RANK-BITMAP storage tier (`sparseDensityThreshold`,
  [`src/bartcore/data.hpp`](../../src/bartcore/data.hpp)), never the
  DENSIFIED tier a bitwise CSC-vs-dense claim is made for elsewhere
  (Context below); codes and tree structure are unchanged (varcount
  bit-identical), but the rank tier's leaf-mean sufficient statistics
  differ from dense in the last bit from the first draw. No EXISTING
  scenario is affected either way - none carries a wide factor - so
  step 13 adds one covering scenario to the equivalence corpus, an
  ADDITION rather than a draw change on any recorded scenario, and it
  owes no oracle under the baselines' rule (MANIFEST P17:
  a re-record that CHANGES a recorded draw must name an oracle for the
  new values; an addition does not - the ordinal/nbinom precedent this
  plan already follows elsewhere). "The equivalence trio" names the
  three bitwise-identical-draws regression harnesses this plan gates
  against throughout: gaussian (50 scenarios, `equivalence.R`), BCF (12,
  `bcf-equivalence.R`), multinomial (11, `multinomial-equivalence.R`),
  each its own script compared against a recorded baseline. S2's own
  gate is stronger still: a Surv-on-formula fit reproduces the
  matrix-interface fit of the same data at the same seed bitwise.
window: S2 and S3 both rewrite [`dbartsData`](../../R/data.R)'s formula
  branch - S2 first, S3 rebases before the predictor-handling block
  below it. S2 and S3 both need front-door.md's S2 slice (its Steps 8
  through 10: family objects, consolidation, `na.action`) landed first
  IN FULL, not just its family-object step - S2 here targets step 8's
  `hazard(breaks = NULL, max.rows = 1e7, link = c("probit",
  "logistic"))` spelling, not the formals it retires, and S3's re-attach
  interacts with the `na.action` formal step 10 adds (see Context).
  S1/S4 touch neither file.
budget: S1 ~230 R + ~150 tests + ~90 Rd; S2 ~180 R + ~200 tests; S3
  ~150 R + ~150 tests; S4 ~10 DESCRIPTION/Rd + ~60 configure output +
  ~30 tests.

Decisions, all in [docs/decisions.md](../decisions.md): dec-B97 (survival),
dec-B99 (posterior removal), dec-B100 (sparse formula columns), dec-B107
(DESCRIPTION wording, configure stubs).

## Goal

Formula LHS takes `Surv` for `aft`/`hazard`, both honour `subset`, and a hazard
fit accepts a `test` set person-period-expanded on the training grid,
`survivalProbabilities` reading stored test draws (no `keepTrees`) when test
data rode the fit call. A sparse `Matrix` column or `sparseFactor` assigned
into a data frame works in `formula`/`.` like any column, no marker; ingestion
lifts such columns out before `model.frame` and re-attaches by row name after.
`posterior` is gone: `summary` computes split-Rhat and bulk/tail ESS itself; a
new `draws()` generic replaces `as_draws_array`/`as_draws_df`, returning a
plain array `posterior::as_draws_array` still accepts from a user who has it.
DESCRIPTION drops the drop-in claim for a compatible-interface one;
`configure.ac` gains three one-release stubs for the removed 0.9-x flags.

## Context

- Survival, matrix interface only today. `aft` response via
  [`extractSurvivalResponse`](../../R/dbarts.R) (wraps
  [`parseSurvivalResponse`](../../R/dbarts.R), `log()`); status rides the
  unsubsetted `bartcore.survival` control attribute, read by the bridge
  ([`applySurvivalAttribute`](../../src/R_interface_bartcore.cpp)) at
  `numObservations` rows - `subset` must subset that channel too. Hazard
  response via [`extractSurvivalTimes`](../../R/dbarts.R), expanded by
  [`expandDiscreteTimeHazard`](../../R/dbarts.R)/
  [`appendHazardPeriodColumn`](../../R/dbarts.R); family then REMAPS to the
  binary token before any family-keyed switch runs, since `node.scale`,
  `control@binary`, `fixedUnitScale` and the weight policy all key on the
  literal token with no hazard arm (docs/design/survival.md). Both refuse
  `subset` in [`dbarts`](../../R/dbarts.R) (the `hazardTokens` block, the
  `directResponse && (family == "aft" || responseIsSurv)` block) and the
  formula interface (same blocks, plus a `Surv` hitting ["survival (Surv)
  responses are not supported by the formula"](../../R/data.R) in
  `dbartsData`). [`hazardSurvivalProbabilities`](../../R/bart.R) ALWAYS
  re-expands and calls `predict(object, bigX, type = "ev")`, requiring
  `object$fit` (["survivalProbabilities on a discrete-time hazard fit requires
  the"](../../R/bart.R)) - training is ragged. `dbarts`'s hazard block sets
  `matchedCall$test <- NULL` unconditionally.
  [`survivalProbabilities.bart`](../../R/bart.R) dispatches on `$periods`,
  never `$family` (load-bearing, docs/design/survival.md).
- [`refuseSparseFormulaColumns`](../../R/mixedMatrix.R) scans
  `all.vars(formula)`, stops on any
  [`isSparseDataFrameColumn`](../../R/mixedMatrix.R) match, called from
  [`dbartsData`](../../R/data.R) before `model.frame` (the MAIN formula site
  this plan rewrites, step 12). `stats::na.pass` is also forced at four OTHER
  `model.frame`/basis-evaluation sites this plan does not touch, each keeping
  its own sparse handling unchanged: `extractMultinomialFormulaData`'s own
  `refuseSparseFormulaColumns` call ([`R/bart.R`](../../R/bart.R), the
  multinomial formula path), the `test =` data-frame recode in
  [`R/data.R`](../../R/data.R) (already bypasses `model.frame` for a sparse
  column, x/y interface only), the `forest()` term's basis evaluation in
  [`R/formulaTerms.R`](../../R/formulaTerms.R), and its prediction-time
  counterpart in [`R/model.R`](../../R/model.R). The x/y interface already
  builds a `dbartsMixedMatrix` via
  [`sparseColumnSlices`](../../R/mixedMatrix.R), reused here. Row subsetting a
  `sparseVector`/`dgCMatrix` this way is exactly what the audit verified (`d$S
  <- M`, docs/plans/sparse-formula-audit.md); a `sparseFactor` was NOT part of
  that audit and has no `[` method at all
  ([`R/sparseFactor.R`](../../R/sparseFactor.R) defines only `show` and
  `length`), so step 12 adds a dedicated row-subset helper for it. The engine
  already stores a CSC column as a rank bitmap AT OR BELOW
  `sparseDensityThreshold` (0.2,
  [`src/bartcore/data.hpp`](../../src/bartcore/data.hpp)) and densifies
  above - a DENSIFIED build already reproduces a dense build of the same
  values bitwise (docs/design/sparse-columns.md's own equivalence-style
  component test, `testSparseEndToEnd`), but that guarantee is NOT the one
  step 12 gets: an indicator column's density is 1/K, so at any practical
  level count the auto-sparse block sits in the RANK tier instead, which
  `testSparseEndToEnd` never covers and which differs from dense in the
  last bit of the leaf-mean sufficient statistics (tree structure and the
  RNG path unmoved; step 13's own forced-sparse-vs-forced-dense test
  states the actual bound). Nothing builds an indicator-expanded factor's
  block as sparse today (`makeCategoricalModelMatrix`/
  `makeModelMatrixFromDataFrame`, [`R/utility.R`](../../R/utility.R) are
  always dense).
- `posterior` (S1 LANDED - this bullet now describes history): was
  Suggests-only. `.onLoad` used to register `as_draws_array`/`as_draws_df`
  for five classes into `posterior`'s table via
  `setHook(packageEvent("posterior", "onLoad"), ...)`. Each in
  [`R/diagnostics.R`](../../R/diagnostics.R) called
  `posterior::as_draws_array`/`as_draws_df` on a base array shaped (iteration,
  chain, variable) by [`toDrawsArray`](../../R/diagnostics.R). `summary.bart`
  called `posterior::summarise_draws` when `posteriorAvailable`, else
  `quantileSummary` (mean/sd/quantiles only); `man/summary.bart.Rd` documented
  the column set (`mean`, `median`, `sd`, `mad`, `q5`, `q95`, `rhat`,
  `ess_bulk`, `ess_tail`) - now produced unconditionally by
  [`summariseDraws`](../../R/diagnostics.R). Five tinytest files called
  `posterior::` directly (named in step 4). dec-B99 already fixes the
  replacement extractor's shape - "a base-R draws extractor returning an
  iterations by chains by variables array with dimnames" that `posterior`'s own
  constructors accept unchanged - so only the NAME is ours to pick; this plan
  calls it `draws()` (a new S3 generic, methods per class) rather than folding
  it into `extract`, since `extract`'s `type` already selects per-observation
  channels under a differently-tuned `sample`/`combineChains` axis.
- DESCRIPTION's Description ends "Also serves as a drop-in replacement for
  package 'BayesTree'."; `man/dbarts-package.Rd` carries the same phrase.
  `configure.ac` has no `AC_ARG_ENABLE`/`AC_ARG_WITH`; the three new stubs
  exist nowhere in the tree.

## Decision

`draws()` naming/shape is not a fork: dec-B99 fixes the shape (Context
above). Both forks were put to VD on 2026-09-08 and are recorded with
the choice; none remains open.

1. `summary`'s column set (VD 2026-09-08, "Use your recommendation"):
   the full set - `mean`, `median`, `sd`, `mad`, `q5`, `q95`, `rhat`,
   `ess_bulk`, `ess_tail` - computed internally, so the documented
   shape is unchanged.
2. `family = "auto"` with a `Surv` left-hand side (VD 2026-09-08, "Use
   your recommendation"): dispatch to `aft`, mirroring the matrix
   interface; the hazard family stays explicit on both interfaces.

## Constraints

- Gates: neutral throughout - the equivalence trio (Context above) IDENTICAL on
  every scenario it already has (50/12/11 today); S2 and S3 each add one
  scenario to the gaussian corpus, a PLAIN ADDITION per MANIFEST P17 (owes no
  oracle), so the count only grows (51, then 52) and no existing scenario may
  move. `lintr::lint_package()` (S1 moves exported names, S3 touches
  NAMESPACE); `R CMD check --as-cran`; `pkgdown::check_pkgdown(".")` for
  `draws.Rd`. S1's rewritten tinytest files must pass with `posterior` absent
  from the library, not merely unused.
- Out of scope: engine/bridge changes (all four items R-layer only);
  left/interval censoring, competing risks, cloglog; a `sparse()` formula term
  (REJECTED by VD - Context above, docs/plans/sparse-formula-audit.md);
  per-column code widths, a streaming range kernel
  (docs/design/sparse-columns.md's open list); `na.action` itself
  (front-door.md's item, honoured passively via `model.frame`'s kept rows). The
  wide-factor auto-sparse cutoff (S3) is a level-count threshold the
  implementer picks and comments, not a VD decision - pending a future
  constants audit, no user-facing argument (dec-B100).

## Steps

S1, posterior removal (R/diagnostics.R, R/hooks.R, tests, Rd, DESCRIPTION;
independent of S2-S4):

1. Rank-normalized split-Rhat and bulk/tail ESS in R/diagnostics.R over one
   (iteration, chain, variable) array (Vehtari, Gelman, Simpson, Carpenter,
   Burkner 2021, "Rank-normalization, folding, and localization"; this MUST
   match `posterior`'s own internals exactly, not merely the paper, since step
   2 pins against `posterior`'s recorded values). Per variable, a shared split
   feeds every leg below: split each chain in half (M chains of N draws -> 2M
   half-chains of N/2). Bulk Rhat: rank-normalize the POOLED SPLIT-CHAIN draws
   (average ranks under ties, `qnorm((rank - 3/8)/(S + 1/4))`, S total draws;
   Blom's constant c = 3/8 gives the denominator S - 2c + 1, i.e. S + 1/4, not
   the S - 1/4 a literal reading of the paper's rounded prose might suggest),
   then the ordinary Gelman-Rubin formula over the 2M half-chains,
   `sqrt(((n-1)/n*W + B/n)/W)`. Folded (tail) Rhat: fold the RAW, PRE-SPLIT,
   PRE-RANK draws first - `abs(x - median(x))` over the whole pooled variable -
   THEN split THAT folded array into 2M half-chains and rank-normalize IT the
   same way, then the same Gelman-Rubin formula; report the max of the two
   Rhats (`posterior::rhat` computes `rhat(z_scale(split_chains(x)))` and
   `rhat(z_scale(split_chains(fold_draws(x))))` and takes their max - the fold
   happens BEFORE the split and BEFORE rank-normalizing, never after). Bulk
   ESS: rank-normalize the POOLED SPLIT-CHAIN draws - the SAME
   split-then-rank-normalize array Bulk Rhat uses, not the unsplit draws - then
   run the autocorrelation-based ESS estimator below on it. Tail ESS: for each
   of the 5% and 95% quantiles, build an indicator on the RAW, UNSPLIT,
   UNRANKED pooled draws (`1(x <= quantile(x, prob))`, no rank-normalization -
   `posterior::ess_quantile` never rank-normalizes an indicator), split THAT
   indicator into 2M half-chains, and run the same ESS estimator directly on
   the 0/1 values; report the smaller of the two quantile ESS values. The
   shared ESS estimator, over 2M half-chains of a (rank-normalized or
   indicator) variable: per half-chain autocovariance via `stats::fft`
   (zero-pad, multiply by its conjugate, inverse-transform, normalize by
   lag-0), averaged across half-chains into `rho_hat_t`; Geyer's initial
   monotone sequence - sum consecutive lag pairs, stop at the first
   non-positive pair, then smooth so consecutive pair sums are non-increasing;
   ESS = (half-chains * half-chain length) / (-1 + 2*sum(kept pairs) + the last
   kept term). ~130 lines.
2. Numeric tests (test-convergence-diagnostics.R): known limits - an iid chain
   gives Rhat within 0.01 of 1, both ESS near n; two concatenated shifted
   normals give Rhat above 1.01 - and against `posterior`'s own values,
   computed once by hand and recorded as fixed literals (posterior is off the
   CI library after this slice).
3. Rewrite the summary path: `summary.bart`, `summary.bartMultinomial`,
   `printSummaryBartBody` drop the `havePosterior`/`quantileSummary` branch,
   always call new `summariseDraws(arr)` (Decision 1);
   `quantileSummary`/`posteriorAvailable` deleted. Ten `as_draws_array.*`/
   `as_draws_df.*` functions deleted; `draws.bart`/`draws.bartMultinomial`/
   `draws.bartOrdinal`/`draws.bartNegbin`/`draws.bartHurdle` added, thin
   wrappers over the existing `bartDrawsArray`/`hurdleDrawsArray`/
   `multinomialDrawsArray` (the `draws()` naming note in Context), unchanged.
   R/hooks.R: delete `registerPosteriorMethods` and the
   `setHook`/`isNamespaceLoaded` block; delete `.onLoad` entirely, leaving
   `.onUnload`. NAMESPACE: add `export(draws)` and one `S3method(draws, x)` per
   class.
4. Rewrite the `posterior::` assertions in test-multinomial-generics.R,
   test-ordinal.R, test-convergence-diagnostics.R, test-nbinom.R and
   test-hurdle.R against `draws()` and the new summary columns;
   test-convergence-diagnostics.R keeps the step-2 numeric tests.
5. man/summary.bart.Rd: rewrite for `draws()`, drop "install posterior"
   language, keep the Rhat > 1.01 note, document `summariseDraws`'s columns as
   unconditional. DESCRIPTION: remove `posterior` from Suggests;
   `Matrix`/`survival` untouched.

S2, survival formula interface (R/dbarts.R, R/data.R's Surv refusal, R/bart.R's
hazard test path; lands before S3):

6. Delete both `!missing(subset)` guards in `dbarts`. Aft block: after
   `survivalStatus <- survival$status`, apply `survivalStatus <-
   survivalStatus[subset]` when given (`subset` is the plain index/logical
   value the x/y branches of `dbartsData` use on `y`; `matchedCall$subset`
   reaches `dbartsData` unchanged and subsets `x`/`y` itself, reading the same
   value). Hazard block: right after `extractSurvivalTimes`, subset the
   covariates, `survival$time`, `survival$status`, `offset`/`weights` BEFORE
   `expandDiscreteTimeHazard` (dec-B97), then set `matchedCall$subset <- NULL`
   (alongside `matchedCall$test <- NULL`) since the original indices no longer
   match the expanded design's N' rows.
7. Hazard's `test` refusal becomes acceptance: person-period expand `test` on
   the SAME grid as training (no event time, every subject expands to all K
   periods, as `hazardSurvivalProbabilities`'s data-frame `newdata` branch does
   today - reuse `appendHazardPeriodColumn`), set `matchedCall$test`;
   `offset.test` replicates likewise. `breaks`/`max.rows`/`link` read off the
   `hazard()` family object (front-door.md step 8), not bare formals.
8. Delete R/data.R's formula-path `Surv` refusal (["survival (Surv) responses
   are not supported by the formula"](../../R/data.R)). A `Surv`
   `model.response(modelFrame)` gets its own short-circuit mirroring `dbarts`'s
   aft/hazard blocks - parse with
   `extractSurvivalResponse`/`extractSurvivalTimes` right after
   `model.response`, before `refuseMultiColumnResponse`/`codeResponse` run;
   status is pulled off that same already-`model.frame`-subsetted response, not
   re-evaluated. Also in this step, `hazardSurvivalProbabilities`: when
   `newdata` is `NULL` and the fit carries `test` (`object[["yhat.test"]]`
   non-NULL), read hazards via `extract(object, type = "ev", sample = "test",
   combineChains = FALSE)` - on `object`, the packaged `bart` fit, NOT
   `object$fit` (the `dbartsSampler`): `extract.dbartsSampler`
   ([`R/generics.R`](../../R/generics.R)) accepts only `type = "predictors"`
   and refuses everything else, while `extract.bart`'s `sample = "test"` arm
   reads `object$yhat.test` directly - instead of re-expanding and calling
   `predict`. No `keepTrees` needed either way
   ([`refuseWithoutTrees`](../../R/generics.R) guards `predict`, not
   `extract`). Re-expand-and-predict stays for any OTHER `newdata`.
9. Tests (test-aft.R, test-hazard.R, test-response-shape.R): formula vs matrix
   bitwise-identical (both families, auto and explicit); `subset` via formula
   bitwise-identical to hand-subsetting; hazard `test =` bitwise-identical to
   re-expansion via `survivalProbabilities(fit, newdata = <same subjects>)` on
   a `keepTrees = TRUE` twin; the Decision-2 auto-dispatch-to-aft assertion.

S3, sparse columns (R/data.R's model.frame site, R/mixedMatrix.R; rebases onto
S2, touches the predictor-handling code below it):

10. First, run the audit's deferred checks - `order`, `split`, `merge`,
    `rbind`, `head`, `str`, a `saveRDS` round trip, `model.frame`'s
    row-dropping under `subset` - against a `d$S <- M` column, append results
    to docs/plans/sparse-formula-audit.md. Then replace
    [`refuseSparseFormulaColumns`](../../R/mixedMatrix.R) with a pull-out: find
    every `isSparseDataFrameColumn` in `data` by class over ALL of `data`, not
    only `all.vars(formula)` (`.` must reach them); remove into a named list,
    leaving `denseData` for `model.frame` unchanged.
11. Expand `.` by hand first: resolve term labels against a placeholder frame
    (dense columns as-is, each sparse name a zero-length stand-in) to expand
    `.`, then rewrite `formula` naming dense terms explicitly - dropping any
    sparse-column term, refusing a sparse name inside a wrapper term (`poly()`,
    `ns()`, `log()`, `offset()`, `:`/`*`), mirroring the existing
    `interactionLabels` refusal. Record which sparse names are used - only
    those are re-attached.
12. Run `model.frame` as today on the dense-only formula against `denseData`
    (`subset`, `weights`, `offset`, whatever `na.action` front-door.md
    installs, unchanged). `rownames(modelFrame)` is CHARACTER, and a
    `dgCMatrix`/`sparseVector` column carries no row names of its own, so
    resolve `pos <- match(rownames(modelFrame), rownames(data))` - integer
    POSITIONS into the original, pre-subset `data` - once, then subset each
    used sparse column by `pos`: ordinary `[` for `sparseVector` (`sv[pos]`)
    and `dgCMatrix` (`M[pos, , drop = FALSE]`, keeping every contributed
    column); for `sparseFactor`, which has NO `[` method
    ([`R/sparseFactor.R`](../../R/sparseFactor.R) defines only `show` and
    `length`, Context above), a new row-subset helper built the same way
    [`remapSparseFactorToTrainingLevels`](../../R/utility.R) already constructs
    a fresh `sparseFactor` via `newValidated` - walk `pos`, look up each row in
    `@i`/`@values` (or the reference code when absent), and re-derive
    `i`/`values`/`length` for the subsetted object. This aligns under `subset`
    and `na.action` together since `pos` is derived from the SAME model frame
    both apply to. Build `x` via `makeModelMatrix(modelFrame[termLabels])` as
    today, append the re-attached sparse columns via `sparseColumnSlices` into
    the existing `dbartsMixedMatrix` assembly - the container the x/y interface
    produces, so no bridge/predict/save-load code learns a second shape. A
    `test` argument with its own sparse columns goes through the same helper.
    Wide-factor auto-sparse, same step: inside
    `makeCategoricalModelMatrix`/`makeModelMatrixFromDataFrame` (R/utility.R),
    a factor under `factors = "indicators"` past a chosen level-count cutoff
    builds its indicator block as a `dgCMatrix` via
    `sparseColumnSlices`/`dbartsMixedMatrix`, automatically, no user-facing
    argument (dec-B100) - but an internal, UNEXPORTED override (a package
    option or a hidden test-only parameter, not documented) lets a test force
    the sparse or dense path independent of level count, since otherwise
    nothing can put the SAME factor through both.
13. Tests (test-data-sparse.R, test-sparse-factor.R, test-bart-formula.R): a
    `sparseFactor`/`dgCMatrix` column via `formula` bitwise-identical to the
    x/y interface; same under `subset`; same via `.`; same with response-NA
    rows dropped; using the step-12 override, the SAME wide-factor data forced
    through the sparse path and forced through the dense path fit the same
    model to EACH OTHER (never two different factors) - NOT bitwise, since an
    indicator column's density (1/K) puts the sparse path in the engine's
    RANK tier rather than the DENSIFIED tier `testSparseEndToEnd`'s bitwise
    claim covers (Context above); the rank tier's own last-bit divergence in
    the leaf-mean sufficient statistics (tree structure and the RNG path
    unmoved) is the bound the test actually checks against. Also add one
    scenario to
    `benchmarks/R/equivalence.R`'s corpus - a `factors = "indicators"` fit with
    a wide factor past the (default) cutoff - recording the usual channels; per
    the rng note above this is an ADDITION, so the existing 50 scenarios still
    reproduce bitwise and the new one needs no oracle, only a clean first
    recording.

S4, DESCRIPTION wording and configure stubs (independent of S1-S3):

14. DESCRIPTION: replace "Also serves as a drop-in replacement for package
    'BayesTree'." with "Also provides a BayesTree-compatible interface."
    `man/dbarts-package.Rd`: same replacement; its `bart` entry keeps "the
    \\pkg{BayesTree}-compatible interface" unchanged. configure.ac: add three
    `AC_ARG_ENABLE`/`AC_ARG_WITH` stubs - `--enable-match-bayes-tree`,
    `--enable-thread-safe-unload`, `--with-xint-size` - each an `AC_MSG_ERROR`
    naming the option, removed for 1.0-0 with no replacement (the exact-ABI
    flag and SIMD kernel selection cover the first two; `XINT_TYPE` is already
    fixed at `uint16_t` engine-wide), firing only when explicitly passed. One
    release only - NEWS states outright removal next release, tracked there and
    in a TODO line (R/tombstones.R is R-level, does not cover configure
    options). Regenerate `configure` with `autoreconf -i` (never hand-edit it);
    commit both. `configure.win` untouched - Windows never ran these checks or
    read the removed flags.
15. Test: a tests/cpp or tinytest check (skipped on Windows) runs `./configure
    --enable-match-bayes-tree` from a clean build directory, asserts non-zero
    exit and the stub's message; repeat for the other two flags.
16. NEWS 1.0-0 UPGRADING: the wording change, the three stubs and their expiry,
    the `posterior` removal and `draws()` rename, the survival
    formula/subset/test surface, the sparse formula-column surface.

## Verification

```
R CMD INSTALL -l <lib> .
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-f0236082.rds
R_LIBS=<lib> Rscript benchmarks/R/bcf-equivalence.R compare benchmarks/baselines/bcf-equivalence-f0236082.rds
R_LIBS=<lib> Rscript benchmarks/R/multinomial-equivalence.R compare benchmarks/baselines/multinomial-equivalence-f0236082.rds
  # 52/12/11 identical, no "max |z|" line (S2 recorded the 51st scenario,
  # aftformula, at 2085cba2; S3 recorded the 52nd, wideFactorIndicators,
  # at deb144d2)
R_LIBS=<lib> Rscript -e 'lintr::lint_package()'
R CMD build <clean copy> && R CMD check --as-cran dbarts_*.tar.gz
R_LIBS=<lib> Rscript -e 'pkgdown::check_pkgdown(".")'
R_LIBS=<lib> Rscript -e 'tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")'
grep -rn posterior R man inst DESCRIPTION NAMESPACE  # only summary.bart.Rd's note
sh -c 'cd $(mktemp -d) && /path/to/configure --enable-match-bayes-tree; echo "exit $?"'
  # non-zero, message names the removal
```

Expected: `bart(Surv(time, status) ~ ., data, family = "hazard")` with
`subset`/`test` matches, draw for draw, the hand-subsetted, hand-expanded
matrix-interface call at the same seed (front-door.md S1 renamed `bart2` to
`bart`; the alias now warns once per session, so the suite calls `bart`
directly); `d$dtm <- M; bart(y ~ ., data = d)`
matches `bart(x = cbind(other, M), y = y)` bitwise; the suite passes with
`posterior` absent; `R CMD check --as-cran` OK.

## Landing note, S1 (2026-09-09)

LANDED at 01634227616afa7e865bbab35bb8a0bc476db940, seven commits:

- 7c9c9c6d68b7461bc1a6babebe2f27c2e78eefc4 Add rank-normalized split-Rhat and bulk/tail ESS to R/diagnostics.R
- 15538877f44d3d74f1d6ebb0011ad46e567bd441 Rewrite the summary/draws path off posterior; add draws() generic
- 3c88d5e40a8a5aa8f44f0c89f9faa30510b48620 Rewrite posterior:: test assertions against draws()
- 62311f95c7a73cf2484e1721d91aae2100ddcb7b Drop posterior from Suggests; document draws(), update CI installs
- 74040f32a36870cbaedf3a38582d8faa87c658d4 Fix the two remaining posterior-removal doc/vignette references
- 68253ca3bd44e4349d82dad982e45ef26af14625 Fix three edge-case bugs in the Rhat/ESS estimator
- 01634227616afa7e865bbab35bb8a0bc476db940 Nits: drop a dead tbl_df branch, fix the plan's rank-normalize formula

[`summariseDraws`](../../R/diagnostics.R) computes rank-normalized split-Rhat
(bulk and folded) and bulk/tail ESS per Vehtari, Gelman, Simpson, Carpenter
and Burkner (2021), matching `posterior`'s own internals exactly - fold
before split and rank-normalize, Blom ranks `(rank - 3/8)/(S + 1/4)`, FFT
autocovariance, Geyer's initial monotone sequence. Against posterior 1.7.0
on 42 arrays and edge cases the max relative differences are rhat 2.9e-16,
ess_bulk 5.1e-15, ess_tail 1.7e-15; NA/NaN/Inf draws return NA as posterior
does, except at 2-3 draws per chain, where posterior's dim-dropping split
returns a number and dbarts's NA is the defensible one. `summary` always
calls `summariseDraws`'s nine columns (retired: [`quantileSummary`](../../R/diagnostics.R),
retired: [`posteriorAvailable`](../../R/diagnostics.R) deleted); the ten
`as_draws_array`/`as_draws_df` methods are deleted and a
[`draws`](../../R/generics.R) generic added (bart, bartMultinomial,
bartOrdinal, bartNegbin, bartHurdle) returning the plain (iteration, chain,
variable) array `posterior::as_draws_array` still accepts.
retired: [`registerPosteriorMethods`](../../R/hooks.R) and `.onLoad` are
deleted, `.onUnload` kept; `posterior` is off Suggests and off both CI
install lists.
[man/summary.bart.Rd](../../man/summary.bart.Rd) rewritten,
[man/draws.Rd](../../man/draws.Rd) new (no `posterior::` example - an
unstated dependency in examples fails check); tests rewritten against
`draws()`, pinned to literals recorded from posterior 1.7.0, passing with
posterior hidden. One design-doc cite and one vignette sentence updated too.

Real diff: 18 files, +631/-362; R 361 changed lines against ~230 budgeted,
tests 245 against ~150, Rd 226 against ~90 (summary.bart.Rd split in two) -
about 1.8x the S1 budget line, no fork.

Gates, run independently: tinytest 8128/0 with posterior visible and hidden
(requireNamespace FALSE confirmed fresh); equivalence 50/12/11 identical, 0
skipped, no "max |z|" line; `R CMD check --as-cran` OK, 0 notes, from a
clean tarball with posterior hidden; `lint_package`, `air format --check`,
`pkgdown::check_pkgdown` clean; doc-freshness and rc-codoc exit 0. Mutation:
dropping the fold step fails the pinned tail-Rhat literal (1.178 to 0.999);
dropping rank-normalization moves the pinned ess_bulk outside tolerance.

Review findings fixed before landing (68253ca3):
[`rhoHatT`](../../R/diagnostics.R) indexed `seq_len(maxT)` where posterior's
`1:max_t` still selects element 1 at `max_t == 0`, inflating ESS on short
chains (a 4-chain, 10-sample summary read 64 against posterior's 20; now
20); [`essQuantile`](../../R/diagnostics.R) and `summariseDraws`'s q5/q95
lacked posterior's NA guard, erroring on an NA draw;
[`rankNormalizeMatrix`](../../R/diagnostics.R) ranked NA as finite, so rhat
read 4.62 where posterior returns NA. 01634227 drops a dead tbl_df test
branch and corrects the plan's formula to `S + 1/4`.

Remaining: S2 (Surv on a formula), S3 (sparse formula columns), S4
(DESCRIPTION wording, configure stubs); front-door S4 owns
[man/bart.Rd](../../man/bart.Rd)'s own-class summary prose and the NEWS
1.0-0 passages still describing as_draws.

## Landing note, S2 (2026-09-09)

LANDED at c38f02d03bfcff8e335cf1fff77dcac20831ff47, fourteen commits:

- df7b554c6a2507c5ec2305a30d97fae29870bc29 Honour subset for aft/hazard, accept a hazard test set
- fdbc107ac04b133631b8636098d40f035a7e60b8 Read a hazard fit's stored test draws in survivalProbabilities
- 4526ebbc8b6df780f26e23dcba8ed365e7ac1cc4 Take a Surv left-hand side on the formula interface
- e24cadd6ee1cc14558bf486293776b963dd23fca Test the Surv formula interface, subset, and the hazard test path
- 37f1ae61b963f2a256ac51c8912adb5fd54f4667 Document the Surv formula interface, subset, and the hazard test path
- 612a5830c3a3b58d6ab2c754af4a60806c051ce2 Add the matrix interface's own hazard/aft subset argument tests
- e084bedca55bd0a15e04282d7b9911d6441049ae Admit a pre-built dbartsData object to the explicit aft/hazard guards
- 1f857d17424604df92c9340ed33bac329809164f Accept a hazard test set on the formula interface too
- fe06cb81f5f8ebe04273877d6b5940b2552b26ea Refuse a Surv formula response in xbart()
- 31513de3c69bb8098e5fda17e63bf634a49a324c Test the pre-built dbartsData object explicit-family fix
- 2085cba2fcfae8fa2e53687b7ecab9ff0c10d430 Add a Surv-on-formula aft scenario to the gaussian equivalence corpus
- caa07105ec7079dcbc22fb495e3c1a845908defb Record equivalence-0ef4c560.rds (the aftformula addition)
- 8c2d3b99fee20c65da8fc74bb2ece6f0d8d99da2 Repoint every live pin of the gaussian equivalence baseline to 0ef4c560
- c38f02d03bfcff8e335cf1fff77dcac20831ff47 Name the re-recorded gaussian baseline after its landed scenario commit

[`dbarts`](../../R/dbarts.R) drops both `!missing(subset)` guards: aft
subsets `survivalStatus` alongside `dbartsData`'s own x/y subsetting;
hazard subsets covariates, time, status, offset and weights BEFORE
`expandDiscreteTimeHazard` (dec-B97) and nulls `matchedCall$subset`, the
expanded rows no longer matching the originals. A hazard fit accepts
`test` on both interfaces - person-period-expanded on the training grid,
`offset.test` replicated, `breaks`/`max.rows`/`link` read off the
`hazard()` family object - and `survivalProbabilities` with `newdata`
`NULL` reads stored test draws via `extract(type = "ev", sample = "test")`
(`hazardSurvivalProbabilities`, [`R/bart.R`](../../R/bart.R)), no
`keepTrees` needed. [`dbartsData`](../../R/data.R)'s formula-path `Surv`
refusal is deleted; a `Surv` `model.response` short-circuits before
`refuseMultiColumnResponse`/`codeResponse`, threading the already-
subsetted time and status back to `dbarts` as attributes; `family =
"auto"` with a `Surv` left-hand side dispatches to `aft` (Decision 2), and
an explicit `family = "aft"`/`"hazard"` on a pre-built `dbartsData`
carrying those attributes (`survivalDataObject`) dispatches too.
[`xbart`](../../R/xbart.R) refuses a `Surv` formula response by name,
closing the silent gaussian-on-log-time gap the deleted refusal opened. A
Surv-on-formula aft scenario (aftformula) joins the gaussian equivalence
corpus, recorded from the reference build as
[equivalence-2085cba2.rds](../../benchmarks/baselines/equivalence-2085cba2.rds)
(recorded as 0ef4c560 before the landing rebase, renamed in c38f02d0 so the
MANIFEST names an ancestor) - an ADDITION under MANIFEST P17, all 50
predecessors reproducing c42b72af bitwise (50 compared / 1 skipped),
c42b72af demoted. Rd sentences on `bart.Rd`, `dbarts.Rd`,
`survivalProbabilities.Rd` and NEWS bullets record the surface.

Real diff: R about 256 against ~180 budgeted (1.4x), tests about 195
against ~200, inside the stop line.

Gates, run independently: tinytest 8194/0; gaussian equivalence 51/51
identical against the new file, 50 identical / 1 skipped against
c42b72af; BCF 12/12, multinomial 11/11 bitwise; the plan's own bitwise
gates as tests (formula equals matrix for aft and hazard, auto and
explicit; subset via formula equals hand-subsetting via the matrix
argument; a hazard test set on either interface equals re-expansion via
`survivalProbabilities` on a `keepTrees` twin; auto dispatch to aft); `R
CMD check --as-cran` OK from a clean tarball; `lint_package`, `air format
--check`, `pkgdown::check_pkgdown` clean; NEWS parses (311); doc-freshness
and rc-codoc exit 0. Mutation: subsetting after expansion fails the
matrix-interface hazard subset test (18750 vs 47700 rows).

Review findings fixed before landing: deleting the `Surv` refusal let
`xbart` silently fit `log(time)` as gaussian, discarding censoring (now
refused by name); explicit `family` on a pre-built Surv `dbartsData` was
refused while auto dispatched (`survivalDataObject`, e084bedc); the
required equivalence-corpus addition had been skipped (added, 2085cba2);
the formula path refused its own hazard test set where the Goal asks for
acceptance, now expanded with `dbartsData`'s own `x.test` machinery
(1f857d17).

Remaining: none.

## Landing note, S4 (2026-09-10)

LANDED at a98c4cad7e0e448e31ace23d365ec8322a0c9dec, two commits:

- 045fad8e1fe9d77355fea048356054171a4cde69 Replace the BayesTree drop-in claim and stub the removed configure flags
- a98c4cad7e0e448e31ace23d365ec8322a0c9dec Make the configure-stub check run under an automated gate

DESCRIPTION and [man/dbarts-package.Rd](../../man/dbarts-package.Rd) say
"provides a BayesTree-compatible interface" (dec-B107); configure.ac
carries `--enable-match-bayes-tree`, `--enable-thread-safe-unload` and
`--with-xint-size` as `AC_MSG_ERROR` stubs that fire only when passed,
`configure` regenerated by autoreconf 2.73 (byte-identical on a reviewer
re-run), `configure.win` untouched. NEWS and TODO record the one-release
expiry. The step-15 check lives twice: test-configure-removed-flags.R
runs from a source checkout (ten assertions via `run_test_file`; it
skips under the installed-package suite, which has no `configure`), and
lint.yaml's `configure-stubs` job runs the same three flags plus a
no-flag configure with `sh` alone on every push.

Gates (reviewer's own libs): tinytest 8167/0; gaussian 51/51 identical,
BCF 12/12, multinomial 11/11; `R CMD check --as-cran` OK from a clean
tarball; `lint_package`, `air format --check`, `pkgdown::check_pkgdown`
clean; NEWS parses (312); doc-freshness, win-drift exit 0. Mutation:
removing one stub and regenerating `configure` fails two of the ten
assertions and the lint.yaml script body.

Review finding fixed before landing: the first cut's test skipped under
every automated gate (no source tree in the installed package), so it
gained the source-checkout fallback and the lint.yaml job.

## Landing note, S3 (2026-09-09)

LANDED at 7165c3521c1afc4f9e3f4b00db72651d1303b7bf, twelve commits:

- aa3c7b4e58d58862d00bcb40dd3d6e81cd285e73 Ingest sparse formula columns by pulling them out and re-attaching
- e950538d0ba018c9dd831f9f9f465a8583678e70 Auto-sparsify a wide factor's indicator expansion
- 0655fe9afe813d08e7375ea1d33fd91147703f05 Fix offset() term detection's indexing; test the sparse formula surface
- deb144d2cdde61fc59bb30d4566d5e475489256b Add a wide-factor auto-sparse scenario to the equivalence corpus
- d7d4baddc540e2857caf6b09d9d340bb8a9e67d6 Record equivalence-3a1db387.rds; repoint every live pin to it
- b6a3db66383e8d350d116c6cf13aba85dc810bb8 Document the sparse formula surface and the auto-sparse wide factor
- dee1d18c196b7d6bc1fbb11cb001e7dc90ee4fda air format R/data.R, R/utility.R, and the two sparse test files
- ea0678b19d425f7364464b45095cdec0725d21bc Fall back to dense past the auto-sparse cutoff when Matrix is absent
- 14df072b02123075ce1b5744e2d5d0aa6a069f84 Test the sparseFactor formula path under subset and response-NA drop
- e0e79cb8597a833fd4859dcb6850507a4026b331 Repoint the gaussian baseline pins in the live plan docs to 3a1db387
- ec4838239764b4c3a7d005111450e415bbf713fb State the storage tier behind the auto-sparse path's last-bit difference
- 7165c3521c1afc4f9e3f4b00db72651d1303b7bf Name the recorded gaussian baseline after its landed scenario commit

[`dbartsData`](../../R/data.R) pulls a formula `data`'s sparse columns
(sparseVector, dgCMatrix, sparseFactor) out by class before `model.frame`
runs and re-attaches them row-subset under `subset` and `na.action`, once
`model.frame` has settled which rows survive
([`subsetSparseFactorRows`](../../R/mixedMatrix.R), sparseFactor having no
`[` method; plain `[` for the other two). `.` is expanded by hand against a
placeholder frame so the rewritten formula names only its dense terms, and
a sparse name wrapped in poly()/ns()/log()/offset()/':'/'*' is refused by
name; a list or environment `data` keeps the old refusal (no row identity
to re-attach by), as does bart2's multinomial formula ingestion. A
`factors = "indicators"` fit auto-sparsifies a wide factor's dummy
expansion past `sparseIndicatorLevelCutoff` (100 levels,
[`R/utility.R`](../../R/utility.R), a memory choice from a level-count
sweep, benchmarks/R/sparse-indicator-cutoff.R), automatic and unmarked
(dec-B100); forcing "sparse" past the cutoff without Matrix installed
falls back to dense, same as never touching the option. A wide-factor
auto-sparse scenario (wideFactorIndicators) joins the gaussian equivalence
corpus, recorded from the reference build as
[equivalence-deb144d2.rds](../../benchmarks/baselines/equivalence-deb144d2.rds)
(recorded as 3a1db387 before the landing rebase, renamed in 7165c352 so the
MANIFEST names an ancestor) - an ADDITION under MANIFEST P17, the 51
predecessors reproducing 2085cba2 bitwise (51 compared / 1 skipped),
2085cba2 demoted. The storage-tier finding: an indicator column's density
(1/K) always lands the auto-sparse block in the engine's rank-bitmap tier
below `sparseDensityThreshold` (0.2), never the densified tier
`testSparseEndToEnd` covers, so its leaf-mean sufficient statistics differ
from dense in the last bit while tree structure and varcount stay
bit-identical; the forced-sparse-vs-forced-dense test therefore compares
with a tolerance, not bitwise identity, and the plan's own bitwise
assumption for this path was corrected in the same slice.

Real diff: about 2x the budget line, judged mechanical fallout by the
reviewer, as S1's 1.8x was.

Gates, run independently (second reader's own libs at the pre-rebase tip
1dfa0663, then the reference build at 7165c352 reproducing the renamed
file 52/52 under --strict-coverage): tinytest shipped 8192/0, reference
8217/0; gaussian equivalence reference strict 52/52, shipped 52/52,
shipped against 2085cba2 51/51 identical + 1 skipped; BCF 12/12,
multinomial 11/11; `R CMD check --as-cran` OK at 9ef3d99f (pre-rebase
equivalent); `lint_package`, `air format --check`, doc-freshness and
rc-codoc clean.

Review findings fixed before landing: no test covered
`subsetSparseFactorRows` under `subset` or a response-NA drop (a mutation
passed the whole suite; two assertions added, mirroring the dgCMatrix
cases, and the mutation now fails exactly them); stale gaussian baseline
pins in four live plan docs (engine-performance.md, front-door.md,
memory-footprint-audit.md, pure-c-header.md); a docs/ path cited from a
shipped `R/data.R` comment; the sparse-vs-dense difference mischaracterized
as a representation-only change explained by floating-point summation
order, corrected to the rank-bitmap storage tier's own last-bit divergence
in leaf-mean sufficient statistics (above).

Remaining: none.
