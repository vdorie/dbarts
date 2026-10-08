# small-rulings-1008-r: eleven small R-surface rulings of 2026-10-08

Status: LANDED 2026-10-08 (6af82179)

agent: sonnet implementer, one, one commit per item.
rng: DRAW-CHANGING on two paths and neutral elsewhere. Item 11 changes the starting sigma of a fit whose indicator design is built sparse (the linear-model estimate where the sparse design used to fall back to the marginal sd), and moves the last bit of the draws of an indicator fit whose columns land in the engine's rank-bitmap tier (5 or more levels at a density of 0.2 or less). Item 7 changes a ppd draw at weight 0. Every other scenario of the equivalence trio is bitwise.
window: pre-release.
budget: ~1200 lines with tests.

## Goal

Eleven maintainer rulings that change the R surface and no draw elsewhere are built, each with help, a test and (where it differs from 0.9-34) a NEWS line.

## Context

Each item is a root TODO item of the same name; the ledger entry (docs/decisions.md) quotes the maintainer and states the rule.

## Constraints

R-surface tier of "Process by risk": one sonnet review, the touched test files, the snapshots, lint, air, doc-freshness and the codoc checks. No engine change under src/bartcore; the bridge only where a refusal lives there. docs/decisions.md and TODO are the orchestrator's.

## Steps

1. binary-one-class-fit (dec-B333). A probit or logistic fit whose response holds one class is fitted with one warning (R/spec.R warnSingleClass, R/xbart.R); a hazard fit keeps its subject-worded refusal. Test: test-family.R, the warning count for every encoding and the latent sampler at one class, both classes, both links.
2. predict-training-rows (dec-B341). predict with no newdata returns each type at the training rows from the stored draws, an offset given replacing the fit's. Every family's predict method. Test: equal to extract/fitted on the fit's own draws, with an offset equal to the draws shifted.
3. survival-offset-training-rows (dec-B340). survivalProbabilities with an offset and no newdata applies it at the training rows off the stored draws. Test: equals a call at newdata = the training rows with the same offset.
4. basis-swap-zero-column (dec-B346). setForestBasis takes a numeric all-zero basis column; creation still refuses one. Test: swap to zeros runs, creation refuses.
5. gp-reanchor-message (dec-B330). The refusal's text says to re-derive during burn-in before draws are saved, or to make a new sampler. Bridge text; test-capi.R and test-leaf-conversions.R.
6. category-test-offset-no-updatestate (dec-B374). setCategoryTestOffset drops updateState; Rd usage, the updateState item's list of methods. Test: the argument is R's own unused-argument error.
7. ppd-weight-zero (dec-B372). A training row at weight 0 draws its ppd at weight 1 in extract and fitted; predict(type = "ppd", weights = 0) is refused by name saying to pass 1. Help for both. Test: finite draws at weight 0, the refusal.
8. predictor-frame-whole (dec-B359, dec-B360). setPredictor(x), setTestPredictor and predict take a data frame, coded by label; a numeric matrix only where every predictor column is numeric, refused by name on a fit with a factor column. Test: each route, both ways.
9. held-value-in-fixed-only (dec-B376). A held sigma or count shape is stored only in fit$fixed; readers (plot, summary, predictive draws, log-likelihood, predict) read it there; fit$fixed carries a forest coefficient held by amplitude = fixed(), per forest. bartCause's reader of a missing sigma is read and reported only.
10. hazard-max-rows-cost (dec-B367). The help of hazard(max.rows = ) and the refusal state the cost of a row; the default stays 1e7. Test: the refusal text.
11. indicator-storage-invisible (dec-B370). An indicator expansion is built sparse or dense per column by sparseDensityThreshold; sigma's starting estimate is the linear-model one on either storage; makeModelMatrixFromDataFrame returns a plain matrix unless the caller supplied a sparse column. Test: sigma and return type identical across storage, with and without Matrix where a test can control it.

## Calls made in the build

- Items 2 and 3 share one commit: the stored fit now carries the training offset (`offset`, absent when none), which both need to replace it.
- predict with no newdata refuses `weights` and `bases` by name; rbart's predict is untouched (it needs a group.by).
- A hazard fit's survivalProbabilities with an offset and no newdata, where the fit stored a test grid, shifts that grid on the latent scale and applies the link, with no trees; where it stored none it replays the trees, with the offset in place of each subject's own. The fit stores `offset.test` for it.
- A one-class hazard fit keeps its subject-worded refusal; only probit and logistic fit.
- A numeric matrix is refused on a design with a factor column at setPredictor, setTestPredictor, setTestPredictorAndOffset and S3 predict; the package's own code that holds codes runs under withCodedPredictors().
- Item 11 takes the starting sigma of a sparse indicator design from Matrix's sparse QR. That QR is wrong on a dependent column (the zero pivot's reflector is arbitrary and later columns are orthogonalised against it: residual sum of squares 248.5 against lm.fit's 245.9 on two factors' full indicator sets), and an indicator design is dependent by construction, so a column whose pivot vanishes is dropped and the QR redone until none does; that gives lm.fit's rank and residual sum of squares. The exported builder always returns a plain matrix.
- A numeric matrix is also refused by the sampler's own `predict` and by `test =` at creation on a design with a factor column.

## Landing note

Equivalence trio at the final tip on the reference build, `--bitwise` where the harness takes it, EQUIVALENCE_CORES=2, baselines as the MANIFEST names them:
- multinomial-equivalence-80b1c8d4: 11 compared / 0 skipped, every channel bitwise.
- bcf-equivalence-1b7d730c: 15 compared / 0 skipped, every channel bitwise.
- equivalence-1b7d730c: 54 of 55 scenarios report identical draws (same RNG stream). The mover is wideFactorIndicators (155 summaries, max |z| = 3.26, one summary above 3, none above 4): a 120-level factor under factors = "indicators", whose sparse design the baseline's build answered with the marginal-sd fallback (1.818) and this one with the linear-model estimate (1.037). It is the ruled change (dec-B370); the main baseline is re-recorded at the landing as equivalence-deb3fe50 (benchmarks/baselines/MANIFEST, with its oracle).
- The four seeded-drift snapshot files pass on the reference build.

Item 3 and the sparse QR of item 11 are checked against the dense and tree-replay results in test-hazard.R and test-indicator-storage.R.

Landed 2026-10-08 as 8b8c2a0f..6af82179 after one sonnet review and two fix rounds, the same reviewer checking each. The review moved the starting sigma from a densified design to the sparse QR the ruling names, and the first fix round's QR then failed at more columns than rows and on a column of a very different scale; both are fixed and pinned, with mutants on the pivot tolerance killed on both sides. Reviewer's gates at the last code tip: the full tinytest suite serially, 19491 results, 0 failures; R CMD check --as-cran --no-manual, one NOTE (CRAN incoming metadata); the hazard-exact, heteroscedastic-exact, hazard-reduction and hurdle-reduction gates in quick mode. Calls: dec-A185. Lockstep: bartCause's dbarts-1.0 branch reads a held sigma from fit$fixed.

## Verification

The touched tinytest files under the slice's library; `lintr::lint_package()`, `air format --check .`, tools/check-rc-codoc.R, tools/check-doc-freshness.R, tools/check-win-drift.R; the four test-reproducibility-*.R files once at the end.
