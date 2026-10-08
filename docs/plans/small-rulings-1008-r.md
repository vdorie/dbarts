# small-rulings-1008-r: eleven small R-surface rulings of 2026-10-08

Status: IN PROGRESS 2026-10-08

agent: sonnet implementer, one, one commit per item.
rng: NEUTRAL, but item 11 moves bartBT's starting sigma for wide-factor data frames and item 7 moves a ppd draw at weight 0, each on its own path only.
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

## Verification

The touched tinytest files under the slice's library; `lintr::lint_package()`, `air format --check .`, tools/check-rc-codoc.R, tools/check-doc-freshness.R, tools/check-win-drift.R; the four test-reproducibility-*.R files once at the end.
