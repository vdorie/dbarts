# basis-formula-terms: a forest's basis formula predicts from the training rows' terms

Status: PLANNED (dec-B257).

agent: sonnet implementer, one; opus reviewer.
rng: NEUTRAL. No fit changes; only what `predict` builds for a forest's basis at new rows.
window: pre-release (dec-B257).
budget: ~400 lines (R ~60, tinytest ~200, manual ~80, records ~60). Plans have run 1.5-2x low.

## Goal

A basis written with a data-dependent transformation, `forest(basis = ~ scale(w))`, `poly(w, 2)`, `ns(w, 3)`,
`bs(w, df = 4)`, is rebuilt at new rows from the training rows' centre, scale and knots, as `lm` rebuilds such
a term. Today the formula is evaluated again on the new rows.

## Context

- The predictors of a fit keep their training terms and predict through them. A forest's basis does not: the
  model stores the basis formula and [`replayForestBasis`](../../R/model.R) evaluates it on the new rows.
- Measured on the branch, `basis = ~ scale(w)`, three new rows whose true fits are 0.95, 1.75, -2.95:
  `predict` returns 0.16, 0.83, -1.76, the rows centred on their own mean; the same rows predicted together
  with the training rows give 0.93, 1.75, -2.82; one new row alone is an error. `~ poly(w, 2)` is further
  off. `~ I(w - 50)`, whose value at a row does not depend on the other rows, is right.
- `stats::makepredictcall` is what `model.frame` uses to turn `scale(w)` into a call that carries the training
  centre and scale; base R and the splines package supply its methods.
- The package does not centre or scale a multiplier (dec-B257), so a caller who wants a standardized column
  writes it, and this is the spelling the help will show.

## Constraints

- A fit's draws are unchanged: the training basis is evaluated as it is now. The seeded snapshot files and
  every equivalence scenario are identical.
- As in `lm`, and no stricter: a transformation with a `makepredictcall` method predicts from the training
  rows; any other expression is evaluated on the new rows, as `lm` evaluates `I(w - mean(w))`. No refusal is
  added for expressions that depend on the other rows: it cannot be decided reliably from the data.
- With `subset`, the training centre and scale are those of the rows the fit used.
- A fit object stored before this change, which carries no rebuilt call, predicts as it does now.
- A level of a factor basis that the fit never saw is refused at predict, as `lm` refuses it; a level seen in
  training and absent from the new rows is not an error.
- What is accepted on the left of `:forest()` does not change here.
- Out of scope, each to TODO if not already there: a basis formula with `+`, `*` or `:` at its top, which is
  evaluated as R code and not as a model formula; `scale(w):forest(x)`.

## Steps

1. Where the basis term of a formula is built ([`R/formulaTerms.R`](../../R/formulaTerms.R)), store beside
   the formula the call `stats::makepredictcall` returns for the evaluated training basis; where the basis is
   rebuilt for new rows ([`replayForestBasis`](../../R/model.R)), evaluate the stored call when there is one.
   The `dbartsData(bases = )` route, which is given a matrix and no formula, is untouched.
2. tinytest: for `scale(w)`, `scale(w, scale = FALSE)`, `poly(w, 2)`, raw `poly`, `ns(w, 3)`, `bs(w, df = 4)`,
   `scale(log(w))` and `scale(cbind(w, v))`, the basis rebuilt for new rows equals what `lm`'s model frame
   gives for the same term and rows, for several rows, for one row, and for rows predicted together with the
   training rows; predictions at new rows do not depend on which other rows are predicted with them; a factor
   basis with a level absent from the new rows predicts, and an unseen level is refused; with `subset` the
   centre is the kept rows'; a fit stored without the rebuilt call predicts as before; `I(w - mean(w))` is
   evaluated on the new rows, pinned as the documented behaviour.
3. On every path that rebuilds a basis: `predict`, `fitted` and `extract` at new rows, `pdbart` and `pd2bart`
   where a basis forest is served, a sampler re-created after a reload, and `xbart` if it takes a basis.
   Each is either covered by a test or listed as not reaching the code.
4. Manual: `forest`'s `basis` and the formula-terms section of `bart` say that the right-hand side is
   evaluated as R code, that `scale()`, `poly()`, `ns()` and `bs()` predict from the training rows as in
   `lm`, that any other expression of the whole column (`w - mean(w)`) is evaluated on the new rows, with a
   pointer to `SafePrediction`, and show standardizing a multiplier.
5. Mutation: with the stored call ignored at predict, the new tests fail.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; the four seeded snapshot files
  unchanged on a reference build; the equivalence compares identical.
- bartCause's suite on a fresh install against this build.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean; `R CMD check --as-cran` shows no new note.
