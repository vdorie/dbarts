# hazard-formula-factors: hazard fits keep the categorical design

agent: sonnet
rng: neutral (every recorded baseline and snapshot unchanged; see RNG class)
budget: ~210 lines (R ~65, tests ~130, man and design ~15)

Status: LANDED 2026-09-28 (2d7b831f, 80e0360a, 0bce4141)

## Goal

A formula fit with `family = "hazard"` (or `"hazard.logistic"`) keeps the
categorical design through the person-period expansion: factors, ordered
factors and sparse columns enter exactly as they do for every other family,
the fit equals the binary fit on the hand-expanded rows, and `predict`,
`survivalProbabilities(newdata = )`, `getTrees(newdata = )` and `extract`
take a data frame, as dec-B97 promised and dec-B78 and dec-B132 require. The
appended `period` column keeps its name, its last position and its ordinal
meaning.

## Context

Root cause. The formula branch of [`dbarts()`](../../R/dbarts.R) hands
`data@x` - the columnar container
[`makeCategoricalModelMatrix`](../../R/utility.R) built, carrying the
`term.labels`, `varTypes` and `factor.levels` attributes - to
[`expandDiscreteTimeHazard`](../../R/dbarts.R), whose non-frame branch calls
`as.matrix` and row-subsets
([R/dbarts.R:168-173](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/dbarts.R#L168-L173)).
That densifies the container to codes and drops every builder attribute;
[`appendHazardPeriodColumn`](../../R/dbarts.R) then `cbind`s the period
column onto a bare matrix, and the branch stores it
([R/dbarts.R:1129-1146](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/dbarts.R#L1129-L1146)).
The aft branch and every ordinary formula fit leave the container and its
attributes in `data@x` untouched.

Measured on the installed tip, private library, in each case against a
hand-expanded frame fit as probit:

- The defect report's "split as if ordered" does not hold. `data@varTypes`
  survives (the branch appends one ordinal entry to it), so the bridge marks
  `g` categorical and its splits are level subsets (direction masks in
  `extract(fit, "trees")`). Training draws equal the by-hand reduction
  bitwise for a factor, an ordered factor and a transformed term.
- The lost `factor.levels` table does change training when a factor
  declares a level no kept subject carries: the bridge's
  [`readDeclaredCategoryCounts`](../../src/R_interface_bartcore.cpp) finds no
  table and infers the count from the codes, so the fit differs from the
  by-hand reduction, which declares the full count.
- Every frame reader is broken, on numeric-only fits too. With no
  `term.labels` there is no term replay and with no `factor.levels` no
  coding, so [`validateXTest`](../../R/data.R) takes the indicator route: a
  factor fails the indicator-count check
  ([R/data.R:939-958](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/data.R#L939-L958)),
  a frame with the response or any extra column fails the column count, a
  transformed term fails, a constant `period` column is dropped by the
  `drop = TRUE` default, and a frame holding any sparse column, used or
  not, is refused
  ([R/data.R:827-833](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/data.R#L827-L833)).
  `factors = "indicators"` loses `drop` and `term.labels` the same way.
- A sparse column fits (densified by `as.matrix`, draws equal to the kept
  sparse design) but nothing can predict from a frame, and `test =` with a
  sparse column fails at fit time with "'x.test' must be numeric": the test
  expansion
  ([R/dbarts.R:1170-1180](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/dbarts.R#L1170-L1180))
  `cbind`s a period column onto the sparse test container.
- The x/y path with a data frame `x` is right: the expander subsets the
  frame and appends `period` before [`dbartsData`](../../R/data.R) builds the
  container, and the fit equals the formula fit bitwise and predicts from an
  extra-column frame. A numeric matrix `x` is unaffected.
- A predictor already named `period` is silently shadowed: the matrix path
  carries two `period` columns, the frame path overwrites the covariate.

Readers of a hazard fit's design:
[`hazardSurvivalProbabilities`](../../R/bart.R) rebuilds the training
subjects from the period-1 rows through `extract(object$fit, "predictors")`,
which densifies a container
([R/bart.R:2861-2879](https://github.com/vdorie/dbarts/blob/479273d329efe4e4ae1070e0389b75223a2853be/R/bart.R#L2861-L2879));
[`hazardPredictRows`](../../R/bart.R) reads `x.train` through
[`validateXTest`](../../R/data.R) and builds its placeholder from row 1 (one
subject, K rows; densifying it is harmless); the formula path's test
expansion above. `data@n.cuts` and the row-name code read only `ncol` and
row counts, which the container answers.

A prototype of Steps 1 and 2, patched over the namespace of the installed
tip, fixed every case above and passed the 16 tinytest files that mention
hazard (1389 pass, 0 fail).

## Constraints

- R only; no engine, bridge or header change.
- The appended column stays named `period`, stays last and stays ordinal.
- Out of scope, each a separate item if wanted: a bare `dgCMatrix` `x` on
  either survival family is refused with a misleading "use the matrix
  interface" message; the formula path places sparse columns after the dense
  terms on every family; the indicator route refuses a frame holding an
  unused sparse column on every family.

## Steps

1. [`expandDiscreteTimeHazard`](../../R/dbarts.R): a container row-subsets
   through [`[.dbartsMixedMatrix`](../../R/mixedMatrix.R), which keeps the
   container, its attributes and its CSC block; a plain matrix keeps its
   builder attributes (`term.labels`, `drop`, `varTypes`, `factor.levels`)
   across the row subset; a data frame is unchanged. Refuse a predictor
   already named `period` on every path ("a hazard fit appends its own
   'period' column; rename the predictor 'period'").
2. [`appendHazardPeriodColumn`](../../R/dbarts.R): a container gains a dense
   double column (the dense list, the map and `columnNames`); on every design
   kind that carries them, `term.labels` gains `"period"`, `drop` gains
   `period = FALSE`, `varTypes` gains `ORDINAL_VARIABLE` and `factor.levels`
   gains `NULL`. This also covers the formula path's sparse test container.
3. [`hazardSurvivalProbabilities`](../../R/bart.R), `newdata = NULL`: take
   the period-1 rows of the stored design by row subset, which keeps a
   container, replicate them K times and overwrite the period column, rather
   than densifying through `extract`. The result must equal the current path
   on a dense fit.
4. Tests, in a new inst/tinytest/test-hazard-factors.R, `n.chains = 1`,
   short runs:
   - `g + z`, ordered `o + z`, a factor with an unused declared level, and a
     transformed term: `yhat.train` identical to probit on the hand-expanded
     frame through `y ~ <rhs> + period`. For a `sparseFactor` the by-hand
     target is the x/y frame in the design's column order (`z, sf, period`),
     since the formula path puts sparse columns last.
   - `attr(fit$fit$data@x, "factor.levels")` holds the training levels and
     `data@varTypes` ends in the period's ordinal code.
   - `survivalProbabilities(newdata = )` on the training frame (response and
     extra columns included) equals the exact-columns frame and
     `newdata = NULL`; the same for the sparse fit. `predict` on the
     hand-expanded frame equals the stored training fit; `getTrees(newdata =
     )` and `extract` run on a frame.
   - `test =` with a factor and with a sparse column at fit time: stored test
     draws equal `predict` on the same rows.
   - `factors = "indicators"` predicts from a frame; an unseen level and an
     NA factor behave as on an aft fit; `period` as a predictor is refused on
     the formula, x/y frame and x/y matrix paths.
   - Numeric-only neutrality is pinned already: test-hazard.R holds the
     formula fit to the matrix fit bitwise, and
     [hazard-reduction.R](../../benchmarks/R/hazard-reduction.R) the matrix
     fit to the reduction.
5. Docs: one sentence in the hazard paragraph of man/bart.Rd (predictors
   enter as for every family; a predictor named `period` is refused); in
   docs/design/survival.md's reduction-gate section, that a formula fit's
   reduction target is the hand-expanded frame through the same formula.
   NEWS: none. Hazard is new in 1.0-0 and this defect never reached main
   (dec-B128); the existing item's claim of a byte-identical reduction
   becomes true for factor designs. At landing, remove the TODO item.

## RNG class

Neutral. With the prototype, `yhat.train` was bitwise unchanged for
numeric-only formula fits, `g + z`, ordered factors, `sparseFactor`,
transformed terms, the x/y frame and matrix paths, `test =` and
`factors = "indicators"`. Draws move only for a formula fit whose factor
declares a level no kept subject carries, and they move to the by-hand
reduction's draws, which is the oracle. No recorded baseline or gate uses
hazard with a factor or a formula: hazard-exact.R, hazard-reduction.R, the
equivalence harness's hazard scenario, composition-matrix.R and
constant-person-period-rows.R all pass numeric matrices through x/y, and no
reproducibility snapshot file fits a hazard model.

## Verification

- Before the change, `saveRDS` the `yhat.train` of a numeric-only formula
  fit and a `g + z` fit; after it, both must be `identical`.
- `R CMD INSTALL -l <lib> .`, then
  `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`, 0 fail.
- `cd tests/cpp && make && ./test_bartcore` (no C++ change; the neutral
  class's gate).
- The lint chain in CLAUDE.local.md, then `R CMD check --as-cran
  --no-manual` on a tarball built outside the tree.
- The fit's `data@x` changes class, so run the exact gates in
  .github/workflows/exact-gates.yaml with `quick`; hazard-exact.R and
  hazard-reduction.R must pass unchanged.

## Agent-made calls

- The report's "split as if ordered" is recorded as not reproduced; the
  training-side defect is the lost declared level count.
- A predictor named `period` is refused rather than renamed or left
  shadowed.
- `newdata = NULL` reconstruction keeps the container (Step 3) rather than
  densifying n x K x p.
- The sparse by-hand target is the x/y frame in design order; the
  sparse-last order of formula designs is left alone.
- Class neutral, not posterior-changing: the unused-level case moves toward
  the documented reduction and no gate or baseline covers it, so there is no
  design note and no re-record.
- No NEWS item (dec-B128). A new test file rather than growing
  test-hazard.R. Sonnet tier (R only).
- The three Constraints exclusions stay out of this item.

## Landing note

Landed as planned, R only. expandDiscreteTimeHazard row-subsets through a new
hazardRowSubset (container kept, matrix builder attributes kept) and refuses a
predictor named period on every path; appendHazardPeriodColumn extends a
container and each builder attribute; the formula path's test expansion uses
the same subset; hazardSurvivalProbabilities with newdata = NULL takes the
period-1 rows of the stored design, replicates them K times and overwrites the
period column. man/bart.Rd and the reduction-gate section of
docs/design/survival.md carry one sentence each; the TODO item is removed.

RNG neutrality was shown against the tip before the change (not pinned in a test, since draws are not bitwise across hosts): the yhat.train sums of a
numeric-only, a g + z and a sf + z formula fit were identical before and after; the g + z fit also equals the x/y frame
fit and the hand-expanded probit fit, and a factor with an unused declared
level now equals the hand-expanded fit. Extra tests beyond the plan: an NA
factor in newdata errors, and under na.pass returns NA.

Gates (R 4.6.1, Darwin 25.6.0 arm64, private library): full tinytest 9657
pass, 0 fail; lintr clean; air format clean; check-rc-codoc, check-win-drift,
check-doc-freshness OK; all 25 exact-gates.yaml gates PASS with quick
(hazard-exact.R and hazard-reduction.R unchanged); R CMD check --as-cran
--no-manual on a tarball: Status 1 NOTE, the stale Date field, none new.

Lines: R 79 added / 21 removed, tests 183, man and design 5, plan note ~30.
The tests ran over the ~130 line estimate for the added cases.

Review fix: the formula path also refuses any term whose variables include
period (transformed, or a factor under factors = "indicators"); the three
pinned sums were dropped from the tests.
