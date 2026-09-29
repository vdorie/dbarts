# sparse-factor-frames: sparse columns behave like data-frame columns

agent: sonnet
rng: neutral (no draw moves; only which rows and columns reach the engine)
budget: ~250 lines (R ~110, tests ~120, man and NEWS ~20)

Status: LANDED-pending 2026-09-28

## Goal

A formula fit with a sparse column (a `sparseFactor`, `Matrix::sparseVector`
or `dgCMatrix` column) predicts from the same kinds of data frame a dense fit
does, and a data frame holding a `sparseFactor` column can be row-subset and
printed (maintainer, 2026-09-28: "Yeah, fix both"; loose ends of dec-A26 under
dec-B100 and dec-B132).

## Context

- Test frames: `validateXTest`'s sparse branch skips the model-frame replay
  the dense branch runs, so it takes `x.test` as the model's columns
  verbatim. `predict(fit, newdata = d)` on the training frame fails with
  "number of columns in 'test' must be equal to that of 'x'" because `d`
  carries the response; `test = d` at fit time fails the same way; a
  transformed term (`y ~ sf + log(z)`) is refused as "'log(z)' present in
  'x' but not in 'test'". The dense-factor equivalents all succeed. Training
  already lifts the S4 columns out of `data`, replays the model frame on the
  dense rest, and re-attaches them (`dbartsData`'s formula path); the test
  path should reuse that, then code the sparse columns over the training
  levels as now.
- Frames: `sparseFactor` has only `show` and `length` methods, so `head(d)`,
  `d[1:5, ]` and `d[d$y > 2, ]` stop with "object of type 'S4' is not
  subsettable" and `print(d)` warns "corrupt data frame". The dgCMatrix
  column case was audited earlier; this plan covers `sparseFactor`, and
  checks the `sparseVector` column case the same way.

## Steps

1. Test path: for a test frame carrying sparse columns, lift them out, run
   the dense remainder through the same term replay (and missing-variable
   message) the dense branch uses, re-attach the sparse columns by name,
   drop any column the model does not use (the response, extras), then code
   as now. Same for `test =` at fit time and `getTrees(newdata = )`,
   `survivalProbabilities(newdata = )` if they share the path.
2. `sparseFactor` methods: `[` for a row index (integer, negative, logical,
   character is not meaningful) returning a `sparseFactor` over the same
   levels and reference; `format` and `as.character` returning the level
   labels, so `print`/`head`/`str` work; `as.factor`/`as.vector` if the
   printing path needs them. Keep the stored representation sparse: a row
   subset maps the stored positions, never densifies.
3. If `sparseVector` columns fail the same frame operations, fix or
   document, whichever is smaller; report which.
4. Tests (inst/tinytest/test-sparse-factor-frames.R): predict on the
   training frame, on a frame with extra columns, with a transformed dense
   term, and `test =` at fit time, each equal to the dense-factor fit's
   column handling (same rows, same values where the fits share trees is not
   required - compare against predicting on the exact-columns frame from the
   same fit); `head`, row subset by index, negative index, logical filter,
   `print` without warning, and a fit on a subset frame equal to a fit on
   the same rows built directly.
5. man/sparseFactor.Rd: the new methods and aliases; NEWS: fold into the
   existing sparseFactor item rather than a new one (the feature is new in
   1.0-0, so no change-against-0.9-x item).

## Verification

Full tinytest; lint gates per CLAUDE.local.md; R CMD check --as-cran
--no-manual on the built tarball; equivalence unaffected (neutral).

## Agent-made calls

Character indexing of a `sparseFactor` is refused (no names are stored).

## Landing note

Landed: validateXTest's sparse branch lifts the sparse columns out, replays
the model terms on the dense remainder through the same replay and
missing-variable message as the dense branch, re-attaches the sparse columns
the model uses by name and drops the rest; predict, `test =` at fit time,
getTrees(newdata = ) and aft survivalProbabilities(newdata = ) share it. Hazard
formula fits densify the design with no level table and refuse the frame
there; that older defect is separate and not addressed here.
`sparseFactor` gained `[` (positive, negative, zero, logical, repeated
indices; character, matrix and out-of-range refused), `as.character`,
`format` and `as.data.frame`, so `data.frame(sf = x)`, head, row subsets,
print and str work. Row subsetting reuses the existing stored-position
mapping, so nothing densifies.

Step 3: a `Matrix::sparseVector` column subsets through Matrix's own method
but prints as a placeholder with a "corrupt data frame" warning, because
Matrix supplies no format method; documented in man/sparseFactor.Rd rather
than registering a method on Matrix's class.

Beyond the plan: `as.data.frame.sparseFactor` (base R's factor has one, and
data.frame(sf = x) failed without it). A frame's row names survive a subset
as base R's do, so a subset fit's yhat.train carries the original row names.

Gates (macOS arm64, private library): tinytest 9431 pass, 0 fail; lintr zero
lints; air format --check clean; check-rc-codoc, check-win-drift,
check-doc-freshness OK; R CMD check --as-cran --no-manual on the built
tarball: Status 1 NOTE (the Date field is over a month old, not new).
Size: R 129 lines added, tests 83, man/NEWS 16 (about 230 total).

Follow-up (review fixes; gates rerun: tinytest 9488 pass 0 fail, lint and air clean, R CMD check --as-cran 1 NOTE, the Date field): as.data.frame delegates to as.data.frame.vector;
`[` honors drop = TRUE; added `[<-`, the double-bracket element read, c, as.factor, as.vector, levels,
is.na, xtfrm, unique, duplicated, Ops (== and !=), summary, str and a refused
`length<-`, so rbind, split-and-rbind, assignment, order, table and na.omit
work on a frame. complete.cases (C code) stays unsupported and is listed in
the Rd with dgCMatrix and sparseVector print limits. Tests cover aft
survivalProbabilities, the sampler's getTrees, na.exclude and a `.` formula.
