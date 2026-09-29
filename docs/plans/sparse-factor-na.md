# sparse-factor-na: a sparseFactor holds missing values

Status: IMPLEMENTED 2026-09-29, pending review (step 2: option B, ruled by the maintainer)

agent: sonnet
rng: neutral (no C or C++ change; every input accepted today codes and routes as before, and the
inputs this adds were refused)
budget: ~400 lines under option A, ~460 under option B (R ~120: sparseFactor.R ~75, A_class.R ~5,
NAMESPACE 1, data.R ~12, mixedMatrix.R ~12, utility.R ~15; tests ~230; man and NEWS ~45)

## Goal

A `sparseFactor` holds a missing value the way a factor and a Matrix sparse vector do: as an explicit
stored entry whose value is `NA_integer_`, never as the reference level. Every method the class has
answers an NA as a base factor does, and a fit reads it as a missing predictor, bitwise as it reads a
dense factor's NA in the same rows (maintainer ruling on dec-A100, 2026-09-29).

## Context

- The class: slots and validity in [`sparseFactor`](../../R/A_class.R); constructor and methods in
  [`sparseFactor`](../../R/sparseFactor.R). Today the constructor, the validity check, `[`
  ([`sparseFactorPositions`](../../R/sparseFactor.R)), `[<-` and `levels<-` refuse NA, and `is.na`
  always answers FALSE. The data-frame methods came in
  [sparse-factor-frames.md](sparse-factor-frames.md).
- Engine hand-off: [`sparseColumnSlices`](../../R/mixedMatrix.R) passes `values - 1` as doubles, so a
  stored NA arrives as a NaN code, which the engine already reads as a missing categorical value. The
  maintainer checked that a fit with a fixed `sigest` is bitwise equal to the dense-factor fit with NA
  in the same rows. No engine, bridge or header change is needed.
- Consumers that already read NA correctly, once the class can hold it: the x/y interface's
  [`applyNaActionToXY`](../../R/data.R) and [`rowsWithMissingPredictors`](../../R/data.R), which run
  on the assembled container; the all-missing column refusal (`sparseAllMissingCheck` in
  [`dbartsData`](../../R/data.R)); the test-set refusal [`refuseTestMissingness`](../../R/data.R)
  and predict's [`resolvePredictRows`](../../R/data.R), both reading the container's stored entries
  through [`unroutableTestColumns`](../../R/data.R) and [`testRowsMissingIn`](../../R/data.R);
  [`subsetSparseColumn`](../../R/mixedMatrix.R) under `subset`; and the sampler's `setPredictor`
  and `setTestPredictor`, which take numeric codes (NaN already means missing there) and check the
  missing policy through [`sourceAnyNA`](../../R/utility.R), which reaches the new `anyNA` method
  for a frame argument.
- Consumers that do not, found by probing the installed tip:
  - [`remapSparseFactorToTrainingLevels`](../../R/utility.R) reads an NA stored code as a level
    unseen in training and refuses the test column. It also refuses a reference level unseen in
    training even when every row is stored, so no row takes it.
  - [`alignContainerFactorLevels`](../../R/utility.R), CSC branch: the same two refusals, for a
    NaN code and for an untaken reference.
  - The formula path in [`dbartsData`](../../R/data.R) lifts sparse columns out
    ([`pullOutSparseFormulaColumns`](../../R/mixedMatrix.R)) before `model.frame` applies
    `na.action`, so the action never sees their missing rows. Probed with a `sparseVector` column
    holding an NA: `na.omit` keeps the row and `na.fail` passes, while the x/y interface drops it.
    This defect already exists for sparse ordinal columns; the fix here covers every lifted column.
  - [`rbart_vi`](../../R/rbart.R) hard-codes `na.omit`. After the fix, a sparse column holding an NA
    drops that row, but the call stops earlier, at the "sparse categorical predictors require
    factors = categorical" refusal, since `rbart_vi` does not take a sparse factor. No change there.
- Base R limits that no method can reach, since `is.atomic` is primitive and the class is not a
  vector: on a data frame holding a `sparseFactor` or a `sparseVector` column, NA or not,
  `complete.cases` and `na.fail` error with "invalid 'type' (unknown) of argument", and
  `stats::na.omit.data.frame` skips the column, so it keeps its NA rows. The fit is unaffected (it
  never hands the S4 column to them). The manual says so and points to `d[!is.na(d$f), ]`,
  replacing its current advice to use `na.omit`.
- `table(sf)` drops declared levels no row takes, where `table(ff)` keeps them: `table` keeps the
  declared levels only for `is.factor`. That is already true on the tip and is not changed here.

## Constraints

- R only. No change under src/, inst/include/ or tests/cpp.
- Rows without NA keep their current representation, so every current fit, test and snapshot is
  unchanged.
- Nothing from base or stats is masked (the frames slice's no-masking test keeps passing).
- R floor 4.2 holds; `rbind`, `match` and `%in%` on the class stay gated on R 4.6.0 (dec-A117).
- Out of scope: `complete.cases`, `na.fail` and `na.omit` on a data frame (above); an NA level
  (`addNA`); `rank`; `relevel`, `rep_len`, `rep.int` and the rest of the frames slice's
  unsupported list.

## Steps

1. Class and constructor (R/A_class.R, R/sparseFactor.R).
   - Validity: `values` may be `NA_integer_`; a non-NA value must be a code in range. `levels` still
     refuses NA and duplicates; `reference` is still a non-NA element of `levels` (unless step 2
     takes option B).
   - [`sparseFactor`](../../R/sparseFactor.R): an NA in `x` (factor, character or integer codes)
     becomes an NA code. A non-NA value absent from `levels` is still refused. With no `levels`
     given, levels come from the non-missing values (`sort(unique(x))` already drops NA). A factor
     `x` that has an NA level (from `addNA`) is refused by name, "'x' has an NA level; a
     sparseFactor stores a missing value, not a missing level", ahead of the validity message.
   - Canonicalization keeps NA entries in both branches: the dense branch selects
     `is.na(codes) | codes != referenceCode` (today's `which(codes != referenceCode)` would turn every
     NA into the reference silently), and the `i` branch's `keep` does the same (today's comparison
     yields NA and indexes garbage rows).
   - Fix the file-header comment that says the formula path refuses the class.
   Budget: ~25 lines.
2. Every row missing (option B, the maintainer's ruling, since a base factor may have no levels). `factor(c(NA, NA))` has no levels, and
   `droplevels` of an all-NA factor leaves none; a `sparseFactor` today needs one level and a
   non-NA reference. The implementer does whichever option is ruled.
   - Option A (not taken): `sparseFactor(x)` with every value NA and no `levels` is refused with
     "'x' has no non-missing values; supply 'levels'"; with `levels`, the reference is `levels[1]`,
     a declared level no row takes. `drop = TRUE` and `droplevels` on an all-NA vector keep the
     reference as the one level, and `levels<-` naming every level NA is refused with "a
     sparseFactor needs at least one level". Cost: ~5 lines in the constructor,
     [`dropSparseFactorLevels`](../../R/sparseFactor.R) and `levels<-`; 3 tests. Each is one
     pinned difference from base.
   - Option B (ruled): allow zero levels with `reference = NA_character_`, only when every row is a
     stored NA (`length(i) == length` and every value NA). Then the all-NA constructor,
     `droplevels` and `levels<-` answer as base does. Cost: ~35 lines over 8 functions in 4
     files. Validity allows that one shape (~8). The constructor builds it (~5).
     [`sparseColumnSlices`](../../R/mixedMatrix.R) and the categorical builder must pass it
     through as K = 0 with an NA reference until the fit's existing all-missing refusal fires (~4).
     The training-level remap and the CSC branch skip the NA reference (~3). `c` takes its
     reference from the first non-NA reference among its arguments, or the first level of the
     union (~5). [`dropSparseFactorLevels`](../../R/sparseFactor.R) and `levels<-` produce the
     shape (~6). `show`, `str` and `summary` print an empty level table (~2). Plus ~25 lines of
     tests. Every other method holds as is, since all rows are stored NA.
3. Reading methods (R/sparseFactor.R, NAMESPACE).
   - `is.na`: TRUE at the stored rows whose value is NA, from the stored entries.
   - New `anyNA` S4 method, `function(x, recursive = FALSE) anyNA(x@values)`, exported in
     `exportMethods` next to `is.na`. This is what `anyNA` on a frame, and so
     [`sourceAnyNA`](../../R/utility.R), reaches.
   - `show`: the stored-entry line adds "(k missing)" when k > 0.
   - `as.character`, `as.integer`, `as.vector`, `xtfrm`, `format`, `summary`, `str`, `unique`,
     `duplicated`, `Ops`, `as.data.frame` and `c` already carry an NA code through
     `sparseFactorCodes`; they change only if a test in step 7 disagrees with base.
   Budget: ~10 lines.
4. Indexing (R/sparseFactor.R, R/mixedMatrix.R).
   - [`sparseFactorPositions`](../../R/sparseFactor.R) returns NA for an NA index and for a positive
     index past the end, as `seq_len(n)[i]` does, instead of refusing. Character, matrix and factor
     indices stay refused (`sf["a"]` is refused where a factor gives NA; there are no names).
   - [`subsetSparseFactorRows`](../../R/mixedMatrix.R) makes an NA position a stored NA entry; today
     `match` would read it as unstored, that is, as the reference.
   - `x[[i]]`: an NA or out-of-range index is refused with base's "subscript out of bounds", since
     `x[i]` now returns an NA element for it.
   - `drop = TRUE` and `droplevels` ([`dropSparseFactorLevels`](../../R/sparseFactor.R)): NA rows
     take no level (`tabulate` already ignores NA); the all-NA case follows step 2.
   Budget: ~15 lines.
5. Assignment (R/sparseFactor.R).
   - `[<-`: an NA value is stored as an NA entry. A non-NA label that is not a level stores NA and
     warns "invalid factor level, NA generated", as `[<-.factor` does. An NA index is dropped when the
     value has length one and refused with "NAs are not allowed in subscripted assignments"
     otherwise. A logical index becomes `which(i)` and the new length is `max(length(x),
     length(i))`, so a logical index longer than the vector extends it, as for a factor. Any index
     past the end fills the gap with NA, so the gap refusal goes.
   - `length<-`: truncates by `x[seq_len(value)]` and pads with stored NA entries, as `length<-` on a
     factor does, replacing the refusal. A data frame holding one still cannot grow by assigning
     past its last row: base's row expansion strips the S4 class before the method is reached, so
     that refusal is base's and stays pinned.
   - `levels<-`: an NA in `value` drops that level and turns its entries NA, as for a factor. If the
     reference's new name is NA, its implicit rows become stored NA entries and the reference moves
     to the most common remaining level (first on ties); no remaining level follows step 2.
     Canonicalization compares with `is.na(v) | v != ref`.
   - `is.na<-` (base default, through `[<-`), `rep` and `c` need no code of their own.
   Budget: ~30 lines.
6. Consumers (R/utility.R, R/mixedMatrix.R, R/data.R).
   - [`remapSparseFactorToTrainingLevels`](../../R/utility.R): check only the non-NA stored codes
     against the training levels and carry NA codes through. Refuse an unseen reference only when a
     row takes it (`length(x@i) < x@length`); otherwise the reference becomes the first training
     level and any stored entry now equal to it is dropped, keeping the canonical form.
   - [`alignContainerFactorLevels`](../../R/utility.R), CSC branch: the same for NaN codes, which
     stay NaN, and for the reference, refused only when the column has fewer stored entries than
     rows; an untaken one is set to code 0 by slot surgery, and stored entries keep their own codes.
   - A new helper in R/mixedMatrix.R beside
     [`pullOutSparseFormulaColumns`](../../R/mixedMatrix.R) gives a lifted column's missing rows
     from its slots, never densifying: `is.na` for a `sparseFactor`, the stored `x` entries for a
     `sparseVector` (its `i` is 1-based; a pattern vector has none), and the stored entries of each
     column for a `dgCMatrix`. [`rowsWithMissingPredictors`](../../R/data.R) is not touched; its
     callers see only containers.
   - [`dbartsData`](../../R/data.R) formula path: when a lifted column the formula uses has a
     missing row, add `dbartsSparseMissing = <numeric, NA on those rows and 0 elsewhere>` to the
     `model.frame` call as a named extra argument, as `weights` rides. It lands as the column
     "(dbartsSparseMissing)", which `modelFrame[termLabels]` never picks up; the name is a prefix of
     no `model.frame.default` formal (formula, data, subset, na.action, drop.unused.levels, xlev) and
     differs from "(weights)" and "(offset)". `subset` then applies to it and the caller's
     `na.action` sees it: `na.omit` and `na.exclude` drop and record the rows, `na.fail` refuses,
     and [`na.keepPredictors`](../../R/data.R) keeps them. Remove the argument from the call right
     after the training frame is evaluated: the call is re-evaluated later against `test`, where a
     training-length extra fails "variable lengths differ" and silently loses `weights.test`.
     Nothing is added when no lifted column has a missing row, so NA-free fits build the same call
     as today.
   Budget: ~40 lines.
7. Tests: new inst/tinytest/test-sparse-factor-na.R, below. Update what the earlier files pin:
   the constructor refusal in [`sparseFactor`](../../inst/tinytest/test-sparse-factor.R), and in
   [`sparseFactor`](../../inst/tinytest/test-sparse-factor-frames.R) the refusals for `[`, `[<-`
   past a gap, an invalid label in a row assignment, `levels<-` and `length<-`, each to its new
   answer. The frame-growth test (base's "cannot set length") stays. Budget: ~230 lines.
8. Docs.
   - [`sparseFactor`](../../man/sparseFactor.Rd): `x` accepts NA; the details paragraph that lists
     NA refusals becomes one stating NA is an explicit stored entry and every method answers it as a
     factor does. It keeps the sentence that a frame holding one cannot grow by assignment past its
     last row. The unsupported list says `complete.cases` and `na.fail` error on a frame holding one
     and `na.omit` keeps its NA rows, pointing to `d[!is.na(d$f), ]`. It adds `rank`, an NA level,
     and `table` dropping unused declared levels (use `table(factor(ff))`'s reading, or
     `levels(sf)`). The `reference` advice says NA rows are always stored, so the most common
     non-missing level is still the best reference. Add the alias `anyNA,sparseFactor-method`.
   - inst/NEWS.Rd: fold into the existing sparseFactor item (new in 1.0-0, so no separate item), and
     make its `complete.cases` exception match the Rd. The formula-path `na.action` fix gets no NEWS
     item: lifted sparse columns never reached main.
   Budget: ~45 lines.
9. Ledger. docs/decisions.md belongs to the orchestrator: do not edit it. The landing report lists
   the calls made below, with the step 2 ruling, for the ledger.

## Tests

In inst/tinytest/test-sparse-factor-na.R, over `ff <- factor(c("a", NA, "b", "a", NA, "c", "a"),
levels = c("a", "b", "c", "d"))` and `sf <- sparseFactor(ff)`. Every comparison below is against
the same call on `ff`, so each test pins base R's own answer, not a transcription of it. Count
warnings with `withCallingHandlers`.

- Storage: NA rows are stored entries with value `NA_integer_`; the reference is "a"; the `i`
  constructor with an NA entry and an integer-code `x` with NA both store it; a non-NA value absent
  from `levels` is still refused; `sparseFactor(addNA(ff))` is refused naming the NA level;
  validity still refuses an out-of-range code. The all-NA cases follow the step 2 ruling.
- Reading: `is.na`, `anyNA` (TRUE, and FALSE for an NA-free sparseFactor), `as.character`,
  `as.integer`, `as.vector`, `xtfrm`, `order`, `sort`, `format`, `summary`, `str` output, `unique`,
  `duplicated` and `factor(sf)` equal the factor's; `table(sf, useNA = "ifany")` equals
  `table(factor(ff), useNA = "ifany")`; on R 4.6.0 or later, `match(sf, "a")` and `sf %in% "a"`
  equal the factor's; `show` names the missing count.
- Indexing: `sf[c(1, NA, 9)]` gives a, NA, NA; `sf[c(TRUE, NA)]` recycles with NA rows; `sf[-10]`
  returns everything; `sf["a"]` stays refused; `sf[[NA]]` and `sf[[9]]` refuse "subscript out of
  bounds"; `sf[, drop = TRUE]` and `droplevels(sf)` keep the used levels.
- Assignment: `sf[2] <- "b"`, `sf[1] <- NA`, `sf[10] <- "b"` (rows 8 and 9 NA),
  `sf[c(rep(FALSE, 8), TRUE)] <- "b"` (length 9, row 8 NA), `sf[c(1, NA)] <- "b"` (row 1 only),
  `is.na(sf) <- 3` and `sf[NA] <- "a"` (no change) each equal the factor's;
  `sf[c(1, NA)] <- c("a", "b")` and `sf[1] <- "z"` give base's error and base's one warning;
  `length(sf) <- 9` pads NA and `length(sf) <- 2` truncates; `levels(sf) <- c("a", NA, "c", "d")`
  turns the b entries NA and drops "b"; renaming the reference to NA leaves the same labels as the
  factor and a non-NA reference; `rep(sf, 2)` and `c(sf, factor(c(NA, "z")))` equal the factor's
  labels and levels.
- Ops: `sf == "a"`, `sf != "a"`, `sf == NA` and `sf == sf` equal the factor's.
- Frames: `print` of a frame holding `sf` has the same output as one holding `ff`, with no
  warning; `head`, `d[is.na(d$f), ]`, `d[c(1, NA), ]` and `is.na(d)` match; `anyNA(d)` is TRUE;
  `na.omit(d)` keeps the NA rows and `complete.cases(d)` and `na.fail(d)` error (the documented
  base limits, pinned so a change surfaces); on R 4.6.0 or later `rbind` gives a factor with NA in
  the same rows.
- Fits: with `sigest = 1`, a fixed seed and `n.threads = 1L`, the x/y fit and the formula fit on a
  frame holding `sf` are identical to the ones on `ff` (training fits and `sigma`), and so is
  `predict` on a test frame with NA in that column. A test `sparseFactor` whose level order differs
  from training and that holds NA predicts identically to the dense test frame. A test
  `sparseFactor` with every row stored and a reference unseen in training predicts, and one whose
  unseen reference some row takes is refused. An NA in a test column that had none in training is
  refused with the column's name, and `na.action = na.pass` answers NA for those rows, as for the
  dense factor. An all-NA column is refused as "predictor columns cannot be entirely missing".
- Formula `na.action`: `na.omit`, `na.exclude` and `na.fail` on a frame holding `sf` keep the same
  rows, record the same `na.action` and refuse the same way as on `ff`; under `na.exclude` the
  fitted values pad back to the frame's row count; the default keeps the rows; `subset` combined
  with `na.omit` keeps the same rows as the dense fit. The same `na.omit` and `na.fail` checks on a
  `sparseVector` column and on a `dgCMatrix` column with an NA, against their dense numeric twins.
  A formula fit with `test` shorter than training, `weights` and `weights.test`, and a sparse
  column holding an NA keeps `weights.test` (equal to the dense fit's).

## Verification

- `R CMD INSTALL .`, then `tinytest::test_package("dbarts")`: 0 failures, the count up by the new
  file's.
- Mutation checks, each reinstalled and then restored with a `touch`: revert the dense-branch
  canonicalization to `which(codes != referenceCode)` and see the storage and fit-identity tests
  fail; remove the `dbartsSparseMissing` argument and see the formula `na.action` tests fail; leave
  it on the call past the training frame and see the `weights.test` test fail.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, and after installing
  `Rscript tools/check-rc-codoc.R .`, `Rscript tools/check-win-drift.R .`,
  `Rscript tools/check-doc-freshness.R .`, each on its own exit status.
- NEWS parses: `tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")` is non-NULL.
- `R CMD check --as-cran --no-manual` on a tarball built from a clean copy outside the tree.
- No equivalence re-record or exact gate: no C change, and no NA-free input moves.

## Calls made

For the maintainer's ledger; each follows base R where base R has an answer. Step 2 (every
row missing) is call 13.

1. `sparseFactor(x)` stores a dense `x`'s NA as explicit NA entries, as `factor(c("a", NA))` keeps
   the NA.
2. The constructor still refuses a non-NA value absent from `levels`, where `factor(x, levels)`
   quietly makes it NA; `[<-` instead follows `[<-.factor` and stores NA with its warning.
3. `levels<-` renaming the reference to NA stores the implicit rows as NA entries and moves the
   reference to the most common remaining level, the rule `drop = TRUE` already uses.
4. An NA level, as `addNA` makes, stays unsupported: an NA value is not an NA level.
   `sparseFactor(addNA(f))` is refused by a message naming the NA level, not by the validity check.
5. `length<-` now pads with NA and truncates, as for a factor, replacing its refusal; a frame still
   cannot grow past its last row, which base refuses before the method is reached.
6. `x[[i]]` refuses an NA or out-of-range index with base's message.
7. `sf["a"]` stays refused where a factor returns NA, since the class stores no names.
8. `show` reports how many stored entries are missing.
9. The most-common-level advice stands: NA rows are always stored, so the best reference is the
   most common non-missing level.
10. The formula-path `na.action` fix covers every lifted sparse column (`sparseVector`,
    `dgCMatrix`), since the defect is shared, and has no NEWS item since it never reached main.
11. A test column's reference level unseen in training is refused only when some row takes it.
12. `complete.cases` and `na.fail` erroring on a frame holding a sparse column, `na.omit` keeping
    its NA rows, `table` dropping unused declared levels and `rank` failing are documented, not
    fixed: no fix exists without masking base or stats.

13. Step 2, option B (ruled): a vector whose every row is a stored NA may have zero levels and
    reference `NA_character_`, so the all-NA constructor, `droplevels` and `levels<-` answer as base R.
14. `c` takes its reference from the first non-NA reference among its arguments, else the first
    level of the union.
15. `na.omit` and `anyDuplicated` have methods on the bare vector, as for a factor.