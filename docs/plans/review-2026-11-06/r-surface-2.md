# Review 4, slice S-F: the R surface, second wave

Status: built, in review (branch review4-r-surface-2).

Each item is a finding of the fourth whole-repo review. RNG class: neutral. No fix changes a sampler
draw and no C++ file is edited; ingest-02's door is S-B's. Most fixes refuse a call that failed late,
gave a wrong answer or read an argument another way, by name and before any compiled call; the rest
make a door read data the way its sibling doors already do. The fourth series adds tests only.

## Fixes as built

### F1: new data and ingestion

**ingest-01, ingest-V-01: a sparse test column against a dense fit.** `densifyUnheldSparseColumns`,
called at the top of `validateXTest`'s data-frame branch, makes each sparse column of a data-frame test
dense (a sparseFactor by `sparseFactorToFactor`, a sparseVector or a one-column sparse matrix a vector,
a wider one a matrix) unless the training design holds that column sparse; a categorical fit made from
a data frame keeps its sparse columns as given. A dgCMatrix fit handed a sparse frame is refused: "a
data-frame 'test' with sparse columns is read against a fit made from a data frame or a formula; give
'test' as a dgCMatrix, or its columns dense". Tests: test-predict-sparse.R, a formula fit predicting on
a frame with a sparseVector and a sparseFactor column, identical to the dense frame through predict,
`bart(test = )`, `$setTestPredictor` and `$predict`; an x/y fit with a sparseVector column; a
`poly(x3, 2)` term; the sparse-trained arms unchanged.

**docs-04, docs-V-01 (Call 5, option (a)): a data frame given to a matrix fit.** `validateXTest`'s
indicators arm calls `makeIndicatorModelMatrix(x.test, if (!is.null(drop)) drop else FALSE, ...)`, so a
matrix fit, which stores no drop record, keeps a test frame's constant column. NEWS, bug fixes. Tests:
test-predict-training-rows.R, `predict(f, dfx[1, ])` and a pinned-column frame equal to the matrix,
through `bart(test = )` and `$setTestPredictor` too; a bartBT data-frame fit unchanged.

**ingest-06: the column-name rule at one predictor.** `validateXTest` matches named against named at
every predictor count, so a renamed one-column test is refused with "column names of 'test' do not
match"; a fit whose own names repeat warns "'x' has repeated column names, so 'test' is matched by
position". NEWS, upgrading. Tests: test-predict-training-rows.R and test-data-compatibility.R.

**rgen-01: a named vector as newdata** is one row matched by name. Test: test-predict-training-rows.R,
on a matrix and a formula fit.

**group B (ingest-08), list newdata.** A plain list of equal-length columns is read with
`as.data.frame(x.test, optional = TRUE)`; any other list is refused naming the forms taken. Test:
test-predict-training-rows.R.

**ingest-05: one NA row name on the matrix route.** `applyNaActionToXY` builds its stand-in frame from
unnamed columns, so rows are read by position. Test: test-na-action.R, a row-named x with a logical
subset holding NA, on `dbartsData` and `bart`.

**ingest-V-02 (dec-B609): an integer subset holding NA.** `refuseOutOfRangeSubset` checks range on the
non-missing entries alone; an NA index pads a row the na.action drops, padded back under na.exclude, as
lm reads it. Tests: test-na-action.R, `c(1:6, NA, 9:14)` on the formula and matrix routes equals
`i[!is.na(i)]` under one seed, na.exclude pads to 13, an index past the data still refused.

**ingest-07: a character matrix column.** `makeIndicatorModelMatrix` and `makeCategoricalModelMatrix`
refuse it by name: "column 'z' is a matrix of characters; give each of its columns as a column of its
own, or as a factor". Logical and numeric matrix columns splice as before. NEWS, bug fixes. Test:
test-data-input-kinds.R, at `bart()`, `dbartsData()` and predict.

**group B (ingest-08), empty designs.** The formula route refuses a design of no term labels and no
sparse columns (`y ~ 1`, `y ~ offset(o)`, a dot less every column), and the matrix route a matrix or
frame of no columns, with the existing "predictors must be specified for regression tree analysis",
before any compiled call; `extractMultinomialFormulaData` the same. A hurdle fit reads its zero part's
offset rows through a formula naming no predictor, so `ingestHurdleData` marks that formula with a
`readsRowsOnly` attribute the check honors. Tests: test-data-errors.R, each spelling; the hurdle case
in test-hurdle-offset.R, which caught the first build's regression.

**install-doors-noted, the formula-name sentence.** Every wholly backticked term label is read by its
bare column name, in `dbartsData` and `extractMultinomialFormulaData`. TODO loses that sentence. NEWS,
bug fixes. Test: test-formula-terms.R, columns `1x`, `x-y`, `x:y`, `if` and `a b` through `y ~ .` and
written terms, predict equal to fitted and draws equal to a renamed twin's.

**na-level-predictor-doors.** `codeLabelsOnTable(value, table)` codes a factor level to level (its NA
level meets the table's), a character vector by `match(..., incomparables = NA)` (a missing label is a
missing value) and a sparseFactor as characters. `mapFactorColumnsToTrainingLevels`,
`mapFactorColumnsToIndicatorLevels` and `codeCategoricalColumnUpdate` take their codes from it;
`alignContainerFactorLevels` too. The TODO entry is closed. Tests: test-unseen-levels.R, an `addNA()`
fit predicting on its rows equal to fitted, through `bart(test = )` and `$setTestPredictor`; a
setPredictor `NA` label coded missing; a "NA" string a label.

**ingest-V-03: a logical or repeated `column` in the setters.** `resolveColumnIndex` refuses a logical
column ("'column' takes column numbers or names, not TRUE/FALSE") and a column given twice ("'column'
names column 'q' more than once"), for setPredictor, setTestPredictor and setCutPoints, the sampler left
as it was. NEWS, upgrading. Tests: test-sampler-predictors.R; three files that pinned the old reading
now pin the refusal (test-mutate-sparse-valued.R, test-cut-points-undo.R and the hazard period block of
test-leaf-standardization-rows.R, where the repeat is now refused before the period is read).

### F2: the fits' methods, pdbart and plots

**docs-02 = rgen-02 (dec-B610).** fitted, residuals and extract on the bart, bartMultinomial,
bartOrdinal, bartNegbin and bartHurdle classes call `warnUnusedDots`, as predict does; extract's trees
arm passes `getTrees`'s own arguments (whole or abbreviated, `charmatch`) and warns of the rest
(`dropUnusedTreesArguments`). Tests: test-generics-errors.R, each class and method.

**rgen-03: ci.level and the "s" draws.** `predict.bart`'s ci.level arm summarizes s(x) by
`posteriorInterval`, padded as the value is (`padPredictedRows(..., first = first)`). Test:
test-heteroscedastic-channels.R. Help: bart.Rd.

**rgen-04: one-draw multinomial plot** keeps its draws a matrix. Test: test-plot-generics.R.

**rgen-06: aft plot labels**: "log observed time" against "posterior interval for E(log T | x)",
censored rows open and events filled under a two-entry legend. Test: test-plot-generics.R.

**rgen-05 (dec-B613): print counts the fit's stored draws.** `fitSynopsis` reads the kept draws from
the stored draws, and the control's n.samples only where no draws are stored; the keepTrees and
keepSampler help say which methods read the kept sampler and which the stored draws, and that running
the sampler further is not a supported use of the fit. Test: test-plot-generics.R. Help: bart.Rd,
bartBT.Rd.

**multiforest-04: pdbart on a basis written as code.** `pdbart.trainingRows` collects the variables
each forest's basis term reads (`get_all_vars` on its terms) beside the predictors'. Test:
test-pdbart.R, pdbart and pd2bart with and without `newdata = d`, a factor and a numeric basis,
gaussian and probit, identical.

**group A (rfit-05, rgen-08, aux-01, aux-02, rminor-V-01): xind.** `pdbart.checkXind(xind, numColumns,
pair, caller)`, from `pdbart.resolveXind` and `plot.pdbart`: the wrong kind ("'xind' takes names, column
numbers, or the bare right-hand side of a formula"), out of range or fractional ("'xind' must be column
numbers between 1 and 3; got 0"), empty, repeated, and for pd2bart not a pair ("pd2bart's 'xind' names a
pair of variables; got 3"). plot also takes the result's variable names. pd2bart's default stays the
first two predictors. NEWS, upgrading. Help: pdbart.Rd. Test: test-pdbart.R, each shape on a matrix and
a formula fit, and `plot(pd, xind = "b")` equal to `xind = 2` by recorded plot.

**aux-03.** `pdbart.setup` refuses a newdata of no rows ("'newdata' has no rows to average over");
`pdbart.defaultLevs` a levquants not numbers in [0, 1]; the grid a `levs` entry holding NA ("'levs'
holds a missing value for 'a'; give the values to evaluate at"). Test: test-pdbart.R.

**group B (rgen-09).** `checkCombineChains` on predict and extract of every fit class and on
`survivalProbabilities.bart` ("'combineChains' must be TRUE/FALSE"); a type that is NA, empty or not a
string gets the type list (`foldTypeAliases`, `resolveHurdleType`, rbart's folds); `summary.bart`'s
vars must be field names ("'vars' must name the fields to summarize"). Tests: test-generics-errors.R,
test-rbart-error.R, test-pdbart.R.

**rfit-08 (dec-B614): sigma() and nobs().** `sigma.bart` is `mean(extract(object, type = "sigma"))`
(1 on a family with no residual scale), refused on a variance-forest fit ("a fit with a variance forest
has no single residual sd; see extract(fit, type = "sigma")"); `nobs.bart` counts the rows of positive
weight and not masked off. Both registered for the five classes. NEWS, new features. Help: bart.Rd.
Test: test-generics-correctValues.R.

### F3: argument checks and messages

**multiforest-01.** `refuseConstantBasis` refuses a basis 0 on every row in the likelihood ("<forest>'s
basis is 0 on every row in the likelihood, so it contributes nothing to a forest; drop the forest")
before offering a basis-free spelling, and builds that spelling with `paste(c(written, sd,
"amplitude = normal()"), collapse = ", ")`. Tests: test-forest-scale-record.R.

**rfit-01 (dec-B612).** `resolveRenamedSeed` refuses `rngSeed` and `seed` both given ("'rngSeed' and
'seed' name the same seed on 'bart'; supply 'seed' alone"), naming bart, bart2 or dbartsControl. NEWS,
upgrading. Test: test-tombstones.R.

**rfit-02.** `refuseMixedArgumentStyles` sends a name that is a prefix of no formal of bart or bartBT
to `refuseForeignFrontDoorArgs`, the ordinary unused-argument refusal. Test: test-front-door.R.

**rfit-06, rfit-07.** `resolveFamily` takes a vector only when it is the door's own default
(`identical(value, tokens)`), else "'family' must name one response family; got "gaussian",
"probit""; xbart reads method and loss as match.arg does (`matchOneToken`); a random subsample's
`n.test` of 1.5 or NA, a k-fold NA `n.test` and `drop = NA` are refused in the package's words. Tests:
test-family.R, test-xbart-error.R (one pin of a raw R error updated).

**rfit-09.** `setControl` gives a fresh control's NA n.samples the sampler's own. Test:
test-sampler-bridge-errors.R.

**rfit-11.** bart drops `factors` from the stored call where the caller did not write it. Test:
test-bart-formula.R, `update(fit)` draws the same under one seed.

**group B, the rest (rfit-03, rfit-04).** `checkFlag(value, name)` in bart for keepCall, samplerOnly
and keepSampler (the last only when supplied); `resolveSigestArg` refuses a non-finite number; bartBT
checks verbose, binaryOffset (numeric, one value or one per row on the matrix route), x.test's column
count and ndpost against keepevery; `resolveSeedArg` refuses a seed of other than one value at every
door; hazard refuses empty breaks; the n.cuts reader refuses half-named values; pdbart and pd2bart
refuse a `pl` that is not one TRUE or FALSE; `resolveProposalProbs` refuses a list entry of another
length. Tests: test-error-quality.R, test-cut-counts-by-name.R, test-pdbart.R.

**docs-10, the R strings.** The `refuseZeroSamples` message and the multinomial refusal ("fields from
the temporary sampler it builds"), and the code comments in R/bart.R, R/utility.R,
R/validateComposition.R and R/family.R, use the standing phrase; the machine meaning beside the BLAS
note is left. TODO's host-wording-sweep keeps only dbarts.h, S-B's share.

**docs-01, option (c) (dec-B611).** The `dbartsData` validity takes only "incorporate" in the
`missing` slot, whose comment says it holds the one mode; `checkMissingPolicy`, `sourceAnyNA`, the
`setData` check and the three setter calls are removed, with inst/common/strictData.R and the tests
that wrote the slot by hand. The mixture vignette says `na.action = na.fail` refuses any missing value,
the response's too. Mutant m20 of the mutation battery, which targeted `checkMissingPolicy`, is retired
with its site. No NEWS: 0.9-34 had no such slot.

**anchor-02.** `refuseLegacyPositionalCall(suppliedCall, callingEnv)` matches the call against
`function(...) NULL` in bart's calling frame before counting, so a wrapper's forwarded dots are counted
by what they hold. Tests: test-front-door.R, the six cases: empty and named dots and a wrapper of
named formals run; four direct positionals, an unnamed fourth through the dots and through two
wrappers refused. The three run cases fail on the build before the fix.

### F4: tests for the mutation lens

| finding | mutant | test | kill shown |
|---|---|---|---|
| mutation-01 | the active-row copy skips the last row | test-leaf-standardization-rows.R, test-forest-basis-scale-update.R: a mask whose only zero is the last row, a far value on it, re-derived standardization (gaussian, probit) and basis scale (probit, logistic) equal to the rows in; cleared, every row | engine mutant, not rebuilt |
| mutation-V-01 | the scan skips the first row | the same tests, the first row | engine mutant, not rebuilt |
| mutation-03 | `extractForest` reads `fits[,, j]` or `bases[[j]]` | test-forest-selection.R: each selection, `c(3, 1)` among them, identical to the full contribution's slices | in process: each mutant 4 failures of 306; 0 on the file before the tests |
| E6 | the zero-multiplier bound over every row | test-bcf-zero-multiplier.R: a weight-0 row with basis 1e9, the fit on the other rows within rmse 0.15 | engine mutant, not rebuilt |
| E32, E42 | the linear-leaf cache kept after setWeights, or after a Student-t sweep | test-linear-leaves.R: live and copy draw alike after a run and after setWeights, gaussian and Student-t | engine mutants, not rebuilt |
| E41 | Student-t setWeights does not mark the rows in | test-update-scale-rows.R: setWeights then `setResponse(updateScale = TRUE)` equals a sampler created with both | engine mutant, not rebuilt |
| E44 | logistic setOffset does not rebuild its working response | test-family-offset.R: sums of squared residuals after setOffset equal a copy's | engine mutant, not rebuilt |

The engine twins read the verifier's tests with the variables renamed for their files; E41's and E44's
data are drawn in a different order, which neither check depends on.

## Calls made

1. ingest-01: densify where the fit cannot hold the column sparse; a categorical fit made from a data
   frame keeps its sparse test columns. Alternative: densify on every fit (loses the categorical fit's
   resident path).
2. List newdata: the refusal names 'test', the function's own word. Alternative: 'newdata'.
3. ingest-V-03: the repeated-column message names the column as written ('q' or 2). Alternative:
   always the position.
4. na-level-predictor-doors: `alignContainerFactorLevels` also takes its codes from
   `codeLabelsOnTable`, beyond the plan's list. Alternative: leave the container path.
5. install-doors-noted: a NEWS line for non-syntactic names, since 0.9-34 fit them under other names.
   Alternative: none.
6. docs-02: extract's trees arm matches getTrees's arguments by `charmatch`, so an abbreviation
   passes. Alternative: whole names only.
7. rgen-03: the "s" band is padded with the value's `first` rows. Alternative: unpadded.
8. rgen-06: open symbols for censored rows, filled for events, under a legend.
9. rgen-05: where no draws are stored, print falls back to the control's n.samples.
10. multiforest-04: the basis variables come from `get_all_vars` on each basis term.
11. Empty designs: the hurdle reader's formula carries a `readsRowsOnly` marker. Alternative: a
    separate reader that builds no design.
12. group A: pd2bart's default is the first two predictors; an empty xind and plot's own words; an
    NA or "" string is not parsed as a formula; a NEWS line; help.
13. aux-03: the levs check sits in the grid builder, where each entry's label is known.
14. rgen-09: summary's vars message says fields, not positions, since vars names fields; rbart's three
    inline folds use `foldTypeAliases`; `survivalProbabilities` checks combineChains too.
15. rfit-08: on Student-t, sigma is the t scale; nobs reads weights and the active mask where their
    lengths match the response; rbart fits get no method; `warnUnusedDots` on both.
16. multiforest-01: the refusal says "in the likelihood" only where rows were left out.
17. rfit-01: bart2's refusal names bart2.
18. rfit-02: a stranger is a name that is a prefix of no formal of either door.
19. rfit-06: a vector is the door's default only when identical to it, match.arg's rule.
20. rfit-07: NA n.test refused in both arms; the (0, 1) wording kept for the fraction; `drop = NA`.
21. group B: keepSampler checked only when supplied (its default reads keepTrees); bartBT's verbose
    may be numeric; binaryOffset's length checked on the matrix route only; the seed length check
    in `resolveSeedArg` serves every door; `pl` may be numeric; a proposal.probs list is unlisted
    after its entry check.
22. Pinned tests: test-xbart-error.R, test-forest-scale-record.R, test-cut-points-undo.R,
    test-cut-counts-by-name.R, test-mutate-sparse-valued.R, test-leaf-standardization-rows.R,
    test-bartcore.R (now an NA thread count for its bridge refusal) and test-predict-na-action.R
    pinned the old behavior and now pin the new.
23. docs-10: R/family.R's hazard comment says "machine", the meaning it had; R/bart.R's comments say
    "temporary sampler" and "sampler-build".
24. docs-01: `sourceAnyNA` goes with its only caller; the mixed-container tests of it are dropped.
25. docs-01: mutant m20 retired with a note, as m19 was. Alternative: delete the entry.
26. anchor-02: the run cases pass `samplerOnly = TRUE`, so no chain runs; the `...` symbol still sits
    after three positionals.
27. F4: the twins sit in the files owning their sites: E6 with the zero-multiplier tests, E32 and E42
    with the linear leaves, E41 with the scale-update rows, E44 with the binary offsets.
28. TODO's host-wording-sweep is narrowed, not closed: dbarts.h is S-B's.
29. The plans index's file count is corrected to the top-level plans it counts.
30. docs/plans/sparse-factor-na.md, a landed plan, keeps `sourceAnyNA` as a name, unlinked, said to
    be removed later. Alternative: leave the stale link (fails doc freshness).

## Gates

Run on arm64 macOS, at most two cores.

- Every new and touched test file passes on the slice's shipped build.
- Full tinytest suite, one process, shipped build: 41840 tests, 0 failed. A first run, 41842 tests,
  failed 3, all pins of behavior this slice changes (a bartcore setControl test of an NA n.samples,
  rfit-09; two predict tests expecting a position warning for a named one-row vector, rgen-01); they
  now pin the new behavior.
- Equivalence trio, reference build, `--bitwise --strict-coverage`, against the current baselines:
  gaussian 55 of 55 scenarios identical draws, BCF 19 of 19, multinomial 11 of 11; none skipped.
- The four seeded-drift snapshot files, reference build: 27 tests, 0 failed, unchanged.
- `R CMD check --as-cran` on a tarball built from a clean copy: 1 NOTE, the incoming-feasibility one
  (the Date field, and an expired certificate on an external URL in bartBT.Rd and pdbart.Rd); the
  tests pass.
- Exact gates in `quick`, the main job's list, shipped build: 44 of 44 passed (33 gates and 11 arms),
  run since a fit's stored call changes (rfit-11).
- `lintr::lint_package()` against the slice's build: no lints. `air format --check .`: clean.
  check-rc-codoc: 53 methods. check-win-drift: OK. check-doc-freshness: OK, after a landed plan's cite
  of the removed `sourceAnyNA` was unlinked. `mutation-battery.R verify-anchors`: every anchor
  resolves. inst/NEWS.Rd parses (252 entries).
- R floor: the added code calls nothing newer than R 4.2.0 (`startsWith`, `lengths`, `anyNA` and
  `match.call`'s `envir` are the newest, all before 3.4).
- Consumers (bartCause, stan4bart, treatSens): no call any fix refuses or newly warns of. treatSens's
  integer setTestPredictor columns and `pdbart(X, Z, levs = , pl = FALSE)` run as before;
  stan4bart renames rngSeed before dbarts sees it; bartCause passes `combineChains` as a logical,
  filters its arguments by bart's and xbart's formals, and forwards only its caller's dots to predict.

## Not verified

- The engine mutants of F4 were not rebuilt; their twins are the verifier's tests, run here on the
  pristine build only.
- Windows and threaded paths beyond the suite's own were not exercised.
