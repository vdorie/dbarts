# predict-na-action (with observation row names)

agent: sonnet implementer per slice (R, Rd, tests); opus diff review
rng: neutral (all three slices)
window: before the merge (TODO predict-na-action; docs/decisions.md dec-B34 as of 31f3bf53)
budget: A1 about +550 plus 20-80 existing test sites and bartCause edits;
A2 about +150; B about +500

## Goal

A1: every observation-indexed output carries the data's row names on its
observation margin, on every fit class and both channels. That covers
predict, extract, fitted, residuals and survivalProbabilities, including
the stored draws arrays those read.

A2: fit-time `na.exclude` padding of fitted and residuals, which today
exists only on class bart, reaches multinomial, ordinal, negbin and
hurdle. This follows VD's consistency ruling; the orchestrator made the
call to do it now.

B: every predict method and survivalProbabilities take `na.action` as
dec-B34 rules; the sampler's `$predict`, `$predictForests` and
`$getTrees(newdata)` keep the default. Order is A1, then A2, then B, run
serially: B's padding and dropping carry A1's names, and A2 shares A1's
fitted and residuals sites.

## Context: today

Missingness. [`na.keepPredictors`](../../R/data.R) keeps predictor-NA
rows, and routes are learned only on columns with training NAs
([MIA missingness](../design/mia-missingness.md#mia-missingness)).
[`bartBT()`](../../R/bart.R) forces `na.omit`.
[`refuseHurdlePositiveMissingness`](../../R/bart.R) gives both hurdle
parts one routable set. The only predict-time refusal is
[`refuseTestMissingness`](../../R/data.R), the tail of
[`validateXTest`](../../R/data.R), which every entrance reaches.
[`predict.bart`](../../R/generics.R) reaches it through the sampler, and
multinomial, ordinal and negbin call it directly.
[`predict.bartHurdle`](../../R/generics.R) goes through predict.bart via
[`hurdleParts`](../../R/generics.R). bartBT fits are class "bart". A
`na.action` passed today is warned and ignored.

0.9-34, checked: dec-B34's account holds. Test data frames went through
`model.frame` under the option default `na.omit`
([R/data.R:38](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/data.R#L38)),
and a matrix NaN binned to cut 0
([src/dbarts/bartFit.cpp:3067-3074](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/src/dbarts/bartFit.cpp#L3067-L3074)).

Row names today:
- No output names anything.
- The training x keeps row names only for a bare double matrix or a
  dgCMatrix. A data frame on the x path, every formula fit and any
  factor design become a dbartsMixedMatrix, whose dimnames are
  list(NULL, cols).
- validateXTest drops test row names for data-frame and integer-matrix
  tests. The integer-matrix branch rebuilds with `matrix()` and loses
  the column names too, a pre-existing bug.
- The only surviving record is the fit's named `na.action`.
- Hazard expansion repeats subject names in the training x
  ("s2 s2 s2"); a data-frame test loses them.
- xbart returns a loss vector, not per-observation output: out of scope.

Memory:
- A name set on a stored array that is shared costs nothing until the
  first C call that needs a writable pointer (colMeans, t, aperm, pnorm,
  [`combineOrUncombineChains`](../../R/generics.R)), which then
  duplicates the whole array (160 MB in the critique's probe).
- extract returns the stored array itself for gaussian ev, every "bart"
  type and matching-layout multinomial ev, and bartCause calls t() and
  aperm() on those.
- Names at n = 1e6: automatic "1".."n" stay a deferred ALTREP string
  (about 0 MB, 76 MB once every element is touched); distinct custom
  names cost 69 MB in memory and 2.5 MB compressed, against 7.6 GB for
  one 1000-draw channel.

Landscape:
- `na.action` at predict: stats, rpart and mgcv take
  `na.action = na.pass` before `...` through `model.frame`, and their
  names identify rows. randomForest drops internally, then pads. gbm,
  xgboost, ranger, BART, bartMachine and stan4bart have no formal.
- Row names: model.frame fitters (lm, glm, rpart) always name rows,
  "1".."n" when the data has none. lm.fit, BART::wbart and ranger
  name nothing for unnamed input.

## Slice A1: row names

Carrier. The names live outside x, in a new dbartsData slot `rowNames`.
It is typed "ANY", like the existing `na.action` slot in
[R/A_class.R](../../R/A_class.R). It holds list(train, test) or NULL and
is filled at dbartsData entry from the raw inputs:
- the model frame's rows after `subset` and the `na.action`;
- `rownames(x)` (a data frame's automatic names count as "1".."n") after
  the keep mask;
- `rownames(test)`, captured before validateXTest.
A matrix without row names gives NULL (precedent: lm.fit, wbart,
ranger); the formula path always names rows, as lm does.

Why outside x. The container alternative must keep the names in step
through:
- `[.dbartsMixedMatrix`, `dimnames` and `as.matrix`;
- every builder, including the sparse lifts;
- setPredictor and xbart fold subsetting;
- the bridge.
It also moves stan4bart, whose `recompute_bart_block` subsets
`bartData@x` containers. The slot's surface is the class definition,
initialize and the two builder paths, and no sampler output needs names.
The slot is an entry-time record, read only when a fit is packaged:
setData brings its own slot, setPredictor keeps row identity, and xbart
reads no names. Old saved data objects lack the slot, so every read goes
through `methods::.hasSlot` and yields NULL.

Packaging. The five packagers ([`packageBartResults`](../../R/bart.R),
[`packageMultinomialResults`](../../R/bart.R),
[`packageOrdinalResults`](../../R/bart.R),
[`packageNegbinResults`](../../R/bart.R),
[`bart2Hurdle`](../../R/bart.R)) copy the slot into fit fields
`row.names.train` and `row.names.test`. They also set dimnames IN PLACE
on the stored draws arrays while those are still unshared: yhat.train,
yhat.test, their means, s.train, s.test and forestFits, so extract's
pass-through returns stay copy-free. The raw fields therefore gain names
too, which VD's list did not name; see the VD item below. Internal
callers compute on unnamed data and name small results last.

Hazard (dec-B34, QA1). Person-period rows are named with `make.unique`
over the subject names, on the training and test channels alike
(s1, s1.1, s1.2). The subject names come from the capture, not from the
expanded x. `make.unique` avoids names already taken: subjects "1" and
"1.1" give "1", "1.2", "1.1", which the tests pin.
survivalProbabilities names its subject margin from the fit fields for
the stored train and test channels (`newdata = NULL`, and AFT), and from
`rownames(newdata)` otherwise.
[`survivalProbabilitiesFromDraws`](../../R/bart.R) must keep names. The
`s` attribute is named too.

Scope. Every extract type with an observation margin (ev, ppd, bart,
loglik, forest, contribution) on every class and sample; fitted;
residuals; every predict type; survivalProbabilities. Out: sigma, k,
varcount, trees, and the sampler's methods. One helper,
`nameObservationMargin(result, names, trailing)`, handles the margin
rule described in slice B, item 4.

Steps A1:
1. [R/A_class.R](../../R/A_class.R) and [R/data.R](../../R/data.R): add
   the slot and the capture on both builder paths. Fix validateXTest's
   integer-matrix branch with `storage.mode<-`, which keeps dimnames.
2. Packagers: fill the fields and name the stored arrays in place.
3. Hazard: apply `make.unique` naming to both channels and derive the
   subject names.
4. Name every return site that computes a fresh result.
5. Tests: a new inst/tinytest/test-row-names.R covering class x channel
   x path (formula, x-path data frame, matrices with and without names,
   dgCMatrix, factor designs, bartBT, `subset`, hazard, AFT).
   - Save and load: a fit and a data object with the new parts stripped
     give unnamed output, plus a one-off probe of an rds saved by the
     pre-slice build.
   - tracemem: extract(type = "bart") and gaussian ev return the stored
     object untraced, and fitted's colMeans makes no copy of the array.
   - Existing sites that now fail on attributes (102 single-line
     accessor comparisons and about 35 data-frame files are candidates):
     compare against a named expectation, or use
     `check.attributes = FALSE` where names are not the point. Never
     re-record a snapshot for names.
6. Harnesses: unname both sides at compare time in equivalence.R,
   bcf-equivalence.R, multinomial-equivalence.R and the two reduction
   gates. Baselines are not re-recorded.
7. Rd: one home paragraph in bart.Rd's Value, with bartBT.Rd and
   dbartsData.Rd (the slot) linking to it. NEWS: a new item, since 0.9-x
   named nothing.
8. Consumers:
   - bartCause dbarts-1.0: extract output flows into mu.hat, and sd.obs
     and sd.cf become named. `as.matrix(responseData@x)` gains no
     names, since the carrier sits outside x. Budget 10-40 lines of
     edits in bartCause, then its full suite.
   - stan4bart bartcore: the container is unchanged and it calls no
     fit accessor, so a smoke run suffices.

Verification A1:
- The full tinytest suite.
- Exact gates in quick mode: the list in
  [.github/workflows/exact-gates.yaml](../../.github/workflows/exact-gates.yaml),
  run as `Rscript benchmarks/R/<gate> quick`.
- On the reference build: the seeded-drift snapshot files and the three
  bitwise compares report "identical draws" on every scenario.
- bartCause's full testthat suite; a stan4bart smoke run.
- The lint chain in CLAUDE.local.md, `R CMD check --as-cran`, and the
  NEWS parse gate.

## Slice A2: fit-time na.exclude parity

The four packagers store the fit's `na.action`, as bart does.
fitted and residuals on multinomial, ordinal, negbin and hurdle pad
through [`padOmittedRows`](../../R/data.R). That covers vectors, factors
and obs x K matrices; draws arrays stay unpadded, as on bart. A1's names
fill the padded positions from the named record. Tests go in
test-row-names.R, one per class under `na.exclude` and `na.omit`.
Also in A2, three fit-time defects the A1 review found (they predate
A1): a formula-path hazard fit under `na.exclude` pads a subject-level
record into person-period outputs (fitted comes back the wrong length,
names out of order); a multinomial fit under `na.omit` with an NA in a
matrix `x` fails on a length mismatch between x and y; an aft fit under
`na.omit` with an NA in a matrix `x` fails on the status length. Each
gets a test.
Verification: the tinytest suite and the lint chain.

## Slice B: na.action

1. Form. A function, or a name through `match.fun`. The default is
   `dbarts::na.keepPredictors`, the fit's own. `na.action = NULL` means
   that default: stated in the Rd, and `getOption("na.action")` is never
   consulted, as in predict.lm.
2. Placement. Immediately before `n.threads` on every predict method,
   keeping [3. D1 after: the signatures](predict-surface.md#3-d1-after-the-signatures)'s
   invariants. On survivalProbabilities.bart it follows `combineChains`,
   the last formal before `...`.
3. Semantics. `na.pass` is special-cased by identity: unroutable rows
   become NA. Every other function is applied through a response-free
   variant of [`applyNaActionToXY`](../../R/data.R), followed by today's
   refusal on the kept rows:
   - `na.keepPredictors` refuses as today;
   - `na.omit` drops;
   - `na.exclude` drops and pads;
   - `na.fail` errors, naming the columns (the column loop already
     exists);
   - a custom function sees a synthetic one-column frame, NA on the
     incomplete rows, on every path including the formula path; the Rd
     says so.
   The refusal message gains a clause naming `na.pass` and `na.omit`.
   The manual says that only the default and `na.pass` look at
   routability.
4. Pad and drop at the end. The row margin is 1 under `ci.level` and for
   `class`. Otherwise it is the last margin, or the one before it when a
   category or forest axis trails, verified in both chain layouts and on
   single-row results. The `s` attribute is padded too. Binary ppd is
   integer, so it pads with NA_integer_. Names come from A1. No
   `na.action` attribute is attached, following predict.lm.
5. No surviving row, and zero-row newdata. Today these hit the engine's
   "requires rows" or validateXTest's "must be numeric". Instead,
   predict one placeholder row:
   - the fit's first training row, which is routable by construction;
     for survivalProbabilities, the first training subject-level row,
     before expansion;
   - its per-row inputs are stubbed: offset NULL (a 1 x K zero matrix
     on a multinomial fit with a category offset), no weights, and each
     basis's first stored row.
   The RNG is protected with `had <- exists(".Random.seed", globalenv())`
   plus the saved value and an `on.exit` that restores it, or removes
   it when it was absent. sampleFromPPD is not a precedent, since it
   never restores a missing seed. Then slice to width zero (`na.omit`,
   zero-row input) or fill NA to nrow(newdata) (`na.exclude`,
   `na.pass`), with names. A single-row result keeps its dimensions
   (probed), so the slice is exact.
6. Per-row channels. A length-1 offset or weight passes through
   untouched, since today it recycles. A length-n offset, an n x K
   offset, length-n weights and bases (bare or listed) are length-checked
   against nrow(newdata), then subset. The raw data frame is subset for
   [`replayForestBasis`](../../R/model.R). Their own NAs are out of
   scope. Today an NA offset gives an NA row silently, and an NA weight
   gives NaN plus an rnorm warning: a TODO candidate.
7. Validate once. predict.bart validates at the top, then calls a
   factored post-validation body of the sampler's `$predict` and
   `$predictForests`, so no warning fires twice. The R5 methods call
   the same body after their own validateXTest, so the sampler surface
   is unchanged.
8. Fix, agent-made. In [`mapFactorColumnsToTrainingLevels`](../../R/utility.R),
   the check `anyNA(refactored) && !anyNA(column)` lets an unseen level
   through as a missing value whenever that column already holds an NA.
   Test `any(is.na(refactored) & !is.na(column))` instead, so an unseen
   level is always refused and never treated as missing. Add a test.
9. RNG. ppd RNG is consumed for kept rows only, so each result is
   bitwise equal to `predict(fit, newdata[kept, ])` under the same seed.
   This is neutral.

Steps B:
1. R/data.R:
   - add `unroutableTestColumns` and `unroutableTestRows` (sparse-aware),
     and give validateXTest an internal `refuseMissing`;
   - add `resolvePredictRows` and `padPredictedRows`.
2. The five predict methods and survivalProbabilities.bart (resolved at
   the subject level before expansion).
3. The placeholder path, the post-validation body (item 7) and the
   factor fix (item 8).
4. Rd: man/na.keepPredictors.Rd is the home for fit-time versus
   predict-time meanings, routability and NULL. man/bartBT.Rd and
   man/bart.Rd hold the usages and link there.
5. NEWS: extend the existing missing-predictors item. At landing, amend
   [Bridge and R surface](../design/mia-missingness.md#bridge-and-r-surface)
   and clear the TODO item.

Verification B:
- A new inst/tinytest/test-predict-na-action.R, per class plus
  survivalProbabilities:
  - the default;
  - `na.pass` and `na.exclude`: positions, names, and identity with
    kept-row predictions across ev, seeded ppd, bart, class, ci.level,
    split chains and the forest arm;
  - `na.omit`;
  - `na.fail` naming the columns;
  - NULL;
  - no surviving row and zero-row input: dimensions and names correct,
    with `.Random.seed` identical whether it was present or absent
    beforehand;
  - a length-1 offset passed through;
  - sparse newdata;
  - the unseen-level fix;
  - a single warning from a positional-match newdata;
  - warnings counted with `withCallingHandlers`.
- The full tinytest suite, the lint chain, `--as-cran`, and
  check-rc-codoc (the sampler signatures are unchanged).
- Mutation proof: make `na.pass` fall through to the refusal, force
  margin 1, and revert the factor fix.

## For VD

1. The raw fit fields (`$yhat.train`, `$yhat.test`, their means, `s.*`,
   `forestFits`) gain names, beyond the four functions VD listed. This
   is what keeps extract copy-free. The alternative is accessor-only
   names, which cost one full copy of a channel on the first C
   operation after each extract. Recommendation: name the raw fields.
2. Pre-existing, and outside this item: at fit time,
   `na.action = NULL` on the formula path falls back to
   `getOption("na.action")` (na.omit; probed, it dropped a
   predictor-NA row), while the matrix path applies no action and
   errors on a missing response. This is a TODO candidate.

## Overlap

The extract scalar-types slice has landed (98524aac, b5ae43e7). Its
sigma, k and varcount branch has no observation margin and is left
alone. A1 and A2 touch every extract, fitted and residuals method in
R/generics.R. B touches the predict methods, predictBlend, hurdleParts,
survivalProbabilities.bart and the sampler's predict body.

## Landing

Slice A1 LANDED 2026-09-27 (e0f80109): observation row names on every
observation-indexed output, recorded in the data object's rowNames slot
and applied in place at packaging; whole-matrix test setters record the
new test set's names; bartCause 2301e4b takes the counterfactual draws'
names from the observed draws. A2 and B remain.
