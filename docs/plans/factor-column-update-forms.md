# factor-column-update-forms: the joint row-by-row update takes a factor column's labels

Status: LANDED 2026-10-07 (cde886b9 to e4afbfbd; dec-B279, dec-A172).

agent: sonnet implementer, one (R, tests, manual); opus reviewer.
rng: by call sequence.
- NEUTRAL, bit for bit, for every call that gives numbers, on any column, and for every call on a column
  no sampler holds as a factor: `updatePredictorPerObservationJointly` installs what it installs today and
  draws its scan order as today.
- A changed result, not a generator matter, for a factor, a character vector, a `sparseFactor` or a
  logical given for a column a sampler holds as a factor: today read through `as.double`, afterwards
  matched to the column's levels by label or refused by name. A refusal draws nothing.
Proved by the full tinytest suite, by a seeded digest of numeric calls on the base and slice builds, and by
the new tests. No C or C++ changes, so no draw path can move.
window: pre-release, before the merge to main (dec-B279). Shares one manual page, in another item, with
[monotone-unforced-refusal.md](monotone-unforced-refusal.md), in review, whose new tests give the joint
form numbers and are not touched by this; shares no file with
[cross-family-state-install.md](cross-family-state-install.md), [leaf-conversions.md](leaf-conversions.md)
or [written-surface.md](written-surface.md). Recommended order: after monotone-unforced-refusal, any time.
budget: ~270 lines (R ~70, C++ none, tinytest ~170, manual and TODO ~30), upper figure 420. The scratch
build under Context came to 69 lines of R.

## Goal

`updatePredictorPerObservationJointly` takes, for a column the samplers hold as a factor or an ordered
factor, what `setPredictor` takes for a named column: labels, matched to the column's levels by name, with
a label the column lacks and a first missing value refused in the same words. It refuses numbers for such a
column, as the column form does. Nothing given for such a column is installed as another level
without a message, and the manual says what each kind of value means.

## Context

All numbers were run on the tip's build (shipped mode): 200 rows, a numeric x1, an unordered factor f
(levels a, b, c, d) and an ordered factor o (lo, mid, hi), 5 trees, one chain; one sampler, and two that
hold the columns at different positions. Each new value is the row's next level.

- How a sampler holds a factor column. `data@x` holds it as codes from 0 in the order of the levels, which
  are in the design's `factor.levels` attribute; ordered and unordered factors alike.
- The joint form today. [`updatePredictorPerObservationJointly`](../../R/updatePredictorPerObservationJointly.R)
  passes its values through `as.double` and the bridge refuses a value that is not a code of the column
  ([`validateColumnValues`](../../src/R_interface_bartcore.cpp)). A factor therefore arrives as R's codes
  from 1. Given for f, with one sampler and with two alike:

  | given | the joint form | by column, and `"partial"` |
  |---|---|---|
  | factor, every level present | E | right |
  | factor, no row at the last level | 200 installed, each as the next level, no message | right |
  | factor, levels declared in another order | E | right |
  | factor, unused levels dropped | 200 installed, each as the next level | right |
  | character labels | 199 installed as missing values, with R's coercion warning | right |
  | `sparseFactor` | an error from the coercion | right |
  | codes from 0, integer or double | right | N |
  | codes from 1 | E; with no row at the last level, 200 installed as the next level | N |
  | a code that is not whole | E | N |
  | logical | 199 installed as codes 0 and 1, 151 of them at another level | N |
  | labels with a first missing value | factor: E; character: 199 installed as missing | M |
  | codes from 0 with a missing one | 200 installed, the missing value with them | N |
  | a label not in the column | factor: E; character: 199 installed as missing | L |

  "Right" is all 200 rows holding the label given. E is the engine's `categorical predictor values must be
  existing category codes`; N, M and L are the column form's `column 'f' is categorical; give its values
  as a factor or character vector of its labels, not numbers`, `column 'f' has missing values, which its
  training values do not` and `column 'f' has label 'z' not among its training levels`.
  The ordered column behaves the same, the engine's message being `ordered factor predictor values must
  be existing level codes`. Every sampler ran three sweeps and copied after each install, right or wrong.
- The code that matches labels. [`codeCategoricalColumnUpdate`](../../R/bartcore.R), called by
  [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) for every update that names columns, the row-by-row
  form included: it reads the column's levels from the design, matches a factor, character vector or
  `sparseFactor` by label, and raises N, M and L. Given the same factor for two samplers that hold the
  column alike it returns identical codes, and those codes handed to the joint form install 200 of 200
  rows at the label given, for a factor in any level order, an ordered factor, a character vector and a
  `sparseFactor`. The joint form can call it as it stands.
- Several samplers whose column differs. Levels in the reverse order in the second sampler: codes install
  in both and mean other labels in the second (0 of 200 rows hold the label the caller coded); the helper
  returns different codes for the two. A column held as a factor in one sampler and as a number in the
  other: codes install in both.
- A numeric column given a factor is read through its integer codes by every form (x1 then holds 1 to 4,
  198 of 200 rows installed, no message), and a character vector becomes missing values with R's
  coercion warning. Not a factor column, and not changed here.
- Callers. The package's seven tinytest files that name the function make 180 calls, 177 with numbers on
  a numeric column and 3 that test a refused argument (traced on the tip's build); the vignette and the
  help example give numbers too. bairrtt calls it with a numeric latent trait; stan4bart, bartCause and
  treatSens do not call it (read only). The tests of monotone-unforced-refusal call it with codes from 0,
  a missing one among them.
- The rule below, tried on a scratch copy of the tip (R only): every label row of the table installs 200
  of 200 rows at the label given, with one sampler and two; a label not in the column, a first missing
  value and a logical are refused in the words below, and after each refusal both samplers' stored
  states and designs are byte for byte the ones before and R's generator has not moved; codes from 0,
  with and without a missing one, install as today; labels through the joint form and through
  `"partial"` on a twin give the same mask, stored state, design and next three draws in 54 of 54
  fixtures, 15 with rows declined; a seeded digest of 18 numeric calls and the sweeps between them is
  equal on the two builds; the tinytest suite passes unchanged (14750 results).

## The rule

By what is given, and by whether the samplers hold the column as a factor (ordered or not):

1. No sampler holds it as a factor: the values go through `as.double`, as today.
2. Numbers (integer or double) for a factor column: refused, as the column form refuses them, with the
   words in 4. For a numeric column they install as they do.
3. Labels - a factor, an ordered factor, a character vector or a `sparseFactor` - for a column every
   sampler holds as a factor: matched to the column's levels by name, by the column form's helper. The
   order and number of the levels of the factor given do not matter. Refused, in the column form's words:
   a label the column does not have (`column 'f' has label 'z' not among its training levels`), and a
   missing value when the column holds none (`column 'f' has missing values, which its training values
   do not`); a column that already holds one takes another.
4. Anything else for a factor column, a logical among them: refused with `column 'f' is categorical;
   give its values as a factor or character vector of its labels, not numbers`, the last two words
   only where the value is a number.
5. One vector goes to every sampler, so labels must code alike in each. Levels that differ between
   samplers: `column 'f' has other levels in sampler 2 than in sampler 1, so its labels cannot be
   installed in both`. A factor in one and a number in another: `column 'f' is categorical
   in sampler 1 and not in sampler 2, so its labels cannot be installed in both; update them in separate calls`.
6. Neither is read as the other. A factor is known by its class before anything is coerced, so its integer
   codes are never read. A number is never matched to a label, even where the labels are numerals: for
   levels "1" to "4", the character "2" is the second level and the number 2 is the third.
7. Every refusal is raised in R before any sampler is touched or the scan order drawn. The value, the
   maintenance of `data@x` and `updateState` are as today.

## Constraints

- Every call that gives numbers installs and returns what it does now, bit for bit.
- The column form and `"partial"` are not changed: they refuse numbers for a factor column, and the
  helper's messages are theirs.
- A first missing value given as a code keeps installing; the tests of monotone-unforced-refusal reach
  their refusal through it.
- No C or C++ changes, no new export, no change to the function's arguments or value.
- No NEWS entry: the function is new in 1.0-0.
- Out of scope: see the end.

## Steps

1. R. A helper beside [`codeCategoricalColumnUpdate`](../../R/bartcore.R) applies the rule to the list of
   samplers and their column indices and returns the one numeric vector: it reads which samplers hold the
   column as a factor from each design's `factor.levels`, returns numbers and columns no sampler holds as
   a factor untouched, raises the refusals of rules 4 and 5, and otherwise calls
   [`codeCategoricalColumnUpdate`](../../R/bartcore.R) once per sampler and compares the results.
   [`updatePredictorPerObservationJointly`](../../R/updatePredictorPerObservationJointly.R) calls it where
   it calls `as.double` today, before any pointer is fetched.
2. tinytest, a new file `test-joint-update-factor.R` on the fixture above; "fails today" names what the
   tip does. After every accepted call the design column is compared with the codes of the labels given.
   - A factor with every level present, one with none at the last level, one with its levels declared in
     reverse, one with unused levels dropped, an ordered factor, a character vector and a `sparseFactor`:
     every row installed holds the label given, one sampler and two, f and o; no warning (counted with
     `withCallingHandlers`); the samplers sweep, copy and reload. Fails today in each (an error, the next
     level, or missing values).
   - Identity with the column form: a sampler given labels through the joint form and a twin given them
     through `setPredictor(x, "f", forceUpdate = "partial")` return the same mask and hold the same
     stored state and design, on a fixture where rows are declined. Fails today.
   - A label not in the column, a first missing value in a factor and in a character vector, a logical:
     the error in the words of the rule; both samplers' stored states and designs identical to before;
     `.Random.seed` unchanged. Fails today in each (installed, or another error).
   - A column that holds a missing value takes labels with another. Fails today: an error.
   - Codes from 0, integer and double, with and without a missing one: installed, the design holding the
     codes given; a code equal to the number of levels and one that is not whole are refused in the
     engine's words. Holds today.
   - Numerals as labels: for levels "1" to "4" the character labels install at their own level and the
     numbers 0 to 3 as codes. The first fails today; the second holds.
   - Two samplers with the levels in another order, and one holding the column as a number: labels are
     refused naming the samplers, numbers install. The refusals fail today.
   - A numeric column: numbers install as today.
3. Mutations (Verification): apply, install, run, report the failing counts, revert, `touch`.
4. Records.
   - Manual, [updatePredictorPerObservationJointly.Rd](../../man/updatePredictorPerObservationJointly.Rd),
     the `x` item, replacing "A numeric vector of new values for the shared column, of length equal to
     the number of observations": "The new values for the shared column, one per observation. For a
     numeric column, numbers. For a column the samplers hold as a factor or an ordered factor, its
     labels: a factor, a character vector or a `sparseFactor`, matched to the column's levels by name,
     as `setPredictor` takes them for a named column. The order of the levels of the factor given does
     not matter; a label the column does not have is refused by name, and so is a missing value when
     the column holds none. Numbers are refused for such a column, as `setPredictor` refuses them for a
     named column, even where the labels are numerals. With several samplers,
     labels need the column to have the same levels in the same order in each; where it does not, or
     where one holds it as a number, update them in separate calls."
   - TODO: the entry `factor-column-update-forms` names this plan for the joint form and keeps its
     sentence on what is not planned.

## Verification

- `R CMD INSTALL` into the slice's own library; the full tinytest suite, 14750 results on the tip, with
  the new file's added; tests/cpp builds and passes as on the base (no file under src/ or tests/cpp
  changes; the diff says so).
- One script on the base and slice builds digesting a seeded run of joint calls with numbers - codes on f
  and o, one with a missing code, values on x1, one sampler and two - with the masks returned, the sweeps
  between them and the final stored states and designs. Equal. With no source under src/ changed, the
  reference-build snapshot files, the three equivalence compares and the exact gates are not re-run; none
  calls the function (checked by search at the tip). If the reviewer runs a compare all the same,
  `EQUIVALENCE_CORES=2`.
- The new tests of monotone-unforced-refusal, once that plan has landed, pass unchanged.
- Mutations, each expected to fail the named test:
  - the helper is not called (values go through `as.double`): "every row installed holds the label given"
    and the identity with the column form;
  - labels are coded by the given factor's own level order: "levels declared in reverse" and "unused
    levels dropped";
  - the check for a label the column lacks is dropped, or the one for a first missing value: the refusal
    tests, each;
  - numbers are matched as labels: "codes from 0" and "numerals as labels";
  - a logical passes: its refusal test;
  - only the first sampler's levels are read: "two samplers with the levels in another order";
  - the samplers' levels are compared after the engine has been called: "stored states and designs
    identical to before";
  - the design keeps the labels instead of the codes: the comparison after every accepted call, and the
    copy and reload.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; `R CMD build` with every vignette rebuilt (one calls the function) and `R CMD check --as-cran`
  on the tarball from a clean copy (R/ and man/ change).

## Out of scope, and where it goes

- The whole-matrix form given a data frame. It fails in `as.double` today, though the manual's `x` item
  speaks of a whole data-frame update. The same lines do not fix it: the helper, given a frame and every
  column, codes it and the forced update then installs 200 of 200 rows right, but only with the frame's
  columns in the design's order (another order is refused as numbers given for a factor), it refuses a
  first missing value that the matrix of codes installs, and a design with no factor column still fails.
  It stays in the TODO entry, not planned.
- The column forms refusing codes, and codes carrying a first missing value where labels may not. Both
  stand; the TODO entry keeps them.
- A numeric column given a factor or a character vector, read through `as.double` by every form. The
  coordinator has the measurement.
- A code vector that was meant as R's codes from 1. Nothing can tell when no row is at the last level;
  the manual states the convention and labels are the way to avoid it.

## Calls made in planning

- Labels are matched by the column form's own helper, not by new code: one rule, one set of messages, and
  a fix to either form reaches both. The helper refuses numbers, so the joint form's helper returns them
  before calling it.
- Numbers stay codes. dec-B279 leaves open matching labels or refusing a factor by name; the first is
  what a caller of `setPredictor` expects, and it costs the same. Refusing numbers too, as the column form
  does, would make the two forms alike; it is not done because `data@x` holds codes, a sampler that moves
  a latent class each sweep carries codes, and the tests of monotone-unforced-refusal give them. It is
  the question in the coordinator's notes.
- A logical is refused for a factor column. Today it is read as codes 0 and 1; no reading of it is a
  label.
- Labels for samplers whose levels differ are refused, not coded per sampler: the engine takes one vector
  for all of them, and a call that installed one label as two different levels would be the defect again.
  Numbers are not checked across samplers: the codes are the caller's then, as today.
- A missing label in a column that holds none is refused, as the column form refuses it, while a missing
  code installs. The asymmetry is the column form's and the TODO entry's "rest"; closing it either way
  changes a form this plan leaves alone.
- The tests get their own file; the existing joint-form file gives only numbers and is not edited.
- Sonnet implements: the change is R and manual text. The reviewer is opus for the identity with the
  column form and the refusals' "nothing touched".
- Order. This and cross-family-state-install share no file but TODO and the plans index and can run in
  parallel worktrees; this one after monotone-unforced-refusal lands, since both edit the function's
  manual page (its Details there, the `x` item here).
- The tip against dec-B279 and the TODO entry. Both say every level is installed as the next one unless
  the last level is present. That holds for a factor in the column's own level order; a character vector
  is installed as missing values with R's coercion warning instead (199 of 200 rows), a `sparseFactor`
  fails in a coercion, and a logical is read as two codes. The entry's "refused as not an existing
  category code" is the engine's message for an unordered column; an ordered one has its own.
  Nothing in this slice was found done already.
- Review changes to rules 4 and 5, and what is known. Samplers that hold the column as a factor with a
  different level table (levels or their order) are refused for numbers as well as labels, before any
  sampler is touched: `the samplers hold column 'f' with different levels (sampler 2 differs from sampler
  1), so one value would be a different level in each; update them in separate calls, or create them with
  the same levels in the same order`. A subset or superset of levels is refused too, so the help's "same
  levels in the same order" is exact. A call that passed numbers to such samplers installed a different
  level in each in silence. Numbers to samplers with one table, and to a factor in one sampler and a
  number in another, are as before.
- A refusal the helper raises for a later sampler (a missing value where that sampler's column holds
  none) ends with `(sampler 3)`, the sampler's position in the list; with one sampler, no suffix.
- A factor or `sparseFactor` for a column no sampler holds as a factor is refused: `column 'x1' is numeric
  and cannot take labels`. Text is refused for it only when some element is not numerals, since numerals
  as text read as numbers today. `setPredictor` by column and `"partial"` read a factor for a numeric
  column through its codes (x1 holds 1 to 4 afterwards) and are not changed here: the helper returns a
  numeric column's values untouched at two early exits, and a whole-frame update goes through them.
- A one-column character matrix gets the joint form's own refusal, not the column form's text about a data
  frame. The six `.Random.seed` assertions are dropped: no accepted call moves R's generator and the
  twin draws test the samplers' own.
- Known: codes given as a `Matrix` `sparseVector`, a `Date` or a `difftime` on a factor column were
  installed on the base build and are refused now; only a base numeric vector is read as codes.

## Landing note

Landed 2026-10-07 as cde886b9 to e4afbfbd on bartcore, 4 commits, 628 lines added over 5 files against a
planned 270 to 420; R, the manual and tests, nothing under src/. One review told to refute: LAND AFTER
FIXES; what the slice installed was right in every case it ran.

What the review changed (dec-A172). The first tests built samplers with no split and read back the
R-side copy of the column alone, so a version that handed the engine the passed factor's own codes
passed all 224; they now run on swept samplers whose trees split on the column and hold every accepted
input to a twin updated through `setPredictor(forceUpdate = "partial")` or through codes, by its
predictions and its next draws, for one to three samplers and for the unordered and the ordered column,
with a copy and a reload. The refusal for samplers whose levels differ told the caller to give numbers,
which installs a different level in each sampler without a message; samplers that hold the column with
different level tables now refuse one vector of any kind, labels or numbers, naming the first sampler
that differs. A refusal raised by a later sampler's levels names that sampler. A factor, a `sparseFactor`
or text that is not numerals given for a numeric column is refused in the joint form. The help says that
`as.integer()` of a factor counts from 1 and is not these codes, and that a missing number is installed
where a missing label is refused.

Left: `setPredictor` by column and `"partial"` still read a factor given for a numeric column through
its integer codes, the fix sitting in a helper that whole-frame updates share; and a data frame as the
whole matrix. Both are in TODO `factor-column-update-forms`. Codes given as a `sparseVector`, a `Date` or
a `difftime` on a factor column were installed before and are refused now.

Gates at landing, on a clean copy of the rebased tree in a library of its own (shipped mode), run in
series: tests/cpp 350 ok, 0 failed; the full tinytest suite 15378 results, 0 failed, 225 files; lintr no
lints; air, rc-codoc, win-drift, doc-freshness and the mutation battery's anchors clean; `R CMD build`
with every vignette rebuilt and `R CMD check --as-cran` with the Date NOTE alone. Neutral for numbers: a
seeded digest of calls that pass numbers, with sweeps between, is equal on the base and the slice (43
items by the implementer, 312 calls by the reviewer). bairrtt, the one consumer that calls the function,
passes numbers for a numeric column; two of its test files and a fit digest were equal on both builds.

Since 2026-10-07 (dec-B286, [small-rulings-batch.md](small-rulings-batch.md)): numbers given for a factor column are
refused, and the help no longer speaks of codes.
