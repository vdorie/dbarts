# leaf-spread-rounding: a leaf covariate whose spread is rounding is one value

agent: opus
rng: neutral for every column above the tolerance (the trio bitwise); a
column under it is refused or kept where it ran before
budget: ~60 engine and bridge lines, ~10 R lines, ~250 test lines, help

Status: BUILT 2026-11-07, reviewed SOUND AFTER NITS, fix round reviewed LAND (2f2457f20, 26940637c)

## Goal

A linear or gp leaf covariate whose spread about its center, over the rows
in the likelihood, is under 64 times machine precision of its size,
`||z - center|| < 64 * DBL_EPSILON * ||z||`, is read as one value: refused
at creation with a message that names the test and advises only what works
for each kind of column (dec-B673), and treated as a column of one value by
every setter that re-derives the center and scale (dec-B654, dec-B655). A
column above the tolerance draws what it drew.

## Context

- `standardizationMomentsForColumn` (src/bartcore/data.hpp) takes the
  center and sd over the observed values of the rows in (`LeafRows`:
  positive weight and kept by the mask), each at its count where the rows
  are counts (frequency weights, a binomial row's trials times its weight,
  nbinom counts), and says where they
  hold fewer than two distinct values; a constant column is centered at its
  value exactly (dec-B358).
- `deriveLeafMoments` (src/bartcore/model.hpp) refuses that case with
  `LeafStandardizationError`, or takes the placeholder (center, sd 1) for a
  gp leaf that states its lengthscales. It serves creation, an xbart fold,
  `setData` that derives (through `Sampler::refuseLeafDerivation` before
  anything moves, then the leaves' `reinitialize`).
- `Sampler::rederiveStandardization` calls the moments function directly
  and keeps the numbers in force for a covariate the rows in hold at one
  value; `warnScaleKept` (R/bartcore.R) warns once per sampler. It serves
  `setPredictor` (both bridge entries) and `setWeights` and
  `setActiveRows` with `updateStandardization = TRUE` (dec-B465).
- The bridge words the refusal in `leafStandardizationMessage`
  (src/R_interface_bartcore.cpp), doors creation, fold and setData.
- The case: half the rows 0.3 and half `0.1 * 3`; scale 3.95e-17 today and
  a linear leaf predicting about -2.7e12 at 0.31 (probe
  scratch/backlog-dossier/probes/c404.R).

## Design

- The test. A new `leafCovariateMoments` in data.hpp returns what
  `standardizationMomentsForColumn` returns, and also true, with sd the
  placeholder 1 and the center left as computed, where `spreadIsRounding`
  finds the observed values of the rows in within the tolerance.
  `deriveLeafMoments` and `Sampler::rederiveStandardization` call it, so
  every door that derives reads one test. The moments function itself is
  unchanged.
- The rows and weights: the rows in, as the center and scale read them, each
  counting once, or, where the rows are counts, at its count, so the test is
  that of the rows written out; a precision weight does not enter, as it
  does not enter the center and scale. Missing values are skipped.
- The norm. Both sums are taken after scaling every value by the power of
  two that brings the largest magnitude among the rows in into [0.5, 1)
  (`std::frexp`, `std::ldexp`): exact, so it cannot overflow or lose bits,
  and the larger sum is at least 0.25. The test compares the squared sums,
  `A < (64 eps)^2 B`.
- The center. The test measures the spread about a refined center, the
  stored center plus the mean of the deviations from it (the corrected
  two-pass mean), not about the stored center: a mean taken through the sum
  misses by up to about n eps of its size, which adds that much to the
  measured spread and would let a rounding column of many rows pass. The
  refined center is good to about eps whatever n, so the test does not grow
  with n. The stored center and sd keep their bits.
- The message. `LeafStandardizationError` gains `rounding`: the rows in hold
  two or more distinct observed values. The bridge words it: the leaf
  covariate "z" holds values whose spread about their center is under 64
  times machine precision of their size (`||z - center|| < 64 *
  .Machine$double.eps * ||z||`) over the rows in ..., so its center and
  scale cannot be taken from them; drop it from linear(), or, if its values
  truly differ, recompute it from its source about an origin near its
  values, or change the rows in (...) (dec-B673: subtracting a constant from
  the stored values restores no lost digit, and on values meant to be equal
  leaves a column that passes the test at a rounding-sized scale). A gp
  leaf's remedies lead with stating its lengthscale, as for one value. Same
  three doors.
- The setters.
  - creation, xbart fold: refused, the message above.
  - `setData` deriving: refused before anything moves, `$setData:` door.
  - `setData` keeping or supplying: derives nothing, untouched.
  - `setPredictor` (whole matrix or columns) with `updateStandardization =
    TRUE`: the covariate keeps its numbers, `dbartsScaleFallbackWarning`
    once per sampler.
  - `setWeights`, `setActiveRows` with `updateStandardization = TRUE`: the
    same.
  - `setForestBasis(updateBasisScale = TRUE)`: a forest basis, not a leaf
    covariate; out of scope.
  - warm start, `setState`, a data object's recorded standardization:
    supplied numbers, no derivation.
- The keep warning says "fewer than two values of the leaf covariate x that
  differ by more than rounding are in the likelihood", true of both cases.
- The gp leaf. A gp leaf without stated lengthscales is refused as the
  linear leaf is. One that states them takes the placeholder, as it does
  for one value: centered at the corrected mean and divided by 1, every
  training row then reading within rounding of zero, and marked as without
  spread, so its record shows the scale NA a column of one value shows.

## Calls (each with its alternative)

1. A wrapper both derivation paths call, the moments function unchanged;
   alternative, the test inside `standardizationMomentsForColumn`.
2. Rows in, read as the center and scale read them: once each, or at their
   counts where the rows are counts (rebase onto frequency weights, which
   gave the center and scale counts), so a frequency fit's test is its
   expansion's; alternatives, every row once whatever its count, which
   would answer differently from the rows written out, or precision weights
   too, which the center and scale do not use.
3. The refined center inside the test; alternatives, the stored center
   (the test grows with n) or a compensated stored center (moves bits on
   columns above the tolerance).
4. Power-of-two scaling; alternative, plain sums of squares, overflowing
   above about 1e154.
5. A rounding column under a gp leaf with stated lengthscales centered at
   the corrected mean the test computes (fix round, the review's F3), so
   every row reads within about a unit in the last place of zero whatever
   n; alternatives, the mean through the sum (20 units off at 200 rows,
   more as n grows) or its first observed value. Its spread mark
   (`standardizedColumnHasSpread`) applies the tolerance to the gathered
   values' distance from the center in force, so the record shows NA as
   for one value; alternative, the mark as first built, reading any nonzero
   value as spread. The distance is from the center, not the values'
   spread about their own mean, so a column moved to one value off its
   center keeps its mark and its conversion by the formula (dec-B328). A gp leaf reads
   the mark only for its record. A linear leaf reaches the placeholder
   scale only from a supplied record; a setPredictor that keeps its numbers
   can then move the column to values at rounding apart, which the next
   conversion (setData, a warm start) now reads as one value, keeping a
   live block as for one value, the one case where a draw moves.
6. The keep warning reworded for both cases; alternative, a separate text
   for the rounding case, carried through the bridge's keep answer.
7. NEWS: none. Linear and gp leaves are new in 1.0-0; 0.9-34 had none.
8. The re-derivation after setPredictor, setWeights and setActiveRows
   (`rederiveOrRelease` in the bridge) catches any C++ exception, releases
   the sampler and raises an R error, the sampler re-created from its
   stored state on next use, as a setData that fails part-way is
   (the review's F2). No path reaches it short of an allocation failure:
   the sampler's test and the leaves' own derivation call one function on
   the same values and rows, so there is no pin; the abort it guards
   against is shown by mutating the sampler's test away. Alternative, the
   abort an uncaught exception gives.
9. The refusal's lead reads "holds values whose spread about their center
   is under ..." rather than "holds values that differ only by rounding",
   which misstated a column of one outlier among equal values (the
   review's F4); alternative, the first words.
10. The orchestrator's call, after the rebase's delta review (F3), from the
   rule that frequency counts read as the rows written out: the spread mark
   at the placeholder scale reads every gathered row at the counts its
   chain reads rows at, a row of count 0 not read, and every row once where
   fewer than two have a positive count (as fallbackCounts lifts them). A
   gp leaf stating its lengthscale over a column at rounding apart at its
   counts then records scale NA, as its rows written out do. A setData
   takes the old mark before its counts move to the new rows, the counts
   then being sized to the old rows. Alternative, the mark once per row,
   which records scale 1 there.

## Constraints

- Draws do not move for a fit whose leaf covariates are above the
  tolerance: the trio bitwise against the MANIFEST's current baselines.
- No change to `inst/include/dbarts/dbarts.h`.

## Steps

1. data.hpp: `spreadIsRounding`, `leafCovariateMoments`; model.hpp:
   `deriveLeafMoments` and the error's `rounding`; sampler.hpp:
   `rederiveStandardization`. tests/cpp coverage of the test.
2. Bridge message; R keep warning; help (dbarts.Rd leaf.prior,
   dbartsSampler-class.Rd updateStandardization, the reference-class
   method docs where they name one value).
3. tinytest pins in a new inst/tinytest/test-leaf-spread-rounding.R.

## Verification

- Pins: the 0.3 against `0.1 * 3` column refused at creation, linear and
  gp, with the message; a timestamp over 100 seconds (relative spread
  1.6e-8) and a column of sd 1e-6 at 0.3 run; a column rounding-close on
  the rows of positive weight and spread on rows of weight 0 refused; each
  re-deriving setter keeps with the warning or refuses as for one value; a
  boundary pin, 1 against 1 + 126 eps (refused) and 1 + 130 eps (runs).
- Mutants: the tolerance at 0; the rows ignored in the test; the test
  skipped in `rederiveStandardization`. Each fails a pin.
- Gates: tests/cpp plain and under ASan/UBSan; the full suite; the trio
  `--bitwise --strict-coverage` on the reference build; the four
  reproducibility files; the exact gates in quick; lint, air, codoc,
  win-drift, doc-freshness; `mutation-battery.R verify-anchors`.

## Build note

Built in 2f2457f20 as planned. The six warning pins in
test-leaf-standardization-rows.R take the reworded keep text. Gates on that
code: tests/cpp plain and under ASan and UBSan; the full suite, 42316
results, 0 failed; the trio `--bitwise --strict-coverage` on the reference
build, 55, 19 and 11 scenarios identical; the four reproducibility files, 27
results; the exact gates in quick, 48 of 48; the three leaf test files under
R-loaded ASan and UBSan; lint, air, codoc, win-drift, doc-freshness,
verify-anchors. Mutants: the tolerance at 0 fails 21 of the 33 R pins and 9
C++ checks; the rows ignored in the test fails the weight-0 creation pin and
the setWeights and setActiveRows keeps; the test skipped in
`Sampler::rederiveStandardization` aborts the process at the first setter
pin, the chain's own derivation throwing where the sampler marked the
covariate derivable; the corrected center removed fails the C++ row-count
check.

Fix round, 26940637c, after the engine review (SOUND AFTER NITS). F1
(dec-B673): the refusal advises dropping the column, recomputing it from
its source about an origin near its values where they truly differ, or
changing the rows in; the help says why subtracting a constant does not
work. F2: calls 8. F3 and F4: calls 5 and 9. Two defects present before
the slice are filed in TODO (leaf-covariate-infinite-value,
leaf-covariate-moments-range). Gates on that code: tests/cpp plain; the
full suite, 42317 results, 0 failed; the trio on the reference build, 55,
19 and 11 scenarios identical, none skipped; the four reproducibility
files, 27 results; lint, air, codoc, win-drift, doc-freshness,
verify-anchors.

Rebase note, onto main past the frequency-weights, state-basis-levels and
setmodel-installs-stated landings, 2026-11-08. The base's moments function
takes counts, so `leafCovariateMoments`, `correctedMean` and
`deviationIsRounding` take them too and read each row as it does: at its
count, a row of count 0 or less not read; the refusal's one-value or
rounding reading skips such a row as well. Without counts every result keeps
its bits. The doors this extends are every door that derives where the rows
are counts: creation and an xbart fold under frequency weights, a binomial
response or nbinom counts; a binomial setData; and the re-derivation after
setPredictor, setWeights and setActiveRows on such a sampler. Pins:
test-leaf-spread-rounding.R's counts section and testLeafSpreadRounding's
counts checks (the counts ignored in the test, and the count-0 skip removed
from the refusal's reading, each fail a check). The refusal pins take the
base's rows wording, "the active-row mask".
Gates on that code, --preclean installs: tests/cpp plain and under ASan and
UBSan; the full suite, 43155 results, 0 failed; the trio `--bitwise
--strict-coverage` on the reference build, 55, 19 and 11 scenarios
identical, none skipped; the four reproducibility files, 27 results; the
exact gates in quick, 49 of 49; the leaf, frequency-weights and binomial
test files under R-loaded ASan and UBSan, 1452 results, 0 failed; lint,
air, codoc, win-drift, doc-freshness, verify-anchors.

Delta fix round, a86bf7f32, after the rebase's delta review (FIX). F1: the
center corrected at the counts is pinned at creation (counts 1000 and 1 on
halves at 1 and 1 + 1e-13), killing the mutant correctedMean ignoring
counts. F2: R pins at counts for the binomial setPredictor and setData,
setActiveRows and the xbart fold, killing the mutant dropping the trials
from the re-derivation's counts. F3: call 10, pinned in tests/cpp and by
the gp record at counts against its rows written out. Gates on that code,
--preclean installs: tests/cpp plain and under ASan and UBSan; the full
suite, 43250 results, 0 failed; the trio `--bitwise --strict-coverage` on
the reference build, 55, 19 and 11 scenarios identical, none skipped; the
four reproducibility files, 27 results; the exact gates in quick, 49 of
49; lint 0, air, codoc, win-drift, doc-freshness, verify-anchors.
