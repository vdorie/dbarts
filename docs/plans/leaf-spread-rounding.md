# leaf-spread-rounding: a leaf covariate whose spread is rounding is one value

agent: opus
rng: neutral for every column above the tolerance (the trio bitwise); a
column under it is refused or kept where it ran before
budget: ~60 engine and bridge lines, ~10 R lines, ~250 test lines, help

Status: BUILT 2026-11-07, not yet reviewed (03f9e1618)

## Goal

A linear or gp leaf covariate whose spread about its center, over the rows
in the likelihood, is under 64 times machine precision of its size,
`||z - center|| < 64 * DBL_EPSILON * ||z||`, is read as one value: refused
at creation with a message that names the test and says to drop the column
or subtract a constant from it, and treated as a column of one value by
every setter that re-derives the center and scale (dec-B654, dec-B655). A
column above the tolerance draws what it drew.

## Context

- `standardizationMomentsForColumn` (src/bartcore/data.hpp) takes the
  center and sd over the observed values of the rows in, membership only
  (`LeafRows`: positive weight and kept by the mask), and says where they
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
  counting once; a row's weight does not enter, as it does not enter the
  center and scale. Missing values are skipped.
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
  covariate "z" holds values that differ only by rounding over the rows in
  ..., its spread about its center under 64 times machine precision of its
  size (`||z - center|| < 64 * .Machine$double.eps * ||z||`), so its center
  and scale cannot be taken from them; drop it from linear() or subtract a
  constant from it, or change the rows in (...). A gp leaf's remedies lead
  with stating its lengthscale, as for one value. Same three doors.
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
  for one value: centered at the computed center and divided by 1, every
  training row then reading within rounding of zero.

## Calls (each with its alternative)

1. A wrapper both derivation paths call, the moments function unchanged;
   alternative, the test inside `standardizationMomentsForColumn`.
2. Rows in, unweighted; alternative, case-weighted sums, which the center
   and scale do not use.
3. The refined center inside the test; alternatives, the stored center
   (the test grows with n) or a compensated stored center (moves bits on
   columns above the tolerance).
4. Power-of-two scaling; alternative, plain sums of squares, overflowing
   above about 1e154.
5. A rounding column under a gp leaf with stated lengthscales centered at
   the computed center; alternative, its first observed value. Its spread
   mark (`standardizedColumnHasSpread`, read for the record's scale) is left
   as built, reading the gathered values, so the record shows scale 1 where
   a column of one value shows NA; alternative, a tolerance-aware mark.
6. The keep warning reworded for both cases; alternative, a separate text
   for the rounding case, carried through the bridge's keep answer.
7. NEWS: none. Linear and gp leaves are new in 1.0-0; 0.9-34 had none.

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

Built in 03f9e1618 as planned. The six warning pins in
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
