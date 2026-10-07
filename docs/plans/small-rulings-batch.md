# small-rulings-batch: monotone words abbreviate, setCutPoints sorts, the joint update refuses numbers for a factor

Status: LANDED 2026-10-07 (3a89f10f to 9d1d6645; dec-B285, dec-B286, dec-B288).

agent: sonnet implementer, one (R, one bridge message, tinytest, manual); opus reviewer.
rng: NEUTRAL, bit for bit, for every call accepted before and still accepted. A call that is refused now
draws nothing; a grid out of order, which was refused, is sorted and draws as the sorted grid does.
window: pre-release, before the merge to main.
budget: ~350 lines.

## Goal

A direction word of `monotone()` matches as base R matches a choice. `setCutPoints` sorts a grid out of
order and refuses a repeated point by name. The joint row update refuses numbers for a factor column.

## Context

- [`parseMonotoneSign`](../../R/model.R) is the one place a direction is read, for `monotone()` and for the
  plain vector alike, through `resolveMonotone`.
- [`bartcoreSamplerSetCutPoints`](../../R/bartcore.R) hands every entry to the bridge, which holds the
  grid to [`cutGridIsValid`](../../src/bartcore/data.hpp) unless it is the grid the column holds.
- [`codeJointColumnUpdate`](../../R/bartcore.R) codes the joint form's values;
  [`codeCategoricalColumnUpdate`](../../R/bartcore.R) holds the words `setPredictor` uses for a number given
  to a factor column.

## The rules

1. dec-B288. A direction word is any unique abbreviation of "increasing" or "decreasing", case-sensitive,
   by `pmatch` with `duplicates.ok = TRUE`, a miss refused with the existing message. The numbers 1, -1 and
   0 and the strings "1", "-1" and "0" are as they were. The help lists the full words and says nothing
   of abbreviation.
2. dec-B285. A numeric grid out of order is sorted and taken. A grid with a repeated point is refused,
   the message saying a point may appear once and to give a denser grid near the value for more splits
   there. The grid the column holds, bit for bit and in whatever order it is given, is taken as it is,
   whether or not it repeats a point. A grid that is not numeric, or holds `NA` or `NaN`, is refused.
   What a restore does with repeated points is not touched here.
3. dec-B286. Numbers given for a factor column of the joint update are refused in `setPredictor`'s words;
   a factor or labels are matched to the levels as before. Differing level tables across samplers refuse,
   and a factor or non-numeral text for a numeric column refuses, as before.

## Constraints

- No engine change: the bridge's grid check is the only C++ touched, and a restore's handling of repeated
  points is another item.
- The numbers 1, -1 and 0 and their strings stay directions, the prior's name stays under `match.arg`, and
  the help says nothing of abbreviation.
- A call accepted before and still accepted draws as before.
- `docs/decisions.md` and the root TODO are the coordinator's.

## Steps

1. `parseMonotoneSign` tries the codes, then `pmatch` on the two words; `MONOTONE_DIRECTION_CODES` no
   longer holds the words.
2. The R method sorts a grid with no `NA`, since `sort()` would drop one; the bridge refuses an `NA` or
   `NaN` and a repeat in two messages. The held-grid test comes after the sort, so the held grid is taken in
   any order it is given.
3. `codeJointColumnUpdate` refuses a number for a categorical column after the level tables are compared; a
   column held as a factor in one sampler and a number in another says so.
4. The help of `monotone`, `setCutPoints` and `updatePredictorPerObservationJointly`, and the method's
   docstring, state the rules. A first missing value is declined by no joint-form call: a factor takes a
   missing label only where the column holds one, and a numeric column's first missing value breaks no
   order.

## Tests

- [test-monotone.R](../../inst/tinytest/test-monotone.R): "i", "inc", "d", "dec" equal the full word; "Inc",
  "increasingly" and "" are refused; a vector mixing "inc" and 0, named and positional.
- [test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R): an unsorted grid gives the sampler
  state and next draws of the sorted one, by column, by list and by data frame; the repeat and `NaN`
  refusals; a constant column's own grid handed back is accepted.
- [test-joint-update-factor.R](../../inst/tinytest/test-joint-update-factor.R): numbers refused for an
  unordered and an ordered column, one sampler and two, with a missing number and a non-code; every label
  case kept; a numeric column moved by numbers (kept).

## Verification

Install into a private library with `--preclean`; the full tinytest suite with `at_home = TRUE`; the lint
set of `docs/plans/README.md`; `R CMD check --as-cran`; the bairrtt and stan4bart suites against the build.

## Landing note

Landed 2026-10-07 as 3a89f10f to 9d1d6645. The independent review, told to refute, found the code of
all three rulings sound under every probe and the faults in the documents: the help still described
the joint form taking level codes, and two landing notes had been rewritten; both corrected before
landing. What the review established about the joint form and a monotone sampler: a factor column's
first missing value cannot be brought by label, R refusing a missing label where the held column has
none, so the engine's row-by-row decline is reached from R only by a whole matrix of codes with
forceUpdate = FALSE and is tested in tests/cpp; a numeric column's first missing value is taken. Left
for the backlog: the refusal of a missing label on a factor column of 70 levels says its training
values had none when they had two (root TODO, missing-label-refusal-text). Gates on a clean copy of
the rebased branch: install, tests/cpp 351, the tinytest suite at home 16445 results and none
failed, lintr, air, rc-codoc, win-drift, anchors, build and check with the Date note alone. On the
same tree bairrtt's suite gave 207 results and stan4bart's 582, none failed.
