# quantile-grid-spread: quantile split points cover the whole column

Status: PLANNED.

agent: opus implementer, one; opus reviewer.
rng: POSTERIOR-CHANGING for a fit with `useQuantiles = TRUE` (or `usequants = TRUE`) on a column with more
distinct values than its cut count plus one, and for any quantile refresh onto more distinct values than the
held count plus one: the set of split points is part of the tree prior. NEUTRAL for every other fit, the default
uniform grid included.
window: pre-release (dec-B235).
budget: ~450 lines (C++ ~50, tests/cpp ~150, tinytest ~80, design note, manual and NEWS ~120, baseline and
records ~50). Plans have run 1.5-2x low.

## Goal

A quantile grid's split points are spread evenly over the midpoints between a column's distinct values, at
creation and at refresh, on every entry point, `bartBT` included, so no part of a column's range is left without
a split point.

## Context

- [`ColumnStore::finishQuantileGrid`](../../src/bartcore/data.hpp) thins U distinct finite values to at most m
  cuts, m the column's cut count, by the integer step `U / m` with offset `step / 2`, keeping the midpoints at
  `k * step + offset`. BayesTree's code has the same arithmetic and 0.9-34 copied it.
- When U is not a multiple of m the points stop short of the top. At m = 100 the share of distinct values above
  the highest cut is 33 percent at U = 150, 50 at 199, 20 at 250, 33 at 299, 10 at 999 and 1.5 at 10099; it is
  near zero at 200, 300 and 1000. No tree can split there.
- A refresh ([`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp), reached by
  `setPredictor(updateCutPoints = TRUE)`) thins at the held count the same way: a column created with 11
  distinct values and refreshed onto 60 keeps the 10 lowest midpoints. 0.9-34 warned there ("ignoring extra
  quantiles"); the branch is silent.
- Measured at the default cut count, 200 trees: on Friedman's function the old rule's test RMSE was 1.97 against
  1.56 for evenly spread midpoints at 150 rows and 2.31 against 1.22 at 199 rows, with no difference from 300
  rows; on the Boston housing data, 20 splits of 405 and 101 rows, 3.60 against 3.12 (paired difference 0.48,
  standard error 0.09).
- `useQuantiles` is taken by `bart`, `dbarts`, `dbartsControl`, `xbart` and `rbart_vi`; `usequants` by `bartBT`,
  by a BayesTree-spelled `bart` call (forwarded to `bartBT`) and by `pdbart` and `pd2bart` (translated to
  `bart`'s name). All build their grid in the one store.

## The rule

For a column with U distinct finite values and cut count m:
- U - 1 <= m: every midpoint, as today.
- otherwise: cut k, for k = 0 to m - 1, is the midpoint between sorted distinct values i and i + 1 with
  i = floor((2k + 1)(U - 1) / (2m)), in integer arithmetic.

The cuts are strictly increasing, are midpoints (never an observed value), and leave at most
floor((U - 1) / (2m)) + 1 distinct values beyond either end cut. A refresh applies the same rule at the held count.

## Constraints

- A column with U - 1 <= m gets the grid it gets today, bit for bit, at creation and at refresh.
- The uniform grid is untouched.
- The same rule on the dense, subset (view) and sparse-backed paths.
- No warning at a refresh: nothing is ignored, the held count is spread.
- The rule that a refresh onto fewer distinct values than the held count is infeasible
  ([`ColumnStore::cutsWouldRemainValid`](../../src/bartcore/data.hpp)) stays.
- Out of scope: raising a held count at a refresh (a column created with few distinct values keeps that count);
  it goes to TODO.

## Steps

1. The rule in [`ColumnStore::finishQuantileGrid`](../../src/bartcore/data.hpp) and whatever
   [`ColumnStore::fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp) and the refresh need to share it.
   tests/cpp, in test_data.cpp: the expectation that hard-codes the old thinning is rewritten; grids at
   U = m + 1, m + 2, 2m - 1, 2m, 3m - 1 and 10m are strictly increasing, are midpoints, reach both ends within
   the bound above, and equal the old grid at U = m + 1; a refresh from 5 onto 50 and from 11 onto 60 distinct
   values spreads the held count, on the dense, subset and sparse-backed paths.
2. tinytest: with `useQuantiles = TRUE` and a column of 150 and of 199 distinct values at the default cut
   count, the highest split point lies above the column's 0.98 quantile; `bartBT(usequants = TRUE)`, a
   BayesTree-spelled `bart` call and `bart(useQuantiles = TRUE)` build the same grid on the same data; a column
   with no more than m + 1 distinct values keeps today's grid.
3. Baselines. The equivalence scenario that uses quantiles is re-recorded on the reference build, with its
   MANIFEST row, as the manifest's earlier re-records were; every other scenario is expected identical and is
   carried. The compare against 0.9-34 is not re-recorded: its quantile row gains an entry in the list of
   explained differences.
4. A design note in docs/design/ (the defect, the rule, the measurements, the lineage, what does not move); the
   `useQuantiles` and `usequants` help; NEWS, under changes against 0.9-34.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes, clean
  under ASan and UBSan.
- The four seeded snapshot files pass unchanged on a reference build (none uses a quantile grid).
- The equivalence compare in statistical (z) mode against the previous baseline: every scenario but the
  quantile one identical; the quantile scenario's z statistics reported.
- Every exact gate in `.github/workflows/exact-gates.yaml` in quick mode: the ones that use quantiles have no
  more distinct values than cuts plus one, so none is expected to move.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean.
