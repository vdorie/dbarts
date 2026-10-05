# The quantile grid covers the whole column

Status: IMPLEMENTED 2026-10-05, landing not yet recorded. Plan: docs/plans/quantile-grid-spread.md. Ruling:
dec-B235 in docs/decisions.md.

With `useQuantiles = TRUE` (`usequants = TRUE` on `bartBT`) a numeric column's split points are midpoints between
its sorted distinct finite values. A column with no more distinct values than its cut count plus one takes every
midpoint. A longer column has to leave some out, and this note is about which.

## The defect

For U distinct values and a cut count m, the rule carried from BayesTree through 0.9-34 kept the midpoints at
positions `k * step + offset` for k from 0 to m - 1, with `step` the integer quotient of U by m and `offset`
half of it. The quotient is rounded down, so the m steps stop short, and everything they do not cover is at the
top of the column. No tree can split there. At the default 100 cuts:

| distinct values | share above the highest split point, old rule | new rule |
|---|---|---|
| 102 | 2.0% | 1.0% |
| 150 | 33.3% | 0.7% |
| 199 | 49.7% | 0.5% |
| 200 | 0.5% | 0.5% |
| 250 | 20.0% | 0.8% |
| 299 | 33.1% | 0.7% |
| 300 | 0.3% | 0.7% |
| 499 | 20.0% | 0.6% |
| 999 | 10.3% | 0.5% |
| 10099 | 1.5% | 0.5% |

A column whose distinct values are a multiple of the cut count loses nothing, which is how fits at 200, 500 and
1000 rows of continuous data never showed it.

A refresh, `setPredictor(updateCutPoints = TRUE)`, keeps the number of cuts a column holds and thinned to it the
same way. A column created with 11 distinct values holds 10 cuts; refreshed onto 60 it kept the 10 lowest
midpoints, the bottom sixth of the column. 0.9-34 warned there ("ignoring extra quantiles"); the rewritten engine
did the same thing silently.

## The rule

In [`ColumnStore::fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp): with M = U - 1 midpoints, midpoint i
lying between sorted distinct values i and i + 1, and c cuts to place,

    cut k = midpoint floor((2k + 1) M / (2c)),   k = 0, ..., c - 1

in integer arithmetic. That is the centre of the k-th of c equal shares of the midpoints. c is the smaller of m and
M at creation, as before ([`ColumnStore::finishQuantileGrid`](../../src/bartcore/data.hpp)), and the held count
at a refresh ([`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp)).

- When c equals M the index is k: every midpoint, the same doubles as before.
- The indices strictly increase, since consecutive ones differ by at least floor(M / c).
- A cut is a midpoint, never an observed value, so no row sits on one.
- At most floor(M / (2c)) + 1 distinct values lie beyond either end cut, and the two ends differ by at most one.
- The product fits 64 bits because a cut count is below 2^16. The floating form, floor((k + 1/2) M / c), gave
  the same indices at every size tried (m from 3 to 1000, U up to 40 m), but the integer form does not depend on
  that.

There is one rule for creation and refresh and for every storage. The dense, sparse-backed and view builds collect
the same set of values and call the same two functions; a row-subset view copies its parent's grid. A refresh
raises no warning, because nothing is ignored. Pinned by
[`testQuantileGridSpread`](../../tests/cpp/test_data.cpp) and, through `dbarts`, `bart`, `bartBT`, a
BayesTree-spelled `bart` call, `rbart_vi`, `pdbart` and `pd2bart`, by
["spreadMidpoints"](../../inst/tinytest/test-quantile-grid.R). `xbart` builds its grids in the same store and
returns no sampler to read one from.

## Measurements

Taken before the rule was in the engine, by installing its grid on a fresh sampler through `setCutPoints`: 200
trees, 100 cuts, 1000 draws after 500, one chain. "BART's rule" is the BART package's, installed the same way.

Friedman's function, ten predictors, noise standard deviation 1; root mean squared error against the true
function at 2000 fresh points, mean of five seeds:

| rows | old rule | new rule | BART's rule | uniform grid |
|---|---|---|---|---|
| 150 | 1.971 | 1.562 | 1.580 | 1.569 |
| 199 | 2.307 | 1.216 | 1.254 | 1.293 |
| 300 | 1.026 | 1.053 | 1.042 | 1.054 |
| 500 | 0.851 | 0.880 | 0.890 | 0.893 |
| 1000 | 0.684 | 0.676 | 0.677 | 0.678 |

At 150 and 199 rows the old rule's error is 1.26 and 1.90 times the new rule's. At 300, 500 and 1000 rows, where
the old rule loses nothing, the four grids are within 0.05 of one another.

Boston housing (MASS), 20 random splits into 405 and 101 rows, holdout root mean squared error:

| | old rule | new rule | BART's rule | uniform grid |
|---|---|---|---|---|
| mean | 3.598 | 3.121 | 3.145 | 3.212 |
| paired difference from the old rule | | -0.477 | -0.452 | -0.385 |
| its standard error | | 0.086 | 0.086 | 0.096 |

On the full data the old rule left 10.5% of the distinct values of `rm`, 12.3% of `lstat`, 16.0% of `age` and
16.2% of `black` above the highest split point.

## Lineage

BayesTree's mbart.cpp has the step-and-offset arithmetic and pgbart shares it; dbarts copied it and shipped it
through 0.9-34. The BART package takes type-7 quantiles of the data at equally spaced probabilities, which weights
by rows rather than by distinct values, can put a split point on an observed value and can repeat one; with fewer
distinct values than cuts it too takes every midpoint. LightGBM splits at midpoints between bins of equal count.
XGBoost and scikit-learn also cover the whole column. None of these leaves the top of a column without a split
point. The new rule keeps what this package has always documented, midpoints between distinct values, and changes
only which are kept; in the measurements above it and BART's rule cannot be told apart.

## What does not move

- The uniform grid, which is the default.
- A quantile column with no more distinct values than its cut count plus one, the single cut of a constant or
  fully missing column included.
- The number of cuts: the smaller of m and U - 1 at creation, fixed afterwards. A refresh does not raise it, so a
  column created with few distinct values keeps that many cuts when refreshed onto many. Raising it is a separate
  question and is not taken here.
- The refusal of a refresh onto fewer distinct values than the held count
  ([`ColumnStore::cutsWouldRemainValid`](../../src/bartcore/data.hpp)).
- Factor columns: an ordered factor's grid is its level midpoints and an unordered one has none.
- A grid installed through `setCutPoints`, and the grid a saved state carries: a state stored under the old rule
  restores with the split points it was stored with.
- The sampler. Given a grid the draws come from the same code, so of the equivalence baseline's 55 scenarios only
  the one that thins a quantile grid moves.

## Against 0.9-34 and BayesTree

A fit that asks for quantiles on a column with more distinct values than cuts plus one now differs from 0.9-34's
and from BayesTree's, `bartBT` included: leaving `bartBT` on the old rule would have kept the flaw for the callers
who write `usequants = TRUE`, and made that spelling behave two ways, since `pdbart` passes it to the modern
function. The 0.9-34 comparison records its quantile row as an explained difference:
[What differs, and why](../plans/classic-compare.md#what-differs-and-why).
