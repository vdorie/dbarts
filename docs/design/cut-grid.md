# A cut grid holds each point once

Status: PLANNED (built on a branch; the plan's status is the record).
Plan: [repeated-cut-restore.md](../plans/repeated-cut-restore.md). Rulings: dec-B285, dec-B297, dec-B298,
dec-B299, dec-B300, dec-B311 and dec-B312 in docs/decisions.md.

A numeric predictor's cut grid is the set of thresholds a tree may split it at. The sampler picks uniformly
among grid positions, in the tree prior and in the proposals, and a stored split names its threshold by value.
Both assume a grid is a set. It was not: the uniform rule placed `n.cuts` points whatever the range, so a
constant column held `n.cuts` copies of its value and a column narrower than the points resolve held a few
values many times; the quantile rule could round two midpoints to one double. On such a grid a threshold was
likelier for being repeated, a constant column looked splittable, and a restore put a split on the first
position holding its value, so the restored chain did not continue as stored.

## The rule

- Every grid a sampler derives, at creation, at `setData` and at a refresh, holds each point once
  ([`ColumnStore::dropRepeatedCuts`](../../src/bartcore/data.hpp)). The points are the ones the rule placed,
  with equal neighbours dropped: the uniform rule's evenly spaced points, the quantile rule's midpoints
  ([`ColumnStore::fillCutsOverRange`](../../src/bartcore/data.hpp),
  [`ColumnStore::fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp)). A constant column, or one with one
  finite value beside missing ones, holds one point; a column with no finite value one point at 0.
- A refresh derives the grid a creation would for the values
  ([`ColumnStore::deriveNumericCuts`](../../src/bartcore/data.hpp)): up to the `n.cuts` asked for at creation,
  whatever number of points the column held, and fewer where the values supply fewer. It never fails for
  want of distinct points.
- Every grid a sampler accepts holds each point once: `setCutPoints` refuses a repeat, and a state or a
  warm-start donor whose grid repeats a point is refused in the bridge, naming the column and the value,
  before anything is installed ([`cutGridIsValid`](../../src/bartcore/data.hpp) is the engine's backstop).
- A stored split's value then names one position, so a restore by value is exact with nothing more stored.
  A stored tree can still hold a split outside the interval its ancestors leave, two splits stacked on one
  value, if no sampler of this version wrote it; the build marks it and the pass that merges an empty side
  merges it, the install reporting itself altered
  ([`Tree::holdsSplitOutsideInterval`](../../src/bartcore/tree.hpp)).

## Where the splits go when a grid is replaced

A refreshed grid may differ from the one it replaces in its points and in their number, and so may a grid given
to `setCutPoints`. The caller chooses ([`SplitPlacement`](../../src/bartcore/data.hpp)):

- by position, `updateCutPoints = "position"` and the default of `setCutPoints(splits = )`: a split keeps its
  position. With the number of points unchanged nothing is touched, which is what 0.9-34 did. With it changed
  position i of n, from 0, goes to floor((2 i + 1) m / (2 n)) of m, the point under the centre of its share of
  the old grid ([`Tree::rescaledSplitIndex`](../../src/bartcore/tree.hpp)).
  At `setCutPoints` 0.9-34 kept the index on a grid of another length and merged a split past its end, so
  that call draws differently from the release (dec-B311).
- by value, `"value"`: a split moves to the new point nearest its old threshold, the move `setData` and a warm
  start from another grid make ([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)).

Either way a split stays inside the interval its ancestors leave. One with no point left there is merged when
the change is forced. An unforced refresh then declines instead, and the grid, its count and every position
moved are put back ([`Tree::tryMapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp),
[`Chain::undoSplitMoves`](../../src/bartcore/chain.hpp)).

## What moves

Measured on a shipped build of the branch before this change, `n.cuts = 100`, 200 rows.

| column | uniform rule, points then | now | quantile rule, then | now |
|---|---|---|---|---|
| five adjacent doubles | 100, 5 distinct | 5 | 4, 3 distinct | 3 |
| constant | 100 copies | 1 | 1 | 1 |
| constant beside missing values | 100 copies | 1 | 1 | 1 |
| 0/1 | 100 | 100 | 1 | 1 |
| six values | 100 | 100 | 5 | 5 |

Fits whose draws change: a column narrower than its points resolve, under either rule, and a constant column
beside missing values under the uniform rule; and any refresh whose grid now differs. A constant column with
no missing value goes from `n.cuts` copies to one point and draws the same, no split on it ever having been
accepted. Every other fit is unchanged bit for bit. Under the uniform rule a column of few values far apart,
a 0/1 column among them, holds `n.cuts` points today; that is what the grid is at this change, and the rule
for such a column is ruled otherwise for a later one (dec-B313: one point per gap between distinct values).

## Checks

Two arms of [bd-balance.R](../../benchmarks/R/bd-balance.R), `narrow` and `constmissing`, put the
birth/death gate on a column of each moved kind against its exact posterior, and each asserts the grid first.
