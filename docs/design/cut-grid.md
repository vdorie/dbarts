# A cut grid holds each point once

Status: LANDED 2026-10-08 (the plan's status is the record).
Plan: [repeated-cut-restore.md](../plans/repeated-cut-restore.md). Rulings: dec-B285, dec-B297, dec-B298,
dec-B299, dec-B300, dec-B311 and dec-B312 in docs/decisions.md. The default rule's one point per gap
([default-rule-per-gap.md](../plans/default-rule-per-gap.md)): dec-B313, dec-B406, dec-B409, dec-B410 and
dec-B411.

A numeric predictor's cut grid is the set of thresholds a tree may split it at. The sampler picks among grid
positions, uniformly or by the grid's weights, in the tree prior and in the proposals, and a stored split names
its threshold by value.
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
- Under the default rule (`useQuantiles = FALSE`) a numeric column with at least two and fewer than `n.cuts`
  distinct finite values holds one point in each gap between neighboring values, the quantile rule's grid for
  the same values, and the tree prior and every proposal choose a gap with probability in proportion to its
  width ([`ColumnStore::weighGapsByWidth`](../../src/bartcore/data.hpp)). Where every width is equal (a 0/1
  column, equally spaced values) the points are equally likely and the column carries no weights. A column
  with `n.cuts` or more distinct values keeps the evenly spaced points. Distinct values are counted as the
  quantile rule counts them: finite values only, -0 and 0 one value, a sparse column's implicit zeros once.
- The point in a gap is its midpoint, or, where the midpoint is not strictly below the upper value (two
  adjacent doubles, or a sum that overflows), the lower value, under both rules; a row at it does not exceed
  it, so the gap stays splittable.
- At a node whose interval holds positions a to b, position k has prior probability w_k / (w_a + ... + w_b),
  which every draw of a position, every count in an acceptance ratio and the rule_gibbs and grow-from-root
  candidates read ([`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp),
  [`drawCutPosition`](../../src/bartcore/model.hpp)). The weights are stored beside the grid as an increasing
  sequence whose differences are the weights ([`ColumnStore::cutMass`](../../src/bartcore/data.hpp)); the
  default rule stores the column's distinct values, so each weight and each sum is one exact difference.
- A grid given to `setCutPoints` is unweighted, except the grid the column holds handed back bit for bit,
  which keeps its weights. A state carries the weights as its `cutMass` attribute, one entry per column, NULL
  where the points are equally likely; a copy and a reload install them with the grid, and a state without
  the attribute installs every grid unweighted. A row-subset view (an xbart fold) copies them with the grid.
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

| column | uniform rule, points then | at this change | per gap | quantile rule, then | at this change | per gap |
|---|---|---|---|---|---|---|
| five adjacent doubles | 100, 5 distinct | 5 | 4 | 4, 3 distinct | 3 | 4 |
| constant | 100 copies | 1 | 1 | 1 | 1 | 1 |
| constant beside missing values | 100 copies | 1 | 1 | 1 | 1 | 1 |
| 0/1 | 100 | 100 | 1 | 1 | 1 | 1 |
| six values | 100 | 100 | 5 | 5 | 5 | 5 |

Fits whose draws change at this change: a column narrower than its points resolve, under either rule, and a
constant column beside missing values under the uniform rule; and any refresh whose grid now differs. A
constant column with no missing value goes from `n.cuts` copies to one point and draws the same, no split on
it ever having been accepted. Every other fit is unchanged bit for bit. The per-gap columns are the default
rule's grid since (below).

## One point per gap, weighted by width (2026-10-09)

Under the default rule every fit with a numeric column of at least two and fewer than `n.cuts` distinct
values draws differently: a 0/1 column, a small count, a sparse or zero-inflated column, any column on fewer
rows than `n.cuts`. Every other fit is unchanged bit for bit, the quantile rule's but on adjacent doubles.
0.9-34 placed 100 evenly spaced points on such a column; the BART package takes the midpoints, unweighted,
under the same boundary.

Weighting by width keeps the gap prior: a gap's share is its width over the node's summed widths, which is the
share of evenly spaced points falling in it to within one point in `n.cuts`. Equal per-gap shares, the ruling's
first reading, moved a zero/nonzero gap's share from about a third to one over the distinct values and cost
5.5% test RMSE on zero-inflated columns; weighted, the planning prototype measured no loss on any design tried
(paired test RMSE: zero-inflated columns +0.0006, se 0.0018; continuous columns on 30 to 90 rows within one
se). What changes is availability: a node may split a column while its interval holds a point, so a node
holding one value of a 0/1 column, which held points before, no longer counts the column as available; in one
tree on a 0/1 column under a noise response the share of draws holding a split falls from 0.177 to 0.131.

## Checks

Two arms of [bd-balance.R](../../benchmarks/R/bd-balance.R), `narrow` and `constmissing`, put the
birth/death gate on a column of each moved kind against its exact posterior, and each asserts the grid first.
Its `weighted` arm puts cells at 0, 1, 3 and 10 under the default `n.cuts`, asserts the grid {0.5, 2, 6.5},
and runs the enumeration with each cut weighted by its gap under five move mixtures (birth/death, change,
swap, perturb, rule_gibbs); the posterior with equally likely cuts is a total variation of 0.23 away.
[monotone-exact-enumeration.R](../../benchmarks/R/monotone-exact-enumeration.R)'s design `u1` puts the
monotone leaf on the same values. heteroscedastic-exact.R's part (b) and multinomial-exact.R's arm 8 weigh a
split on a 0/1 column by `base` alone, which holds only while the children hold no point; each asserts the
column's one point.

## Checks

Two arms of [bd-balance.R](../../benchmarks/R/bd-balance.R), `narrow` and `constmissing`, put the
birth/death gate on a column of each moved kind against its exact posterior, and each asserts the grid first.
