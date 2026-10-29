# Forests with no basis share one prior total

Status: IN PROGRESS 2026-10-29: slice B, `blocks()` per group, and slice 1, the shared total's
engine and bridge, built on its branch; slice 2, the R surface that reaches it, not built. The plan is [forests-shared-total.md](../plans/forests-shared-total.md), which carries the steps,
the gates and the calls made in planning; this note is the standing statement of the rule. Rulings: dec-B292,
dec-B503, dec-B505 and dec-B506, in docs/decisions.md.

Several forests with no basis share one prior total, each an equal part of its variance whatever its tree count
(dec-B505); a forest partitioned by `blocks()` gives each of its groups an equal part of its own prior the same
way (dec-B506), so the formula of several such forests and the matching `blocks()` call are one model at every
count. This note states the `blocks()` half, which is built; the forests' half lands with slice 2.

## `blocks()`: each group an equal part

A forest of m trees has prior sd s at every point: each leaf is N(0, s^2 / m), and a point's fit sums one leaf
from each tree. Partitioned into G groups of m_g trees, the forest gives each group 1 / G of the variance: a
leaf of a tree in group g is N(0, s^2 / (G m_g)), so group g's function has sd s / sqrt(G) and the forest's
sum keeps sd s. A tree's leaf scale is the forest's times f_g = sqrt(m / (G m_g)). `trees.per.group` says how
complex each group's function may be, not how much of the response it explains, the reason the maintainer
gave for the forests (dec-B505: "That doesn't mean that each forest shouldn't attempt to capture the same
amount of residual variation."). Before this, every tree took s^2 / m, so a group's part followed its count.

With equal counts f_g is 1 and nothing changes. At the default split the change is small, 75 trees over two
groups of 38 and 37 moving each group's prior sd by under 1 percent; at stated uneven counts it is the
prior asked for, 5 and 15 trees of 20 moving the groups from 0.50 and 0.87 of s to 0.71 each.

How the engine holds it, by symbol:

- [`installBlockMasks`](../../src/bartcore/chain.hpp) stores the factors (`Forest::blockScaleFactor`, one per
  group) only where the groups' counts differ; a forest with equal groups or none stores nothing and runs the
  code it ran before, operand for operand. A copy, a reload and a warm start rebuild the factors from the
  counts at creation and place trees by slot, the slot's group giving the scale; no state carries them.
- A tree's update reads its group's scale: [`TreeLeafScale`](../../src/bartcore/combiner.hpp) puts it on the
  forest's leaf for the length of one tree's moves and draw, in the sweep, a grow-from-root sweep and a prior
  leaf draw, so the leaf models and the rule scan see it with no change of their own.
- The reads outside a tree's update take the tree's factor by name: the k statistic sums (mu / f_g)^2, so the
  k draw on the forest's scale is exact; the level step ([`Chain::drawLevelShift`](../../src/bartcore/chain.hpp))
  weighs each leaf by its own tree's sd; the multinomial grand-level shift
  ([`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp)) weighs each leaf's count and
  value by 1 / f_g^2. The probit rescaling step reads no leaf scale.

Gated by benchmarks/R/blocks-exact.R: three trees in groups of two and one, the 8 structures enumerated and
the leaves integrated in closed form, at a fixed k and under a drawn k; the per-tree prior misses its fits by
over 300 standard errors.
