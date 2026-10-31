# Forests with no basis share one prior total

Status: IN PROGRESS 2026-10-30: slice B, `blocks()` per group, slice 1, the shared total's engine and
bridge, slice 2, the R surface that reaches it, and slice 3, the constant column, built on its branch; slice F, the
default flip, not built. The plan is
[forests-shared-total.md](../plans/forests-shared-total.md), which carries the steps, the gates and the calls
made in planning; this note is the standing statement of the rule. Rulings: dec-B292, dec-B503, dec-B505,
dec-B506, dec-B509 and dec-B513, in docs/decisions.md.

Several forests with no basis share one prior total, each an equal part of its variance whatever its tree count
(dec-B505); a forest partitioned by `blocks()` gives each of its groups an equal part of its own prior the same
way (dec-B506), so the formula of several such forests and the matching `blocks()` call are one model at every
count.

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

## Forests: each an equal part of the fitting function's prior

The sharers are the forests with no basis that state no `sd` and no `amplitude = student()`
([`leafPriorForests`](../../R/model.R) in R, [`forestSharers`](../../src/bartcore/combiner.hpp) in the engine,
each reading what the forests state; no flag passes between them). The total is the fitting function's leaf
prior, as `bart()`'s one forest has it (dec-B509): k 2 under gaussian and k drawn under chi(1.5, 2) under probit
and logistic by default, a stated `k`, `leaf.prior = normal(k = )` or a `chi()` on k, and `normal(sd = )` naming
the total's size. Each sharer is held at 1, and with T the total's sd, sharer j of K over m_j trees has sd
T / sqrt(K) and leaves N(0, T^2 / (K m_j)): [`forestMapEntry`](../../src/bartcore/combiner.hpp) gives it the
single-forest chain's node scale times sqrt(1 / K), so the sum has sd T at every point, whatever the counts.
One sharer alone beside other forests has the whole of the total, the k path with share 1 (dec-B509), its
coefficient held and its reader as one forest's: `prior.sd.of` is "forest", no group line prints. At a held k of
2 stating no sd it keeps the map's held entry operand for operand (dec-B290), which is the same prior, so a held
gaussian fit keeps its draws.

A forest with no basis is held unless it states a df ([`forestMapEntry`](../../src/bartcore/combiner.hpp) reads
the statement, isnan and not a test of finiteness). With `student()` and no `sd`, its median is the size the leaf
prior gives ([`Chain::leafPriorSize`](../../src/bartcore/chain.hpp)): a named sd, the node scale over a held k, or
over the chi prior's median k, 1.574 under probit and 2.855 under logistic at chi(1.5, 2). A flat prior on k gives
no size and is refused by name. bcf's mu states df 1 ([`expandForestSpecs`](../../src/bartcore/combiner.hpp)) so
the K = 2 fixtures keep their half-Cauchy. A leaf prior given to the fitting function where every forest with no
basis states its size has no forest and is refused by name, the default restated taken. The check behind the
flip (dec-B509) found a held prognostic forest not clearly worse for the effects, and f's coverage about half a
point lower than under the drawn coefficient.

- One k. A drawn k is drawn once for the group ([`Chain::drawSharedK`](../../src/bartcore/chain.hpp)), from
  every sharer's leaves each standardized by its own leaf scale, and reported on each sharer: `getK` and a
  fit's `k` give the group's.
- A forest that states `sd` has that size and is outside the share; the others share
  among themselves. `amplitude = student()` draws one forest's coefficient, so it is refused on a sharer and
  taken with `sd` stated or on a forest with no basis that shares with none.
- The tree count (dec-B503). The sharers that state no `n.trees` split the fitting function's as `blocks()`
  splits a forest, the first n mod K taking one more (75 over two is 38 and 37); a stated count is outside the
  split and decides how finely that forest divides its equal part, not its size
  ([`plainTreeCounts`](../../R/model.R)). The fitting function's `n.trees`, `power` and `base` are refused as
  given twice only where every forest with no basis states its own. `interactions` and `blocks` at the fitting
  function, a partition of one forest's predictors, are refused beside two or more forests with no basis.
- The probit rescaling step ([`Chain::drawForestRescale`](../../src/bartcore/chain.hpp)) scales the latents,
  every sharer's leaves and k by one factor over the group's leaf count, and reads every other forest's fit,
  a forest with a basis or one stating sd, as part of the offset, which its conditional already leaves
  unscaled; so it stays exact in every model whose sharers draw k.
- The reader. `getLeafPrior` gives a sharer the fitting function's leaf prior as `leaf.prior`, its part in
  force as `sd` (NA where k is drawn), `sd.stated` FALSE, amplitude `fixed()` and `prior.sd.of` "forest,
  shared". Print shows each forest's tree count, the total, and which forests share at what k. A fit records
  the group's k wherever the sharers stand, each sharer's held coefficient in `fixed$amplitude` whether or not
  `fixed()` is written; where k is drawn `summary` tabulates it and `extract(type = "leaf.prior.sd")` gives each
  sharer's `k.scale` over each draw.
- On a live sampler, with two or more sharers, a sharer's `sd` through `$setLeafPrior(forests = )` and the
  sampler's `leaf.prior` are refused by name; a size or total is stated at creation, and changing either live
  waits on the forest door. A lone sharer takes the `sd` and leaves the leaf prior, k pinned at 1, as creation
  stating it does (dec-B525's end state); the reader's own `leaf.prior` written back moves nothing, and
  `$setLeafPrior(leaf.prior)` takes only a restatement of the one in force.

Without a forest with a basis, `y ~ forest(x1) + forest(x2) + forest(x3)` is then the one-forest model
`y ~ forest(x1 + x2 + x3, blocks = blocks(list("x1", "x2", "x3")))` with the same trees per group, at every
count. Gated by benchmarks/R/forest-prior-exact.R's shared-total arms (held sharers integrated in closed form,
beside a held factor, beside a forest stating sd, at k 2 and 3 and at unequal counts, where an oracle at the
full total or at T / K misses its fits by over 100 standard errors), probit-k-scale-exact.R's
shared arms (two frozen sharers under one drawn k, alone, beside a held factor read as offset and with a user
offset) and forests-blocks-equivalence.R (the two spellings by Welch z over seeds under gaussian and probit, at
equal counts, at 38 and 37 and beside a stated count, and under gaussian with k drawn).

A numeric basis constant over the rows in the likelihood is refused at creation, sd stated or not, naming the
spelling with no basis (dec-B528; [`refuseConstantBasis`](../../R/forestBasis.R)). A forest with no basis stating
`sd` is that model; a later swap to such a column stays taken (dec-B346).
