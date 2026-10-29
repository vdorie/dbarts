# A single declared forest keeps to its columns

Status: LANDED 2026-10-06 (9e691f35). Plan: [single-forest-vars.md](../plans/single-forest-vars.md). Rulings: dec-B241
(one declared forest is the ordinary way to state a single forest) and dec-B254 (a forest's columns are structure)
in docs/decisions.md.

## The defect

`forest(vars = )` restricts a forest to some of the predictor columns, and the manual says any forest may be
restricted. On a single declared forest the argument was dropped unread: `dbarts(x, y, forests =
list(forest(vars = c("a", "b"))))` on three columns drew exactly what the call without `vars` draws. Measured on
150 rows whose response leans on the third column, 20 trees: 4463 splits on the third column in 300 sweeps, in
every sweep, with or without `vars`. A name that is no column and an empty selection were accepted too.

Every other forest already had the mechanism. A forest of a multi-forest fit and the variance forest each
carry a 0/1 byte per column, and every tree reads its forest's bytes when it lists the columns a node may split
on ([`Tree::collectAvailableVariables`](../../src/bartcore/tree.hpp)). The single mean forest was built from
[`SamplerOptions`](../../src/bartcore/chain.hpp), which had no such list, and
[`resolveSamplerSpec`](../../R/spec.R) resolved `vars` only where there were several forests.

## The rule

A single declared forest restricted to columns C is the fit on `x[, C]`: the same trees and the same draws under
the same seed and the same residual scale estimate, with zero split counts reported for the other columns. The
residual scale estimate still comes from a linear fit on every column, since it belongs to the data and not to
a forest, so the identity is stated at an equal estimate.

The rule has no exception by tree prior. Under DART the Dirichlet is laid over C: its dimension, the default
`rho` and the concentration grid use the number of allowed columns, only those columns are drawn, and every
other column reports probability 0 ([`DartPrior::initialize`](../../src/bartcore/model.hpp),
[`DartPrior::update`](../../src/bartcore/model.hpp)). Under `split.probs` the caller's ratios hold among the
allowed columns; a vector that gives none of them a positive probability states no ratios and is refused by
name, at creation and by `setModel`.

One gap is left open on purpose. A state is installed as it is given, so a state that carries a DART
probability on an excluded column - an edited one, or an unrestricted DART donor's whose live trees happen to
lie within the list - reports that probability until the next Dirichlet draw puts it back at 0. Under the
default delay that draw is at most half the burn-in away: the count of updates already skipped rides the
state. Zeroing at install was not built: the split weights are state and the list is model, a state in the
sampler's own units installs bit for bit, and no draw reads the weight of an excluded column, the trees
listing allowed columns only and the concentration update summing over them.

Three things follow from the columns being part of what the model is, not a setting of it.

- The list is stored on the model object, as the `forest.columns` attribute beside the interaction and block
  constraints, so a copy, a reload and a sampler built from `dbartsSpec`'s pieces rebuild it. Naming every
  column restricts nothing and stores nothing. A multi-forest fit keeps every forest's columns where it kept
  them, on the control's forests attribute.
- [`dbartsSampler$setModel`](../../R/dbarts.R) refuses a model whose list differs from the sampler's, a model
  with none included, before anything is stored.
- A warm-start donor or a state holding a tree that splits outside the list is declined by `installTrees`
  and by `setState` and, forced, cut back at the split ([install-surface.md](install-surface.md)).

`blocks` beside `vars` partitions the allowed columns, on a single forest and on the first forest of several;
the first forest of two used to be refused there for not naming the columns it may not split on. The category
forests of a multinomial fit all take the one list
([`MultinomialForestSpec`](../../src/bartcore/combiner.hpp)).

A hazard fit's design carries a `period` column the caller did not supply. It stays allowed in every forest
whatever that forest's `vars` names (dec-B365; until 2026-10-08 a fit of several forests took each `vars` as
written): `vars` restricts the caller's columns, and a discrete-time hazard that cannot vary over periods is
not the model asked for. A forest constant over time has no spelling.

## Why not zeros in the split probabilities

Giving the excluded columns a split probability of zero needs no engine change. It was not taken.

- It makes structure ride on a parameter. Split probabilities are a parameter `setModel` may change, so a
  plain model handed to a forest confined by zeros is accepted and the confinement is gone: in the 300 sweeps
  after such a call 4099 of 8267 splits were on the excluded column.
- When this was decided, a zero confined a forest only while some allowed column could still be split on.
  Where none could, the split-variable draw proposed the first column available, whatever its probability,
  and the split was accepted. With two 0/1 columns held to one cut each beside a continuous column given
  probability zero, that column was split on 44425 times in 2000 sweeps of 20 trees, in 1995 of the sweeps;
  the column list gives 0 on the same fixture. A zero now excludes its column on every path
  ([A zero split probability](#a-zero-split-probability)).
- DART has no split probabilities to zero, and a multi-forest fit refuses `split.probs` and DART outright, so
  one argument would have had two mechanisms and a gap.
- The caller would read split probabilities they did not write on the stored tree prior.

## A zero split probability

A column of split probability zero is never split on: it is outside the forest's allowed set, as a column
`vars` leaves out is, and a node whose remaining columns all have probability zero cannot grow. The rule
rides the forest's column mask, so every site that asks which columns a node may split on - the prior's
growth probability and draws, the grow-from-root builder, the rule Gibbs and change moves, a state install -
agrees without a check of its own. The mask is the `vars` set intersected with the columns of positive
probability ([`Chain::buildColumnMask`](../../src/bartcore/chain.hpp)), built at creation and rebuilt when
`setModel` installs new probabilities, the block rows re-intersected with it
([`Chain::pointTreesAtColumnMask`](../../src/bartcore/chain.hpp)). A forest whose probabilities are all
positive, or that has none, keeps the `vars` mask alone, or none, and draws as before. Under DART nothing
changes: its weights are drawn and never zero, and DART and `split.probs` are exclusive.

[`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp) returns no column where none available has
positive probability, a case the mask keeps every caller from reaching.

What follows:

- `split.probs = c(0.5, 0.5, 0)` draws bit for bit what the same probabilities draw beside a `vars` naming
  the first two columns.
- [`dbartsSampler$setModel`](../../R/dbarts.R) refuses a zero on a column that a live tree of any chain
  splits on, naming the column, before any chain is written
  ([`Sampler::liveSplitOnZeroProbability`](../../src/bartcore/sampler.hpp)). It is refused rather than cut
  back: `setModel` has no forced form, and the trees would change under a call that states a prior. Kept
  draws are not read; they are replayed for prediction, not priced.
- A state, a copy or a warm-start donor that splits on such a column meets the mask: declined unforced, cut
  back forced, as for `vars`. The reported split probabilities are the ones stated.
- A fit that states a zero draws from the prior it states, which differs from earlier releases wherever the
  columns of positive probability run out along a path. On 150 rows with two 0/1 columns held to one cut
  each and a continuous third that carries the signal, 2000 sweeps of 20 trees under probabilities 0.5, 0.5
  and 0 give 27726, 27375 and 0 splits; 0.9-34 gives 40009, 39757 and 49623.

## Measurements

On the fixture above, 150 rows and three columns with `vars` naming the first two:

| fit | splits on the excluded column | splits on the allowed columns |
|---|---|---|
| before, 300 sweeps of 20 trees | 4463, in every sweep | 4337 |
| after, 300 sweeps of 20 trees | 0 | 8049 |
| after, 1000 sweeps of 75 trees | 0 | 102666 |
| after, allowed columns 0/1 held to one cut each, 2000 sweeps of 20 trees | 0 (30171 before; 44425 under zero split probabilities and no list) | 55069 |

The restricted fit is draw for draw the fit on the two-column matrix - fitted values, sigma and split counts
identical - for a gaussian and a probit response, with an interaction limit, with `split.probs` of 0.1, 0.3
and 0.6 against 0.25 and 0.75, and under DART, whose reported probabilities agree as well. No split falls
outside the list under the probit, logistic, ordinal, negative-binomial, aft and Student-t families, linear
and Gaussian-process leaves, a monotone column, a factor column, every category forest of a multinomial fit,
or after `copy`, a reload, `setPredictor`, `setData`, `sampleTreesFromPrior` and `growFromRoot`. Twenty-four
kinds of fit that state no `vars` on a single forest draw bit for bit what they drew before.

## Left out

When the control migration gives every model a forest record, the `forest.columns`
attribute moves into the first forest's record with the rest.
