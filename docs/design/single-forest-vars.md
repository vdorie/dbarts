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
- A state or a warm-start donor holding a tree that splits outside the list is refused with the messages a
  restricted variance forest gives; the same predicate judges both
  ([`Chain::columnMaskStateFeasible`](../../src/bartcore/chain.hpp)).

`blocks` beside `vars` partitions the allowed columns, on a single forest and on the first forest of several;
the first forest of two used to be refused there for not naming the columns it may not split on. The category
forests of a multinomial fit all take the one list
([`MultinomialForestSpec`](../../src/bartcore/combiner.hpp)).

A hazard fit's design carries a `period` column the caller did not supply. It stays allowed whatever `vars`
names: `vars` restricts the caller's columns, and a discrete-time hazard that cannot vary over periods is not
the model asked for.

## Why not zeros in the split probabilities

Giving the excluded columns a split probability of zero needs no engine change. It was not taken.

- It makes structure ride on a parameter. Split probabilities are a parameter `setModel` may change, so a
  plain model handed to a forest confined by zeros is accepted and the confinement is gone: in the 300 sweeps
  after such a call 4099 of 8267 splits were on the excluded column.
- It confines a forest only while some allowed column can still be split on. Where none can,
  [`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp) proposes the first column available,
  whatever its probability, and the split is accepted. With two 0/1 columns held to one cut each beside a
  continuous column given probability zero, that column is split on 44425 times in 2000 sweeps of 20 trees,
  in 1995 of the sweeps; the column list gives 0 on the same fixture. On the default grid of 100 cuts a 0/1
  column rarely runs out and the zeros mostly hold: no such split on five fixtures of six, and 3514 of them,
  in 1335 of the 2000 sweeps, on the sixth.
- DART has no split probabilities to zero, and a multi-forest fit refuses `split.probs` and DART outright, so
  one argument would have had two mechanisms and a gap.
- The caller would read split probabilities they did not write on the stored tree prior.

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

A zero split probability still does not exclude a column once the columns with positive probability run out,
with `split.probs` alone and no `vars`; the fix belongs with `setModel`'s treatment of a zero probability on a
column the trees use. When the control migration gives every model a forest record, the `forest.columns`
attribute moves into the first forest's record with the rest.
