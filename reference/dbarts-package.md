# Discrete Bayesian Additive Regression Trees Sampler

Fits Bayesian additive regression trees (BART; Chipman, George, and
McCulloch 2010), a Bayesian “sum-of-trees” model in which each tree is
held to be a weak learner by its prior. The package offers two ways in,
described below, and provides a BayesTree-compatible interface.

What distinguishes dbarts from other BART implementations is that its
sampler is *mutable*: predictors, response, offset, and weights can be
swapped *between* MCMC draws, so that BART can be a conditional model
inside a larger Gibbs or Metropolis-Hastings sampler rather than only a
standalone fit.

## Fitting a model

For an ordinary fit, use one of the fitting functions. Each runs the
sampler to completion and returns posterior draws.

- [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) - the
  current interface, taking a formula or matrices and reaching the full
  feature set, including every response `family`. This is the one to
  reach for first.

- [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) - the
  BayesTree-compatible interface, kept at its historical defaults (200
  trees, one chain, factors expanded to indicator columns, binary
  responses probit, missing data rejected).

- [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) -
  crossvalidates over `k`, `power`, `base`, and the tree count.

- [`pdbart`](https://vdorie.github.io/dbarts/reference/pdbart.md) -
  partial dependence plots for one or two variables.

Fits support
[`predict`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`fitted`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`residuals`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`plot`](https://vdorie.github.io/dbarts/reference/bartBT.md), and
[`summary`](https://vdorie.github.io/dbarts/reference/summary.bart.md);
see [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
‘Generics’ section for what each returns and on which scale.
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
turns a survival fit into survival-curve draws.

## Driving the sampler directly

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) builds a
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
without running it. The sampler is a reference class whose methods -
`run`, `setResponse`, `setPredictor`, `setOffset`, `setWeights`, and the
rest - map onto the underlying C++ engine, so an outer sampler can
alternate its own draws with BART's.
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
and
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
build its data and control objects,
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
its priors, and
[`samplePriorPredictive`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
draws from the prior for calibration before any fitting.
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
resolves the same specification without constructing a sampler, for
packages that embed dbarts and supply their own design matrix.

## Model features

Reached through
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) unless
noted.

- Response families (`family`): `"gaussian"`, `"probit"`, `"logistic"`,
  accelerated failure time (`"aft"`), discrete-time survival hazard
  (`"hazard"`, `"hazard.logistic"`), multinomial (`"multinomial"`,
  `bart` only), ordered categorical (`"ordinal"`), negative-binomial
  counts (`"nbinom"`), and semicontinuous two-part
  (`"hurdle.lognormal"`, `bart` only). `"auto"`, the default, resolves
  the family from the response.

- Outlier-robust Student-t errors (`family = student(...)`) and a
  heteroscedastic variance forest (`variance`).

- Structural constraints: monotonicity (`monotone`), interaction limits
  ([`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)),
  and block-additivity
  ([`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)).

- Leaf models: constant, linear, or Gaussian-process (`node.prior`; see
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)).

- Variable selection by the DART prior (`tree.prior = dart()`).

- Predictors: missing values incorporated in place unconditionally,
  unordered factors split on level subsets (`factors`), and sparse or
  mixed dense/sparse input (`Matrix::dgCMatrix`,
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)).

- Warm starts from a previous fit (`warm.start`) or by XBART-style
  grow-from-root (`n.grow.sweeps`).

## Reproducibility

Every chain runs its own random number generator, so results never
depend on the thread count and sampling does not advance R's stream.
Calling [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand
suffices; passing `seed` instead makes a fit reproducible without
touching R's stream at all. See
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
‘Reproducibility’ section.

## Memory

Two arrays of 4 bytes per observation per tree per chain dominate what
the sampler itself holds, so its footprint is linear in the number of
observations, the number of trees, and the number of chains, and
everything else is a rounding error beside them. The predictor store
costs 2 bytes per cell once, whatever the chain count, because every
chain reads the same store. A test set costs its own store, 2 bytes per
test cell, plus 16 bytes per test row per chain.

Worked examples, for the sampler alone (the engine and its predictor
store, not the returned draws):

|                  |                |           |            |             |
|------------------|----------------|-----------|------------|-------------|
| **observations** | **predictors** | **trees** | **chains** | **sampler** |
| 10,000           | 10             | 75        | 1          | 7 MB        |
| 100,000          | 20             | 200       | 1          | 170 MB      |
| 100,000          | 20             | 200       | 4          | 660 MB      |
| 1,000,000        | 50             | 200       | 1          | 1.8 GB      |

A large-\\n\\ fit is usually limited not by the sampler but by the
training predictions it returns. `yhat.train` is 8 bytes per observation
per draw per chain, and two copies of it are live while the fit is
packaged - the one the sampler wrote and the one the fit returns - so
the 100,000-observation four-chain fit above at the default 500 draws
peaks near 3.2 GB of returned array against 660 MB of sampler. The lever
is `keepTrainingFits = FALSE` - an argument of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and of
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md),
and `keeptrainfits` on
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) - which
drops the array and the fitted means taken from it;
[`predict`](https://vdorie.github.io/dbarts/reference/bartBT.md) on a
`keepTrees` fit recovers them for whichever rows are wanted. A per-draw
`callback` (see
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) goes
further alongside `keepFits = FALSE` by never allocating the array at
all, taking the 100,000-observation four-chain fit above from about 3891
MB to about 698 MB, and a 1,000,000-observation one-chain fit at 50
predictors from about 10560 MB to about 2576 MB.

Two options cost more than their names suggest. `keepTrees` keeps every
draw's trees, about 24 bytes per node per tree per draw per chain, which
is 23 MB per chain at 200 trees and 500 draws and grows with the tree
count and the draw count but not with \\n\\. A designated-covariate leaf
(`node.prior = linear(...)` or `gp(...)`) caches sufficient statistics
keyed on leaf membership; measured, that cache runs about 13 to 14 bytes
per observation per tree per chain, and at 100,000 observations, 200
trees and one chain it was 276 MB - larger than everything else the fit
holds put together. It also replaces one of the two 4-byte arrays above
with an 8-byte one.

## Using dbarts from another package

Calling the R interface needs no special linkage. A package that drives
the sampler from C should declare `LinkingTo: dbarts` and use the flat C
API declared in the installed header `dbarts/dbarts.h`, whose entry
points are reached through `R_GetCCallable` and versioned by the
two-component handshake `DBARTS_C_API_MAJOR`/`DBARTS_C_API_MINOR`, with
`dbarts_apiHash()` as an opt-in exact-ABI check, which moves on additive
releases as well. That header creates no sampler and names no R type:
the sampler is built in R -
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md) -
and the handle the entry points take is the address in that object's
external pointer, read with `R_ExternalPtrAddr`. See
[`dbarts-embedding`](https://vdorie.github.io/dbarts/reference/dbarts-embedding.md).
The C++ interface of releases before 1.0-0, `dbarts/R_C_interface.hpp`,
has been removed.

## References

Chipman, H., George, E., and McCulloch, R. (2010) BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, **4**(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285) .

## Author

Vincent Dorie <vdorie@gmail.com>, with Hugh Chipman and Robert
McCulloch. See `citation("dbarts")` for how to cite the package.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) to fit a
model, [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
to build a sampler, and the package vignettes:
[`vignette("gibbs_sampler_mixture_model", package = "dbarts")`](https://vdorie.github.io/dbarts/articles/gibbs_sampler_mixture_model.md)
and
[`vignette("working_with_saved_trees", package = "dbarts")`](https://vdorie.github.io/dbarts/articles/working_with_saved_trees.md).
