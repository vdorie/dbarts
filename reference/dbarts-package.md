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
