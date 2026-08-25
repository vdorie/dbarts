# Crossvalidation For Bayesian Additive Regression Trees

Fits the BART model against varying `k`, `power`, `base`, and `n.trees`
parameters using \\K\\-fold or repeated random subsampling
crossvalidation, sharing burn-in between parameter settings. Results are
returned as an array of evaluations of a loss function on the held-out
sets.

## Usage

``` r
xbart(
    formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
    method = c("k-fold", "random subsample"), n.test = c(5, 0.2),
    n.reps = 40L, n.burn = c(200L, 150L),
    loss = c("rmse", "log", "mcr"), n.threads = dbarts::guessNumCores(), n.trees = 75L,
    k = NULL, power = 2, base = 0.95,
    split.probs = NULL, dart = FALSE, drop = TRUE,
    resid.prior = chisq, sigest = NA_real_,
    seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    missing = c("incorporate", "error"),
    node.prior = NULL, n.cuts = 100L, useQuantiles = FALSE, n.thin = 1L,
    storage = c("double", "single"), tree.prior = NULL)
```

## Arguments

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  following an analogous model description syntax as
  [`lm`](https://rdrr.io/r/stats/lm.html). For backwards compatibility,
  can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) matrix
  `x.train`. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).

- data:

  An optional data frame, list, or environment containing predictors to
  be used with the model. For backwards compatibility, can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) vector
  `y.train`.

- subset:

  An optional vector specifying a subset of observations to be used in
  the fitting process.

- weights:

  An optional vector of weights to be used in the fitting process. For a
  gaussian response, BART fits a model with observations \\y \mid x \sim
  N(f(x), \sigma^2 / w)\\, where \\f(x)\\ is the unknown function.
  Binary responses differ: a `"probit"` model does not support weights
  (a weighted probit has no tractable latent-variable form), except that
  weights identically 1 are treated as absent; a `"logistic"` model
  treats them as observation counts and so requires positive integers
  (its Polya-Gamma latent for a count \\w\\ is a sum of \\w\\ unit
  draws).

- offset:

  An optional vector specifying an offset from 0 for the relationship
  between the underlying function, \\f(x)\\, and the response \\y\\.
  Useful for all three response families, though each applies it
  differently: for a gaussian response, \\y = f(x) + \mathrm{offset} +
  \epsilon\\, a fixed component of the mean (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for
  the interaction with BART's internal range-scaling); for binary
  responses it enters the link directly, \\P(Y = 1 \mid X = x) =
  \Phi(f(x) + \mathrm{offset})\\ for a `"probit"` fit and the analogous
  model on the logistic link for `"logistic"`.

- verbose:

  A logical determining if additional output is printed to the console.

- n.samples:

  A positive integer, setting the number of posterior samples drawn for
  each fit of training data and used by the loss function.

- method:

  Character string, either `"k-fold"` or `"random subsample"`.

- n.test:

  For each fit, the test sample size or proportion. For method
  `"k-fold"`, is expected to be the number of folds, and in \\\[2,
  n\]\\. For method `"random subsample"`, can be a real number in \\(0,
  1)\\ or a positive integer in \\(1, n)\\. When a given as proportion,
  the number of test observations used is the proportion times the
  sample size rounded to the nearest integer.

- n.reps:

  A positive integer setting the number of cross validation steps that
  will be taken. For `"k-fold"`, each replication corresponds to fitting
  each of the \\K\\ folds in turn, while for `"random subsample"` a
  replication is a single fit.

- n.burn:

  Unlike the single-scalar `n.burn` of
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
  `bart2`, and `rbart_vi`, here it is one or two non-negative integers,
  specifying 1) the burn-in when a chain is freshly started against a
  data split and 2) the burn-in when moving from one parameter setting
  to another over the same split. Chains are never carried between data
  splits or folds - the held-out observations of one were training
  observations of the previous, so continuing a chain lets slowly-mixing
  settings score against data they have effectively seen.

- loss:

  Either one of the pre-set loss functions as character-strings (`mcr` -
  misclassification rate for binary responses, `rmse` -
  root-mean-squared-error for continuous response), `log` - negative
  log-loss for binary response (`rmse` serves this purpose for
  continuous responses), a function, or a function-evaluation
  environment list-pair. Functions should have prototypes of the form
  `function(y.test, y.test.hat, weights)`, where `y.test` is the held
  out test subsample, `y.test.hat` is a matrix of dimension
  `length(y.test)` \\\times\\ `n.samples`, and `weights` are an optional
  vector of user-supplied weights. See examples.

- n.threads:

  Replications are independent, and for `n.threads > 1` they are divided
  into approximately equal chunks and executed on that many parallel
  workers (a
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html) cluster).
  The default uses
  [`guessNumCores`](https://vdorie.github.io/dbarts/reference/guessNumCores.md),
  which should work across the most common operating system/hardware
  pairs.

- n.trees:

  A vector of positive integers setting the BART hyperparameter for the
  number of trees in the sum-of-trees formulation. See
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md).

- k:

  A vector of positive real numbers, setting the BART hyperparameter for
  the node-mean prior standard deviation. If `NULL`, the grid default of
  2 is used for every response family. Binary responses do not inherit
  `bart2`'s Chi hyperprior default: a hyperprior is not a grid, so
  taking it here would leave the `k` axis a single cell. Hyperprior
  crossvalidation not possible at this time. A hyperprior `k` is held,
  not swept, and is DRAWN every sweep in every cell, so the reported
  loss is computed under a shrinkage that moves within each fit rather
  than under the named value.

- power:

  A vector of real numbers greater than one, setting the BART
  hyperparameter for the tree prior's growth probability, given by
  \\{base} / (1 + depth)^{power}\\.

- base:

  A vector of real numbers in \\(0, 1)\\, setting the BART
  hyperparameter for the tree prior's growth probability.

- split.probs:

  Prior probabilities that a variable is used in a splitting rule, as in
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md). A
  single value or `NULL` yields the uniform default; a named or unnamed
  vector assigns per-column probabilities. Fixed for the whole
  crossvalidation, not part of the swept grid. Cannot be combined with
  `dart`.

- dart:

  Use the DART sparsity-inducing prior of Linero (2018) instead of a
  fixed `split.probs`. `TRUE` enables it with default hyperparameters; a
  prior created by `dbartsPriors$dart` supplies its own. The split
  probabilities are then sampled rather than fixed, while `power` and
  `base` are still swept.

- drop:

  Logical, determining if dimensions with a single value are dropped
  from the result.

- resid.prior:

  An expression of the form `chisq` or `chisq(df, quant)` that sets the
  prior used on the residual/error variance.

- sigest:

  A positive numeric estimate of the residual standard deviation. If
  `NA`, a linear model is used with all of the predictors to obtain one.
  Fitting functions (`xbart`,
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)/`bart2`,
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md))
  spell this `sigest`; sampler constructors
  ([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
  `dbartsSpec`) spell the same concept `sigma`.

- seed:

  Optional integer specifying the desired pRNG
  [seed](https://rdrr.io/r/base/Random.html). Results are reproducible
  for a fixed pair of `seed` and `n.threads`; the caller's random stream
  is left untouched when a seed is given. Without one,
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand suffices.

- factors:

  How factor columns in a data frame enter the model; as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).

- family:

  The response model, resolved as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md):
  `"auto"` fits gaussian models to continuous responses and probit
  models to those coded 0/1, `"gaussian"` forces a continuous fit, and
  `"probit"` and `"logistic"` require a 0/1 response. A two-level
  factor, logical, or two-level character response is detected and fit
  as probit; a factor with three or more levels is an error, as `xbart`
  does not cross-validate the multinomial model. The built-in binary
  losses transform test predictions through the family's link. This
  vocabulary is narrower than
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s by
  design - `xbart` cross-validates a single scalar loss per fold, which
  the own-class families' K-forest or two-part fits have no single
  counterpart of; the wider family set lives on `bart2`.

- missing:

  How missing values in the predictors enter the model; as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).

- node.prior:

  An optional expression of the form `normal(k)`, `linear(columns, k)`,
  or `gp(columns, k, ...)` selecting the leaf model, as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). The
  default fits constant leaves. A `k` given inside the prior stands in
  for a missing `k` argument; the `k` argument otherwise drives the
  crossvalidation grid as usual.

  A named leaf calibration (`normal(scale = )`, and likewise for
  `linear` and `gp`) is NOT a grid axis: it is held fixed across every
  cell, so the `k` grid sweeps the prior standard deviation `scale / k`
  about a fixed anchor rather than about whatever the response range
  happened to imply. It is honored in cells that re-use a sampler as
  well as in cells that create one, so the loss surface does not depend
  on the order the cells run in. The `sd` spelling meets the same
  refusal here as everywhere when `k` is a hyperprior.

- n.cuts:

  A positive integer giving the number of decision rules used for each
  predictor, as in
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).

- useQuantiles:

  Logical; as in
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md),
  determines whether decision rules are placed at the empirical
  quantiles of each predictor's values rather than spaced uniformly
  through its range.

- n.thin:

  A positive integer; as in
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md),
  thins each cell's chain against serial correlation. `n.samples` are
  still returned regardless.

- storage:

  A character string selecting the precision of the internal running
  residual, spelled and defaulted as in
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).
  Every cell's sampler is created over a shared per-fold data handle, a
  path the engine keeps in double precision regardless of family or leaf
  model; `"single"` is refused rather than silently ignored.

- tree.prior:

  An expression of the form `cgm` or `cgm(power, base, split.probs)`
  that sets the tree structure prior, or a prior object built with
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md) -
  `cgm(...)` or `dbartsPriors$dart(...)`. `power` and `base` are xbart's
  grid axes, so a supplied object's own `power`/`base` are replaced by
  the swept grid values every cell exactly as the `k` argument replaces
  a supplied `node.prior`'s `k`; the object's other content - a `cgm`
  object's `split.probs`, a DART object's Dirichlet hyperparameters -
  rides every cell unchanged. Because `power`/`base` are grid axes here
  rather than ordinary scalars, they may be supplied alongside
  `tree.prior`; `dart` and `split.probs` would only duplicate what a
  supplied `tree.prior` already specifies, so combining either with it
  is an error naming both.

## Details

Crossvalidates `n.reps` replications against the crossproduct of given
hyperparameter vectors `n.trees` \\\times\\ `k` \\\times\\ `power`
\\\times\\ `base`. For each fit, either one fold is withheld as test
data and `n.test - 1` folds are used as training data or `n * n.test`
observations are withheld as test data and `n * (1 - n.test)` used as
training. A replication corresponds to fitting all \\K\\ folds in
`"k-fold"` crossvalidation or a single fit with `"random subsample"`.
The training data is used to fit a model and make predictions on the
test data which are used together with the test data itself to evaluate
the `loss` function.

`loss` functions are either the default of average negative log-loss for
binary outcomes and root-mean-squared error for continuous outcomes,
misclassification rates for binary outcomes, or a `function` with
arguments `y.test` and `y.test.hat`. `y.test.hat` is of dimensions equal
to `length(y.test)` \\\times\\ `n.samples`, on the latent scale for
binary outcomes. A third option is to pass a list of
`list(function, evaluationEnvironment)`, so as to provide default
bindings. The binary losses `"log"` and `"mcr"` cannot apply to
continuous responses.

## Value

An array of dimensions `n.reps` \\\times\\ `length(n.trees)` \\\times\\
`length(k)` \\\times\\ `length(power)` \\\times\\ `length(base)`. If
`drop` is `TRUE`, dimensions of length 1 are omitted. If all
hyperparameters are of length 1, then the result will be a vector of
length `n.reps`. When the result is an array, the `dimnames` of the
result shall be set to the corresponding hyperparameters.

For method `"k-fold"`, each element is an average across the \\K\\ fits.
For `"random subsample"`, each element represents a single fit.

The result is a bare array with no class, so the fit generics -
`predict`, `extract`, `fitted`, `residuals` - do not apply to it; it is
a table of losses, not a fit.

## Author

Vincent Dorie: <vdorie@gmail.com>

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)

## Examples

``` r
f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

mad <- function(y.test, y.test.hat, weights) {
    # note, weights are ignored
    mean(abs(y.test - apply(y.test.hat, 1L, mean)))
}



## low iteration numbers to to run quickly
xval <- xbart(x, y, n.samples = 15L, n.reps = 4L, n.burn = c(10L, 3L),
              n.trees = c(5L, 7L),
              k = c(1, 2, 4),
              power = c(1.5, 2),
              base = c(0.75, 0.8, 0.95), n.threads = 1L,
              loss = mad)
```
