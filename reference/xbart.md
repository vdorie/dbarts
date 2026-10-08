# Crossvalidation For Bayesian Additive Regression Trees

Fits the BART model against varying `k` (or `sd`), `power`, `base`, and
`n.trees` parameters using \\K\\-fold or repeated random subsampling
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
    k = NULL, sd = NULL, power = 2, base = 0.95,
    split.probs = NULL, drop = TRUE,
    sigest = NULL,
    seed = NULL,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    leaf.prior = NULL, n.cuts = 100L, useQuantiles = FALSE, n.thin = 1L,
    storage = c("double", "single"), tree.prior = NULL,
    parallel = getOption("dbarts.parallel", "auto"), cl = NULL,
    control = dbarts::dbartsControl(), sigma = NULL, ...)
```

## Arguments

- sigma:

  The 0.9-x spelling of `sigest`, accepted for one release with a
  once-per-session warning and removed in dbarts 1.1-0. Supplying both
  is an error.

- ...:

  Not used for new code: the channel that lets a retired argument
  spelling (`dart`, now `tree.prior = dart()`; `resid.prior`, now
  `family = gaussian(sigma = )`, its one home, so writing it both ways
  is refused where the two disagree) reach a message naming its
  successor, instead of R's own “unused argument” error. Any other name
  is refused. Removed in dbarts 1.1-0.

- control:

  A
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object carrying the sampler and engine settings, including the ones
  `xbart` spells no flat name for: `categoricalExhaustiveCap`,
  `testFitParallelCutoff`, `predictParallelCutoff`,
  `sparseDensityThreshold`, `treeShift`, and the tree-move mixture
  `proposal.probs`. Precedence, one rule: a flat argument named in the
  call wins over the control's slot of the same name, and a slot the
  control speaks for - one its own
  [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  call named, or one edited afterwards to differ from a fresh
  control's - wins over this function's default, while a slot it never
  spoke for leaves that default standing. The fields a sweep forces -
  `n.chains`, the control's own `n.threads`, `keepTrees`,
  `keepTrainingFits`, `updateState`, `verbose` - are applied after both,
  and `n.burn` and `n.threads` here are the sweep's own grid axis and
  worker count rather than the control fields of the same name. A
  control taken from a fitted sampler carries that fit's model
  configuration and is refused by name - pass a fresh
  [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  following an analogous model description syntax as
  [`lm`](https://rdrr.io/r/stats/lm.html). For backwards compatibility,
  can also be the
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) matrix
  `x.train`. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). A
  `dbartsData` object carrying `bases` (a multi-forest model) is
  refused: `xbart` cross-validates single-forest models.

- data:

  An optional data frame, list, or environment containing predictors to
  be used with the model. For backwards compatibility, can also be the
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) vector
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
  weights identically 1 are treated as absent - the all-0/1 vector
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  installs as an active-row mask is refused here, cross-validation
  partitioning the rows itself; a `"logistic"` model treats them as
  observation counts and so requires positive integers (its Polya-Gamma
  latent for a count \\w\\ is a sum of \\w\\ unit draws, so a sweep's
  time grows with the total count, and a count above \\10^6\\ is
  refused).

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

  A logical determining if additional output is printed to the console,
  including the one-line message `family = "auto"` prints, once per
  call, naming the family it resolves to.

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
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and
  `bart`, here it is one or two non-negative integers, specifying 1) the
  burn-in when a chain is freshly started against a data split and 2)
  the burn-in when moving from one parameter setting to another over the
  same split. A longer vector is an error naming the argument rather
  than being silently truncated to its first two entries; dbarts 0.9-x
  read a third element as a per-replication burn-in, which no longer
  exists. Chains are never carried between data splits or folds - the
  held-out observations of one were training observations of the
  previous, so continuing a chain lets slowly-mixing settings score
  against data they have effectively seen.

- loss:

  Either one of the pre-set loss functions as character-strings (the
  default is `"rmse"` for a continuous response and `"log"` for a binary
  one; `mcr` - misclassification rate for binary responses, `rmse` -
  root-mean-squared-error for continuous response), `log` - negative
  log-loss for binary response (`rmse` serves this purpose for
  continuous responses), a function, or a function-evaluation
  environment list-pair. A function keeps its own environment, so a
  closure reads what it captured; the list form calls the function from
  the given environment. Functions should have prototypes of the form
  `function(y.test, y.test.hat, weights)`, where `y.test` is the held
  out test subsample, `y.test.hat` is a matrix of dimension
  `length(y.test)` \\\times\\ `n.samples`, and `weights` are an optional
  vector of user-supplied weights. See examples.

- n.threads:

  Every (replication, fold) pair is an independent unit of work, and for
  `n.threads > 1` the units are divided into approximately equal chunks
  and executed on that many parallel workers (see `parallel`). A
  `k`-fold run of a single replication therefore uses up to `k` workers.
  The default uses
  [`guessNumCores`](https://vdorie.github.io/dbarts/reference/guessNumCores.md),
  which should work across the most common operating system/hardware
  pairs, and is one where the cores cannot be counted; a stated `NA` is
  refused. Warnings raised while fitting or scoring, including by a
  supplied `loss` function, are collected from every unit and signalled
  in the calling session once all units have finished, in unit order and
  each distinct warning once, so the same warnings are seen at any
  `n.threads`.

- parallel:

  One of `"auto"`, `"fork"` or `"socket"`, and by default
  `getOption("dbarts.parallel", "auto")`. `"fork"` runs the workers as
  forked copies of the calling session, which is fast to start and
  shares the data; `"socket"` starts separate R sessions (a
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html) cluster)
  that each load dbarts and receive their own copy of the data, which
  costs roughly 80 MB per worker beyond the data itself. `"auto"` forks
  unless the platform is Windows, the session runs under RStudio or
  Positron, or R.app, in which case it uses sockets. `"fork"` is an
  error on Windows. Forking can hang or crash when the calling session
  has already used a library that is not fork-safe, such as a
  multithreaded BLAS; if so, set `parallel = "socket"` or
  `options(dbarts.parallel = "socket")`. Results do not depend on the
  kind of worker. Ignored when `cl` is given or only one worker is used.

- cl:

  `NULL` or a cluster from the `parallel` package that the caller has
  already started, used instead of starting workers whatever `parallel`
  says. The units are divided into `min(n.threads, length(cl))` chunks,
  and `xbart` does not stop the cluster. The workers must load the same
  dbarts version as the calling session, since the work function is sent
  by reference and each worker runs its own installed dbarts. A one-node
  `cl`, like `n.threads = 1`, runs the units in the calling session.

- n.trees:

  A vector of positive integers setting the BART hyperparameter for the
  number of trees in the sum-of-trees formulation. See
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md).

- k:

  The grid for the BART hyperparameter setting the leaf-mean prior
  standard deviation: a vector of positive real numbers, each a cell of
  its own, or a [`list`](https://rdrr.io/r/base/list.html) whose entries
  are positive numbers and `k` hyperpriors
  ([`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)),
  so that `k = list(1, 2, chi())` sweeps two fixed values against the
  modelled one. A modelled cell holds its hyperprior and DRAWS `k` every
  sweep, so its loss is computed under a shrinkage that moves within the
  fit rather than under a named value; it is labelled in the result by
  the constructor call that rebuilds it. If `NULL`, one cell is run at
  the default
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) would fit
  for the response type: fixed 2 for a continuous response and
  `chi(1.5, 2)` for a binary one. A `k` carried by a supplied
  `leaf.prior` stands in for a missing argument. Fixed cells are always
  swept largest to smallest, with a modelled cell last, and the reported
  `k` axis is un-permuted back to the order given, so results do not
  depend on the order `k` is listed in; cells still warm-start off the
  previous one, so this is order-invariance, not an unbiased estimate
  for each cell taken alone.

- sd:

  The same grid axis stated in absolute spreads instead, exclusive with
  `k`: a vector of leaf-prior sds on the family's scale (see
  `normal(sd = )` in
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)),
  or a [`list`](https://rdrr.io/r/base/list.html) mixing them with
  `invchi` laws. Each spread is held fixed across folds, where a `k` is
  relative to each fold's own training range; on a binary family, whose
  `k.scale` is a constant of the latent scale, a grid of fixed cells is
  the same sweep either way (`sd = k.scale / k`); a modelled cell is
  not, since `invchi(df, c)` starts its chain at the spread `c` where
  `chi(df, k.scale / c)` starts at `k.scale / 2`. A `k` inside
  `leaf.prior` beside an `sd` grid is refused, as an `sd` there is
  beside either grid. The cells are swept most shrunk first - the
  smallest sd - with the warm starts and unit parallelism the `k` grid
  has, and the result's axis is labelled `sd`.

- power:

  A vector of real numbers greater than one, setting the BART
  hyperparameter for the tree prior's growth probability, given by
  \\{base} / (1 + depth)^{power}\\.

- base:

  A vector of real numbers in \\(0, 1)\\, setting the BART
  hyperparameter for the tree prior's growth probability.

- split.probs:

  Prior probabilities that a variable is used in a splitting rule, as in
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md). A single
  value or `NULL` yields the uniform default; a named or unnamed vector
  assigns per-column probabilities. Fixed for the whole crossvalidation,
  not part of the swept grid. Cannot be combined with a DART
  `tree.prior`.

- drop:

  Logical, determining if dimensions with a single value are dropped
  from the result.

- sigest:

  A positive numeric estimate of the residual standard deviation. If
  `NULL` (the default), a linear model is used with all of the
  predictors to obtain one, fit separately for each fold on that fold's
  training rows, so no fold's prior reads its held-out responses (a fit
  on all rows runs once beforehand to raise any refusal or fallback
  once; where that fit falls back to the marginal standard deviation, as
  described below, each fold takes its own training rows' marginal
  standard deviation without repeating the linear model); an explicit
  `NA` is a missing value and is refused. Every entry point spells this
  `sigest` - `xbart`,
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)/`bart`
  and the sampler constructors
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)/`dbartsSpec`
  alike; `sigma` is the retired 0.9-x spelling on `dbarts` and on
  `xbart` itself (`dbartsSpec` refuses it); a family object's own
  `sigma` is the prior this estimate calibrates, not the estimate. That
  estimate falls back to the marginal standard deviation of the response
  when the linear model's residual standard error comes out non-finite,
  warning as it does so (class `dbartsSigmaFallbackWarning`); a design
  with sparse-backed predictor columns skips the linear model altogether
  and falls back the same way (class `dbartsSparseSigmaFallbackWarning`,
  a `dbartsSigmaFallbackWarning`). It is the estimate a `chisq` residual
  prior's quantile is calibrated against, so it stands beside
  `family = gaussian(sigma = chisq(df, quant))` and is refused beside
  `family = gaussian(sigma = fixed(value))`, which fixes the residual
  scale outright and would overwrite the estimate with its square root.

- seed:

  Optional integer specifying the desired pRNG
  [seed](https://rdrr.io/r/base/Random.html). `NULL` (the default) means
  not given here and defers to a seed already sitting in `control`, if
  any; an `NA` is a missing value, as in
  [`set.seed`](https://rdrr.io/r/base/Random.html): it reads as `NULL`
  for one release, with a once-per-session warning. From the seed in
  force, `xbart` draws a split seed for each replication and a seed for
  each sampler its units create - one per (replication, fold) unit at
  each distinct tree count - with
  [`sample.int`](https://rdrr.io/r/base/sample.html), in one pass, under
  the caller's own [`RNGkind`](https://rdrr.io/r/base/Random.html): a
  given seed therefore draws different values, and so gives different
  results, under a different kind. Each unit's sampler then takes its
  seed through `control`'s seed slot, exactly as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `seed` does, driving a dedicated generator that never reads R's
  stream, so results are reproducible at any `n.threads` for a given
  seed and `RNGkind`. A supplied `loss` function that draws random
  numbers draws on the worker running it, from that process's own
  default generator; at `n.threads` greater than 1 this makes the run
  not reproducible, whatever seed `xbart` itself is given. The caller's
  random stream is left untouched when a seed is given, and advanced
  only by that derivation when one is not; without a seed,
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand suffices.
  See the Reproducibility section of
  [bart](https://vdorie.github.io/dbarts/reference/bart.md).

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
  does not cross-validate the multinomial model, and neither is a matrix
  response (per-category counts or a (time, status) pair), which is
  refused naming `bart`. A survival response - a `Surv` object or
  two-column `(time, status)` pair, on `formula` or as `data` directly -
  is likewise an error naming `bart`/`dbarts` with
  `family = "aft"`/`"hazard"`: `xbart` does not cross-validate a
  survival model. The built-in binary losses transform test predictions
  through the family's link. This vocabulary is narrower than
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s by
  design - `xbart` cross-validates a single scalar loss per fold, which
  the own-class families' K-forest or two-part fits have no single
  counterpart of; the wider family set lives on `bart`. Base R family
  objects map as [`glm`](https://rdrr.io/r/stats/glm.html) takes them
  (`binomial` is the logit link, so `family = binomial` is `"logistic"`,
  not probit) and the rest are refused; see
  [`dbartsFamilies`](https://vdorie.github.io/dbarts/reference/dbartsFamilies.md).

- leaf.prior:

  An optional expression of the form `normal(k)`, `linear(columns, k)`,
  or `gp(columns, k, ...)` selecting the leaf model, as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). The
  default fits constant leaves. A `k` given inside the prior stands in
  for a missing `k` argument; the `k` argument otherwise drives the
  crossvalidation grid as usual.

  A named `sd` inside the prior (`normal(sd = )`, and likewise for
  `linear` and `gp`) stands in for a missing `sd` argument, as a
  one-cell `sd` axis. Beside a `k` or `sd` grid it is refused, since
  both would state the spread: name the spreads in the `sd` grid. Two
  calls with the same `seed` use the same folds whatever their spreads,
  so one call per `sd` reproduces each cell a grid's fresh start would.

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
  a supplied `leaf.prior`'s `k`; the object's other content - a `cgm`
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

An array with up to six dimensions, in order `rep` (length `n.reps`),
`n.trees`, `k` (or `sd`), `power`, `base`, and `loss`. `rep` is always
present. `n.trees`, `k`, `power`, and `base` are each omitted when
`drop` is `TRUE` and the corresponding grid has length 1; with
`drop = FALSE` they are always present, an absent `k` included, whose
single cell is then named for the default it ran at. The trailing `loss`
dimension, sized to however many values a single call to `loss` returns,
is present only when that count is greater than one - independent of
`drop` entirely; the default losses and an ordinary scalar-returning
custom `loss` never contribute it. When none of the above survive, the
result collapses to a plain vector of length `n.reps`. When the result
remains an array, its `dimnames` name the swept values on each surviving
hyperparameter axis - exact integer labels for `n.trees`, values rounded
to 2 significant digits for the double-valued `k`/`power`/`base` axes,
and the constructor call for a modelled `k` cell; the `loss` axis, when
present, carries no per-slot names, and neither does `rep`.

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
