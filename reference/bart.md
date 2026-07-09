# Bayesian Additive Regression Trees

BART is a Bayesian “sum-of-trees” model in which each tree is
constrained by a prior to be a weak learner.

- For numeric response \\y = f(x) + \epsilon\\, where \\\epsilon \sim
  N(0, \sigma^2)\\.

- For binary response \\y\\, \\P(Y = 1 \mid x) = \Phi(f(x))\\, where
  \\\Phi\\ denotes the standard normal cdf (probit link).

## Usage

``` r
bart(
    x.train, y.train, x.test = matrix(0.0, 0, 0),
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2.0,
    power = 2.0, base = 0.95, splitprobs = 1 / numvars,
    binaryOffset = 0.0, weights = NULL,
    ntree = 200,
    ndpost = 1000, nskip = 100,
    printevery = 100, keepevery = 1, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100, printcutoffs = 0,
    verbose = TRUE, nchain = 1, nthread = 1, combinechains = TRUE,
    keeptrees = FALSE, keepcall = TRUE, sampleronly = FALSE,
    seed = NA_integer_,
    proposalprobs = NULL,
    keepsampler = keeptrees)

bart2(
    formula, data, test, subset, weights, offset, offset.test = offset,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL,
    power = 2.0, base = 0.95, split.probs = 1 / num.vars,
    dart = FALSE,
    n.trees = 75L,
    n.samples = 500L, n.burn = 500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = TRUE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 1L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE, keepTrees = FALSE, 
    keepCall = TRUE, samplerOnly = FALSE,
    seed = NA_integer_,
    proposal.probs = NULL,
    keepSampler = keepTrees,
    warm.start = NULL,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    missing = c("incorporate", "error"),
    ...)

# S3 method for class 'bart'
plot(
    x,
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    ...)

# S3 method for class 'bart'
predict(
    object, newdata, offset, weights,
    type = c("ev", "ppd", "bart"),
    combineChains = TRUE,
    n.threads,
    ci.level = NULL,
    ...)

extract(object, ...)
# S3 method for class 'bart'
extract(
    object,
    type = c("ev", "ppd", "bart", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE, ...)

# S3 method for class 'bart'
fitted(
    object,
    type = c("ev", "ppd", "bart"),
    sample = c("train", "test"),
    ci.level = NULL,
    ...)

# S3 method for class 'bart'
residuals(object, type = "ev", ...)
```

## Arguments

- x.train:

  Explanatory variables for training (in sample) data. May be a matrix
  or a data frame, with rows corresponding to observations and columns
  to variables. If a variable is a factor in a data frame, it is
  replaced with dummies. Note that \\q\\ dummies are created if \\q \>
  2\\ and one dummy is created if \\q = 2\\, where \\q\\ is the number
  of levels of the factor. A sparse `Matrix::dgCMatrix` of ordinal
  predictors is also accepted, as is a data frame mixing ordinary
  columns with
  [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  or `dgCMatrix` columns; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).

- y.train:

  Dependent variable for training (in sample) data. If `y.train` is
  numeric a continous response model is fit (normal errors). If
  `y.train` is a binary factor or has only values 0 and 1, then a binary
  response model with a probit link is fit.

- x.test:

  Explanatory variables for test (out of sample) data. Should have same
  column structure as `x.train`. `bart` will generate draws of \\f(x)\\
  for each \\x\\ which is a row of `x.test`.

- sigest:

  For continuous response models, an estimate of the residual standard
  deviation (residual standard error), \\\sigma\\, used to calibrate an
  inverse-chi-squared prior on the error variance. If not supplied, the
  least-squares estimate is derived instead. See `sigquant` for more
  information. Not applicable when \\y\\ is binary. Same concept as
  `sigma` in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md).

- sigdf:

  Degrees of freedom for error variance prior. Not applicable when \\y\\
  is binary. Default 3, Chipman, George, and McCulloch's calibration
  (see References); an aggressive choice rather than a derived one.

- sigquant:

  The quantile of the error variance prior that the rough estimate
  (`sigest`) is placed at. The closer the quantile is to 1, the more
  aggressive the fit will be as you are putting more prior weight on
  error standard deviations (\\\sigma\\) less than the rough estimate.
  Not applicable when \\y\\ is binary. Default 0.90, from the same
  source as `sigdf`.

- k:

  For numeric \\y\\, `k` is the number of prior standard deviations
  \\E(Y\|x) = f(x)\\ is away from \\\pm 0.5\\. The response (`y.train`)
  is internally scaled to range from \\-0.5\\ to \\0.5\\. For binary
  \\y\\ fit with the default probit link, `k` is the number of prior
  standard deviations \\f(x)\\ is away from \\\pm 3\\, the probit
  reference scale; `family = "logistic"` (`bart2` only) widens this to
  \\\pm \pi \sqrt{3}\\, the standard deviation of the standard logistic
  latent variable. In both cases, the bigger \\k\\ is, the more
  conservative the fitting will be. The value can be either a fixed
  number, or the a *hyperprior* of the form
  `chi(degreesOfFreedom = 1.5, scale = Inf)`. For `bart2`, the default
  of `NULL` uses the value 2 for continuous reponses and a `chi`
  hyperprior for binary ones. The default `chi` hyperprior is improper,
  and slightly penalizes small values of `k`. The default of 2 for
  continuous responses follows Chipman, George, and McCulloch's argument
  (see References) that with node prior standard deviation \\\sigma\_\mu
  = 0.5 / (k \sqrt{m})\\ for \\m\\ trees, \\k\\ prior standard
  deviations of \\f(x)\\ span the whole coded response range regardless
  of \\m\\ – so \\k = 2\\ places that range at roughly a 95% prior
  interval. The range scaling itself is an outlier-sensitive convention
  shared with BayesTree/bartMachine: extreme \\y\\ values compress the
  effective prior on everything else, for which the published workaround
  is to log-transform or winsorize such values before fitting, or to let
  `k` adapt via the `chi` hyperprior.

- power:

  Power parameter for tree prior. Default 2, Chipman, George, and
  McCulloch's empirical recommendation (see References) for the
  split-probability decay \\base (1 + depth)^{-power}\\, not derived
  from a formula.

- base:

  Base parameter for tree prior. Default 0.95, from the same source as
  `power`.

- splitprobs, split.probs:

  Prior and transition probabilities of variables used to generate
  splits. Can be missing/empty/`NULL` for equiprobability, a numeric
  vector of length equal to the number variables, or a named numeric
  vector with only a subset of the variables specified and a `.default`
  named value. Values given for factor variables are replicated for each
  resulting column in the generated model matrix. `numvars` and
  `num.vars` symbols will be rebound before execution to the number of
  columns in the model matrix.

- dart:

  `TRUE` places a DART prior (a Dirichlet prior over the split-variable
  probabilities, inducing variable selection; Linero 2018) on the trees
  with default settings, using `power` and `base` for the growth
  probabilities; a prior created by
  [`dbartsPriors$dart`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
  supplies its own settings, overriding `power` and `base`. Cannot be
  combined with `split.probs`. Per-sample probabilities are returned in
  the `varprobs` element of the fit.

- binaryOffset:

  Used for binary \\y\\. When present, the model is \\P(Y = 1 \mid x) =
  \Phi(f(x) + \mathrm{binaryOffset})\\, allowing fits with probabilities
  shrunk towards values other than \\0.5\\. \\\Phi\\ is the link for the
  default probit family; `bart2` fits with `family = "logistic"` use
  `plogis` instead, as noted in the value section below.

- weights:

  An optional vector of weights to be used in the fitting process. For a
  gaussian response, BART fits a model with observations \\y \mid x \sim
  N(f(x), \sigma^2 / w)\\, where \\f(x)\\ is the unknown function. A
  probit fit (`bart`, or `bart2` with the default binary family) does
  not support weights, except that weights identically 1 are treated as
  absent; `bart2` with `family = "logistic"` treats them as observation
  counts and requires positive integers.

- ntree, n.trees:

  The number of trees in the sum-of-trees formulation.

- ndpost, n.samples:

  The number of posterior draws after burn in, `ndpost / keepevery` will
  actually be returned.

- nskip, n.burn:

  Number of MCMC iterations to be treated as burn in.

- printevery, printEvery:

  As the MCMC runs, a message is printed every `printevery` draws.

- keepevery, n.thin:

  Every `keepevery` draw is kept to be returned to the user. Useful for
  “thinning” samples.

- keeptrainfits, keepTrainingFits:

  If `TRUE` the draws of \\f(x)\\ for \\x\\ corresponding to the rows of
  `x.train` are returned.

- usequants, useQuantiles:

  When `TRUE`, determine tree decision rules using estimated quantiles
  derived from the `x.train` variables. When `FALSE`, splits are
  determined using values equally spaced across the range of a variable.
  See details for more information.

- numcut, n.cuts:

  The maximum number of possible values used in decision rules (see
  `usequants`, details). If a single number, it is recycled for all
  variables; otherwise must be a vector of length equal to
  `ncol(x.train)`. Fewer rules may be used if a covariate lacks enough
  unique values.

- printcutoffs, printCutoffs:

  The number of cutoff rules to printed to screen before the MCMC is
  run. Given a single integer, the same value will be used for all
  variables. If 0, nothing is printed.

- verbose:

  Logical; if `FALSE` suppress printing.

- nchain, n.chains:

  Integer specifying how many independent tree sets and fits should be
  calculated. Default 1 for `bart`, BayesTree's historical single-chain
  convention; `bart2`'s (and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)'s)
  default is 4.

- nthread, n.threads:

  Integer specifying how many threads to use. Depending on the CPU
  architecture, using more than the number of chains can degrade
  performance for small/medium data sets. As such some calculations may
  be executed single threaded regardless.

- combinechains, combineChains:

  Logical; if `TRUE`, samples will be returned in arrays of dimensions
  equal to `nchain` \\\times\\ `ndpost` \\\times\\ number of
  observations. Default `TRUE` across `bart`, `bart2`, and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md); the
  number of chains is tracked regardless on the returned object's
  `n.chains` component (see ‘Value’).

- keeptrees, keepTrees:

  Logical; must be `TRUE` in order to use `predict` with the result of a
  `bart` fit. Note that for models with a large number of observations
  or a large number of trees, keeping the trees can be very memory
  intensive.

- keepcall, keepCall:

  Logical; if `FALSE`, returned object will have `call` set to
  `call("NULL")`, otherwise the call used to instantiate BART.

- seed:

  Optional integer specifying the desired pRNG
  [seed](https://rdrr.io/r/base/Random.html). A
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand suffices
  for reproducibility; supplying `seed` instead gives reproducible
  results without touching R's stream. See Reproducibility section
  below.

- proposalprobs, proposal.probs:

  Named numeric vector or `NULL`, optionally specifying the proposal
  rules and their probabilities. Elements should be `"birth_death"`,
  `"change"`, and `"swap"` to control tree change proposals, and
  `"birth"` to give the relative frequency of birth/death in the
  `"birth_death"` step. Defaults are 0.5, 0.1, 0.4, and 0.5
  respectively.

- keepsampler, keepSampler:

  Logical that can be used to save the underlying
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  object even if `keepTrees` is false.

- warm.start:

  `bart2` only. A previous fit whose forests seed the initial state
  instead of drawing trees from the prior: either a `bart` object fit
  with `keepSampler = TRUE` or a
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
  The donor must share this fit's predictors, tree count, cut grid, and
  DART setting; only its trees, `sigma`, and `k` carry over, and each
  chain starts from a different donor sample so multiple chains stay
  overdispersed. A warm start biases early draws toward the donor, so it
  shortens burn-in rather than removing it; keep a non-zero `n.burn`.
  See
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
  `installTrees` method for finer control.

- factors, family, missing:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md);
  `bart2` only. `bart` keeps the historical behavior: factors always
  expand into indicator columns (`bart2`'s default, `"categorical"`,
  instead keeps each unordered factor as a single predictor - a
  different model-matrix representation, so a factor-predictor fit
  changes if moved from one interface to the other), binary responses
  are probit, and missing data are rejected.

- formula:

  The same as `x.train`, the name reflecting that a formula object can
  be used instead.

- data:

  The same as `y.train`, the name reflecting that a data frame can be
  specified when a formula is given instead.

- test:

  The same as `x.train`. Can be missing.

- subset:

  A vector of logicals or indices used to subset of the data. Can be
  missing.

- offset:

  The same as `binaryOffset`. Can be missing.

- offset.test:

  A vector of offsets to be used with test data, in case it is different
  than the training offset. If `offset.test` is missing, defaults to
  `NULL`.

- object:

  An object of class `bart`, returned from either the function `bart` or
  `bart2`.

- newdata:

  Test data for prediction. Obeys all the same rules as `x.train` but
  cannot be missing.

- sampleronly, samplerOnly:

  Builds the sampler from its arguments and returns it without running
  it. Useful to use the `bart2` interface in more complicated models.

- x:

  Object of class `bart`, returned by function `bart`, which contains
  the information to be plotted.

- plquants:

  In the plots, beliefs about \\f(x)\\ are indicated by plotting the
  posterior median and a lower and upper quantile. `plquants` is a
  double vector of length two giving the lower and upper quantiles.

- cols:

  Vector of two colors. First color is used to plot the median of
  \\f(x)\\ and the second color is used to plot the lower and upper
  quantiles.

- type:

  The quantity to be returned by generic functions. Options are `"ev"` -
  samples from the posterior of the individual level expected value,
  `"bart"` - the sum of trees component; same as `"ev"` for linear
  models but on the latent scale (probit or logistic, matching the fit's
  family) for binary ones, `"ppd"` - samples from the posterior
  predictive distribution, and `"trees"` - a data frame with tree
  information for when model was fit with `keepTrees` equal to `TRUE`.
  To synergize with
  [`predict.glm`](https://rdrr.io/r/stats/predict.glm.html),
  `"response"` can be used as a synonym for `"ev"` and `"link"` can be
  used as a synonym for `"bart"`. For information on extracting trees,
  see the subsection below.

- sample:

  Either `"train"` or `"test"`.

- ci.level:

  For `fitted` and `predict`, an optional single number in \\(0, 1)\\.
  When `NULL` (the default) the posterior mean is returned as before.
  When given, the result is instead a matrix with one row per
  observation and columns `est` (the posterior mean), `ci.lower`, and
  `ci.upper` giving a symmetric credible band at that level, pooled over
  all samples and chains. The kind of interval follows `type`: `"ev"`
  yields a credible interval for the conditional mean (a probability for
  binary responses), `"ppd"` a prediction interval that additionally
  carries the residual noise (and so is wider), and `"bart"` a credible
  interval on the latent scale.

- ...:

  Additional arguments passed on to `plot`, `dbartsControl`, or
  `extract` when `type` is `"trees"`. Not used in `predict`.

## Details

BART is an Bayesian MCMC method. At each MCMC interation, we produce a
draw from the joint posterior \\(f, \sigma) \mid (x, y)\\ in the numeric
\\y\\ case and just \\f\\ in the binary \\y\\ case.

Thus, unlike a lot of other modeling methods in R, `bart` does not
produce a single model object from which fits and summaries may be
extracted. The output consists of values \\f^\*(x)\\ (and \\\sigma^\*\\
in the numeric case) where \* denotes a particular draw. The \\x\\ is
either a row from the training data (`x.train`) or the test data
(`x.test`).

### Decision Rules

Decision rules for any tree are of the form \\x \le c\\ vs. \\x \> c\\
for each ‘\\x\\’ corresponding to a column of `x.train`. `usequants`
determines the means by which the set of possible \\c\\ is determined.
If `usequants` is `TRUE`, then the \\c\\ are a subset of the values
interpolated half-way between the unique, sorted values obtained from
the corresponding column of `x.train`. If `usequants` is `FALSE`, the
cutoffs are equally spaced across the range of values taken on by the
corresponding column of `x.train`.

The number of possible values of \\c\\ is determined by `numcut`. If
`usequants` is `FALSE`, `numcut` equally spaced cutoffs are used
covering the range of values in the corresponding column of `x.train`.
If `usequants` is `TRUE`, then for a variable the minimum of `numcut`
and one less than the number of unique elements for that variable are
used.

### End-node prior parameter `k`

The amount of shrinkage of the node parameters is controlled by `k`. `k`
can be given as either a fixed, positive number, or as any value that
can be used to build a supported hyperprior. At present, only
\\\chi\_\nu s\\ priors are supported, where \\\nu\\ is a degrees of
freedom and \\s\\ is a scale. Both values must be positive, however the
scale can be infinite which yields an improper prior, which is
interpretted as just the polynomial part of the distribution. If \\nu\\
is 1 and \\s\\ is \\\infty\\, the prior is “flat”.

For BART on binary outcomes, the degree of overfitting can be highly
sensitive to `k` so it is encouraged to consider a number of values. The
default hyperprior for binary BART, `chi(1.5, Inf)`, has been shown to
work well in a large number of datasets, however crossvalidation may be
helpful. Running for a short time with a flat prior may be helpful to
see the range of values of `k` that are consistent with the data.

### Generics

`bart` and
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) support
[`fitted`](https://rdrr.io/r/stats/fitted.values.html) to return the
posterior mean of a predicted quantity, as well as
[`predict`](https://rdrr.io/r/stats/predict.html) to return a set of
posterior samples for a different sample. In addition, the `extract`
generic can be used to obtain the posterior samples for the training
data or test data supplied during the initial fit.

Using `predict` with a `bart` object requires that it be fitted with the
option `keeptrees`/`keepTrees` as `TRUE`. Keeping the trees for a fit
can require a sizeable amount of memory and is off by default.

All generics return values on the scale of expected value of the
response by default. This means that `predict`, `extract`, and `fitted`
for binary outcomes return probabilities unless specifically the
sum-of-trees component is requested (`type = "bart"`). This is in
contrast to `yhat.train`/`yhat.test` that are returned with the fitted
model.

`residuals` returns \\y\\ minus
`fitted(object, type = type, sample = "train")`; the default
`type = "ev"` gives ordinary response-scale residuals, so for a binary
fit this is \\y - P(Y = 1 \mid x)\\. Binary fits do not have
`yhat.train.mean`/`yhat.test.mean` components (see ‘Value’); use
`fitted` to obtain the equivalent posterior-mean probabilities directly.

### Saving

[`save`](https://rdrr.io/r/base/save.html)ing and
[`load`](https://rdrr.io/r/base/load.html)ing fitted BART objects for
use with `predict` requires that R's serialization mechanism be able to
access the underlying trees, in addition to being fit with
`keeptrees`/`keepTrees` as `TRUE`. For memory purposes, the trees are
not stored as R objects unless specifically requested. To do this, one
must “touch” the sampler's state object before saving, e.g. for a fitted
object `bartFit`, execute `invisible(bartFit$fit$state)`.

### Reproducibility

Every chain runs its own pseudo random number generator, so worker
threads never touch R's generator and results do not depend on the
thread count. Without a `seed`, chain generators are seeded from R's
stream when the sampler is created, so calling
[`set.seed`](https://rdrr.io/r/base/Random.html) beforehand makes
results reproducible. Sampling itself never advances R's stream: after a
fit, [`.Random.seed`](https://rdrr.io/r/base/Random.html) has moved only
by the draws taken at creation.

Setting `seed` makes results reproducible without involving R's stream
at all: the seed drives a dedicated generator that hands each chain its
own seed. A single-chain run with a given seed reproduces the first
chain of a multi-chain run with the same seed.

### Extracting Trees

When a model is fit with `keeptrees` (`bart`) or `keepTrees` (`bart2`)
equal to `TRUE`, the generic `extract` can be used to retrieve a data
frame containing the tree fit information. In this case, `extract` will
accept the additional, optional arguments: `chainNums`, `sampleNums`,
and `treeNums`. Each should be an integer vector detailing the desired
trees to be returned. A further optional argument `newdata` routes a new
set of predictors (in the same form accepted by `predict`) through the
frozen trees so that the `n` column counts those observations instead of
the training data.

The result of `extract` will be a data frame with columns:

- `chain`, `sample`, `tree` - index variables

- `n` - number of observations in node, drawn from the training data or
  from `newdata` when supplied

- `var` - either the index of the variable used for splitting or -1 if
  the node is a leaf

- `value` - either the value such that observations less than or equal
  to it are sent down the left path of the tree or the predicted value
  for a leaf node

The order of nodes in the result corresponds to a depth-first traversal,
going down the left-side first. The names of variables used in splitting
can be recovered by examining the column names of the `fit$data@x`
element of a fitted `bart` or `bart2` model. See the package vignette
“Working with dbarts Saved Trees”.

## Value

`bart` and `bart2` return lists assigned the class `bart`. For
applicable quantities, `ndpost / keepevery` samples are returned. In the
numeric \\y\\ case, the list has components:

- `yhat.train`:

  A array/matrix of posterior samples. The \\(i, j, k)\\ value is the
  \\j\\th draw of the posterior of \\f\\ evaluated at the \\k\\th row of
  `x.train` (i.e. \\f^\*(x_k)\\) corresponding to chain \\i\\. When
  `nchain` is one or `combinechains` is `TRUE`, the result is a
  collapsed down to a matrix.

- `yhat.test`:

  Same as `yhat.train` but now the \\x\\s are the rows of the test data.

- `yhat.train.mean`:

  Vector of means of `yhat.train` across columns and chains, with length
  equal to the number of training observations.

- `yhat.test.mean`:

  Vector of means of `yhat.test` across columns and chains.

- `sigma`:

  Matrix of posterior samples of `sigma`, the residual/error standard
  deviation. Dimensions are equal to the number of chains times the
  numbers of samples unless `nchain` is one or `combinechains` is
  `TRUE`.

- `first.sigma`:

  Burn-in draws of `sigma`.

- `varcount`:

  A matrix with number of rows equal to the number of kept draws and
  each column corresponding to a training variable. Contains the total
  count of the number of times that variable is used in a tree decision
  rule (over all trees).

- `varprobs`:

  Present only under a DART prior: the sampled split-variable
  probabilities, in the same layout as `varcount`. Each draw sums to one
  across variables.

- `sigest`:

  The rough error standard deviation (\\\sigma\\) used in the prior.

- `y`:

  The input dependent vector of values for the dependent variable. This
  is used in `plot.bart`.

- `fit`:

  Optional sampler object which stores the values of the tree splits.
  Required for using `predict` and only stored if `keeptrees` or
  `keepsampler` is `TRUE`.

- `n.chains`:

  Information that can be lost if `combinechains` is `TRUE` is tracked
  here.

- `k`:

  Optional matrix of posterior samples of `k`. Only present when `k` is
  modeled, i.e. there is a hyperprior.

- `first.k`:

  Burn-in draws of `k`, if modeled.

- `family`:

  The resolved response family (`"gaussian"`, `"probit"`, or
  `"logistic"`). `predict`, `extract`, `fitted`, and `plot` use it to
  transform latent draws to probabilities.

In the binary \\y\\ case, the returned list has the components
`yhat.train`, `yhat.test`, and `varcount` as above, but not
`yhat.train.mean`/`yhat.test.mean` - use `fitted` to get the posterior
mean of \\P(Y = 1 \mid x)\\ instead. In addition the list has a
`binaryOffset` component giving the value used.

Note that in the binary \\y\\, case `yhat.train` and `yhat.test` are
\\f(x) + \mathrm{binaryOffset}\\. For draws of the probability \\P(Y = 1
\| x)\\, apply the family's inverse link to these values: the normal cdf
(`pnorm`) for `"probit"` fits, `plogis` for `"logistic"` ones.

For continuous response fits, the `plot` method sets `mfrow` to
`c(1, 2)` and makes two plots. The first plot is the sequence of kept
draws of \\\sigma\\ including the burn-in draws. Initially these draws
will decline as BART finds a good fit and then level off when the MCMC
has burnt in. The second plot has \\y\\ on the horizontal axis and
posterior intervals for the corresponding \\f(x)\\ on the vertical axis.
For binary response fits, only this second kind of plot is drawn, with
the posterior median of \\P(Y = 1 \mid x)\\ on the horizontal axis and
its posterior interval on the vertical axis.

## References

Chipman, H., George, E., and McCulloch, R. (2009) BART: Bayesian
Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble
Learning. Advances in Neural Information Processing Systems 19,
Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

both of the above at: <https://www.rob-mcculloch.org>

Friedman, J.H. (1991) Multivariate adaptive regression splines. *The
Annals of Statistics*, **19**, 1–67.

Linero, A.R. (2018) Bayesian regression trees for high-dimensional
prediction and variable selection. *Journal of the American Statistical
Association*, **113**(522), 626–636.

## Author

Hugh Chipman: <hugh.chipman@gmail.com>, Robert McCulloch:
<robert.mcculloch1@gmail.com>, Vincent Dorie: <vdorie@gmail.com>.

## See also

[`pdbart`](https://vdorie.github.io/dbarts/reference/pdbart.md)

## Examples

``` r
## simulate data (example from Friedman MARS paper)
## y = f(x) + epsilon , epsilon ~ N(0, sigma)
## x consists of 10 variables, only first 5 matter

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

## run BART
set.seed(99)
bartFit <- bart(x, y)
#> 
#> Running BART with numeric y
#> 
#> number of trees: 200
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 1
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.002181
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  proposal probabilities: birth/death 0.50, swap 0.10, change 0.40; birth 0.50
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 10
#>  init sigma: 2.756573, curr sigma: 2.756573
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) (3: 100) (4: 100) (5: 100) 
#> (6: 100) (7: 100) (8: 100) (9: 100) (10: 100) 
#> 
#> Running mcmc loop:
#> iteration: 100 (of 1000)
#> iteration: 200 (of 1000)
#> iteration: 300 (of 1000)
#> iteration: 400 (of 1000)
#> iteration: 500 (of 1000)
#> iteration: 600 (of 1000)
#> iteration: 700 (of 1000)
#> iteration: 800 (of 1000)
#> iteration: 900 (of 1000)
#> iteration: 1000 (of 1000)
#> total seconds in loop: 0.225102
#> 
#> Tree sizes, last iteration:
#> [1] 2 3 3 2 2 2 2 2 4 2 3 3 3 1 2 1 2 3 
#> 2 2 4 2 2 3 3 2 2 2 2 3 2 2 2 1 3 3 2 2 
#> 2 2 3 4 2 2 2 4 3 2 2 3 1 2 3 2 2 3 3 2 
#> 3 2 2 2 2 3 2 2 3 2 2 2 2 2 2 2 2 3 2 2 
#> 2 2 2 2 3 5 2 2 3 2 2 2 1 3 2 2 2 3 2 2 
#> 1 2 5 1 3 3 3 4 2 2 2 2 3 2 2 2 2 2 2 1 
#> 2 4 2 2 2 2 3 2 2 2 2 4 2 2 3 2 2 2 2 3 
#> 2 3 2 2 2 2 2 1 2 4 4 3 2 4 4 3 2 1 2 3 
#> 3 4 3 2 3 2 2 2 2 2 2 4 2 3 2 3 3 2 2 2 
#> 3 2 2 3 2 2 2 4 4 2 3 1 3 2 2 2 1 3 2 4 
#> 2 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 34) (2: 29) (3: 29) (4: 30) (5: 25) 
#> (6: 25) (7: 36) (8: 21) (9: 28) (10: 16) 
#> 
#> DONE BART
#> 

plot(bartFit)


## compare BART fit to linear matter and truth = Ey
lmFit <- lm(y ~ ., data.frame(x, y))

fitmat <- cbind(y, Ey, lmFit$fitted, bartFit$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
#>              y        Ey        lm      bart
#> y    1.0000000 0.9847984 0.8841787 0.9984931
#> Ey   0.9847984 1.0000000 0.9009389 0.9886903
#> lm   0.8841787 0.9009389 1.0000000 0.8975062
#> bart 0.9984931 0.9886903 0.8975062 1.0000000

## binary response with a logistic link, via bart2
set.seed(0)
n <- 60L
x.bin <- matrix(runif(n * 2), n, 2)
f.bin <- 3 * x.bin[,1] - 1.5
y.bin <- rbinom(n, 1L, plogis(f.bin))
fit.logit <- bart2(y.bin ~ x.bin, family = "logistic",
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L)
#> 
#> Running BART with binary y
#> 
#> number of trees: 20
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 1
#> Prior:
#>  prior on k: chi with 1.500000 degrees of freedom and inf scale
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  proposal probabilities: birth/death 0.50, swap 0.10, change 0.40; birth 0.50
#> data:
#>  number of training observations: 60
#>  number of test observations: 0
#>  number of explanatory variables: 2
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) 
#> Running mcmc loop:
#> total seconds in loop: 0.001358
#> 
#> Tree sizes, last iteration:
#> [1] 3 4 3 2 2 2 2 3 3 2 3 3 3 1 2 3 2 1 
#> 2 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 16) (2: 12) 
#> DONE BART
#> 

## fit with missing predictor values, using missing = "incorporate"
## (the default): every split rule learns a direction for NAs
set.seed(1)
x.na <- x.bin
x.na[1:5, 1] <- NA
y.na <- x.na[,2] + rnorm(n, 0, 0.3)
y.na[is.na(x.na[,1])] <- y.na[is.na(x.na[,1])] + 2
fit.na <- dbarts(y.na ~ x.na, control = dbartsControl(
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L))
samples.na <- fit.na$run()
```
