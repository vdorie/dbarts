# Discrete Bayesian Additive Regression Trees Sampler

Creates a sampler object for a given problem which fits a Bayesian
Additive Regreesion Trees model. Internally stores state in such a way
as to be mutable.

## Usage

``` r
dbarts(
    formula, data, test, subset, weights, offset, offset.test = offset,
    verbose = FALSE, n.samples = 800L,
    tree.prior = cgm, node.prior = normal, resid.prior = chisq,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    control = dbarts::dbartsControl(), sigma = NA_real_, seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic", "aft"),
    missing = c("incorporate", "error"))
```

## Arguments

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  following an analogous model description syntax as
  [`lm`](https://rdrr.io/r/stats/lm.html). For backwards compatibility,
  can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) matrix
  `x.train`, including a sparse `Matrix::dgCMatrix`: its columns enter
  as ordinal predictors, sufficiently sparse columns are stored in a
  compact rank-bitmap layout instead of being expanded, and the
  predictor-mutation surface (`setPredictor` and relatives, `setData`)
  is fixed at creation. Sparse inputs are not supported by
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) or
  the `node.prior = linear()` and `gp()` leaf models; test matrices
  remain dense. A data frame may mix ordinary columns with
  [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  or `dgCMatrix` columns (assign them into the frame; they do not
  survive `data.frame(...)` or
  [`I()`](https://rdrr.io/r/base/AsIs.html)): sparse columns behave as
  in an all-sparse design, while dense columns keep categorical splits
  and linear-leaf designation. A
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  column is recognized but not yet supported; construction refuses it.
  See
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  for the full data-frame-to-predictor mapping.

- data:

  An optional data frame, list, or environment containing predictors to
  be used with the model. For backwards compatibility, can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) vector
  `y.train`.

- test:

  An optional matrix or data frame with the same number of predictors as
  `data`, or `formula` in backwards compatibility mode. If column names
  are present, a matching algorithm is used.

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
  Useful for both response families, though the two apply it
  differently. For a gaussian response, \\y = f(x) + \mathrm{offset} +
  \epsilon\\: the offset is netted out of \\y\\ before BART's internal
  range-scaling (see ‘Details’), so it acts as a fixed component of the
  mean rather than a free parameter. For binary responses it enters the
  link directly: a `"probit"` fit assumes \\P(Y = 1 \mid X = x) =
  \Phi(f(x) + \mathrm{offset})\\, where \\\Phi\\ is the standard normal
  cumulative distribution function, and a `"logistic"` fit the analogous
  model on the logistic link.

- offset.test:

  The equivalent of `offset` for test observations. Will attempt to use
  `offset` when applicable.

- verbose:

  A logical determining if additional output is printed to the console.
  See
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).

- n.samples:

  A positive integer setting the default number of posterior samples to
  be returned for each run of the sampler. Can be overriden at run-time.
  See
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).

- tree.prior:

  An expression of the form `cgm` or `cgm(power, base, split.probs)`
  setting the tree prior used in fitting, or a prior object built with
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md).

- node.prior:

  An expression of the form `normal` or `normal(k)` that sets the prior
  used on the averages within nodes, or a prior object built with
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md).
  `linear(columns, k)` instead fits each leaf with an intercept plus a
  linear term in the designated continuous predictor columns (character
  names or numeric indices into the model matrix), standardized
  internally, all coefficients sharing the `normal(k)` prior; factor
  columns cannot be designated.
  `gp(columns, k, lengthscale, max.leaf.size)` fits each leaf with a
  smooth Gaussian-process function of the designated columns under a
  squared-exponential kernel; leaves larger than `max.leaf.size` fall
  back to constant fits.
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) accepts
  the same specifications through its own `node.prior` argument. See
  “Response scaling” below for how `k` interacts with the response's
  internal scaling.

- resid.prior:

  An expression of the form `chisq` or `chisq(df, quant)` that sets the
  prior used on the residual/error variance, or a prior object built
  with
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md).

- proposal.probs:

  Named numeric vector or `NULL`, optionally specifying the proposal
  rules and their probabilities. Elements should be `"birth_death"`,
  `"change"`, and `"swap"` to control tree change proposals, and
  `"birth"` to give the relative frequency of birth/death in the
  `"birth_death"` step.

- control:

  An object inheriting from `dbartsControl`, created by the
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  function.

- sigma:

  A positive numeric estimate of the residual standard deviation. If
  `NA`, a linear model is used with all of the predictors to obtain one.
  Same concept as `sigest` in
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)/`bart2`
  and [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md);
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) uses
  this function's `sigma` name directly.

- seed:

  Optional integer seed for the random number generator, a convenience
  mirror of `dbartsControl(rngSeed = )`. When not `NA` it overrides the
  seed in `control`; the fitting-function wrappers
  ([`bart2`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)) accept
  the same argument.

- factors:

  How factor columns in a data frame enter the model. The default
  `"categorical"` keeps each unordered factor as a single predictor
  whose splits send a subset of its levels (at most 65535) down each
  branch, and codes ordered factors as ordinal; level tables are
  retained so that test data are coded identically. `"indicators"`
  expands each factor into binary indicator columns, as previous
  versions always did and as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) still
  does.

- family:

  The response model. `"auto"` fits gaussian models to continuous
  responses and probit models to those coded 0/1, as always.
  `"gaussian"` forces a continuous fit even for a 0/1 response;
  `"probit"` and `"logistic"` require a 0/1 response and fit
  latent-variable models, with fits and predictions on the latent scale.
  `"logistic"` uses Polya-Gamma augmentation. `"aft"` fits an
  accelerated failure time (log-normal) survival model: the response is
  a `Surv` object (from the survival package) or a two-column
  `(time, status)` matrix or data frame (status 1 an event, 0
  right-censored; logical status also works, factor status does not),
  and
  [`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
  produces survival-probability draws. Fits and predictions are on the
  log-time scale (\\E\[\log T \mid x\]\\), not the time scale. A `Surv`
  response selects `"aft"` automatically when `family` is `"auto"`; an
  explicitly conflicting family (e.g. `"gaussian"`) with a `Surv`
  response is an error. A `Surv`-like object carrying no `type`
  attribute is treated as right-censored; `type`s other than `"right"`
  are rejected. Survival fits currently use the matrix
  (`x.train`/`y.train`) interface and do not support `subset` or case
  weights. See `weights` for how the families differ in their support
  for weights.

- missing:

  How missing values in the predictors enter the model. The default
  `"incorporate"` keeps them: every split rule learns a direction for
  missing values on its variable, so an observation whose split value is
  `NA` follows that rule's chosen branch (“Missingness Incorporated in
  Attributes”, Twala et al. 2008), in training and test data alike.
  `"error"` rejects predictors containing `NA`. The response, `weights`,
  and `offset` must always be complete; note that previous versions
  silently dropped incomplete rows for formula inputs.

## Details

“Discrete sampler” refers to that `dbarts` is implemented using
ReferenceClasses, so that there exists a mutable object constructed in
C++ that is largely obscured from R. The `dbarts` function is the
primary way of creating a
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
for which a variety of methods exist.

### Data frame ingestion

Data frame input is stored columnar: no \\n \times p\\ predictor matrix
is retained. Code that previously reached into the sampler's stored data
expecting a plain matrix should call
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)`(sampler, "predictors")` -
or `as.matrix` on it - to obtain the numeric predictor matrix, factor
columns as their integer codes. See
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
for the column-type mapping (numeric, factor, sparse, and
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)).

### Test offset synchronization

When `offset.test` is left at its default of tracking `offset`, the
sampler links the two: the sampler's `setOffset` method re-derives
`offset.test` from each new `offset` it is given. Calling the sampler's
`setTestOffset` or `setTestPredictorAndOffset` methods breaks this link
– `offset.test` is set independently from then on, and the link never
re-forms, even if the two are later set to equal values.

### Response scaling

Continuous responses are range-scaled internally: `y` (net of `offset`)
is mapped to \\\[-0.5, 0.5\]\\ by its observed minimum and maximum, the
convention of the entire BART software lineage (BayesTree, BART,
bartMachine), which is what lets `k` (see `node.prior` above) and its
defaults transfer across packages and papers. The known caveat is
outlier sensitivity: extreme `y` values stretch the range and compress
the effective leaf prior for everything else; the published workaround
is to log-transform or winsorize such values before fitting, or to let
`k` adapt via the `chi` hyperprior
([`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)).
The scale is fixed when the sampler is created; the sampler's
`setResponse` and `setOffset` methods only re-anchor it when called with
`updateScale = TRUE`, intended for burn-in only, since re-anchoring
mid-run makes fits across iterations no longer comparable (see
[dbartsSampler-class](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).

## Value

A reference object of
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

## See also

[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md),
[`extract.dbartsSampler`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)

## References

Twala, B.E.T.H., Jones, M.C., and Hand, D.J. (2008) Good methods for
coping with missing data in decision trees. *Pattern Recognition
Letters*, **29**(7), 950–956.
