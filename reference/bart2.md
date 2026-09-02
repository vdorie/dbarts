# Bayesian Additive Regression Trees with Extended Response Families

BART is a Bayesian “sum-of-trees” model in which each tree is
constrained by a prior to be a weak learner; see
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) for the base
gaussian/probit model definition. `bart2` is dbarts' own formula-first
fitting interface: besides the default gaussian/probit response it
reaches further families through `family` - logistic, accelerated
failure time (`family = "aft"`), discrete-time survival hazard
(`family = "hazard"`, `"hazard.logistic"`), K-category multinomial
(`family = "multinomial"`), ordered categorical (`family = "ordinal"`, a
cumulative probit), negative-binomial counts (`family = "nbinom"`), and
semicontinuous two-part responses (`family = "hurdle.lognormal"`) - and,
on its formula interface, an additional per-observation-modulated forest
declared inline with a
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
(see ‘Formula Terms’ below). These remaining own-class families are not
reachable from
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md): naming one
of them to `bart`'s own `family` is refused, pointing back here.

## Usage

``` r
bart2(
    formula, data, test, subset, weights, offset, offset.test = offset,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = NULL, prior.scale = NA_real_,
    power = 2.0, base = 0.95, split.probs = NULL,
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
    proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL,
    interactions = NULL,
    blocks = NULL,
    variance = NULL,
    keepSampler = keepTrees,
    warm.start = NULL,
    n.grow.sweeps = 0L,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic", "aft", "multinomial",
               "ordinal", "nbinom", "hazard", "hazard.probit",
               "hazard.logistic", "hurdle.lognormal", "twopart"),
    missing = c("incorporate", "error"),
    resid.dist = gaussian,
    dispersion = NA_real_,
    breaks = NULL, max.rows = 1e7,
    tree.prior = NULL, node.prior = NULL, resid.prior = NULL,
    storage = c("double", "single"), updateState = TRUE)

# S3 method for class 'bartMultinomial'
extract(
    object,
    type = c("ev", "ppd", "bart", "forest", "loglik"),
    sample = c("train", "test"),
    combineChains = TRUE, ...)

# S3 method for class 'bartMultinomial'
fitted(
    object, type = c("ev", "class", "bart"),
    ci.level = NULL, ...)

# S3 method for class 'bartMultinomial'
predict(
    object, newdata,
    type = c("ev", "ppd", "bart", "forest", "class"),
    offset = NULL,
    combineChains = TRUE, ci.level = NULL, n.threads, ...)

# S3 method for class 'bartMultinomial'
print(x, ...)

# S3 method for class 'bartMultinomial'
residuals(object, ...)

# S3 method for class 'bartOrdinal'
extract(
    object, type = c("ev", "ppd", "bart", "loglik"),
    sample = c("train", "test"),
    combineChains = TRUE, ...)

# S3 method for class 'bartOrdinal'
fitted(
    object, type = c("ev", "class", "bart"),
    ci.level = NULL, ...)

# S3 method for class 'bartOrdinal'
predict(
    object, newdata,
    type = c("ev", "ppd", "bart", "class"),
    offset = NULL,
    combineChains = TRUE, ci.level = NULL, n.threads, ...)

# S3 method for class 'bartOrdinal'
print(x, ...)

# S3 method for class 'bartOrdinal'
residuals(object, ...)

# S3 method for class 'bartNegbin'
extract(
    object, type = c("ev", "ppd", "bart", "loglik"),
    sample = c("train", "test"),
    combineChains = TRUE, ...)

# S3 method for class 'bartNegbin'
fitted(
    object, type = c("ev", "ppd", "bart"),
    ci.level = NULL, ...)

# S3 method for class 'bartNegbin'
predict(
    object, newdata,
    type = c("ev", "ppd", "bart"),
    offset = NULL,
    combineChains = TRUE, ci.level = NULL, n.threads, ...)

# S3 method for class 'bartNegbin'
print(x, ...)

# S3 method for class 'bartNegbin'
residuals(object, ...)

# S3 method for class 'bartHurdle'
extract(
    object, type = c("ev", "ppd", "prob", "bart", "loglik"),
    sample = c("train", "test"),
    combineChains = TRUE, ...)

# S3 method for class 'bartHurdle'
fitted(
    object, type = c("ev", "ppd", "prob", "bart"),
    ci.level = NULL, ...)

# S3 method for class 'bartHurdle'
predict(
    object, newdata,
    type = c("ev", "ppd", "prob", "bart"),
    offset = NULL,
    combineChains = TRUE, ci.level = NULL, n.threads, ...)

# S3 method for class 'bartHurdle'
print(x, ...)

# S3 method for class 'bartHurdle'
residuals(object, type = "ev", ...)

# S3 method for class 'bartMultinomial'
plot(x, plquants = c(0.05, 0.95), cols = NULL, ...)

# S3 method for class 'bartOrdinal'
plot(x, plquants = c(0.05, 0.95), cols = NULL, ...)

# S3 method for class 'bartNegbin'
plot(x, plquants = c(0.05, 0.95), cols = c("blue", "black"), ...)

# S3 method for class 'bartHurdle'
plot(x, plquants = c(0.05, 0.95), cols = c("blue", "black"), ...)

# S3 method for class 'bartMultinomial'
summary(object, ...)

# S3 method for class 'bartOrdinal'
summary(object, vars = c("thresholds", "sigma", "k", "tau"), ...)

# S3 method for class 'bartNegbin'
summary(object, vars = c("dispersion", "sigma", "k", "tau"), ...)

# S3 method for class 'bartHurdle'
summary(object, vars = c("sigma", "k", "tau"), ...)

# S3 method for class 'summary.bartHurdle'
print(x, ...)
```

## Arguments

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  following an analogous model description syntax as
  [`lm`](https://rdrr.io/r/stats/lm.html), optionally carrying one or
  more [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)
  terms that declare additional per-observation-modulated forests (see
  ‘Formula Terms’ below). For backward compatibility, can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) matrix
  `x.train` (`data` then plays the role of `y.train`); see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `formula` item for the full data-frame-to-predictor mapping, including
  sparse `Matrix::dgCMatrix` and
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  columns.

- data:

  An optional data frame, list, or environment containing predictors,
  used together with a formula `formula`. For backward compatibility,
  can also be the
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) vector
  `y.train`.

- test:

  An optional matrix or data frame with the same number of predictors as
  `data`, or `formula` in backward-compatibility mode; if column names
  are present, a matching algorithm is used. Refused, by name, when
  `formula` carries a
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  (‘Formula Terms’): an amplitude-coupled fit has no per-observation
  test replay in this version.

- subset:

  An optional vector of logicals or indices used to subset the data. Can
  be missing. A
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)
  term's basis is evaluated inside the same subsetted model frame as
  every other predictor (‘Formula Terms’); a basis given directly
  through
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `forests =` instead follows that argument's own subsetting rule (see
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)'s
  `basis` item).

- weights:

  An optional vector of weights to be used in the fitting process. For a
  gaussian response, BART fits a model with observations \\y \mid x \sim
  N(f(x), \sigma^2 / w)\\, where \\f(x)\\ is the unknown function. A
  probit fit (the default binary family here, and
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s) does
  not support weights, except that weights identically 1 are treated as
  absent; a `family = "logistic"` fit treats them as observation counts
  and requires positive integers. For a weighted logistic fit, the
  `"ppd"` draw at an observation with weight \\w\\ is the number of
  successes among \\w\\ trials, \\\mathrm{Binomial}(w, p)\\ with \\p\\
  the fitted probability.

- offset:

  An additive offset entering the family's mean or link, applied to
  training data; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `offset` item for the exact per-family formula. The same concept as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  `binaryOffset` (which only ever applies to that function's
  probit/logistic default family). Can be missing.

  For `predict` on a `bartMultinomial` fit, the same name carries a
  different shape: the per-category shift at the PREDICTED rows, an
  `nrow(newdata)` x K numeric matrix entering the K raw fits before the
  softmax, in the same layout as `bart2`'s own train-side `offset`
  above. It is never taken from the fit - these rows are not the fit's
  rows - so a fit trained with an `offset` REQUIRES one here, refusing
  by name rather than reporting the offset-free surface, which an
  explicit all-zero matrix asks for. `NULL` (the default) is the only
  value a fit trained without one takes by default, and passing the
  training offset back at the training rows reproduces `yhat.train`. A
  large negative entry is a shift like any other: the probabilities it
  drives are the exact softmax, so one far enough below its row's others
  underflows to zero (and its log to `-Inf`) as it would anywhere else.

  For `predict` on a `bartNegbin` fit, the same name is the log-exposure
  shift at the PREDICTED rows, entering the replayed log-odds latent
  \\\psi\\ additively before \\r e^{\psi}\\. `NULL` (the default)
  applies none.

  For `predict` on a `bartOrdinal` or `bartHurdle` fit the name is a
  formal only so that the fourth position means the same thing on all
  six `predict` methods: neither family has an out-of-sample offset
  channel, so any value but `NULL` is refused by name.

- offset.test:

  The equivalent of `offset` for test observations. Defaults to tracking
  `offset`: a scalar, or a vector already matching the test row count,
  is reused directly, and any other length is refused by name rather
  than silently recycled. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  ‘Test offset synchronization’ details for how a live sampler keeps the
  two linked after fitting.

- sigest:

  For continuous response models, an estimate of the residual standard
  deviation (residual standard error), \\\sigma\\, used to calibrate an
  inverse-chi-squared prior on the error variance. If not supplied, the
  least-squares estimate is derived instead. See `sigquant` for more
  information. Not applicable when \\y\\ is binary. Same concept as
  `sigma` in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md);
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) spells
  it `sigest` too.

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
  \\E(Y\|x) = f(x)\\ is away from \\\pm 0.5\\. The response is
  internally scaled to range from \\-0.5\\ to \\0.5\\. For binary \\y\\
  fit with the default probit link, `k` is the number of prior standard
  deviations \\f(x)\\ is away from \\\pm 3\\, the probit reference
  scale; `family = "logistic"` widens this to \\\pm \pi \sqrt{3}\\,
  three standard deviations of the standard logistic latent variable.
  The value can be either a fixed number, or a *hyperprior* of the form
  `chi(degreesOfFreedom = 1.5, scale = 2)`. The default, `NULL`, uses
  the value 2 for continuous responses and the `chi(1.5, 2)` hyperprior
  for binary ones, which centers the sampled `k` near the field-standard
  fixed value of 2 (prior median 1.9) while adapting to the data; pass
  `k = 2` for the fixed BART-package default, or `chi(1.5, Inf)` for the
  old improper prior. See
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s `k`
  item, and the ‘End-node prior parameter `k`’ details there, for the
  full calibration argument and its outlier-sensitivity caveat.

- prior.scale:

  Names the leaf calibration in response units instead of inheriting it
  from the range of the response: the prior standard deviation of the
  forest total \\f\\ at `k = 1`, so the prior standard deviation in
  force is `prior.scale / k` and the prior mean is the response
  transform's shift. `NA` (the default) leaves the family's own
  calibration in place and nothing changes. Useful when the fit is one
  block of a larger sampler and the response it is handed varies in
  scale between sweeps. For the leaf-model qualifications on what the
  named quantity bounds, and for the `sd` spelling and its refusal under
  a sampled `k`, see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  Details and
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md).

- power:

  Power parameter for tree prior. Default 2, Chipman, George, and
  McCulloch's empirical recommendation (see References) for the
  split-probability decay \\base (1 + depth)^{-power}\\, not derived
  from a formula.

- base:

  Base parameter for tree prior. Default 0.95, from the same source as
  `power`.

- split.probs:

  Prior and transition probabilities of variables used to generate
  splits. `NULL` (the default) requests equiprobability; can also be a
  numeric vector of length equal to the number of variables, or a named
  numeric vector with only a subset of the variables specified and a
  `.default` named value. Values given for factor variables are
  replicated for each resulting column in the generated model matrix.
  The symbol `num.vars` is rebound before execution to the number of
  columns in the model matrix. Collides with a supplied `tree.prior`
  (see below).

- dart:

  `TRUE` places a DART prior (a Dirichlet prior over the split-variable
  probabilities, inducing variable selection; Linero 2018) on the trees
  with default settings, using `power` and `base` for the growth
  probabilities; a prior created by
  [`dbartsPriors$dart`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
  supplies its own settings, overriding `power` and `base`. Cannot be
  combined with `split.probs`. Per-sample probabilities are returned in
  the `varprobs` element of the fit.

- n.trees:

  The number of trees in the sum-of-trees formulation.

- n.samples:

  The number of posterior samples requested; `n.samples %/% n.thin`
  draws are actually kept and returned. This is a SWEEP budget in the
  BayesTree tradition - the same role
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s `ndpost`
  plays (`ndpost / keepevery` are what is actually returned) and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)'s own
  `n.samples` plays too. It differs from
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)'s
  (and the lower-level
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s)
  `n.samples`, which already counts the draws one `run()` call returns
  and is NOT divided by thinning: that layer's sampler can be `run()`
  repeatedly with a different count each time, so its `n.samples` is a
  per-call return count rather than a one-shot total budget.

- n.burn:

  Number of MCMC iterations to be treated as burn in.

- n.chains:

  Integer specifying how many independent tree sets and fits should be
  calculated. Default 4;
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s default
  is 1, BayesTree's historical single-chain convention.

- n.threads:

  Integer specifying how many threads to use. Depending on the CPU
  architecture, using more than the number of chains can degrade
  performance for small/medium data sets. As such some calculations may
  be executed single threaded regardless.

  On `predict`, `n.threads` is a per-call worker count for the
  saved-tree replay, defaulting to the fit's own: the replay is
  partitioned by (chain, posterior draw), each partition writing its own
  rows and nothing being reduced across workers, so the answer is
  identical bit for bit at every value and only the time taken changes.

- combineChains:

  Logical; if `TRUE`, samples are returned collapsed across chains
  rather than in a chains-by-samples array. Default `TRUE`, as for
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md); the
  number of chains is tracked regardless on the returned object's
  `n.chains` component (see ‘Value’).

- n.cuts:

  The maximum number of possible values used in decision rules (see
  `useQuantiles`, and
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  ‘Decision Rules’ details). If a single number, it is recycled for all
  variables; otherwise must be a vector of length equal to the number of
  predictor columns. Fewer rules may be used if a covariate lacks enough
  unique values. Factor predictors take no rules from it, of either
  kind: a factor's grid follows its level table.

- useQuantiles:

  When `TRUE`, determine tree decision rules using estimated quantiles
  derived from the predictor variables. When `FALSE`, splits are
  determined using values equally spaced across the range of a variable.
  See [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  ‘Decision Rules’ details.

- n.thin:

  Every `n.thin` draw is kept to be returned to the user. Useful for
  “thinning” samples against serial correlation; see `n.samples` above
  for how the two combine.

- keepTrainingFits:

  If `TRUE` the draws of \\f(x)\\ for \\x\\ corresponding to the rows of
  the training data are returned.

- printEvery:

  As the MCMC runs, a message is printed every `printEvery` draws.

- printCutoffs:

  The number of cutoff rules printed to screen before the MCMC is run.
  Given a single integer, the same value is used for all variables. If
  0, nothing is printed.

- verbose:

  Logical; if `FALSE` suppress printing.

- keepTrees:

  Logical; must be `TRUE` in order to use `predict` with the result of a
  fit. Note that for models with a large number of observations or a
  large number of trees, keeping the trees can be very memory intensive.

- keepCall:

  Logical; if `FALSE`, the returned object has `call` set to
  `call("NULL")`, otherwise the call used to instantiate BART.

- samplerOnly:

  Builds the sampler from its arguments and returns it without running
  it. Useful to use the `bart2` interface in more complicated models.

- seed:

  Optional integer specifying the desired pRNG
  [seed](https://rdrr.io/r/base/Random.html). A
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand suffices
  for reproducibility; supplying `seed` instead gives reproducible
  results without touching R's stream. See
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  Reproducibility section.

- proposal.probs:

  Named numeric vector, optionally specifying the proposal rules and
  their probabilities. Elements should be `"birth_death"`, `"change"`,
  and `"swap"` to control tree change proposals, and `"birth"` to give
  the relative frequency of birth/death in the `"birth_death"` step. The
  default is the named vector
  `c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)`,
  identical to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s own
  default.

- monotone:

  Optional per-predictor monotonicity constraints, passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md); see
  its `monotone` argument. A named vector selects predictors by
  model-matrix column name, each direction one of
  `"+"`/`"increasing"`/`1` or `"-"`/`"decreasing"`/`-1` - matching is
  case-insensitive, so `"Increasing"` is accepted; only numeric and
  ordered columns are eligible. A constraint forces birth/death-only
  proposals and a fixed `k = 2`. `NULL` (the default) fits the
  unconstrained model.

- interactions:

  Optional interaction constraints, passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md); see
  its `interactions` argument and
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md).
  Built with `interactions(max.order = , groups = , forbid = )`: cap the
  distinct predictors per root-to-leaf path, confine interactions to
  declared groups, and/or forbid named predictors from co-occurring.
  `NULL` (the default) fits the unconstrained model.

- blocks:

  Optional block-additive constraint, passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md); see
  its `blocks` argument and
  [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md). Built
  with `blocks(groups = , trees.per.group = )`: confine each whole tree
  to one declared group of predictors so the ensemble is exactly a sum
  of per-group functions. The groups must form a total, disjoint
  partition of the predictors, validated at fit time. `NULL` (the
  default) fits the ordinary model.

- variance:

  Heteroscedastic BART (Pratola et al. 2020): passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Declares a second forest modeling the residual variance surface
  \\s^2(x)\\ as a product of scaled-inverse-chi-squared leaves, coupled
  to the mean forest through the per-observation precision. Either the
  plain selector naming which predictors drive the variance - a
  one-sided formula or column selector (`~ x1 + x2`), `TRUE` for every
  predictor, or `NULL` (the default) for a homoscedastic fit - or a
  [`varianceForest`](https://vdorie.github.io/dbarts/reference/varianceForest.md)`(vars = , n.trees = , base = , power = )`
  object, which additionally sets the variance forest's own tree count
  and tree-structure prior (defaulting to `40` trees and the mean
  forest's `tree.prior`, exactly as the plain selector does). Gaussian
  responses only; `resid.dist = student()` residuals and a grouped
  ([`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)) fit
  are refused with it too, unadjudicated rather than unsupported by
  design (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)). The
  fit gains `s.train`/`s.test` (posterior draws of \\s(x)\\), and
  `predict` attaches an `"s"` attribute carrying \\s(x)\\ for new data.
  \\s(x)\\ is the residual scale the reporting channels use:
  `extract(type = "loglik")`, the `type = "ppd"` draws, and
  [`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md)'s
  `mean.s` row all read it instead of the fixed `sigma` such a fit
  stores (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)).

- keepSampler:

  Logical that can be used to save the underlying
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  object even if `keepTrees` is false.

- warm.start:

  A previous fit whose forests seed the initial state instead of drawing
  trees from the prior: either a `bart` object fit with
  `keepSampler = TRUE` or a
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
  The donor must share this fit's predictors, tree count, and DART
  setting; a donor fit on a different cut grid is supported, its splits
  remapped onto this fit's grid (collapsing any the grid starves), the
  same way a data replacement remaps existing splits. Only its trees,
  `sigma`, and `k` carry over, and each chain starts from a different
  donor sample so multiple chains stay overdispersed. A warm start
  biases early draws toward the donor, so it shortens burn-in rather
  than removing it; keep a non-zero `n.burn`. See
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
  `installTrees` method for finer control.

- n.grow.sweeps:

  A non-negative integer; when positive, the initial forest is built by
  `n.grow.sweeps` sweeps of XBART-style grow-from-root (He, Yalov and
  Hahn 2019) instead of a draw from the prior, an init-then-refine
  workflow that reaches a good fit in far fewer sweeps before the exact
  sampler takes over. The default `0L` draws trees from the prior,
  exactly as before. Like `warm.start`, a grow-from-root start biases
  the early draws toward the grown fit, so it shortens burn-in rather
  than removing it; keep a non-zero (if smaller) `n.burn`. The posterior
  the sampler targets is unchanged - only the starting point moves.
  Constant-leaf models only, and mutually exclusive with `warm.start`
  (both request an initialization). See
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
  `growFromRoot` method.

- factors:

  How factor columns in a data frame enter the model, exactly as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). The
  default `"categorical"` keeps each unordered factor as a single
  predictor whose splits send a subset of its levels down each branch -
  a different model-matrix representation than
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s, which
  always expands every factor into indicator columns; `"indicators"`
  instead expands each factor into binary indicator columns, matching
  `bart`. A factor-predictor fit therefore changes if moved between the
  two interfaces.

- family:

  The response family. `"auto"` (the default) fits gaussian models to
  continuous responses and probit models to those coded 0/1, exactly as
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md); a
  two-level factor, logical, or two-level character response is also
  detected and fit as probit, an unordered factor (or character)
  response with three or more levels is detected and fit as multinomial,
  and an ordered factor with three or more levels is detected and fit as
  ordinal - each reporting the choice in a one-line message.
  `"gaussian"`, `"probit"`, and `"logistic"` force those fits directly;
  `"aft"`, `"hazard"`, `"hazard.probit"` (an accepted alias for
  `"hazard"`), `"hazard.logistic"`, `"multinomial"`, `"ordinal"`,
  `"nbinom"`, and `"hurdle.lognormal"` (alias `"twopart"`) reach the
  extended families described below.
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) fits only
  a narrower, appended `c("auto", "logistic", "aft")` of this token set
  directly; naming one of the remaining tokens to `bart`'s `family` is
  refused, pointing here.

  `family = "multinomial"` fits a K-category softmax classifier: K
  forests, one per category, coupled through an interleaved Polya-Gamma
  one-vs-rest augmentation. When `formula` is a formula, `data` is
  required (not optional as for every other family):
  `bart2(y ~ x1 + x2, family = "multinomial", ...)` with no `data`
  errors by name rather than resolving `y`/`x1`/`x2` from the calling
  environment. `y.train` must be a factor (or character, coerced with
  [`factor`](https://rdrr.io/r/base/factor.html)) with at least two
  levels and no `NA`s or `NA` levels; `K = nlevels(y.train)` is
  inferred, never given explicitly, and `levels(y.train)` is captured at
  fit time and carried on every K-shaped output (`yhat.train`'s trailing
  dimension, `fitted`'s columns, the `fitted(type = "class")` factor).
  Alternatively, `y.train` may be an n x K numeric matrix of nonnegative
  integer trial counts, one row per observation and one column per
  category, with row sum (trials) \\n_i \ge 1\\: category k's count is
  then a binomial(\\n_i\\, p_k) draw under the same one-vs-rest
  augmentation, `K = ncol(y.train)` is inferred, and the levels carried
  on every K-shaped output are `colnames(y.train)` when present,
  otherwise `as.character(seq_len(K))`. A one-hot count matrix (every
  \\n_i = 1\\) is the same model as the factor response above,
  reproduced bit for bit. The formula interface is supported alongside
  the matrix one: `bart2(formula, data, family = "multinomial")` with a
  factor (or character) left-hand side routes to the factor response
  above, and a `cbind(c1, ..., cK) ~ x` left-hand side (the same idiom
  `glm`'s binomial family uses) routes to the count-matrix response
  above, `K` and the levels taken from the `cbind` column names; the
  right-hand side is coded exactly as it is for every other family's
  formula fit, including `predict` on a data frame `newdata`. Under the
  default `family = "auto"`, a factor (or character) response with three
  or more levels is detected and fit as multinomial, reporting the
  choice in a one-line message; a two-level factor (or logical) response
  resolves to probit, and a numeric response is unchanged. A
  count-matrix (`cbind(c1, ..., cK) ~ x`) response is only ever
  multinomial when `family = "multinomial"` is given explicitly - it is
  never inferred. `weights`, `subset`, `samplerOnly`, `warm.start`,
  `n.grow.sweeps`, `dart` (or a DART `tree.prior`), `split.probs`,
  `monotone`, and `variance` are all refused with an error naming the
  limitation (an integer weight is already expressible as row-wise count
  replication in the response, a non-integer one has no exact
  augmentation sampler, and the K-forest engine copies only
  power/base/proposal-probability fields from the host sampler it
  briefly builds, so DART, fixed split probabilities, monotone direction
  constraints, and a variance forest never reach it). `offset` is
  accepted only as an n x K numeric matrix, entering the K forests' raw
  fits before the softmax, one column per category, in the same layout
  as a count-matrix `y.train`; a flat (length-n) offset is refused by
  name, since a common per-observation shift is the softmax's own null
  direction and is identically inert. `offset` is a TRAIN-side argument
  only: `offset.test` is refused by name too, and `yhat.test` is always
  computed WITHOUT any category offset, even when `offset` was supplied
  for training - a caller comparing an offset-fitted `yhat.train`
  against `yhat.test` should keep this asymmetry in mind. A category
  test offset on the fit-time `test` rows is a sampler-level capability
  only (a
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
  own `$setCategoryTestOffset` method, reached through
  `keepSampler = TRUE`, or the internal creators' own `offset.test`
  argument), not reachable from `bart2`; `predict`'s own `offset` is the
  supported route to an offset test surface, taking a matrix for the
  rows it is given. `test` is supported: an `x.test` of the same column
  structure as `x.train` reports the K-category softmax probabilities on
  the held-out rows as `yhat.test`, shaped and levels-named exactly like
  `yhat.train` (see ‘Value’). `keepTrees` is supported too: it retains
  every one of the K forests' trees so `predict` can replay them at new
  predictors afterward, reproducing `yhat.test` bitwise when `newdata`
  matches the fit-time `test`; without `keepTrees`, `predict` errors. A
  fit trained with an `offset` replays as well, given `predict`'s
  `offset` at the new rows, and refuses by name without it: no resident
  offset describes rows the fit never saw. The per-forest leaf scale
  follows its own K-dependent calibration (the K = 2 anchor is the
  logistic scale \\\pi\sqrt{3}\\ divided by \\\sqrt{2}\\, for the
  identified pairwise log-odds); `k` is read from the usual node prior
  exactly as for any other family, but the node prior's `node.scale`
  itself is NOT consulted - the multinomial engine calibrates its own.
  The fit's class is `"bartMultinomial"`, not `"bart"`: see ‘Value’
  below and the `extract`/`fitted`/`predict` methods for
  `bartMultinomial` objects.

  `family = "ordinal"` fits an ordered categorical response by a
  cumulative probit (a single forest, unlike multinomial's K): a latent
  \\z = f(x) + \epsilon\\ with \\\epsilon \sim N(0, 1)\\ is cut at
  ordered thresholds \\-\infty = \gamma_0 \< \gamma_1 = 0 \< \gamma_2 \<
  \ldots \< \gamma\_{K-1} \< \gamma_K = \infty\\, so \\P(Y = k \mid x) =
  \Phi(\gamma_k - f(x)) - \Phi(\gamma\_{k-1} - f(x))\\. The free
  thresholds are sampled by a one-at-a-time marginal Metropolis update
  (latents integrated out) under a log-gap prior \\N(0, 1.5^2)\\; at \\K
  = 2\\ the model is exactly the probit family. `y.train` should be an
  ordered factor ([`is.ordered`](https://rdrr.io/r/base/factor.html));
  its level order defines the category order, and `levels(y.train)` is
  captured at fit time and carried on every K-shaped output. Under the
  default `family = "auto"`, an ordered factor with three or more levels
  is detected and fit as ordinal, reporting the choice in a one-line
  message (an unordered factor stays multinomial - the two never
  overlap, the dispatch key is `is.ordered`); a two-level ordered factor
  is binary and resolves to probit. Given explicitly,
  `family = "ordinal"` also accepts an unordered factor or character
  response (with a message that the category order is taken from the
  level order - usually alphabetical, so order the levels yourself) and
  a numeric response, whose ordered levels are `sort(unique(y.train))`.
  `weights` are not supported (a weighted truncated-normal latent
  likelihood is not a coherent model); `warm.start` and `n.grow.sweeps`
  are refused with an error naming the limitation. `samplerOnly` is
  supported, unlike multinomial: it returns the fitted sampler unrun,
  the same object `fit` is under `keepSampler`. `test` and `keepTrees`
  are supported as for multinomial, and `predict` replays the saved
  trees at new predictors, differencing the cumulative probit at the
  stored per-draw thresholds. `rbart_vi` and `xbart` do not fit ordinal
  responses (grouped ordinal models and ordinal-scale cross-validation
  losses are recorded follow-ups). The fit's class is `"bartOrdinal"`,
  not `"bart"`: see ‘Value’ below.

  `family = "nbinom"` fits a non-negative integer (count) response by a
  negative-binomial model with the Polya-Gamma augmentation (a single
  forest, like logistic): the forest fits a log-odds latent \\\psi =
  f(x) + o\\, the count law is \\y_i \sim \mathrm{NB}(r,
  \mathrm{plogis}(\psi_i))\\ with dispersion \\r\\, and the mean is
  \\E\[y_i \mid x\] = r e^{\psi_i} = r e^{f(x_i) + o_i}\\. The offset
  \\o_i\\ therefore enters the mean multiplicatively and is a
  log-exposure (\\o_i = \log(\mathrm{exposure}\_i)\\). `y.train` must be
  a non-negative integer count, and `family = "nbinom"` is always
  explicit - a count carries no unambiguous class, so it is never
  inferred under `family = "auto"` (a numeric response there stays
  gaussian). The dispersion `r` is estimated by default on a capped
  positive-integer grid (\\\\1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30,
  50\\\\) under a renormalized \\\mathrm{gamma}(2, 0.1)\\ prior by a
  closed-form discrete full conditional; passing `dispersion` as a
  positive integer fixes \\r\\ instead (v1 ships the exact integer
  envelope, so a real fixed dispersion is refused - continuous
  dispersion is a recorded follow-up). Like probit, the latent scale is
  fixed at 1 and `weights` are not supported (exposure belongs in the
  offset, not in observation replication); `warm.start` and
  `n.grow.sweeps` are refused with an error naming the limitation.
  `samplerOnly` is supported, unlike multinomial: it returns the fitted
  sampler unrun, the same object `fit` is under `keepSampler`. `test`
  and `keepTrees` are supported as for ordinal, and `predict` replays
  the saved trees at new predictors, reporting mean counts \\r
  e^{\psi}\\ at the stored per-draw dispersion (a log-exposure
  `offset.test` enters \\\psi\\ additively). `rbart_vi` and `xbart` do
  not fit count responses (grouped negative-binomial models, real
  dispersion, and a Poisson family are recorded follow-ups). `y.train`
  is additionally capped at \\10^6\\: the dispersion grid's count
  histogram is sized from the largest count, so a larger one allocates
  without bound, and a count above it is refused at creation and at
  every response swap alike. The fit's class is `"bartNegbin"`, not
  `"bart"`: see ‘Value’ below.

  `family = "hazard"` and `family = "hazard.logistic"` fit a
  discrete-time survival hazard model as *ingestion sugar* over the
  binary families - the model adds no engine code (`"hazard.probit"` is
  an accepted alias for `"hazard"`). The response is the same `Surv`
  object or two-column `(time, status)` the `"aft"` family takes; each
  subject observed to time \\t_i\\ is person-period-expanded into its
  at-risk rows (one per period \\k = 1, \ldots, t_i\\), each carrying
  the subject's covariates, an ordinal period column (appended last, so
  trees split on contiguous periods and borrow strength across time),
  and the binary indicator \\y\_{ik} = \mathrm{status}\_i \cdot 1\\k =
  t_i\\\\. The discrete hazard is \\h(k \mid x) = g(f(x, k) + o)\\ for
  the chosen binary link \\g\\ (probit for `"hazard"`, logistic for
  `"hazard.logistic"`), and the fit is EXACTLY an ordinary
  `"probit"`/`"logistic"` fit on the expanded rows - a hazard fit is
  byte-identical, draw for draw, to the binary fit on the hand-expanded
  design with the same seed. The time grid defaults to the sorted
  distinct observed times (the BART `surv.bart` convention); `breaks`
  coarsens it (a single integer bins at that many quantiles, a boundary
  vector gives explicit right-closed intervals), and `max.rows` guards
  the expansion size. Offsets are on the link scale and replicate per
  subject; `weights` replicate and then follow the chosen binary
  family's policy (probit refuses non-unit weights, logistic requires
  integer counts). A `Surv` response requires the family to be requested
  explicitly (`family = "auto"` selects `"aft"`); the matrix interface
  is required (no formula, `subset`, or `test` - expand held-out
  subjects with `survivalProbabilities(fit, times, newdata = )`, which
  needs `keepTrees`). The returned object is an ordinary `"bart"` fit
  whose `$family` records the binary link (so every prediction generic
  transforms correctly with no special case) and whose `$periods`
  element carries the grid; `predict`/`extract` return the per-(subject,
  period) hazards on the expanded rows, and
  [`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
  produces survival-curve draws \\S(t \mid x) = \prod\_{k :
  \mathrm{periods}\[k\] \le t} (1 - h(k \mid x))\\. `rbart_vi` and
  `xbart` do not fit hazard responses (grouped/frailty hazard and a
  `cloglog` link are recorded follow-ups).

  `family = "hurdle.lognormal"` (alias `"twopart"`, which resolves and
  prints as `"hurdle.lognormal"`) fits a semicontinuous two-part
  (hurdle) model for a non-negative response with exact zeros: an
  OCCUPANCY probit fit of \\z = 1\\y \> 0\\\\ over all n observations,
  glued at report time to a POSITIVE-PART gaussian fit of \\\log y\\
  over the subset \\\\i : y_i \> 0\\\\. The two component fits share no
  parameters and are composed entirely R-side from two ordinary `bart2`
  fits at independently derived seeds - no engine code, and no coupling
  between the parts - so a shared variable-selection prior across the
  occupancy and positive parts is foreclosed by this composition (a
  recorded limitation, alongside a Duan-smearing retransformation and a
  heteroscedastic positive part, both follow-ups). `y.train` must be
  non-negative and finite and must carry at least one exact zero and one
  positive value (a response with no zeros, or none positive, is refused
  by name). By default `predict`/`fitted`/`extract` report the NATURAL
  (response) scale via posterior-predictive Monte Carlo: \\E\[y \mid x\]
  = P(y \> 0 \mid x)\\e^{f(x) + \sigma^2 / 2}\\ computed per posterior
  draw from the positive part's single \\\sigma\\ per draw, recycled
  across observations (the positive part is always homoscedastic; a
  heteroscedastic positive part is not reachable, see the recorded
  limitation above). `type = "prob"` returns the occupancy probability
  \\\pi(x)\\ through the probit link; `type = "link"`/`"log"` the
  positive part's log-scale linear predictor \\f(x)\\; `type = "ppd"`
  draws the proper BIMODAL predictive - a Bernoulli(\\\pi\\) spike at
  zero, else a lognormal draw - which the plain gaussian ppd path cannot
  produce. Fits currently use the matrix interface only
  (`bart2(x.train, y.train, family = "hurdle.lognormal")`); `weights`,
  `subset`, `offset`/`offset.test`, and `test` are all refused with an
  error naming the limitation - the positive-part fit is instead given
  the full training `x` as its own `x.test`, so its in-sample fitted
  values cover the zero rows it never trained on. `predict` requires
  `keepTrees = TRUE` and replays both saved forests at `newdata`.
  [`dbarts()`](https://vdorie.github.io/dbarts/reference/dbarts.md) does
  not fit this family - it composes two samplers, which only `bart2()`
  builds, so requesting it there is an error directing here. `rbart_vi`
  and `xbart` do not fit it either. The fit's class is `"bartHurdle"`,
  holding both component fits (`$occupancy`, a `"bart"` probit fit of
  the occupancy indicator; `$positive`, a `"bart"` gaussian fit of
  \\\log y\\ on the positive subset) under their own
  `extract`/`fitted`/`predict`/`residuals`/`print` methods.

- missing:

  How missing predictor values enter the model, exactly as in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). The
  default `"incorporate"` keeps them, every split rule learning a
  direction for `NA`; `"error"` rejects predictors containing `NA` - the
  historical behavior
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) always
  uses. The learned route is per column: an `NA` in `test` or `newdata`
  is refused, naming the column, wherever that column was complete in
  training.

- resid.dist:

  The residual error *law* for a continuous (gaussian) response, passed
  through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). The
  default [`gaussian()`](https://rdrr.io/r/stats/family.html) gives
  normal errors; `student(df)` fits outlier-robust Student-t errors by a
  Gaussian scale-mixture augmentation, with `student(df = nu)` fixing
  the degrees of freedom and `student()` estimating them on a capped
  grid. Only continuous gaussian responses carry it (an error
  otherwise). Two caveats: with a flexible tree-forest mean, tail
  inference is partially confounded with fit flexibility, so an
  estimated \\\nu\\ should be posterior-checked; and \\\sigma\\ is then
  the *conditional* scale, the marginal residual variance being
  \\\sigma^2\\\nu/(\nu-2)\\. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for
  the full description.

- dispersion:

  The negative-binomial dispersion \\r\\ (`family = "nbinom"` only;
  ignored otherwise), passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). `NA`
  (the default) estimates \\r\\ on a capped positive-integer grid under
  a renormalized \\\mathrm{gamma}(2, 0.1)\\ prior; a supplied value
  fixes it and must be a positive integer (v1 ships the exact integer
  envelope, so a real fixed dispersion is refused). Larger \\r\\
  approaches the Poisson limit; smaller \\r\\ is more overdispersed.

- breaks, max.rows:

  Discrete-time hazard controls (`family = "hazard"` only; ignored
  otherwise), passed through to
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  `breaks` sets the period grid: `NULL` (the default) uses the sorted
  distinct observed times, a single positive integer bins at that many
  quantiles, and a numeric boundary vector gives explicit right-closed
  intervals. `max.rows` (default \\10^7\\) refuses an over-large
  person-period expansion, naming the coarsening levers.

- tree.prior, node.prior, resid.prior:

  The full
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) prior
  objects, forwarded unevaluated so a bare
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
  vocabulary name inside them (`cgm`, `dart`, `normal`, `linear`, `gp`,
  `chisq`, `fixed`) resolves exactly as it does for `dbarts`:
  `node.prior = linear(columns = 1:3)`, `gp(...)`, and
  `resid.prior = fixed(1)` are reachable this way. `NULL` (the default
  for all three) instead builds the tree/node/residual priors from
  `power`/`base`/`split.probs`/`dart`, `k`/`prior.scale`, and
  `sigdf`/`sigquant` respectively, exactly as before this argument
  existed. Supplying an object alongside a shorthand that would
  otherwise help build the same prior is an error naming both:
  `tree.prior` collides with any of `power`/`base`/`split.probs`/`dart`;
  `node.prior` with `k`/`prior.scale`; `resid.prior` with
  `sigdf`/`sigquant`/`sigest`. A supplied `resid.prior` under a
  fixed-unit-scale family (`"probit"`, `"logistic"`, `"ordinal"`,
  `"nbinom"`, the hazard families, `"multinomial"`) is diagnosed the
  same way `sigest`/`sigdf`/`sigquant` already are - a
  `dbartsFamilyGatedWarning`, since the family overwrites it with
  `fixed(1)` regardless. `tree.prior`/`node.prior` are honored on every
  family, including both component fits of
  `family = "hurdle.lognormal"`.

- storage:

  Passed through to
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md);
  a formal on both `bart` and `bart2`. A character string selecting the
  precision of the internal running residual; see `dbartsControl`'s
  `storage` for the full description. The default `"double"` reproduces
  existing draws bitwise; `"single"` changes the numbers a fit returns -
  the sampled values shift slightly - in exchange for a bandwidth-bound
  speedup, and is currently supported only for continuous (gaussian)
  responses with constant leaves.

- updateState:

  Passed through to
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).
  A logical setting the default behavior for many
  [`sampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  methods with regards to the immediate updating of the cached state of
  the object; see `dbartsControl`'s `updateState`.

- ...:

  Present on the
  `bartMultinomial`/`bartOrdinal`/`bartNegbin`/`bartHurdle` methods
  below for S3 generic compatibility, but not silently discarded: a name
  foreign to the method called is refused by name, and any other
  unrecognized name warns (class `dbartsUnusedArgsWarning`). See
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `...` item.

- object:

  An object of class `bartMultinomial`, `bartOrdinal`, `bartNegbin`, or
  `bartHurdle`, returned by `bart2` under the corresponding `family`.

- newdata:

  Test data for prediction. Obeys the same rules as `data`/`test` but
  cannot be missing.

- type:

  The quantity to return; see ‘Value’ below for what each family's
  generics accept and produce. As on a `"bart"` fit
  ([`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s `type`
  item), `"response"` is a synonym for `"ev"` and `"link"` for `"bart"`,
  on every method here; a value a method does not offer is an error
  naming the ones it does, and one it names only to refuse is an error
  naming the reason.

- sample:

  Either `"train"` or `"test"`. It is `extract`'s own argument (and
  `fitted`'s, on the families that still carry it - fitted is always the
  training rows); refused by name on `predict`, whose stored train and
  test channels are `extract`'s `sample` instead.

- vars:

  Character vector of fields to gather, as
  [`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md)'s
  own `vars`; requested fields absent from `object` are silently
  dropped. `summary.bartOrdinal`'s `"thresholds"` contributes one draws
  variable per threshold \\\gamma_1, \ldots, \gamma\_{K-1}\\, the
  ordinal analog of `sigma`; `summary.bartNegbin`'s `"dispersion"`
  contributes the per-draw dispersion \\r\\, the count analog of
  `sigma`. `summary.bartHurdle` applies `vars` separately to its
  `$occupancy` and `$positive` component fits. `summary.bartMultinomial`
  has no `vars` formal at all: it always pools the per-category
  mean-probability channel (its only scalar convergence instrument), so
  a supplied `vars` is refused by name rather than silently ignored.

- x:

  An object of one of the four own-class families (`bartMultinomial`,
  `bartOrdinal`, `bartNegbin`, `bartHurdle`), as returned by `bart2`
  under the matching `family`, for `print`/`plot`; an object of class
  `summary.bartHurdle`, as returned by `summary` on a `bartHurdle` fit,
  for `print.summary.bartHurdle`.

- cols:

  For `plot.bartMultinomial`/`plot.bartOrdinal`: a vector of colors, one
  per category, tracing each category's training-mean predicted
  probability (`bartOrdinal`'s threshold panel, when drawn) or each
  category's own P2 panel point; `NULL` (the default) chooses colors
  automatically. For `plot.bartNegbin`/`plot.bartHurdle`: two colors, as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `cols` item - the first for a posterior interval's bars, the second
  for the reference line.

- plquants:

  For
  `plot.bartOrdinal`/`plot.bartNegbin`/`plot.bartHurdle`/`plot.bartMultinomial`,
  as [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `plquants` item: the lower and upper quantiles bounding each
  posterior-interval panel's bars.

- ci.level:

  For `fitted`/`predict` on a `bartMultinomial`, `bartOrdinal`, or
  `bartNegbin` fit, as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `ci.level` item. On a K-margined `type` (`"ev"`/`"ppd"` for the first
  two), the result keeps the category margin: an (observation \\\times\\
  K \\\times\\ 3) array of `est`/`ci.lower`/`ci.upper` rather than a
  3-column matrix, since categories are not poolable. On a
  non-K-margined `type` (`"bart"`; any type for `bartNegbin`) the result
  is the ordinary 3-column matrix. `bartHurdle`'s `ci.level` is
  documented under
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) directly,
  since it never carries a K margin. Refused, by name, on `residuals`
  (as [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
  `ci.level` item); on `bartMultinomial` and `bartOrdinal` also refused,
  by name, with `type = "class"` - a class prediction is a label rather
  than a quantity with a credible band.

## Formula Terms

A [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) call
written inside `bart2`'s `formula` - bare, or as one operand of `:` -
declares a second forest whose amplitudes multiply a basis, without a
separate [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
`forests =` list (the two routes are mutually exclusive: a term-bearing
formula given together with `forests =` is refused by name). The
canonical spelling uses colon sugar:


    y ~ x1 + x2 + z:forest(x1 + x2)             # numeric z: one amplitude
    y ~ x1 + x2 + factor(z):forest(x1 + x2)     # factor z: one amplitude per level
    y ~ x1 + x2 + (a + b):forest(x1 + x2)       # two-column basis, one forest

Each desugars to the general named form, `forest(x1 + x2, basis = ~ z)`
(respectively `~ factor(z)`, `~ cbind(a, b)`) - textually identical to a
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md) object
built directly for
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
`forests =`, so the term route and the constructor route share one
vocabulary and one validator; write `basis =` explicitly for anything
the colon sugar cannot say (an arbitrary formula, or an
already-evaluated vector/matrix). `z * forest(...)` is refused, naming
both explicit spellings (`z:forest(...)`, or `z + z:forest(...)` to also
keep `z` as a main effect) instead of guessing which one was meant.

In this position
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)'s
unnamed first slot is a SYMBOLIC predictor set, not `basis`:
`forest(x1 + x2)` means the variable set `c("x1", "x2")`, read from the
call and never evaluated - a bare symbol, or symbols joined by `+`, and
nothing else; a transformation, a removal, `.`, or `:`/`*` inside the
slot is refused by name, pointing at `vars = ` for what the grammar
cannot say. A factor name in that slot expands to every column the
factor produces (all of its indicator columns under
`factors = "indicators"`). Every NAMED argument inside a term-context
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) -
including `basis` itself, for the general form
`forest(x1 + x2, basis = ~ z)` - evaluates normally, exactly as it does
for [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
own `forests =`.

A term's basis is evaluated inside the fit's own model frame, so
`subset`, missing-data handling, and `data` scoping apply to it exactly
as they do to any predictor - unlike a basis given directly through
`forests =`, which is evaluated against the raw `data` and then
restricted to the subsetted rows (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)'s
`basis` item for that route's rule; the two reach the same rows whenever
`subset` keeps every row unaltered).

A handful of shapes are refused outright rather than silently
misinterpreted, each error naming the offending expression: a
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
anywhere other than a top-level additive term or a `:` operand (inside
[`I()`](https://rdrr.io/r/base/AsIs.html), a removal, or the formula's
left-hand side); a compound `:` operand containing anything other than
bare numeric/logical symbols joined by `+` (a factor or character
member, a transformation, a multi-way `:` chain, or a second
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) on the
other side); `test` given together with a term (an amplitude-coupled fit
has no per-observation test replay in this version); a family the
multiplier model does not support (`"aft"`, `"ordinal"`, `"nbinom"`, the
hazard families, `"multinomial"`, and `"hurdle.lognormal"`/`"twopart"` -
gaussian, probit, and logistic are the only families a term can join);
and a formula whose only right-hand-side content is the term, which
would leave nothing to split on. See ‘Value’ below, and
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)'s
`type = "forest"`, for reading the resulting per-forest fits back out.

## Value

Under the default and
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)-compatible
families (`"gaussian"`, `"probit"`, `"logistic"`, `"aft"`, `"hazard"`
and its aliases), `bart2` returns a list of class `"bart"`, documented
in full under
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s ‘Value’
section - the same `yhat.train`/`sigma`/`varcount`/... components,
including `forestFits`, `glue`, `bases`, `n.forests`, and a
`"forest.labels"` attribute when the formula carries a
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
(‘Formula Terms’ above): that term route is the only way `bart2` itself
reaches the per-forest channel, since it has no `forests =` formal of
its own. A heteroscedastic fit (`variance` supplied) additionally
carries `s.train`/`s.test`, posterior draws of the residual scale
\\s(x)\\ laid out as `yhat.train`/`yhat.test`; see the `variance`
argument above. A discrete-time hazard fit
(`family = "hazard"`/`"hazard.logistic"`) additionally carries
`periods`, the ordered period grid the person-period expansion used;
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
dispatches its survival-curve branch on it.

The remaining families - `"multinomial"`, `"ordinal"`, `"nbinom"`, and
`"hurdle.lognormal"`/`"twopart"` - return their own list class instead,
never `"bart"`, each with its own generics, described below: precisely
so a K-category or two-part fit cannot silently flow through a method
that assumes the single-forest shape.

`bart2(family = "multinomial")` returns a *different* list, of class
`"bartMultinomial"`. Components: `call`, `family` (`"multinomial"`),
`levels` (the resolved category names, length K: `levels(y.train)` for a
factor response, or `colnames(y.train)` - falling back to
`as.character(seq_len(K))` when unnamed - for a count-matrix response),
`K`, `n.chains`, `n.trees`, `y` (the original response: the factor
`y.train`, or the validated n x K count matrix when a count response was
supplied), and `yhat.train` - the posterior draws of the K softmax
probabilities, an array of dimension (`n.chains` \\\times\\, when
`n.chains > 1` and `combineChains = FALSE`) `n.samples` \\\times\\
number of training observations \\\times\\ K, with the resolved category
levels as the trailing dimension's names; `combineChains = TRUE` (the
default) folds the chain margin into the samples margin as usual.
`yhat.train` already holds PROBABILITIES, not a latent score - there is
no `"bart"`-scale component to convert, unlike the binary families. When
`test` was supplied, `yhat.test` is the same K softmax probabilities on
the held-out rows, dimensioned like `yhat.train` with the
training-observation margin replaced by the test one. `varcount` is the
per-sample per-category split-usage channel: each category forest's
per-draw variable counts, dimensioned like `yhat.train` but with the
number of training predictors in place of the number of observations,
and `colnames(x.train)` (when present) named on that margin. `fit` is
present whenever `keepTrees` is `TRUE` *or* `keepSampler` is set,
independent of `keepTrees`: it is the `dbartsSampler` whose K-forest
engine actually ran (one `bartcore_create`, not a discarded host), fully
mutable and readable on the channels the softmax gives meaning to -
`$setCounts`, `$setCategoryOffset`, `$setPredictor`, and the rest - and
refused by name on the ones it does not (`$setResponse`, `$setOffset`,
`$setSigma`, `$setCalibration`, `$setForestWeights`); `fit$storeState()`
followed by
[`save`](https://rdrr.io/r/base/save.html)/[`load`](https://rdrr.io/r/base/load.html)
restores a sampler `predict.bartMultinomial` can replay through.
`predict` still requires `keepTrees = TRUE` - a kept `fit` alone carries
no saved trees to replay - and errors otherwise.

Generics for a `"bartMultinomial"` fit: `fitted(object)` (`type = "ev"`,
the default) returns the posterior-mean n \\\times\\ K probability
matrix, columns named by the fit's captured `levels`;
`fitted(object, type = "class")` returns the argmax of that mean matrix
as a factor over the original levels - a classification convenience.
`extract(object, type = "ev", sample = "train")` returns the full
`yhat.train` draws, and `sample = "test"` returns `yhat.test` (an error
if the fit carries no test channel); `extract(object, type = "ppd")`
draws one category per posterior draw from its probability vector,
returned as an integer code (1-based, indexing the fit's captured
`levels`) in an array shaped like `"ev"` minus the K margin.
`extract`/`fitted`/`predict` with `type = "bart"`, and
`extract`/`predict` with `type = "forest"`, error naming the reason: the
run records only the identified softmax probabilities, and a category's
forest is a latent whose level is reproducibly structured yet not
identified (in the same sense as BCF's \\a\\), so a raw replay would
read as signal - the identified content is the log-ratio, which the logs
of the reported probabilities carry. `predict(object, newdata)` requires
`object` fit with `keepTrees = TRUE` (error otherwise); it codes
`newdata` to the training columns exactly as `predict.bart` does and
returns a levels-named (`n.chains` \\\times\\) `n.samples` \\\times\\
number of new observations \\\times\\ K probability array, the
`yhat.test`/`yhat.train` convention; `type = "ppd"` draws one category
per posterior draw the same way `extract(object, type = "ppd")` does;
`type = "bart"` and `type = "forest"` are named only to be refused, for
the non-identification reason above; `type = "class"` returns the same
argmax factor `fitted(object, type = "class")` does, over the
probability draws `predict` itself replays; paired with `ci.level` it is
refused by name instead - a class prediction is a label rather than a
quantity with a credible band. `offset` supplies the per-category shift
at the new rows, required on a fit trained with an `offset` and
reproducing `yhat.train` when the training offset is passed back at the
training rows. `residuals(object)` returns an n \\\times\\ K matrix, the
observed proportion (an indicator for a factor response,
`y / rowSums(y)` for a count-matrix one) minus the fitted probability in
`fitted(object)`, per category. `extract`'s `combineChains` (default
`TRUE`) reshapes `"ev"`/`"ppd"`/`"loglik"` to a combined or per-chain
layout, exactly as `extract.bart`'s does; refused by name on `fitted`
and `residuals`, which always summarize the combined draws.
`fitted`/`predict`'s `ci.level` returns an (observation \\\times\\ K
\\\times\\ 3) array of `est`/`ci.lower`/`ci.upper` instead of a 3-column
matrix, since a category margin survives alongside the observation one;
refused by name on `residuals` (as
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
`ci.level` item). `extract(object, type = "loglik")` is the multinomial
log-density of the observed row (its coefficient included, so it reduces
to \\\log p\_{i, y_i}\\ for a one-trial row), extract-only,
`sample = "test"` refused by name, shaped like `"ppd"` (the K margin
dropped) - loo/WAIC on it is leave-one-*row*-out, since the likelihood
unit is the whole count row. `forest`/`contribution` on
`extract`/`predict`, `sample` on `fitted`, `type` on `residuals`, and
`vars` on `summary` are all refused by name rather than silently
ignored: none is meaningful on this K-widened, non-identified-latent
shape. `plot(object)` adds a second panel to the per-category trace: the
posterior median and interval of the predicted probability of each
observation's OWN observed category against that median (a count
response with multi-trial rows instead plots the observed proportion
\\y\_{ik} / n_i\\ against the interval of the drawn \\p\_{ik}\\, one
point per row-category cell). `plotTree` and
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
are refused by name (a multinomial fit's trees live on its sampler; it
has no hazard channel). `as_draws_array`/`as_draws_df` (posterior, if
installed) expose the same per-category mean-probability channel
`summary(object)` pools, named `meanProb[<level>]` - this family has no
other scalar posterior parameter. `summary(object)` pools that same
per-category mean-probability channel into posterior mean/sd/quantiles
(R-hat/ESS when the posterior package is installed), the multinomial
analog of `summary.bart`'s \\\sigma\\/k/\\\tau\\ summary.

`bart2(family = "ordinal")` likewise returns its own list, of class
`"bartOrdinal"`. Components: `call`, `family` (`"ordinal"`), `levels`
(the ordered category names, length K), `K`, `n.chains`, `n.trees`, `y`
(the observed categories as an ordered factor over `levels`),
`yhat.train` (and `yhat.test` when `test` was supplied) - the posterior
draws of the K category probabilities, shaped and levels-named exactly
as the multinomial arrays above - `latent.train` (and `latent.test`) -
the corresponding draws of the latent \\f(x)\\ on the unit-variance
probit scale, shaped like a binary family's `yhat.train` -
`thresholds` - the posterior draws of the K - 1 finite thresholds
\\(\gamma_1 = 0, \gamma_2, \ldots)\\, dimensioned draws \\\times\\ (K -
1), the ordinal analog of gaussian's `sigma`, from which probabilities
at any latent value can be reconstructed - and `varcount`.
`thresholds.raw` (the per-draw thresholds in the internal layout
`predict` consumes) is present only under `keepTrees`. `fit` is present
whenever `keepTrees` is `TRUE` *or* `keepSampler` is set, independent of
`keepTrees`: it is the `dbartsSampler` whose engine actually ran, fully
mutable and readable like any sampler
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) returns,
and `fit$storeState()` followed by
[`save`](https://rdrr.io/r/base/save.html)/[`load`](https://rdrr.io/r/base/load.html)
restores a sampler `predict.bartOrdinal` can replay through.
`summary(object)` reports the thresholds alongside whatever
mean-function scale `vars` finds present, pooled into posterior
mean/sd/quantiles (R-hat/ESS when posterior is installed), the ordinal
analog of `summary.bart`'s \\\sigma\\/k/\\\tau\\ summary.

Generics for a `"bartOrdinal"` fit: `fitted(object)` returns the
posterior-mean n \\\times\\ K probability matrix (columns named by
`levels`); `fitted(object, type = "class")` the argmax category of that
mean as an ordered factor; `fitted(object, type = "bart")` the
posterior-mean latent. `extract(object, type = "ev")` returns the
probability draws, `type = "bart"` (alias `"link"`) the latent draws,
and `type = "ppd"` draws one category per posterior draw (1-based codes
indexing `levels`). `predict(object, newdata)` requires
`keepTrees = TRUE` and returns the probability draws at `newdata`
(`type = "bart"` the replayed latent; `type = "ppd"` category draws;
`type = "class"` the argmax ordered factor
`fitted(object, type = "class")` returns, over the replayed draws;
paired with `ci.level` it is refused by name instead, for the same
reason `bartMultinomial` refuses it); when `newdata` matches the
fit-time `test` it reproduces `yhat.test`. `residuals(object)` returns
the n \\\times\\ K matrix of observed-indicator minus fitted
probability. `extract`'s `combineChains` (default `TRUE`) and
`fitted`/`predict`'s `ci.level` work as for `bartMultinomial` above,
including the (observation \\\times\\ K \\\times\\ 3) interval shape for
a K-margined `type`, the refusal of both on `residuals`, and the refusal
of `combineChains` on `fitted`/`residuals`.
`extract(object, type = "loglik")` is \\\log P(y_i \mid \eta,
\gamma)\\ - exactly the stored category probability at the observed
level, since `yhat.train` already holds the cumulative-probit
difference - extract-only, `sample = "test"` refused by name, shaped
like `type = "ppd"`. `forest`/`contribution` on `extract`/`predict` and
`sample` on `fitted` are refused by name (this family has a single
forest and `fitted` is always the training rows). `plot(object)` traces
the free thresholds \\\gamma_2, \ldots, \gamma\_{K-1}\\ (\\\gamma_1 =
0\\ is pinned and carries no information; at K = 2 there is none, and
the plot degrades to one full-device panel) alongside the
per-observation latent interval, observations ordered by median
\\\eta\\, coloured by observed level, with dashed reference lines at the
posterior-median thresholds. `plotTree` and
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
are refused by name. `as_draws_array`/`as_draws_df` default to
`vars = c("thresholds", "sigma", "k", "tau")`, matching `summary`'s own
default; `thresholds` contributes `threshold[1]` (the pinned 0) through
`threshold[K - 1]`.

`bart2(family = "nbinom")` likewise returns its own list, of class
`"bartNegbin"`. Components: `call`, `family` (`"nbinom"`), `n.chains`,
`n.trees`, `y` (the observed counts), `yhat.train` (and `yhat.test` when
`test` was supplied) - the posterior draws of the mean counts \\\mu = r
e^{f(x) + o}\\, shaped like a binary family's `yhat.train` -
`latent.train` (and `latent.test`) - the corresponding draws of the
log-odds latent \\\psi = f(x) + o\\ - `dispersion` - the per-draw
dispersion \\r\\, the count analog of gaussian's `sigma` - and
`varcount`. `dispersion.raw` (the per-draw \\r\\ in the internal layout
`predict` consumes) is present only under `keepTrees`. `fit` is present
whenever `keepTrees` is `TRUE` *or* `keepSampler` is set, independent of
`keepTrees`: it is the `dbartsSampler` whose engine actually ran, fully
mutable and readable - `$getDispersion()` on it answers with the fit's
own \\r\\ - and `fit$storeState()` followed by
[`save`](https://rdrr.io/r/base/save.html)/[`load`](https://rdrr.io/r/base/load.html)
restores a sampler `predict.bartNegbin` can replay through.
`summary(object)` reports the dispersion draws alongside whatever
mean-function scale `vars` finds present, pooled into posterior
mean/sd/quantiles (R-hat/ESS when posterior is installed), the count
analog of `summary.bart`'s \\\sigma\\/k/\\\tau\\ summary.

Generics for a `"bartNegbin"` fit: `fitted(object)` returns the
posterior-mean count per observation; `fitted(object, type = "bart")`
the posterior-mean log-odds latent; `fitted(object, type = "ppd")` a
Monte Carlo mean over `extract(object, type = "ppd", sample = "train")`
draws. `extract(object, type = "ev")` returns the mean-count draws,
`type = "bart"` the latent draws, and `type = "ppd"` draws one count per
posterior draw from \\\mathrm{NB}(r, \mathrm{plogis}(\psi))\\;
`sample = "test"` selects the test channel. `predict(object, newdata)`
requires `keepTrees = TRUE` and returns the mean-count draws at
`newdata` (`type = "bart"` the replayed latent; `type = "ppd"` count
draws), with an optional log-exposure `offset`; when `newdata` matches
the fit-time `test` and no offset is given it reproduces `yhat.test`.
`residuals(object)` returns the observed count minus the posterior-mean
count per observation. `extract`'s `combineChains` and
`fitted`/`predict`'s `ci.level` (a plain 3-column matrix here - this
family has no K margin) work as for
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md), including
`combineChains`'s refusal on `fitted`/`residuals` and `ci.level`'s on
`residuals`. `extract(object, type = "loglik")` is
\\\mathrm{dnbinom}(y_i, \mathrm{size} = r_s, \mathrm{mu} =
\mu\_{s,i})\\, extract-only, `sample = "test"` refused by name.
`forest`/`contribution` on `extract`/`predict` and `sample` on `fitted`
are refused by name. `plot(object)` traces the integer-gridded
dispersion \\r\\ (a step plot; there is no burn-in channel to bridge
from, since `bart2` drives one `n.burn`/`n.samples` sweep for this
family) alongside
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own
gaussian observed-vs-fitted panel, applied to the counts. `plotTree` and
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
are refused by name. `as_draws_array`/`as_draws_df` default to
`vars = c("dispersion", "sigma", "k", "tau")`, matching `summary`'s own
default.

`bart2(family = "hurdle.lognormal")` returns its own list, of class
`"bartHurdle"`, holding both component fits (`$occupancy`, a `"bart"`
probit fit of the occupancy indicator; `$positive`, a `"bart"` gaussian
fit of \\\log y\\ on the positive subset) plus `call` and `family`
(`"hurdle.lognormal"`). `extract`/`fitted`/`predict`/`residuals`/`print`
methods for `"bartHurdle"` compose the two components: `type = "ev"`
(the default) is the natural-scale posterior mean or draws described
under the `family` argument above, `type = "prob"` the occupancy
probability, `type = "link"`/`"log"` the positive part's log-scale
linear predictor, and `type = "ppd"` the bimodal predictive draw.
`ci.level` on `fitted`/`predict` works as for
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md), including
its refusal on `residuals`; `combineChains` is likewise refused by name
on `fitted`/`residuals`, and `fitted`'s own `sample` formal is gone
outright - a hurdle fit has no separate test channel, so `fitted` is
always the training rows and a supplied `sample` is refused by name, as
it already was on `extract`. `residuals(object)` returns the observed
response minus `fitted(object)` on the natural scale.
`extract(object, type = "loglik")` is the density of \\y_i\\ at its own
natural scale (a \\-\log y_i\\ Jacobian against the components'
log-scale channels, so it is comparable to any other model OF y, and
differs from a log-y computation by the constant \\-\sum\_{y_i \> 0}
\log y_i\\): \\\log(1 - \pi\_{s,i})\\ where \\y_i = 0\\, else \\\log
\pi\_{s,i} + \mathrm{dnorm}(\log y_i, f\_{s,i}, \sigma_s, \log =
\mathrm{TRUE}) - \log y_i\\, with NO truncation (the lognormal positive
part's support is already \\(0, \infty)\\, unlike a count hurdle's
positive part, which would have to be truncated at zero - the two rules
must not be conflated); it is NOT the sum of the two components' own
`"loglik"` channels, extract-only, `sample = "test"` refused by name.
\\\pi\\ here is \\P(Y \> 0 \mid x)\\, the complement of the zero
probability other implementations report (brms's `hu`), so the zero
branch is \\\log(1 - \pi)\\ rather than \\\log \mathrm{hu}\\; a formula
ported from one of those must be inverted. `forest`/`contribution` on
`extract`/`predict` are refused by name (each component is a single
forest). `plot(object)` draws four panels: the positive component's
sigma trace; the occupancy probability \\\pi(x)\\; the positive part on
the scale it fit (\\\log y\\ over the \\y \> 0\\ rows); and the composed
natural-scale mean \\E\[y \mid x\] = \pi e^{f + \sigma^2/2}\\ over all n
rows, the only panel showing the model this family exists for.
`plotTree` and
[`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
are refused by name, naming `object$occupancy$fit`/`object$positive$fit`
as the route to the trees. `as_draws_array`/`as_draws_df` return the
union of both components' present scalar fields, labelled
`occupancy.<field>`/`positive.<field>` (a dot, not a bracket, since
posterior parses a bracket as an index) - the same two blocks
`summary(object)` prints under. `summary(object)` reports both
components - `$occupancy`'s and `$positive`'s own `summary.bart` tables,
under their own headers - rather than pooling them into one, since the
two component fits share no parameters.

## References

Chipman, H., George, E., and McCulloch, R. (2010) BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, **4**(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285) .

He, J., Yalov, S., and Hahn, P.R. (2019) XBART: accelerated Bayesian
additive regression trees. *Proceedings of the 22nd International
Conference on Artificial Intelligence and Statistics (AISTATS)*, PMLR
**89**, 1130–1138.

Linero, A.R. (2018) Bayesian regression trees for high-dimensional
prediction and variable selection. *Journal of the American Statistical
Association*, **113**(522), 626–636.

Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch, R.E. (2020)
Heteroscedastic BART via multiplicative regression trees. *Journal of
Computational and Graphical Statistics*, **29**(2), 405–417.

## Author

Hugh Chipman: <hugh.chipman@gmail.com>, Robert McCulloch:
<robert.mcculloch1@gmail.com>, Vincent Dorie: <vdorie@gmail.com>.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) for the
BayesTree-compatible interface and the shared `"bart"`-class
Value/Details (Decision Rules, end-node `k`, Generics, Saving,
Reproducibility, Extracting Trees).

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for the
mutable sampler bart2 builds on,
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md) for the
constructor
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)
formula terms desugar to, and
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)/[`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)
for the other formula-first entry points.

## Examples

``` r
## binary response with a logistic link
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
#>  prior on k: chi with 1.500000 degrees of freedom and 2.000000 scale
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
#> total seconds in loop: 0.001514
#> 
#> Tree sizes, last iteration:
#> [1] 2 2 2 2 4 3 2 2 3 1 2 2 2 2 3 2 2 2 
#> 3 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 13) (2: 12) 
#> DONE BART
#> 

## a forest() formula term: a second forest whose amplitude is modulated
## by a binary z, alongside the fit's own (prognostic) forest - the
## Bayesian causal forest shape, declared inline rather than through
## dbarts()'s forests =
set.seed(1)
n <- 60L
x1 <- runif(n); x2 <- runif(n)
z  <- rbinom(n, 1L, 0.5)
y  <- x1 + z * (1 + x2) + rnorm(n, 0, 0.2)
fit.bcf <- bart2(y ~ x1 + x2 + z:forest(x1 + x2),
                 n.samples = 20L, n.burn = 20L, n.chains = 1L,
                 n.trees = 10L, n.threads = 1L)
#> 
#> Running BART with numeric y
#> 
#> number of trees: 10
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 1
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.011208
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  proposal probabilities: birth/death 0.50, swap 0.10, change 0.40; birth 0.50
#> data:
#>  number of training observations: 60
#>  number of test observations: 0
#>  number of explanatory variables: 2
#>  init sigma: 0.774053, curr sigma: 0.774053
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) 
#> Running mcmc loop:
#> total seconds in loop: 0.001658
#> 
#> Tree sizes, last iteration:
#> [1] 2 2 2 3 3 2 2 1 2 3 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 4) (2: 8) 
#> DONE BART
#> 
fit.bcf$n.forests
#> [1] 2
modulatingFit <- extract(fit.bcf, type = "forest", forest = 2L)
```
