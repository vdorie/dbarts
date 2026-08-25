# Bayesian Additive Regression Trees

BART is a Bayesian “sum-of-trees” model in which each tree is
constrained by a prior to be a weak learner.

- For numeric response \\y = f(x) + \epsilon\\, where \\\epsilon \sim
  N(0, \sigma^2)\\.

- For binary response \\y\\, \\P(Y = 1 \mid x) = \Phi(f(x))\\, where
  \\\Phi\\ denotes the standard normal cdf (probit link).

- `bart` is the BayesTree-compatible surface: besides the default
  gaussian/probit response it reaches `family = "logistic"` and
  `family = "aft"` directly (below), each packaging as an ordinary
  `"bart"` object.
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) reaches
  further response families - logistic, accelerated failure time,
  discrete-time survival hazard, K-category multinomial, ordered
  categorical, negative-binomial counts, and semicontinuous two-part
  responses - plus its own formula-term multi-forest interface; naming
  one of `bart2`'s remaining own-class families to `bart`'s `family` is
  refused, pointing there.

## Usage

``` r
bart(
    x.train, y.train, x.test = matrix(0.0, 0, 0),
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2.0, prior.scale = NA_real_,
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
    keepsampler = keeptrees,
    resid.dist = gaussian,
    subset = NULL,
    storage = c("double", "single"),
    family = c("auto", "logistic", "aft"))

# S3 method for class 'bart'
plot(
    x,
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    ...)

# S3 method for class 'bart'
predict(
    object, newdata, offset, weights,
    type = c("ev", "ppd", "bart", "forest"),
    combineChains = TRUE,
    n.threads,
    ci.level = NULL,
    forest = NULL,
    bases = NULL,
    ...)

extract(object, ...)
# S3 method for class 'bart'
extract(
    object,
    type = c("ev", "ppd", "bart", "loglik", "trees", "forest"),
    sample = c("train", "test"),
    combineChains = TRUE,
    forest = NULL,
    contribution = FALSE, ...)

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
  numeric a continuous response model is fit (normal errors). If
  `y.train` is a two-level factor (or logical, or two-level character)
  or has only values 0 and 1, then a binary response model with a probit
  link is fit. A factor with three or more levels is an error directing
  to [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s
  `family = "multinomial"` (unordered) or `family = "ordinal"`
  (ordered), neither of which `bart` fits.

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

- k:

  For numeric \\y\\, `k` is the number of prior standard deviations
  \\E(Y\|x) = f(x)\\ is away from \\\pm 0.5\\. The response (`y.train`)
  is internally scaled to range from \\-0.5\\ to \\0.5\\. For binary
  \\y\\ fit with the default probit link, `k` is the number of prior
  standard deviations \\f(x)\\ is away from \\\pm 3\\, the probit
  reference scale; `family = "logistic"` widens this to \\\pm \pi
  \sqrt{3}\\, which is three standard deviations of the standard
  logistic latent variable (\\\pi / \sqrt{3}\\ each), the same three-sd
  span probit's \\\pm 3\\ covers. In both cases, the bigger \\k\\ is,
  the more conservative the fitting will be. The value can be either a
  fixed number, or a *hyperprior* of the form
  `chi(degreesOfFreedom = 1.5, scale = 2)`. `bart`'s own default is the
  fixed value 2; see
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s `k`
  item for its `NULL` default, a `chi(1.5, 2)` hyperprior on binary
  responses. The default of 2 for continuous responses follows Chipman,
  George, and McCulloch's argument (see References) that with node prior
  standard deviation \\\sigma\_\mu = 0.5 / (k \sqrt{m})\\ for \\m\\
  trees, \\k\\ prior standard deviations of \\f(x)\\ span the whole
  coded response range regardless of \\m\\ – so \\k = 2\\ places that
  range at roughly a 95% prior interval. The range scaling itself is an
  outlier-sensitive convention shared with BayesTree/bartMachine:
  extreme \\y\\ values compress the effective prior on everything else,
  for which the published workaround is to log-transform or winsorize
  such values before fitting, or to let `k` adapt via the `chi`
  hyperprior.

- prior.scale:

  Names the leaf calibration in response units instead of inheriting it
  from the range of `y.train`: the prior standard deviation of the
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

- splitprobs:

  Prior and transition probabilities of variables used to generate
  splits. Can be missing/empty/`NULL` for equiprobability, a numeric
  vector of length equal to the number variables, or a named numeric
  vector with only a subset of the variables specified and a `.default`
  named value. Values given for factor variables are replicated for each
  resulting column in the generated model matrix. The symbol `numvars`
  is rebound before execution to the number of columns in the model
  matrix.

- binaryOffset:

  Used for binary \\y\\. When present, the model is \\P(Y = 1 \mid x) =
  \Phi(f(x) + \mathrm{binaryOffset})\\, allowing fits with probabilities
  shrunk towards values other than \\0.5\\. \\\Phi\\ is the link for the
  default probit family; a `family = "logistic"` fit uses `plogis`
  instead, as noted in the value section below.

- weights:

  An optional vector of weights to be used in the fitting process. For a
  gaussian response, BART fits a model with observations \\y \mid x \sim
  N(f(x), \sigma^2 / w)\\, where \\f(x)\\ is the unknown function. A
  probit fit (`bart`, or
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) with the
  default binary family) does not support weights, except that weights
  identically 1 are treated as absent; a `family = "logistic"` fit
  treats them as observation counts and requires positive integers. For
  a weighted logistic fit, the `"ppd"` draw at an observation with
  weight \\w\\ is the number of successes among \\w\\ trials,
  \\\mathrm{Binomial}(w, p)\\ with \\p\\ the fitted probability.

- ntree:

  The number of trees in the sum-of-trees formulation.

- ndpost:

  The number of posterior draws after burn in, `ndpost / keepevery` will
  actually be returned.

- nskip:

  Number of MCMC iterations to be treated as burn in.

- printevery:

  As the MCMC runs, a message is printed every `printevery` draws.

- keepevery:

  Every `keepevery` draw is kept to be returned to the user. Useful for
  “thinning” samples.

- keeptrainfits:

  If `TRUE` the draws of \\f(x)\\ for \\x\\ corresponding to the rows of
  `x.train` are returned.

- usequants:

  When `TRUE`, determine tree decision rules using estimated quantiles
  derived from the `x.train` variables. When `FALSE`, splits are
  determined using values equally spaced across the range of a variable.
  See details for more information.

- numcut:

  The maximum number of possible values used in decision rules (see
  `usequants`, details). If a single number, it is recycled for all
  variables; otherwise must be a vector of length equal to
  `ncol(x.train)`. Fewer rules may be used if a covariate lacks enough
  unique values.

- printcutoffs:

  The number of cutoff rules to printed to screen before the MCMC is
  run. Given a single integer, the same value will be used for all
  variables. If 0, nothing is printed.

- verbose:

  Logical; if `FALSE` suppress printing.

- nchain:

  Integer specifying how many independent tree sets and fits should be
  calculated. Default 1 for `bart`, BayesTree's historical single-chain
  convention;
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s (and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)'s)
  default is 4.

- nthread, n.threads:

  Integer specifying how many threads to use. Depending on the CPU
  architecture, using more than the number of chains can degrade
  performance for small/medium data sets. As such some calculations may
  be executed single threaded regardless.

- combinechains, combineChains:

  Logical; if `TRUE`, samples from every chain are collapsed into a
  single matrix/vector of `nchain` \\\times\\ `ndpost` rows (see ‘Value’
  for the row order); if `FALSE`, they are kept in arrays of dimensions
  equal to `nchain` \\\times\\ `ndpost` \\\times\\ number of
  observations. Default `TRUE` across `bart`,
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md), and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md); the
  number of chains is tracked regardless on the returned object's
  `n.chains` component (see ‘Value’).

- keeptrees:

  Logical; must be `TRUE` in order to use `predict` with the result of a
  `bart` fit. Note that for models with a large number of observations
  or a large number of trees, keeping the trees can be very memory
  intensive.

- keepcall:

  Logical; if `FALSE`, returned object will have `call` set to
  `call("NULL")`, otherwise the call used to instantiate BART.

- seed:

  Optional integer specifying the desired pRNG
  [seed](https://rdrr.io/r/base/Random.html). A
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand suffices
  for reproducibility; supplying `seed` instead gives reproducible
  results without touching R's stream. See Reproducibility section
  below.

- proposalprobs:

  Named numeric vector or `NULL`, optionally specifying the proposal
  rules and their probabilities. Elements should be `"birth_death"`,
  `"change"`, and `"swap"` to control tree change proposals, and
  `"birth"` to give the relative frequency of birth/death in the
  `"birth_death"` step. The default is
  `c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)`.

- keepsampler:

  Logical that can be used to save the underlying
  [`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  object even if `keeptrees` is false.

- storage:

  Passed through to
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md);
  a formal on both `bart` and
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md). A
  character string selecting the precision of the internal running
  residual; see `dbartsControl`'s `storage` for the full description.
  The default `"double"` reproduces existing draws bitwise; `"single"`
  changes the numbers a fit returns - the sampled values shift
  slightly - in exchange for a bandwidth-bound speedup, and is currently
  supported only for continuous (gaussian) responses with constant
  leaves.

- family:

  `bart`'s own `family` is a narrower, appended formal,
  `c("auto", "logistic", "aft")`: the first-token default `"auto"`
  reproduces today's behavior exactly (a binary `y.train` still resolves
  to probit), and `"logistic"`/`"aft"` are newly reachable directly,
  both packaging as an ordinary `"bart"` object (`"aft"` needs a
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) or two-column
  `(time, status)` `y.train`, as below). The other ten
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) tokens
  are refused BY NAME rather than falling through to
  [`match.arg`](https://rdrr.io/r/base/match.arg.html)'s generic “should
  be one of” message, which names neither the token nor `bart2`:
  `"multinomial"`, `"ordinal"`, `"nbinom"`, `"hurdle.lognormal"`, and
  `"twopart"` (the `"hurdle.lognormal"` alias) package as their own S3
  class or compose two samplers, so
  `bart2(x.train, y.train, family = ...)` is where they are fit;
  `"gaussian"`, `"probit"`, and `"hazard.probit"` are already reachable
  through `family = "auto"` plus the response, so they add no capability
  as separate `bart` tokens; and `"hazard"`/`"hazard.logistic"` need
  `breaks`/`max.rows`, which `bart` does not have.

- subset:

  A vector of logicals or indices used to subset the data; an appended
  formal defaulting to `NULL` (today's behavior, unchanged), forwarded
  to [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) the
  same way. See
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s
  `subset` item for how a
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  composes with it on that interface.

- object:

  An object of class `bart`, returned from either the function `bart` or
  `bart2`.

- newdata:

  Test data for prediction. Obeys all the same rules as `x.train` but
  cannot be missing.

- offset:

  For `predict`: an optional offset for `newdata`, applied the same way
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `offset.test` is at fit time. `NULL` (the default, when missing)
  applies none. Refused, by name, with `type = "forest"`: an offset
  shifts the recombination of the forests, never any one forest's own
  total, so it belongs to the arithmetic under ‘Value’ rather than to
  what comes back. On an amplitude-coupled fit's combined arms
  (`type = "ev"`, `"ppd"`, `"bart"`) it is that recombination, and
  enters the sum of trees before the family's link, exactly as it does
  for a single forest.

- sampleronly:

  Builds the sampler from its arguments and returns it without running
  it. Useful to use the `bart` interface in more complicated models.

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
  predictive distribution, `"loglik"` - for `extract` only, the
  log-likelihood of each training observation at each posterior draw,
  `"trees"` - a data frame with tree information for when model was fit
  with `keepTrees` equal to `TRUE`, and `"forest"` - for `extract` and
  `predict`, the per-forest channels of an amplitude-coupled
  multi-forest fit (see `forest`, `contribution`, and the
  `forestFits`/`glue`/`bases` components under ‘Value’); an error naming
  the reason on any other fit, a
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)
  multinomial one included, whose K forests are per-category latents
  rather than additive components of one location.
  `extract(type = "forest")` reads the stored in-sample channel, while
  `predict(type = "forest")` replays each forest at `newdata`; both
  report the raw per-forest total, leaving the recombination to the
  caller, where the combined arms (`"ev"`, `"ppd"`, `"bart"`) perform it
  given the bases at those rows (see `bases`). For `"ppd"`, a weighted
  logistic fit draws the number of successes among the observation-count
  weight, \\\mathrm{Binomial}(w_i, p_i)\\ (see `weights`), an aft
  (survival) fit draws on the log-time scale (the fit models \\\log
  T\\), and a heteroscedastic fit
  ([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `variance`) draws its noise at that observation's own
  \\s(x_i)/\sqrt{w_i}\\; the scalar `sigma` such a fit reports is a
  fixed unit residual times the range of the response and is not its
  residual scale. For `"loglik"`, gaussian fits evaluate \\y_i \mid x_i
  \sim N(\hat{f}(x_i), \sigma^2 / w_i)\\ in logs at each draw of \\f\\
  and \\\sigma\\ - or, for a heteroscedastic fit, \\y_i \mid x_i \sim
  N(\hat{f}(x_i), s^2(x_i) / w_i)\\ at each draw of \\f\\ and \\s\\, a
  `resid.dist = student()` fit the corresponding marginal \\t\_\nu\\
  density at that draw's \\\nu\\ (the fit's `$resid.df`) rather than the
  normal one, binary fits evaluate the Bernoulli log-likelihood of the
  fitted probability, multiplied for a weighted logistic fit by the
  observation-count weight \\w_i\\, and an aft fit contributes the log
  density for an event and the log survival tail \\\log P(T \> C)\\ for
  a right-censored observation, both on the log-time scale. When chains
  are combined the result is a samples-by-observations matrix directly
  consumable by WAIC/PSIS-LOO implementations such as those in the loo
  package; the chains-first convention is kept, so a per-chain array
  (`combineChains = FALSE`, dimension chains-by-samples-by-observations)
  is reordered to the draws-by-chains-by-observations that
  `loo::relative_eff` expects with `aperm(x, c(2, 1, 3))`. To synergize
  with [`predict.glm`](https://rdrr.io/r/stats/predict.glm.html),
  `"response"` can be used as a synonym for `"ev"` and `"link"` can be
  used as a synonym for `"bart"`; the
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)
  extended-family methods take the same two synonyms, each against its
  own set of types. For information on extracting trees, see the
  subsection below.

- sample:

  Either `"train"` or `"test"`. `"test"` is refused, by name, when
  `type` is `"forest"`: an amplitude-coupled fit has no test fits (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `forests =`).

- forest:

  For `extract(type = "forest")` and `predict(type = "forest")`: which
  forest(s) to return, by 1-based index or by margin name (`"forest1"`,
  `"forest2"`, ...); `NULL` (the default) returns every forest. The
  returned array always keeps the trailing forest margin, subset to the
  requested forests, even when only one is selected.

- bases:

  For `predict` on an amplitude-coupled multi-forest fit: the bases each
  forest's amplitudes multiply AT the predicted rows, which off the
  training rows only the caller knows. A length-K list, entry k either
  `NULL` (for a forest that declared no basis, whose amplitude
  multiplies an implicit all-ones column) or anything
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)'s
  `basis` accepts as a value - a numeric vector or matrix, a factor, a
  character or a logical vector - expanded by the same rule and required
  to have `nrow(newdata)` rows and the width that forest's amplitudes
  take. When exactly one forest carries a basis the bare value may be
  given on its own, without the list: for a Bayesian causal forest,
  `bases = cbind(1 - zstar, zstar)` is the counterfactual arm being
  predicted under. Unlike at fit time, a column of all zeros and a
  single-level factor are accepted, since a constant arm is the point of
  a counterfactual. `NULL` (the default) re-derives each basis from
  `newdata` when the fit declared it as a
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  (which then needs the variables that term names, coded against the
  fit's own factor levels), and is otherwise an error naming the forest.
  Refused on a single-forest fit, and with `type = "forest"`, which
  reports each forest's total before any basis.

- contribution:

  For `extract(type = "forest")` only: `FALSE` (the default) returns the
  raw per-forest total `forestFits` already carries (see ‘Value’).
  `TRUE` instead returns each selected forest's per-observation
  CONTRIBUTION to the fit, \\(\mathrm{basis}\_k \\ \mathrm{glue}\_k)
  \times \mathrm{forestFits}\_k\\, computed on demand from the stored
  `forestFits`, `glue`, and `bases` rather than stored itself.

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

  Additional arguments passed on to `plot` (via `plot.bart`), or to the
  sampler's `getTrees` method via `extract` when `type` is `"trees"`
  (`chainNums`, `sampleNums`, `treeNums`, `newdata`; see ‘Extracting
  Trees’ below). Not used in `predict`.

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
interpreted as just the polynomial part of the distribution. If \\\nu\\
is 1 and \\s\\ is \\\infty\\, the prior is “flat”.

For BART on binary outcomes, the degree of overfitting can be highly
sensitive to `k` so it is encouraged to consider a number of values. The
default hyperprior for binary `bart2`/`dbarts` fits, `chi(1.5, 2)`,
centers the sampled `k` near the field-standard fixed value of 2 (prior
median 1.9) while adapting to the data; pass `k = 2` for the fixed
BART-package default. Crossvalidation may still be helpful, and running
for a short time with a flat prior (`scale = Inf`) can show the range of
`k` values the data are consistent with.

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
[`load`](https://rdrr.io/r/base/load.html)ing a fitted BART object for
use with `predict` requires that it be fit with `keeptrees`/`keepTrees`
as `TRUE`, and that the sampler's tree state be captured before saving
by calling `storeState()` on the sampler: `bartFit$fit$storeState()`
(for `rbart_vi`, on each element of `fit`:
`lapply(rbartFit$fit, function(f) f$storeState())`). The state is not
captured automatically because it duplicates the trees on the R side and
can exceed the fit's own sample blocks; it is materialized only on
request. A fit saved without it reloads, but `predict`,
`extract(type = "trees")`, and `plotTree` all stop identically, each
with an error naming `storeState()`. The same convention covers
[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s
own-class families - `family = "multinomial"`, `"ordinal"`, and
`"nbinom"` - whose `$fit` is the sampler (K-forest or single-forest)
that actually ran; see
[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s Value
section for each. A `family = "hurdle.lognormal"` fit has no `$fit` of
its own: it is the pair `$occupancy` and `$positive`, each a `bart2` fit
in its own right, so both of their samplers need the call before saving:
`fH$occupancy$fit$storeState(); fH$positive$fit$storeState()`.

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

- `chain`, `sample`, `tree` - index variables; `chain` is omitted on a
  single-chain fit

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

A fit made with `keeptrees`/`keepTrees` `FALSE` but a sampler kept
anyway (`keepSampler = TRUE`, or a `bart2` fit with `fit` present per
‘Value’ below) still answers `extract(fit, "trees")`, the same fallback
[`plotTree`](https://vdorie.github.io/dbarts/reference/plotTree.md)
uses: the frame holds the sampler's CURRENT trees rather than a saved
history, so it carries no `chain`/`sample` column and reflects whatever
predictors and structure the sampler holds at the time of the call, not
a fixed draw.

## Value

`bart` and [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)
return lists assigned the class `bart` for the shared
gaussian/probit/logistic/aft/hazard\* families (`bart2` reaches further
response families under their own list classes and generics; see
[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s Value
section). For applicable quantities, `ndpost / keepevery` samples are
returned. In the numeric \\y\\ case, the list has components:

- `call`:

  The matched call, or `call("NULL")` when the fit was made with
  `keepcall`/`keepCall` equal to `FALSE`.

- `yhat.train`:

  A array/matrix of posterior samples. The \\(i, j, k)\\ value is the
  \\j\\th draw of the posterior of \\f\\ evaluated at the \\k\\th row of
  `x.train` (i.e. \\f^\*(x_k)\\) corresponding to chain \\i\\. When
  `nchain` is one or `combinechains` is `TRUE`, the result is collapsed
  down to a matrix: with \\m = \\`ndpost %/% keepevery` draws kept per
  chain, row \\r\\ is chain \\1 + (r - 1) \\/\\ m\\'s draw \\1 + (r - 1)
  \\\\ m\\ - chain 1's whole run first, then chain 2's, and so on
  (chain-major), not interleaved draw by draw. Every collapsed field a
  fit returns (`sigma`, `k`, `tau`, `varcount`, `ranef`, ...) uses this
  same row order, so row \\r\\ of one pairs with row \\r\\ of another -
  e.g. `sigma[r]` is the residual scale that produced `yhat.train[r, ]`.

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
  number of samples unless `nchain` is one or `combinechains` is `TRUE`,
  in which case it collapses to a vector in the same chain-major order
  as `yhat.train` (see above), so combined `sigma[r]` pairs with
  combined `yhat.train[r, ]`.

- `first.sigma`:

  Burn-in draws of `sigma`.

- `varcount`:

  A matrix with number of rows equal to the number of kept draws and
  each column corresponding to a training variable. Contains the total
  count of the number of times that variable is used in a tree decision
  rule (over all trees). Row order (when combined) is `yhat.train`'s
  chain-major order (see above).

- `varprobs`:

  Present only under a DART prior: the sampled split-variable
  probabilities, in the same layout as `varcount`. Each draw sums to one
  across variables.

- `forestFits`, `glue`, `bases`, `n.forests`:

  Present only for an amplitude-coupled multi-forest fit: built through
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `forests =`, or through
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)'s
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)
  formula term (see its ‘Formula Terms’ section); not reachable from
  `bart` directly, which has no multi-forest front door. `n.forests` is
  the forest count K. `forestFits` is a (`n.chains` \\\times\\, when
  uncombined) `n.samples` \\\times\\ number of training observations
  \\\times\\ K array: forest k's own RESPONSE-scale raw total at each
  draw, i.e. \\\mathrm{response.scale} \times f_k(x)\\, with NO
  amplitude (`glue`) folded in - the same quantity the sampler's
  `$getForestFits` reports, up to that one scalar. The trailing margin
  is named `forest1`, ..., `forestK`; a declaration's own names
  (`names(forests)` on the
  [`dbarts()`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  route; a `bart2`
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  supplies none) ride separately as a `"forest.labels"` attribute on the
  fit object when given. `glue` is the ragged (`n.chains` \\\times\\,
  when uncombined) `n.samples` \\\times\\ sum(q_k) matrix of amplitudes
  multiplying each forest's basis columns, forest-major (the layout
  `$getForestAmplitudes` documents); its `"forest"` attribute names each
  column's forest. `bases` is a length-K list of each forest's expanded
  basis matrix (no draw axis; `NULL` for a forest with no declared
  basis, whose amplitude then scales an implicit all-ones column).
  Together they reconstruct the training fit exactly: \\\hat y =
  \mathrm{response.shift} + \sum_k (\mathrm{bases}\_k \\
  \mathrm{glue}\_k) \times \mathrm{forestFits}\_k\\, with
  `response.scale`/`response.shift` read from the sampler's
  `$getCalibration`. See `extract`'s `type = "forest"`.
  `predict(object, newdata, type = "forest")` reports the same quantity
  at NEW rows - an (`n.chains` \\\times\\, when uncombined) `n.samples`
  \\\times\\ `nrow(newdata)` \\\times\\ K array, same trailing margin,
  same `forest` selection - replayed from the saved trees, and so
  requiring `keeptrees`/`keepTrees`. It carries no `glue`, no
  `response.shift` and no offset either, because off the training rows
  the bases are the caller's, so an offset is refused there rather than
  folded in and there is no `contribution` argument for the same reason.
  The identity itself is what `predict(object, newdata, type = "ev")`
  (or `"ppd"`, or `"bart"`) performs at those rows in one call, taking
  each `bases` entry's replacement - the basis that forest's amplitudes
  multiply at the new rows, for a Bayesian causal forest the \\(1 -
  z^\*, z^\*)\\ indicator pair of whichever assignment is being
  predicted under - through its own `bases` argument, or re-derived from
  `newdata` on the
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  route. Written out by hand it remains the fallback for a caller
  holding only these components and no sampler. A fit built through a
  [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term
  additionally carries `basis.terms`, the length-K list of declaring
  formulas and fit-time factor levels that re-derivation reads; it is
  absent when the bases arrived as values.

- `sigest`:

  The rough error standard deviation (\\\sigma\\) used in the prior.

- `resid.dist`, `resid.df`:

  The residual error law the fit was made under, `"gaussian"` or
  `"student"` (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `resid.dist`). `resid.df` is present only for a `student()` fit: the
  degrees of freedom \\\nu\\ each draw was conditioned on, in `sigma`'s
  layout - a fixed \\\nu\\ repeats its value, an estimated one gives
  that draw's grid value. `extract(type = "loglik")` reads both.

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
  modeled, i.e. there is a hyperprior. Collapses (when combined) in
  `sigma`'s chain-major order.

- `first.k`:

  Burn-in draws of `k`, if modeled.

- `binaryOffset`:

  Present only for a binary fit: the offset value used.

- `family`:

  The resolved response family (`"gaussian"`, `"probit"`, `"logistic"`,
  or `"aft"`). `predict`, `extract`, `fitted`, and `plot` use it to
  transform latent draws to probabilities. For an `"aft"` fit,
  predictions and fitted values are on the LOG-TIME scale: these
  functions return the linear predictor \\E\[\log T \mid x\]\\, exactly
  as for a gaussian fit of \\\log T\\, and never the time scale that
  `survreg` or `flexsurv` users might expect. The `type` aliases
  `"response"` and `"link"` are inert for `"aft"` - both return the same
  log-scale linear predictor; neither exponentiates. The
  posterior-median survival time is `exp` of the linear predictor, and
  [`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
  gives survival-probability draws. A discrete-time hazard fit
  (`family = "hazard"`/`"hazard.logistic"`, reachable only through
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md)) records
  the binary link here (`"probit"`/`"logistic"`) - the fit IS that
  binary fit on the expanded rows - and carries a separate `$periods`
  element (the ordered period grid); `predict`/`extract` return the
  per-(subject, period) hazard through the link, and
  [`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
  dispatches its survival-curve branch on the `$periods` marker.

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

Chipman, H., George, E., and McCulloch, R. (2010) BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, **4**(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285) .

Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble
Learning. Advances in Neural Information Processing Systems 19,
Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
<https://www.rob-mcculloch.org>

Friedman, J.H. (1991) Multivariate adaptive regression splines. *The
Annals of Statistics*, **19**, 1–67.

## Author

Hugh Chipman: <hugh.chipman@gmail.com>, Robert McCulloch:
<robert.mcculloch1@gmail.com>, Vincent Dorie: <vdorie@gmail.com>.

## See also

[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) for the
formula-first interface, further response families, and formula-term
multi-forest models.

[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
`$predictForests` for the same per-forest replay from a mutable sampler,
which is what `predict(type = "forest")` reads and what
`extract(type = "forest")` reports in sample.

[`pdbart`](https://vdorie.github.io/dbarts/reference/pdbart.md)

[`dbarts-embedding`](https://vdorie.github.io/dbarts/reference/dbarts-embedding.md)
and
[`vignette("dbarts-as-a-component", package = "dbarts")`](https://vdorie.github.io/dbarts/articles/dbarts-as-a-component.md)
for using BART as one block of a sampler you write, rather than as the
whole fit: these functions own their own MCMC loop, while
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) returns
a
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
that runs one sweep at a time.

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
#> total seconds in loop: 0.218644
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

## fit with missing predictor values, using missing = "incorporate"
## (the default): every split rule learns a direction for NAs
set.seed(1)
n <- 60L
x.na <- matrix(runif(n * 2), n, 2)
x.na[1:5, 1] <- NA
y.na <- x.na[,2] + rnorm(n, 0, 0.3)
y.na[is.na(x.na[,1])] <- y.na[is.na(x.na[,1])] + 2
fit.na <- dbarts(y.na ~ x.na, control = dbartsControl(
                   n.samples = 20L, n.burn = 20L, n.chains = 1L,
                   n.trees = 20L, n.threads = 1L))
samples.na <- fit.na$run()

## recipe: between-draws row subsetting on a latent family (e.g. principal
## stratification or mediation, where an outer Gibbs step redraws which rows
## belong to a stratum every iteration). $setActiveRows changes which rows
## are "in" the model from one run() to the next, without rebuilding the
## sampler and without losing any row's fitted value; see the "active" entry
## of \link{dbartsSampler-class} for the full semantics. Unlike zero case
## weights (gaussian only), this reaches probit, ordinal, logistic, nbinom
## and aft samplers - and an emptied stratum is a legal all-zeros mask
## rather than a case the caller must special-case and skip.
set.seed(3)
n <- 120L
x.strat <- matrix(runif(n * 2), n, 2)
y.strat <- rbinom(n, 1L, plogis(2 * x.strat[,1] - 1))
sampler.strat <- dbarts(x.strat, y.strat, family = "probit",
                        control = dbartsControl(n.chains = 1L, n.threads = 1L,
                          n.trees = 20L, updateState = FALSE))

stratum <- rbinom(n, 1L, 0.6)  # this iteration's stratum membership
sampler.strat$setActiveRows(stratum)
draws.strat <- sampler.strat$run(20L, 20L)

## the next outer iteration redraws the stratum; even a stratum that empties
## out is accepted and runs, rather than needing to be skipped
stratum <- rep(0, n)
sampler.strat$setActiveRows(stratum)
draws.empty <- sampler.strat$run(0L, 5L)
```
