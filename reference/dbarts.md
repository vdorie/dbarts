# Discrete Bayesian Additive Regression Trees Sampler

Creates a sampler object for a given problem which fits a Bayesian
Additive Regression Trees model. Internally stores state in such a way
as to be mutable.

## Usage

``` r
dbarts(
    formula, data, test, subset, weights, offset, offset.test = offset,
    verbose = FALSE, n.samples = 800L,
    tree.prior = cgm, node.prior = normal, resid.prior = chisq,
    resid.dist = gaussian,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL,
    interactions = NULL,
    blocks = NULL,
    variance = NULL, n.trees.variance = 40L,
    power.variance = NULL, base.variance = NULL,
    forests = NULL,
    control = dbarts::dbartsControl(), sigma = NA_real_, seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic", "aft", "ordinal",
               "nbinom", "hazard", "hazard.probit", "hazard.logistic",
               "hurdle.lognormal", "twopart"),
    missing = c("incorporate", "error"), dispersion = NA_real_,
    breaks = NULL, max.rows = 1e7)
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
  predictor-mutation surface accepts both whole-matrix replacement,
  `setPredictor(x)`, and the column-granular
  `setPredictor(x, column = j)`, which replaces a sparse column whole;
  replacing a sparse-backed column densifies its storage permanently.
  Per-observation replacement of a sparse-backed column, and `setData`,
  are fixed at creation. Sparse inputs are not supported by
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) or
  the `node.prior = linear()` and `gp()` leaf models (a sparse-backed
  test column cannot serve a designated leaf covariate either). A sparse
  or data-frame test set stays resident (rank-bitmap or dense per
  column, the same rule as training) rather than densifying at
  ingestion; only `predict` and `getTrees(newdata = )` materialize it to
  a dense matrix, since those entry points evaluate already-built trees
  on resolved values. A data frame may mix ordinary columns with
  [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  or `dgCMatrix` columns (assign them into the frame; they do not
  survive `data.frame(...)` or
  [`I()`](https://rdrr.io/r/base/AsIs.html)): sparse columns behave as
  in an all-sparse design, while dense columns keep categorical splits
  and linear-leaf designation. A
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  column is accepted here (the x/y, data-frame path) exactly as the
  other sparse columns are; a bare `formula` object cannot carry it
  through `model.frame` and refuses it, with a message to use the x/y
  interface instead. See
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

- resid.dist:

  The residual error *law* for a continuous (gaussian) response,
  orthogonal to `resid.prior` (which is the prior on the error
  variance). The default
  [`gaussian()`](https://rdrr.io/r/stats/family.html) gives the usual
  normal errors. `student(df)` instead fits outlier-robust Student-t
  errors by the classic Gaussian scale-mixture augmentation: each
  observation carries a latent precision \\\lambda_i \sim
  \mathrm{Gamma}(\nu/2, \nu/2)\\ so that \\\sqrt{w_i}\\\epsilon_i /
  \sigma \sim t\_\nu\\, and a gross outlier draws a small \\\lambda_i\\,
  downweighting its leverage on the fit and on \\\sigma\\.
  `student(df = nu)` fixes the degrees of freedom at `nu` (a positive
  number; `nu = 4` is a conventional robust default), while `student()`
  (equivalently `student(df = NULL)`) estimates them on a capped grid
  under a proper tail-bounding prior. Only a continuous gaussian
  response carries this: `student()` with a `"probit"`, `"logistic"`, or
  `"aft"` family is an error, as those have their own fixed latent
  scale. Two caveats: (a) with a flexible tree-forest mean, tail
  inference is partially confounded with fit flexibility, so an
  estimated \\\nu\\ should be posterior-checked rather than read as a
  pure noise property; (b) \\\sigma\\ is now the *conditional* scale, so
  the marginal residual variance is \\\sigma^2\\\nu/(\nu-2)\\, not
  \\\sigma^2\\. See “Response scaling” below - robust errors mitigate an
  outlier's leverage but not the range compression that outlier causes.

- proposal.probs:

  Named numeric vector or `NULL`, optionally specifying the proposal
  rules and their probabilities. Elements should be `"birth_death"`,
  `"change"`, and `"swap"` to control tree change proposals, and
  `"birth"` to give the relative frequency of birth/death in the
  `"birth_death"` step.

- monotone:

  Optional per-predictor monotonicity constraints (monotone BART;
  Chipman, George, McCulloch, and Shively 2022). The forest is
  constrained so that the fitted function is monotone increasing or
  decreasing in each named predictor. Supply a named vector selecting
  predictors by model-matrix column name, each element one of
  `"+"`/`"increasing"`/`1` (increasing) or `"-"`/`"decreasing"`/`-1`
  (decreasing); an unnamed integer or character vector of length equal
  to the number of columns assigns directions positionally, with `0` for
  unconstrained. Only numeric and ordered columns are eligible - a
  direction on a categorical (unordered factor) predictor is an error. A
  constraint forces birth/death-only tree proposals (an explicit
  non-default `proposal.probs` is then an error) and a fixed `k = 2` (an
  explicit `k` hyperprior is an error); linear and Gaussian-process
  leaves are not supported under the constraint. `NULL` (the default) or
  an all-zero vector fits the ordinary unconstrained model.

- interactions:

  Optional interaction constraints, built with
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md):
  `interactions(max.order = 2)` caps the number of distinct predictors
  that may appear together on any single root-to-leaf path,
  `interactions(forbid = list(c("a", "b")))` bars named predictors from
  ever sharing a path, and
  `interactions(groups = list(c("a", "b"), "c"))` confines interactions
  to declared groups (named predictors may co-occur only with their
  group-mates). The three may be combined. Named or indexed predictors
  resolve against the model matrix at fit time, where an unknown name,
  an empty group, a `max.order` below 1, or a reference to a dropped
  column is an error. `NULL` (the default) fits the ordinary
  unconstrained model.

- blocks:

  Optional block-additive constraint, built with
  [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md):
  `blocks(groups = list(c("a", "b"), "c"))` confines each whole tree to
  one declared group of predictors, so the ensemble is exactly a sum of
  per-group functions (a functional-ANOVA / grouped-GAMI decomposition).
  Unlike `interactions(groups = )`, which is a per-path allow-list, this
  is a static per-tree mask, so the groups must form a total, disjoint
  partition of the predictors - every predictor named exactly once.
  `trees.per.group` optionally fixes how many trees each group gets;
  otherwise the trees are spread as evenly as possible. Validated at fit
  time (an unpartitioned or doubly-named predictor, or a
  `trees.per.group` that does not sum to `n.trees`, is an error). `NULL`
  (the default) fits the ordinary model. May be combined with
  `interactions`. See
  [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md).

- variance:

  Optional heteroscedastic variance forest (HBART; Pratola, Chipman,
  George, and McCulloch 2020). When supplied, a second ensemble models
  the residual variance surface \\s^2(x)\\ as a product of
  scaled-inverse-chi-squared (multiplicative) leaves, so \\y_i =
  f(x_i) + s(x_i)\epsilon_i\\, coupled to the mean forest through the
  per-observation precision. `variance` selects which predictors drive
  the variance: a one-sided formula or character/integer column selector
  (`~ x1 + x2`), `TRUE` for every predictor, or `NULL` (the default) for
  the ordinary homoscedastic fit. Gaussian responses and constant leaves
  only; monotone constraints and the latent families are not supported.
  The per-tree leaf prior is calibrated from the residual
  (`resid.prior`) hyperparameters so that a constant variance surface
  reproduces the homoscedastic `sigma` posterior. The fit gains
  posterior draws `s.train`/`s.test` of \\s(x)\\, and `predict` attaches
  an `"s"` attribute with \\s(x)\\ at new predictors (requires
  `keepTrees`).

- forests:

  Optional list of
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  specifications declaring the ensembles the mean is a weighted sum of,
  instead of the single ensemble `NULL` (the default) fits. Each forest
  carries its own tree count, structure prior, leaf scale, column
  restriction and `interactions`/`blocks` constraints; a forest given a
  `basis` enters the mean multiplied by the amplitudes on that basis's
  columns. Any number of forests may be declared, each carrying a basis
  of any width; the two-forest case with a two-level factor basis is the
  Bayesian causal forest \\y = a \mu(x) + b_z \tau(x) + \epsilon\\: a
  prognostic forest \\\mu\\ over every predictor and a modulating forest
  \\\tau\\ over the columns its `vars` allows. Every forest past the
  first needs a `basis`, the amplitudes multiplying it being what
  distinguish it from the first. Gaussian, `"probit"` and `"logistic"`
  responses: under a latent family the forests combine into the index
  rather than into the mean, on the link's own fixed scale, so every
  forest's `sd` is stated in latent standard deviations, `sigma` is
  pinned and there is no response transform to re-anchor. `"aft"`,
  `"ordinal"` and `"nbinom"` are refused at creation, each naming what
  it is missing. The result is an ordinary `dbartsSampler`; per-forest
  fits and the per-forest amplitudes are read off the sampler rather
  than from the run's combined `train` channel. The columns a basis
  expands to are stored on the data object, so they are subset by
  `subset` and survive a sampler's re-creation. A basis declared here
  has nowhere to ride when `formula` is already a `dbartsData` object:
  `dbarts()` refuses that combination by name rather than silently
  discarding the declaration and fitting a single-forest model - build
  the data object with the bases already on it, `dbartsData(bases = )`,
  or declare them through
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
  which always takes a pre-built data object and installs the
  declaration, replacing whatever bases it carried. Options a two-forest
  model does not read - `monotone`, `variance`, a DART tree prior,
  `split.probs`, a linear or Gaussian-process node prior, a `k`
  hyperprior or non-default `k`, a non-default `proposal.probs`,
  Student-t residuals, grouped random effects, `storage = "single"`,
  per-column cut counts, and a `test` set - are refused at creation
  rather than ignored, as is any declaration the engine cannot honour.

- n.trees.variance, power.variance, base.variance:

  The number of variance trees (default `40`) and their tree-structure
  prior. `power.variance` and `base.variance` default to the mean
  forest's `tree.prior` values. Read only when `variance` is supplied.

- control:

  An object inheriting from `dbartsControl`, created by the
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  function.

- sigma:

  A positive numeric estimate of the residual standard deviation. If
  `NA`, a linear model is used with all of the predictors to obtain one.
  Same concept as `sigest` in
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)/`bart2`,
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md), and
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md).

- seed:

  Optional integer seed for the random number generator, a convenience
  mirror of `dbartsControl(seed = )`. When not `NA` it overrides the
  seed in `control`; the fitting-function wrappers
  ([`bart2`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)) accept
  the same argument.

- factors:

  How factor columns in a data frame enter the model. The default
  `"categorical"` keeps each unordered factor as a single predictor
  whose splits send a subset of its levels (at most 65535) down each
  branch, and codes ordered factors as ordinal; level tables are
  retained so that test data are coded identically. The number of
  categories is the declared level table, not the levels the training
  rows happen to take, so a declared but unobserved level keeps its own
  bin and test or replacement data carrying it are accepted
  ([`droplevels`](https://rdrr.io/r/base/droplevels.html) beforehand to
  model only the observed levels). `"indicators"` expands each factor
  into binary indicator columns, as previous versions always did and as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) still
  does.

- family:

  The response model. `"auto"` fits gaussian models to continuous
  responses and probit models to those coded 0/1, as always; a two-level
  factor, logical, or two-level character response is also detected and
  fit as probit, reporting the choice in a one-line message, while a
  factor (or character) response with three or more levels is an error
  directing to
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  `family = "multinomial"` (which `dbarts` does not fit). An explicit
  family that a factor response cannot support (e.g. `"gaussian"`) is
  also an error rather than a silent fit of the integer level codes.
  `"gaussian"` forces a continuous fit even for a 0/1 numeric response;
  `"probit"` and `"logistic"` require a 0/1 response and fit
  latent-variable models, with fits and predictions on the latent scale.
  `"logistic"` uses Polya-Gamma augmentation.

  `"aft"` fits an accelerated failure time (log-normal) survival model:
  the response is a `Surv` object (from the survival package) or a
  two-column `(time, status)` matrix or data frame (status 1 an event, 0
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

  `"ordinal"` fits an ordered categorical response by a cumulative
  probit: a latent \\z = f(x) + \epsilon\\, \\\epsilon \sim N(0, 1)\\,
  is cut at ordered thresholds \\\gamma_1 = 0 \< \gamma_2 \< \ldots \<
  \gamma\_{K-1}\\ into the \\K\\ ordered categories, the free cutpoints
  sampled by a marginal Metropolis update with the latents integrated
  out. The response should be an ordered factor
  ([`is.ordered`](https://rdrr.io/r/base/factor.html)), whose level
  order defines the category order; under `family = "auto"` an ordered
  factor with three or more levels is detected and fit as ordinal,
  reporting the choice in a one-line message (a two-level ordered factor
  is binary and resolves to probit). `family = "ordinal"` given
  explicitly also accepts an unordered factor or character response -
  with a message that the category order is taken from the level order -
  and a numeric response, whose ordered levels are `sort(unique(y))`.
  Like probit, the latent scale is fixed at 1, fits are on the latent
  scale, and weights are not supported. `bart2` reports ordinal fits as
  \\n \times K\\ category probabilities; see
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md).

  `"nbinom"` fits a non-negative integer (count) response by a
  negative-binomial model with the Polya-Gamma augmentation: the forest
  fits a log-odds latent \\\psi = f(x) + o\\ and \\y \sim \mathrm{NB}(r,
  \mathrm{plogis}(\psi))\\, with mean \\E\[y \mid x\] = r e^{\psi}\\, so
  the offset enters multiplicatively as a log-exposure. The response
  must be a non-negative integer no larger than \\10^6\\ (the dispersion
  grid's count histogram is sized from the largest count), and
  `"nbinom"` is never inferred - a count carries no unambiguous class,
  so it must be requested explicitly. The dispersion `r` is estimated by
  default (see `dispersion`); like probit, the latent scale is fixed at
  1, fits are on the latent (log-odds) scale, and weights are not
  supported (exposure belongs in the offset). `bart2` reports mean
  counts; see
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md).

  `"hazard"` and `"hazard.logistic"` fit a discrete-time survival hazard
  model by person-period expansion (`"hazard.probit"` is an accepted
  alias for `"hazard"`). The response is the same `Surv` object or
  two-column `(time, status)` the `"aft"` family takes, but each subject
  observed to time \\t_i\\ is expanded into its at-risk rows (one per
  period \\k = 1, \ldots, t_i\\), each carrying the subject's
  covariates, an ordinal period column appended last, and the binary
  indicator \\y\_{ik} = \mathrm{status}\_i \cdot 1\\k = t_i\\\\; the
  discrete hazard is then \\h(k \mid x) = g(f(x, k) + o)\\ for the
  chosen binary link \\g\\ (probit for `"hazard"`, logistic for
  `"hazard.logistic"`), and the fit IS an ordinary
  `"probit"`/`"logistic"` fit on the expanded rows - it adds no engine
  code. The time grid is the sorted distinct observed times by default;
  `breaks` coarsens it. The offset is on the link scale and replicates
  per subject; weights replicate and follow the chosen binary family's
  policy. A `Surv` response requires the family to be requested
  explicitly (`"auto"` selects `"aft"`). Survival curves come from
  [`survivalProbabilities`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md),
  which requires `keepTrees`; `bart2` reports the fit with its
  `$periods` grid, and its `$family` records the binary link. Like
  `"aft"`, hazard fits use the matrix interface and do not support
  `subset` or a `test` set (expand test subjects with
  `survivalProbabilities(..., newdata = )`).

  `"hurdle.lognormal"` (alias `"twopart"`, which resolves and prints as
  `"hurdle.lognormal"`) fits a semicontinuous two-part (hurdle) model
  for a non-negative response with exact zeros: an occupancy probit fit
  of \\z = 1\\y \> 0\\\\ over all n observations, glued at report time
  to a lognormal positive-part fit - an ordinary gaussian fit of \\\log
  y\\ - over the subset \\\\i : y_i \> 0\\\\; the two parts share no
  parameters and are composed from two ordinary fits at independently
  derived seeds, so no engine code is added and a shared
  variable-selection prior across the parts is not available (a recorded
  limitation). `y.train` must be non-negative and finite and must carry
  at least one exact zero and one positive value. By default,
  predictions report the natural (response) scale via
  posterior-predictive Monte Carlo, \\E\[y \mid x\] = P(y \> 0 \mid
  x)\\e^{f(x) + \sigma^2 / 2}\\, heteroscedasticity-aware when the
  positive part carries a `variance = ~x` surface; see
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md) for the
  full `type` options (`"prob"`, `"link"`/`"log"`, the bimodal `"ppd"`)
  and the refusal list (`weights`, `subset`, `offset`, and `test` are
  all unsupported this arc). `dbarts()` does not fit this family: it
  composes two samplers, which only
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md) builds,
  so requesting it here is an error directing to
  [`bart2()`](https://vdorie.github.io/dbarts/reference/bart.md).

- dispersion:

  The negative-binomial dispersion \\r\\ (family `"nbinom"` only;
  ignored otherwise). `NA` (the default) estimates \\r\\ on a capped
  positive-integer grid under a renormalized \\\mathrm{gamma}(2, 0.1)\\
  prior. A supplied value fixes \\r\\ and must be a positive integer: v1
  ships the exact integer envelope, so a real fixed dispersion is
  refused. Larger \\r\\ approaches the Poisson limit (variance \\\to\\
  mean); smaller \\r\\ is more overdispersed.

- missing:

  How missing values in the predictors enter the model. The default
  `"incorporate"` keeps them: every split rule learns a direction for
  missing values on its variable, so an observation whose split value is
  `NA` follows that rule's chosen branch (“Missingness Incorporated in
  Attributes”, Twala et al. 2008), in training and test data alike.
  `"error"` rejects predictors containing `NA`. The response, `weights`,
  and `offset` must always be complete; note that previous versions
  silently dropped incomplete rows for formula inputs.

- breaks:

  The discrete-time grid for a `"hazard"` fit (ignored otherwise).
  `NULL` (the default) uses the sorted distinct observed times, one
  period per distinct time (the BART `surv.bart` convention). A single
  positive integer bins the times at the \\(1{:}K)/K\\ quantiles, giving
  \\K\\ periods. A numeric vector of length two or more gives explicit
  strictly-increasing interval boundaries \\b_0 \< \ldots \< b_K\\,
  defining right-closed periods \\(b\_{k-1}, b_k\]\\; every time must
  lie in \\(b_1, b_K\]\\.

- max.rows:

  A guard on the person-period expansion for a `"hazard"` fit (ignored
  otherwise). If the expanded design would exceed `max.rows` rows
  (\\\sum_i t_i\\), the fit is refused with a message naming the
  coarsening levers. The default \\10^7\\ catches an over-fine grid on
  heavily continuous times; coarsen with `breaks` or raise `max.rows`.

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

### Naming the leaf calibration

Because the leaf prior is calibrated off the range of the response the
sampler was CONSTRUCTED on, an R program that drives the sampler between
sweeps – feeding it latents, residuals, or offsets from an outer model –
inherits whatever calibration its construction vector happened to imply.
`node.prior = normal(scale = )` names that calibration instead, in
response units: it is the prior standard deviation of the forest total
\\f\\ at `k = 1`, so the prior standard deviation in force is
`scale / k` and the prior mean is the response transform's shift
(`(max(y) + min(y)) / 2`, net of `offset`, for a continuous response; 0
for the latent-scale families). `normal(sd = )` names the same quantity
at the resolved `k`, and is refused under a `k` hyperprior.
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md), and
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) take it
as a `prior.scale` argument directly. Unset, nothing changes.

The named quantity is the LEAF-PARAMETER scale of the forest total. It
equals the prior standard deviation of \\f(x)\\ at every \\x\\ for the
constant leaf only; under the other leaf models the prior of \\f(x)\\ is
x-dependent and `scale / k` bounds it in a leaf-model-specific
direction. Under a `linear` node prior it is a LOWER bound, attained at
the standardized covariate origin, with \\sd(f(x))\\ equal to
`scale / k` times \\\sqrt{1 + \\z(x)\\^2}\\ for \\z\\ the internally
standardized leaf covariates (a missing value maps to \\z_j = 0\\, and a
constant column contributes 0); the prior mean is exact. Under a `gp`
node prior it is an UPPER bound over \\x\\, attained at rows reproducing
a leaf member and on leaves past `max.leaf.size` (which draw as constant
leaves); elsewhere the prior variance is \\(scale / k)^2 c(x)' C^{-1}
c(x)\\ and decays to 0 as \\x\\ leaves the leaf's data cloud, at which
point every prior draw equals the prior mean exactly. Under a `monotone`
constraint it is a LOWER bound in the interior – the realized standard
deviation runs a few per cent to about 20% above it – and the prior mean
is NOT the prior mean of \\f(x)\\: the constrained marginal is skew,
with an x-dependent mean that tracks the constraint direction and spans
several prior standard deviations along the constrained axis (see
`monotone` above).

A named calibration records an INTENT. It is applied at creation and
re-derived at every `setModel` against the response transform then in
force, and the engine never writes it back, so a channel that re-anchors
the transform (`setResponse` or `setOffset` with `updateScale = TRUE`,
or `setData`) moves the calibration actually in force while leaving the
recorded intent alone; re-issue the model to restate it. Two-forest
(`forests`) and multinomial models take both forests' leaf scales from
their own calibration maps and refuse a named calibration rather than
drop it.

## Value

A reference object of
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

## See also

[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md),
[`extract.dbartsSampler`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)

## References

Chipman, H.A., George, E.I., McCulloch, R.E., and Shively, T.S. (2022)
mBART: multidimensional monotone BART. *Bayesian Analysis*, **17**(2),
515–544.

Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch, R.E. (2020)
Heteroscedastic BART via multiplicative regression trees. *Journal of
Computational and Graphical Statistics*, **29**(2), 405–417.

Twala, B.E.T.H., Jones, M.C., and Hand, D.J. (2008) Good methods for
coping with missing data in decision trees. *Pattern Recognition
Letters*, **29**(7), 950–956.

## Examples

``` r
set.seed(8)
n <- 100L
x <- matrix(runif(n * 3L), n, 3L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.2)

control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.burn = 0L,
                         n.samples = 1L, n.trees = 25L, updateState = FALSE)

## The sampler is created once and then run repeatedly, so that BART can act
## as a conditional model inside a larger Gibbs/Metropolis-Hastings sampler.
sampler <- dbarts(y ~ x, control = control)
samples <- sampler$run(numBurnIn = 50L, numSamples = 1L)
str(samples$train)
#>  num [1:100, 1] 0.107 0.148 1.545 0.975 -0.282 ...

## an outer step revises the response; the sampler picks it up on the next
## run() without being rebuilt. The prior scale set at creation stays locked
## unless updateScale = TRUE is passed (burn-in only).
sampler$setResponse(y + rnorm(n, 0, 0.05))
samples <- sampler$run()

## priors are given as expressions in dbarts's own vocabulary
sampler <- dbarts(y ~ x, control = control,
                  tree.prior = cgm(power = 1.5),
                  node.prior = normal(k = chi(1.5, 2)),
                  resid.prior = chisq(df = 5))

## an additive fit that is monotone increasing in the first predictor
sampler <- dbarts(y ~ x, control = control,
                  interactions = interactions(max.order = 1L),
                  monotone = c(1, 0, 0))
samples <- sampler$run(numBurnIn = 20L, numSamples = 5L)

## a non-gaussian family: counts by negative binomial, offset a log-exposure
counts <- rnbinom(n, size = 4, mu = exp(0.5 * x[, 1L]) * 4)
sampler <- dbarts(x, counts, family = "nbinom", control = control)
samples <- sampler$run(numBurnIn = 20L, numSamples = 5L)
```
