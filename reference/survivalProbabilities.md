# Survival Probability Draws from a Survival Fit

Posterior draws of the survival probability \\S(t \mid x)\\ from a
survival fit: an accelerated failure time (AFT) log-normal fit produced
by [`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) (or
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) with
`family = "aft"`, the grouped (random-intercept) AFT fit produced by
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) with
`family = "aft"`, or a discrete-time hazard fit produced by
[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) with
`family = "hazard"` (or `"hazard.logistic"`).

## Usage

``` r
survivalProbabilities(object, ...)

# S3 method for class 'bart'
survivalProbabilities(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...
)

# S3 method for class 'rbart'
survivalProbabilities(
  object,
  times,
  newdata = NULL,
  group.by,
  combineChains = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted `bart` or `rbart` object from a survival model: an AFT model
  (its `family` element equals `"aft"`) or a discrete-time hazard model
  (it carries a `$periods` grid; its `family` records the binary link).
  A hazard fit must have been made with `keepTrees = TRUE`.

- times:

  A numeric vector of times (on the original, non-logarithmic scale) at
  which to evaluate the survival function. For a discrete-time hazard
  fit the times are horizons on the period grid, and the default (when
  `times` is missing) is the training grid, `object$periods`.

- newdata:

  Optional predictors at which to evaluate. When `NULL` the training
  observations are used; otherwise `object` must have been fit with
  `keepTrees = TRUE`.

- group.by:

  For the `rbart` method with `newdata`, the grouping factor for the new
  observations, as in
  [`predict`](https://vdorie.github.io/dbarts/reference/rbart.md). A
  group not seen in training draws its intercept from \\N(0, \tau)\\.
  Ignored (and unnecessary) when `newdata` is `NULL`.

- combineChains:

  A logical determining whether the chain dimension is collapsed into
  the draw dimension, as elsewhere in the package.

- ...:

  Not used; for compatibility with the generic.

## Details

Under the log-normal AFT model \\\log T = f(x) + \sigma \epsilon\\ with
\\\epsilon \sim N(0, 1)\\, the survival function is \$\$S(t \mid x) =
1 - \Phi\left(\frac{\log t - f(x)}{\sigma}\right),\$\$ where \\f(x)\\ is
the linear predictor \\E\[\log T \mid x\]\\ - the log-time-scale
quantity that
[`predict`](https://vdorie.github.io/dbarts/reference/bart.md) and
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md) return
for an `"aft"` fit - and \\\sigma\\ the residual standard deviation. The
probability is evaluated at every posterior draw of \\f(x)\\ and
\\\sigma\\, following the package's convention that draw-level functions
return draws: take means or quantiles over the draw margin for point
estimates and credible bands.

For a grouped fit from
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) with
`family = "aft"`, the `rbart` method uses the same formula, with
\\f(x)\\ replaced by the linear predictor \\E\[\log T \mid x, g\] =
f(x) + \alpha\_{g}\\ that includes the drawn random intercept for the
observation's group. It is sourced from the expected-value (`"ev"`)
channel, so dropping the intercepts (which would misplace every grouped
curve) is not possible.

For a discrete-time hazard fit, the survival function is the cumulative
product \$\$S(t \mid x) = \prod\_{k \\:\\ \mathrm{periods}\[k\] \le t}
(1 - h(k \mid x)),\$\$ where \\h(k \mid x) = g(f(x, k) + o)\\ is the
per-period hazard through the fit's binary link \\g\\. Because the
training design is ragged (each subject carries only its at-risk rows),
the method ALWAYS re-expands its subjects onto the full grid and replays
the trees - so it requires `keepTrees = TRUE` even for the training data
(`newdata = NULL`). With `newdata`, each new subject is expanded to one
row per period and its curve evaluated the same way.

## Value

An array of survival probabilities. With `combineChains = TRUE` (or a
single chain), the dimensions are draws \\\times\\ times \\\times\\
observations; with `combineChains = FALSE` on a multi-chain fit, a chain
dimension precedes them (chains \\\times\\ draws \\\times\\ times
\\\times\\ observations). Observations are those of `newdata`, or the
training data when `newdata` is `NULL`.

## See also

[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md),
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md),
[`predict`](https://vdorie.github.io/dbarts/reference/bart.md),
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)

## Author

Vincent Dorie: <vdorie@gmail.com>.

## Examples

``` r
set.seed(42)
n <- 100L
x <- matrix(runif(n * 2L), n, 2L)
eventTime <- exp(x[, 1L] + 0.5 * rnorm(n))
censorTime <- exp(x[, 1L] + 0.25 + 0.5 * rnorm(n))
status <- as.numeric(eventTime <= censorTime)
observedTime <- pmin(eventTime, censorTime)

fit <- bart2(
  x,
  cbind(observedTime, status),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 100L,
  n.chains = 1L,
  verbose = FALSE
)

surv <- survivalProbabilities(fit, times = c(0.5, 1, 2))
dim(surv) # draws x times x observations
#> [1] 100   3 100

## posterior-mean survival curve and 90% band for the first observation
apply(surv[, , 1L], 2L, mean)
#> [1] 0.9959171 0.9259428 0.5725221
apply(surv[, , 1L], 2L, quantile, probs = c(0.05, 0.95))
#>          [,1]      [,2]      [,3]
#> 5%  0.9874998 0.8314811 0.3862876
#> 95% 0.9998579 0.9851896 0.7723411
```
