# Survival Probability Draws from an AFT Fit

Posterior draws of the survival probability \\S(t \mid x)\\ from an
accelerated failure time (AFT) log-normal fit produced by
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md) (or
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) with
`family = "aft"`.

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
survivalProbabilities(object, ...)
```

## Arguments

- object:

  A fitted `bart` object from an AFT model (its `family` element equals
  `"aft"`).

- times:

  A numeric vector of times (on the original, non-logarithmic scale) at
  which to evaluate the survival function.

- newdata:

  Optional predictors at which to evaluate. When `NULL` the training
  observations are used; otherwise `object` must have been fit with
  `keepTrees = TRUE`.

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

There is no method for `rbart` fits:
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) has no
`family` argument, so a grouped (random-intercept) AFT model is not yet
fittable from R, and a survival curve that ignored the fitted intercepts
would be wrong. The `rbart` method fails with an error saying so.

## Value

An array of survival probabilities. With `combineChains = TRUE` (or a
single chain), the dimensions are draws \\\times\\ times \\\times\\
observations; with `combineChains = FALSE` on a multi-chain fit, a chain
dimension precedes them (chains \\\times\\ draws \\\times\\ times
\\\times\\ observations). Observations are those of `newdata`, or the
training data when `newdata` is `NULL`.

## See also

[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md),
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
#> [1] 0.9978512 0.9363851 0.5622351
apply(surv[, , 1L], 2L, quantile, probs = c(0.05, 0.95))
#>          [,1]      [,2]      [,3]
#> 5%  0.9941061 0.8735136 0.4033542
#> 95% 0.9998750 0.9804595 0.6923473
```
