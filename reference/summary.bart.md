# Convergence Diagnostics for BART Fits

Reports a per-variable posterior summary of the scalar parameters
(`sigma` and `k`) of a `bartBT`/`bart` fit, along with split-\\\hat{R}\\
and bulk/tail effective sample size, computed by dbarts itself (no
posterior dependency). A heteroscedastic fit
([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
`variance`) has no scalar residual scale and reports `mean.s` in place
of `sigma`; see `vars`.

`bart`'s four own-class fits (`"bartMultinomial"`, `"bartOrdinal"`,
`"bartNegbin"`, `"bartHurdle"`; see
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) summarize
through this same method, each exposing the scalar posterior parameters
that family carries rather than its per-observation channels - never
`yhat.train` itself.

[`draws`](https://vdorie.github.io/dbarts/reference/draws.md) returns
the same fit's chain-dimensioned draws as a plain array, the shape this
summary is computed from.

## Usage

``` r
# S3 method for class 'bart'
summary(object, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'summary.bart'
print(x, ...)
```

## Arguments

- object, x:

  An object of class `bart`, as returned by
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) or
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md); for the
  four own-class fits, a `bart` fit of the matching class.

- vars:

  Character vector of fields to summarize. Requested fields absent from
  `object` (e.g. `k` when unmodeled, or `tau`, which no shipped family
  carries) are silently dropped. `sigma`, `k`, and `tau` contribute one
  summarized variable each; any other field (`varcount`, `varprobs`,
  `yhat.train`, `yhat.test`) contributes one variable per column, named
  `"field[column]"`.

  On a heteroscedastic fit the `sigma` token resolves to `mean.s`: the
  mean over the training observations of that draw's variance surface
  \\s(x)\\, one value per draw; see
  [`draws`](https://vdorie.github.io/dbarts/reference/draws.md) for the
  full rule.

  A `"bartMultinomial"` fit has no `vars` argument: its only scalar
  posterior parameter is the per-category mean predicted probability,
  reported as `meanProb[<level>]`.

- ...:

  Unused.

## Value

An object of class `summary.bart` with elements `call`, `stats` (a data
frame with one row per requested scalar, or `NULL` if none of `vars` are
present on the fit), and `vars`. `stats` has columns `mean`, `median`,
`sd`, `mad`, `q5`, `q95`, `rhat`, `ess_bulk`, and `ess_tail` - the same
columns
[`posterior::summarise_draws`](https://mc-stan.org/posterior/reference/draws_summary.html)
reports for a plain array, computed without that package. Its print
method notes when any R-hat exceeds 1.01; this does not withhold the
rest of the summary or error, as dbarts does not refuse to summarize a
non-converged fit.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`draws`](https://vdorie.github.io/dbarts/reference/draws.md)

## Examples

``` r
# \donttest{
fit <- bart(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)
#> 
#> Call:
#> bart(formula = y ~ x, data = data.frame(y = rnorm(100), x = rnorm(100)), 
#>     n.trees = 5L, n.samples = 20L, n.burn = 20L, n.chains = 2L, 
#>     n.threads = 1L, verbose = FALSE, factors = "categorical", 
#>     proposal.probs = c(birth_death = 0.6, swap = 0, change = 0.4, 
#>     perturb = 0, rule_gibbs = 0, birth = 0.5))
#> 
#>   variable      mean    median        sd        mad       q5      q95     rhat
#> 1    sigma 0.9767962 0.9624961 0.0614717 0.05889664 0.894719 1.084578 1.104571
#>   ess_bulk ess_tail
#> 1 49.57946 49.57265
#> 
#> Note: some R-hat values exceed 1.01; chains may not have converged.
# }
```
