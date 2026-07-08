# Convergence Diagnostics and Posterior-Package Draws for BART Fits

`summary` reports a per-variable posterior summary of the scalar
parameters (`sigma`, `k`, and, for
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) fits,
`tau`) of a `bart`/`bart2`/`rbart_vi` fit, along with split-\\\hat{R}\\
and effective sample size when the posterior package is installed.

`as_draws_array` and `as_draws_df` convert a fit's chain-dimensioned
draws to posterior's array/data frame conventions. Both are methods for
posterior's generics and so require it to be loaded
([`library(posterior)`](https://mc-stan.org/posterior/) or a qualified
[`posterior::as_draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
call); dbarts registers them only when posterior is available, as it is
a `Suggests`, not a hard dependency.

## Usage

``` r
# S3 method for class 'bart'
summary(object, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'rbart'
summary(object, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'summary.bart'
print(x, ...)

# S3 method for class 'bart'
as_draws_array(x, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'rbart'
as_draws_array(x, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'bart'
as_draws_df(x, vars = c("sigma", "k", "tau"), ...)
# S3 method for class 'rbart'
as_draws_df(x, vars = c("sigma", "k", "tau"), ...)
```

## Arguments

- object, x:

  An object of class `bart` or `rbart`, as returned by
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md), or
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md).

- vars:

  Character vector of fields to gather. Requested fields absent from
  `object` (e.g. `k` when unmodeled, `tau` for a plain `bart`/`bart2`
  fit) are silently dropped. `sigma`, `k`, and `tau` contribute one
  draws variable each; any other field (`varcount`, `varprobs`,
  `yhat.train`, `yhat.test`, `ranef`) contributes one variable per
  column, named `"field[column]"`. \\f(x)\\ draws (`yhat.*`) are
  reachable this way but are not summarized automatically, as they carry
  one variable per observation.

- ...:

  Unused.

## Value

`summary` returns an object of class `summary.bart` with elements
`call`, `stats` (a data frame with one row per requested scalar, or
`NULL` if none of `vars` are present on the fit; `rhat`, `ess_bulk`, and
`ess_tail` columns are added when posterior is installed), and
`posterior` (whether those columns were computed). Its print method
notes when posterior is missing, and separately when any R-hat exceeds
1.01; neither withholds the rest of the summary or errors.

`as_draws_array` and `as_draws_df` return a posterior
`draws_array`/`draws_df`: the fit's native (chain, sample\[, variable\])
storage, transposed to posterior's (iteration, chain, variable)
convention.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)

## Examples

``` r
# \donttest{
fit <- bart2(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)
#> 
#> Call:
#> bart2(formula = y ~ x, data = data.frame(y = rnorm(100), x = rnorm(100)), 
#>     n.trees = 5L, n.samples = 20L, n.burn = 20L, n.chains = 2L, 
#>     n.threads = 1L, verbose = FALSE)
#> 
#> # A tibble: 1 × 10
#>   variable  mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
#>   <chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 sigma     1.05   1.04 0.0725 0.0815 0.942  1.16 0.983     39.4     56.9

if (requireNamespace("posterior", quietly = TRUE)) {
  library(posterior)
  as_draws_array(fit)
}
#> This is posterior version 1.7.0
#> 
#> Attaching package: ‘posterior’
#> The following objects are masked from ‘package:stats’:
#> 
#>     mad, sd, var
#> The following objects are masked from ‘package:base’:
#> 
#>     %in%, match
#> # A draws_array: 20 iterations, 2 chains, and 1 variables
#> , , variable = sigma
#> 
#>          chain
#> iteration    1    2
#>         1 1.01 1.01
#>         2 0.92 1.16
#>         3 1.11 0.99
#>         4 1.16 1.02
#>         5 0.99 1.06
#> 
#> # ... with 15 more iterations
# }
```
