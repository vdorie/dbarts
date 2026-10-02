# Convergence Diagnostics for BART Fits

Reports a per-variable posterior summary of the scalar parameters
(`sigma` and the leaf prior's `k` or `leaf.prior.sd`) of a
`bartBT`/`bart` fit, along with split-\\\hat{R}\\ and bulk/tail
effective sample size, computed by dbarts itself (no posterior
dependency). A heteroscedastic fit
([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
`variance`) has no scalar residual scale and reports `mean.s` in place
of `sigma`; see `vars`.

`bart`'s four own-class fits (`"bartMultinomial"`, `"bartOrdinal"`,
`"bartNegbin"`, `"bartHurdle"`; see
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) summarize
through this same method, each exposing the scalar posterior parameters
that family carries rather than its per-observation channels - never
`yhat.train` itself. A chain- separated (rather than pooled) view of
these same draws is
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md), whose
`"sigma"`, `"k"`, `"leaf.prior.sd"`, and other scalar types return a
vector or, with `combineChains = FALSE`, a chains-by-samples matrix. The
table holds the parameters the fit sampled; one line under it names each
requested parameter the fit held fixed, with its value, in the manner of
`summary.glm`'s dispersion line.

## Usage

``` r
# S3 method for class 'bart'
summary(object, vars = c("sigma", "k", "leaf.prior.sd", "resid.df"), ...)
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

  Character vector of fields to summarize. By default the leaf-scale
  field follows how the leaf prior was named: `k` for a fit that named
  it by `k`, `leaf.prior.sd` (the anchor over `k`; see
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)) for
  one that named it by `sd`, and either can be requested by name, as can
  `first.k`, the burn-in draws. A requested parameter the fit held fixed
  (`sigma` under `fixed()`, a fixed `k` or `leaf.prior.sd`, a fixed
  `shape` or `resid.df`) is not tabulated but named on the line under
  the table; with the default `vars` that includes a fixed Student-t
  `resid.df`, and a sampled one is tabulated; other requested fields
  absent from `object` are silently dropped. `sigma`, `k` and
  `leaf.prior.sd` contribute one summarized variable each; any other
  field (`varcount`, `varprobs`, `yhat.train`, `yhat.test`) contributes
  one variable per column, named `"field[column]"`.

  On a heteroscedastic fit the `sigma` token resolves to `mean.s`: the
  mean over the training observations of that draw's variance surface
  \\s(x)\\, one value per draw, since such a fit has no scalar `sigma`;
  `s.train` itself is reachable by name, contributing one variable per
  observation.

  For the four own-class fits, `vars` is scoped to that family's own
  vocabulary (see
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)): a
  `"bartOrdinal"` fit's `"thresholds"` contributes `threshold[2]`
  through `threshold[K - 1]`, the first being pinned at 0 and named
  under the table; a `"bartNegbin"` fit's `"shape"` contributes the
  per-draw shape \\r\\; a `"bartMultinomial"` fit has no `vars` argument
  at all - its only scalar posterior parameter is the per-category mean
  predicted probability, reported as `prob[<level>]`; a `"bartHurdle"`
  fit applies `vars` to both components separately.

- ...:

  Unused.

## Value

An object of class `summary.bart` with elements `call`, `stats` (a data
frame with one row per requested sampled scalar, or `NULL` if none of
`vars` are present on the fit), `vars`, and `fixed`, a named list of the
requested parameters the fit held fixed and their values (a named
vector, by forest, on a fit with several forests). `stats` has columns
`mean`, `median`, `sd`, `mad`, `q5`, `q95`, `rhat`, `ess_bulk`, and
`ess_tail` - the same columns `posterior::summarise_draws` reports for a
plain array, computed without that package. Its print method notes when
any R-hat exceeds 1.01; this does not withhold the rest of the summary
or error, as dbarts does not refuse to summarize a non-converged fit.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md),
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md)

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
#>     n.threads = 1L, verbose = FALSE, factors = "categorical")
#> 
#>   variable     mean   median         sd        mad      q5      q95     rhat
#> 1    sigma 1.062273 1.056108 0.05847884 0.05357023 0.96293 1.145456 1.046318
#>   ess_bulk ess_tail
#> 1 51.85035 36.29032
#> (Fixed, not sampled: k = 2)
#> 
#> Note: some R-hat values exceed 1.01; chains may not have converged.
# }
```
