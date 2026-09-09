# Chain-Dimensioned Draws From a BART Fit

Converts a fit's chain-dimensioned draws to a plain array with
dimensions `(iteration, chain, variable)` and dimnames on the variable
margin - the shape
[`posterior::as_draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
accepts unchanged from a caller who has that package installed
(`as_draws_array(draws(fit))`).
[`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
is built on the same array.

[`bart2`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md)'s
four own-class fits (`"bartMultinomial"`, `"bartOrdinal"`,
`"bartNegbin"`, `"bartHurdle"`; see
[`bart2`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md))
have their own `draws` methods, each exposing the scalar posterior
parameters that family carries rather than its per-observation
channels - never `yhat.train` itself.

## Usage

``` r
# S3 method for class 'bart'
draws(x, vars = c("sigma", "k", "tau"), ...)

# S3 method for class 'bartMultinomial'
draws(x, vars = "meanProb", ...)
# S3 method for class 'bartOrdinal'
draws(x, vars = c("thresholds", "sigma", "k", "tau"), ...)
# S3 method for class 'bartNegbin'
draws(x, vars = c("dispersion", "sigma", "k", "tau"), ...)
# S3 method for class 'bartHurdle'
draws(x, vars = c("sigma", "k", "tau"), ...)
```

## Arguments

- x:

  An object of class `bart`, as returned by
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
  [`bart2`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md);
  for the four own-class methods, a `bart2` fit of the matching class.

- vars:

  Character vector of fields to gather. Requested fields absent from `x`
  (e.g. `k` when unmodeled, or `tau`, which no shipped family carries)
  are silently dropped. `sigma`, `k`, and `tau` contribute one draws
  variable each; any other field (`varcount`, `varprobs`, `yhat.train`,
  `yhat.test`) contributes one variable per column, named
  `"field[column]"`. \\f(x)\\ draws (`yhat.*`) are reachable this way
  but are not summarized automatically by
  [`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md),
  as they carry one variable per observation.

  On a heteroscedastic fit the `sigma` token resolves to `mean.s`: the
  mean over the training observations of that draw's variance surface
  \\s(x)\\, one value per draw. The `sigma` such a fit stores is the
  variance forest parameterization's fixed unit residual times the range
  of the response, a constant with no posterior content, so it is not
  reported as a parameter. `s.train` itself is reachable by name,
  contributing one variable per observation.

  For the four own-class methods, `vars` is scoped to that family's own
  vocabulary (see
  [`bart2`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md)):
  a `"bartOrdinal"` fit's `"thresholds"` contributes `threshold[1]`
  (pinned at 0) through `threshold[K - 1]`; a `"bartNegbin"` fit's
  `"dispersion"` contributes the per-draw dispersion \\r\\; a
  `"bartMultinomial"` fit has a single channel, so its `vars` only ever
  means `"meanProb"` and it always reports `meanProb[<level>]`, its only
  scalar posterior parameter; a `"bartHurdle"` fit applies `vars` to
  both components and labels the result
  `occupancy.<field>`/`positive.<field>` (a dot, not a bracket, since
  posterior parses a bracket as an index).

- ...:

  Unused.

## Value

A numeric array with dimensions `(iteration, chain, variable)` and
dimnames on the variable margin, matching the fit's native (chain,
sample\[, variable\]) storage transposed to that convention. The four
own-class methods return the same convention over their own
family-scoped `vars`.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md),
[`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md)

## Examples

``` r
# \donttest{
fit <- bart2(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
d <- draws(fit)
dim(d)
#> [1] 20  2  1
# }
```
