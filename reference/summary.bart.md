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
`summary.glm`'s dispersion line. The summary of a fit of several forests
adds the full block of what each forest resolved to, after the table
(see the `sd` argument of
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), and
the line under the table then leaves out each forest's `leaf.prior.sd`,
which the block gives. Under the block, each forest with a basis has a
line for the effect of one unit of its basis (for a factor, of its
second level against its first): its *average* over the training rows in
the likelihood, and its *sd across rows*, the standard deviation of the
effect across those rows, each with its posterior median and the 5% to
95% interval `stats` uses. Both are computed for each draw and the
median taken last, on the scale the forests add on (the response's
units, or the probit or logit scale of a binary response), from each
forest's own fits times the coefficient on its basis, which is what
`predict` gives with the basis one unit higher less the basis as it
stands; a fit that kept no per-forest fits (`keepFits = FALSE`) says so
in one line instead. On a binary (probit or logistic), binomial or
discrete-time hazard response, a second line under each effect gives it
as a probability: for a factor the change in probability from its first
level to its second, for a numeric basis the slope of the probability
per unit of the basis (as `marginaleffects::avg_slopes` reports it),
each an average and an sd across rows computed per draw and row from the
kept fits and averaged last, so the average is the posterior of the mean
change in probability. A 0/1 numeric basis has a slope that is not the 0
to 1 contrast, which `factor(z)` as the basis gives exactly. A gaussian
response has no such line. Where the weights are counts of identical
rows (the count families' default, and gaussian with
`weights.type = "frequency"`), the average and the sd across rows, as
the forest line's column mean and sd, count each row at its count, as
the data written out row by row would. A binomial response, written as
two columns, counts each row at its trials times any such weight, and
the header then reads “over 371 trials (18 training rows)”.

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
  it by `k`, `leaf.prior.sd` (the `k.scale` over `k`; see
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md), where
  on a fit with several forests it is each forest's `sd` in the
  response's units, printed in the block after the table) for one that
  named it by `sd`, and either can be requested by name, as can
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
requested parameters the fit held fixed and their values (on a fit with
several forests a vector with one entry for each forest, `k` and
`leaf.prior.sd` named by the forests' labels), and, on a fit of several
forests that kept its fits, `forests`, the lines of the full block that
the print method gives after the table, and `effects`, where a forest
has a basis, a data frame with a row for each effect (columns `forest`,
`label`, `effect`, `average`, `average.q5`, `average.q95`, `spread`,
`spread.q5` and `spread.q95`, and, for a binary, binomial or hazard
response, `probability.kind` (`"change"` or `"slope"`) and the same six
summaries prefixed `probability.`, with the number of rows averaged over
and the scale as attributes `"rows"` and `"scale"`) or the string
`"dropped"` for a fit that kept no per-forest fits. `stats` has columns
`variable`, `mean`, `median`, `sd`, `mad`, `q5`, `q95`, `rhat`,
`ess_bulk`, and `ess_tail` - the same columns
`posterior::summarise_draws` reports for a plain array, computed without
that package. Its print method shows each R-hat as a value, with no bar
and no note, and dbarts does not refuse to summarize a fit whose chains
disagree.

## Reading R-hat

R-hat compares each chain's spread with the spread of the chains pooled,
so a value above 1 says the chains disagree; pooling them does not cure
it. At the defaults BART's chains routinely disagree on `sigma`, on `k`
and on a single row's prediction, so values above 1 there are expected,
and no single bar says whether a fit is usable. Averaged estimands, such
as the mean prediction over a subgroup or a contrast between two
subgroups, mix far better, so read R-hat on the estimand you report.
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md) with
`combineChains = FALSE` gives each chain's draws; the example computes
R-hat of a contrast this way, with
[`rhat`](https://vdorie.github.io/dbarts/reference/rhat.md).

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
#>     n.threads = 1L, verbose = FALSE)
#> 
#>   variable     mean   median         sd        mad        q5      q95      rhat
#> 1    sigma 1.071629 1.056783 0.08315055 0.08252557 0.9504061 1.195286 0.9996105
#>   ess_bulk ess_tail
#> 1 44.60904 43.93939
#> (Fixed, not sampled: k = 2)
# }

# \donttest{
# R-hat of the estimand reported: the difference between the mean prediction
# for two subgroups, one value per draw and chain
set.seed(1)
n <- 200L
x <- matrix(runif(n * 5L), n, 5L)
y <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
  10 * x[, 4L] + 5 * x[, 5L] + rnorm(n)
fit <- bart(x, y, n.chains = 4L, n.samples = 500L, n.burn = 100L,
            n.trees = 50L, n.threads = 1L, verbose = FALSE)
# chains by samples by rows
yhat <- dbarts::extract(fit, "ev", combineChains = FALSE)
group <- x[, 1L] > 0.5
contrast <- apply(yhat[, , group, drop = FALSE], 1:2, mean) -
  apply(yhat[, , !group, drop = FALSE], 1:2, mean)
dbarts::rhat(t(contrast)) # rhat() takes iterations by chains
#> [1] 1.022329
# }
```
