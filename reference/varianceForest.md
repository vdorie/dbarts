# Variance Forest Specification for Heteroscedastic BART

Build the specification of the heteroscedastic variance forest (HBART;
Pratola, Chipman, George, and McCulloch 2020). Pass the result as the
`variance` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
or [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) - the
same argument that already accepts the plain selector (`NULL`/`FALSE`
for none, `TRUE`/a one-sided formula/a character or integer column
selector for the predictors driving the variance). This constructor is
the other accepted type of that one argument, adding the variance
forest's own tree count and structure prior alongside the column
selection.

The constructor is not exported: it resolves by bare name inside the
arguments that take it, and elsewhere is written
`dbartsForests$varianceForest(...)`; see
[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md).

## Usage

``` r
varianceForest(
    vars = NULL, n.trees = NULL, base = NULL, power = NULL,
    proposal.probs = NULL, n.perturb.cuts = NULL)
```

## Arguments

- vars:

  Which predictors drive the variance surface, written as
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)'s
  `vars` is: bare terms (`x1 + x2`), a one-sided formula in place or
  held in a variable (`~ x1 + x2`, the tilde accepted and never
  required), or a character or integer column selector (a bare factor
  term name expands to its indicator columns). `NULL` (the default)
  reads *every* predictor - the object's own “unrestricted” reading,
  distinct from `variance = NULL`, which declares no variance forest at
  all.

- n.trees:

  The number of variance trees. `NULL` (the default) takes the variance
  forest's own default of `40`, the same value `variance = TRUE` and the
  other plain-selector forms use.

- base, power:

  The variance forest's tree-structure prior. `NULL` (the default)
  inherits the mean forest's `base` and `power`, those of the fitting
  function or the lone mean forest's own.

- proposal.probs, n.perturb.cuts:

  The variance forest's proposal mixture and perturb window, as
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  describes each. Each one stated here governs the variance forest and
  each one left out takes the control's: `$setControl` changes it only
  where none is stated here, and one stated here is fixed for the
  sampler's life. A stated `proposal.probs` is a whole mixture, filled
  against the package defaults as
  [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)'s
  is, and a positive `"rule_gibbs"` share is refused: that move does
  nothing on the variance forest's leaves. `n =` and `p =` each
  abbreviate more than one argument, so write `n.trees` and `power` in
  full.

## Details

A second ensemble models the residual variance surface \\s^2(x)\\ as a
product of scaled-inverse-chi-squared (multiplicative) leaves, so \\y_i
= f(x_i) + s(x_i)\epsilon_i\\, coupled to the mean forest through the
per-observation precision. Gaussian and `"aft"` (survival) responses and
constant leaves only - the latent families route their own precisions
through the same channel; monotone constraints are not supported. Under
`"aft"` the surface is the dispersion of log survival time, and a
right-censored observation's latent log-time is redrawn at its own
\\s(x_i)\\. Every value is validated when a sampler is built.

## Value

A `dbartsVarianceForest` specification object, resolved when a sampler
is built.

## References

Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch, R.E. (2020)
Heteroscedastic BART via multiplicative regression trees. *Journal of
Computational and Graphical Statistics*, **29**(2), 405–417.

## See also

[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)

## Examples

``` r
set.seed(0)
n <- 200L
x <- matrix(runif(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
s <- ifelse(x[, 1] > 0.5, 0.3, 1.2)
y <- 2 * x[, 1] + s * rnorm(n)

## every predictor drives the variance, 20 variance trees
fit <- bart(x, y,
             variance = varianceForest(n.trees = 20L),
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## restricted to x1, with its own tree-structure prior
fit.restricted <- bart(x, y,
                        variance = varianceForest(vars = ~x1, n.trees = 20L,
                                                  base = 0.9, power = 1.5),
                        n.trees = 25L, n.samples = 20L, n.burn = 20L,
                        n.chains = 1L, verbose = FALSE)

dbartsForests$varianceForest(vars = ~x1, n.trees = 20L, base = 0.9, power = 1.5)
#> dbarts variance forest specification
#>   vars           = x1
#>   n.trees        = 20L
#>   base           = 0.9
#>   power          = 1.5
#>   proposal.probs = <control's>
#>   n.perturb.cuts = <control's>
```
