# Monotone Constraints for BART

Build a monotone-constraint specification making a BART fit monotone
increasing or decreasing in chosen predictors (monotone BART; Chipman,
George, McCulloch, and Shively 2022), and name the prior it is read
under. Pass the result as the `monotone` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md);
a plain direction vector passed there is shorthand for
`monotone(directions)` at the default prior.

The constructor is not exported: it resolves by bare name inside the
arguments that take it, and elsewhere is written
`dbartsForests$monotone(...)`; see
[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md).

## Usage

``` r
monotone(directions, prior = c("leaf", "joint"))
```

## Arguments

- directions:

  The per-predictor directions. A named vector selects predictors by
  model-matrix column name, each element `"increasing"` or `1`,
  `"decreasing"` or `-1`, or `0` for unconstrained. An unnamed vector as
  long as the model matrix is wide assigns the directions by position.
  Matching is case-sensitive and takes no abbreviations. A predictor
  named `prior` is written like any other,
  `monotone(c(prior = "increasing"))`.

- prior:

  The prior the constraint is read under. `"leaf"` restricts the
  leaf-value prior to the leaf values that are monotone given the tree,
  leaving the tree prior unchanged. `"joint"` conditions the tree
  structure and its leaf values on monotonicity together; it is not
  available yet, and a fit naming it is an error.

## Details

Only numeric and ordered columns are eligible: a direction on a
categorical (unordered factor) predictor is an error. Names, positions
and every direction are validated against the model matrix at fit time;
`prior` is matched when the specification is built. An all-zero
specification fits the unconstrained model.

A constraint forces birth/death-only tree proposals (a `control` naming
a non-default `proposal.probs` is then an error) and a fixed `k = 2` (an
explicit `k` hyperprior is an error); linear and Gaussian-process
leaves, variance forests and multi-forest models are not supported under
it. It may be combined with
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)
and [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md).

## Value

A `dbartsMonotone` specification object, a list of `directions` and
`prior`, resolved when a sampler is built.

## References

Chipman, H. A., George, E. I., McCulloch, R. E., and Shively, T. S.
(2022) mBART: multidimensional monotone BART. *Bayesian Analysis*
**17**(2), 515–544.

## See also

[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md),
[`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)

## Examples

``` r
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] - x[, 2] + sin(2 * pi * x[, 3]) + rnorm(n, 0, 0.2)
df <- data.frame(y, x)

## increasing in x1, decreasing in x2, free in x3
fit <- bart(y ~ x1 + x2 + x3, df,
             monotone = monotone(c(x1 = "increasing", x2 = "decreasing"),
                                 prior = "leaf"),
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## the plain vector, positionally, at the default prior
fit.plain <- bart(y ~ x1 + x2 + x3, df,
                   monotone = c(1, -1, 0),
                   n.trees = 25L, n.samples = 20L, n.burn = 20L,
                   n.chains = 1L, verbose = FALSE)
```
