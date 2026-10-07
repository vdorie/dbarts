# Forest Specification Constructors

A list of the constructor functions building the forest specifications
that the `monotone`, `interactions`, `blocks`, `variance` and `forests`
arguments of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
and [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) accept.
The constructors are not exported, so that names like `forest` and
`blocks` stay out of the search path, where another package could mask
them or be masked depending on attach order.

Inside those arguments, and inside a formula's
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) term,
the constructors resolve by bare name: a call is always the constructor,
and a bare name the caller has assigned a value to is that value.
Anywhere else - a specification built ahead of the call, or an argument
a wrapper forces before passing it on - write
`dbartsForests$interactions(...)`.
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md)'s
first argument, the predictors a forest splits on, is read as terms when
it names a predictor of the fit and otherwise takes the value it has
where and when the constructor is called. Its `basis` is read the same
way: a name of a column of the fit's data is that column, and anything
else is what it was where and when the constructor is called; see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md).

## Format

A list of functions:

- `interactions`:

  Interaction constraints; see
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md).

- `blocks`:

  Block-additive constraints; see
  [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md).

- `monotone`:

  Monotone constraints and their prior; see
  [`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md).

- `forest`:

  One forest of a model, in a `forests` list or as a term of a formula;
  see [`forest`](https://vdorie.github.io/dbarts/reference/forest.md).

- `varianceForest`:

  The heteroscedastic variance forest; see
  [`varianceForest`](https://vdorie.github.io/dbarts/reference/varianceForest.md).

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md),
[`dbartsFamilies`](https://vdorie.github.io/dbarts/reference/dbartsFamilies.md)

## Examples

``` r
## built ahead of the call, then passed as a value
constraint <- dbartsForests$interactions(max.order = 1)

set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] + ifelse(x[, 2] > 0.5, 1, -1) + rnorm(n, 0, 0.2)

fit <- bart(x, y,
             interactions = constraint,
             n.trees = 25L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## the same constraint, written by bare name inside the argument
fit.bare <- bart(x, y,
                  interactions = interactions(max.order = 1),
                  n.trees = 25L, n.samples = 20L, n.burn = 20L,
                  n.chains = 1L, verbose = FALSE)
```
