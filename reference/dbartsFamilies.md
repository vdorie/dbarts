# Response Family Constructors

A list of constructor functions building the response families that the
`family` argument of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
and [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)
accepts. Bundling them keeps generic names like `gaussian` and `probit`
out of the search path, where another package -
[`stats::gaussian`](https://rdrr.io/r/stats/family.html) first among
them - could mask them or be masked depending on attach order.

`family` takes either a token string (`family = "gaussian"`), meaning
that family at its defaults, or one of these calls
(`family = student(df = 4)`), meaning that family with settings. Inside
the `family` argument the constructors resolve by bare name whatever the
caller has attached, exactly as the prior vocabulary does inside
`tree.prior`; elsewhere reach them through this list. That vocabulary
shadows the caller's own frame for these ten names only, so a variable
named `gaussian` holding `student(3)` resolves there to the gaussian
family, not to its own value; name such a variable something else, or
write `dbartsFamilies$student(3)`.

Every setting that only one family reads rides its family object rather
than a formal of the fitting function. The retired spellings
`resid.dist`, `dispersion`, `breaks` and `max.rows` are accepted for one
release with a once-per-session warning and are removed in dbarts 1.1-0.

## Usage

``` r
# S4 method for class 'dbartsFamily'
show(object)
```

## Arguments

- object:

  A family object, as one of the constructors below returns. Printing
  one names the call that would build it.

## Format

A list of functions:

- [`gaussian()`](https://rdrr.io/r/stats/family.html):

  A continuous response with normal errors; the default for a numeric
  response under `family = "auto"`.

- `student(df = NULL)`:

  A continuous response with outlier-robust Student-t errors, drawn by
  the Gaussian scale-mixture augmentation. `df = NULL` estimates the
  degrees of freedom on a capped grid; a positive number fixes them. The
  former `resid.dist = student(df)` spelling.

- `probit()`, `logistic()`:

  Binary latent-variable models for a 0/1 response.

- `multinomial()`, `ordinal()`:

  Unordered K-category softmax and ordered cumulative-probit models.

- `nbinom(dispersion = NA)`:

  Negative-binomial counts. `NA` estimates the dispersion \\r\\ on a
  capped positive-integer grid; a positive integer fixes it. The former
  `dispersion` argument.

- `aft()`:

  Accelerated failure time (log-normal) survival.

- `hazard(breaks = NULL, max.rows = 1e7, link = c("probit", "logistic"))`:

  Discrete-time survival hazard by person-period expansion. `breaks`
  sets the period grid: `NULL` uses the sorted distinct observed times,
  a single positive integer bins at that many quantiles, and a numeric
  boundary vector gives explicit right-closed intervals. `max.rows`
  refuses an over-large expansion, naming the coarsening levers. `link`
  selects the binary link the expanded rows are fit under, and is what
  the tokens `"hazard.probit"` and `"hazard.logistic"` name. The former
  `breaks` and `max.rows` arguments.

- `hurdle.lognormal()`:

  A semicontinuous two-part response: an occupancy probit glued to a
  lognormal positive part. Composed from two samplers, so only
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) fits it.

## Details

Which entry points accept which family:

|                    |          |            |                |           |
|--------------------|----------|------------|----------------|-----------|
| **family**         | **bart** | **dbarts** | **dbartsSpec** | **xbart** |
| `gaussian`         | yes      | yes        | yes            | yes       |
| `student`          | yes      | yes        | yes            | no        |
| `probit`           | yes      | yes        | yes            | yes       |
| `logistic`         | yes      | yes        | yes            | yes       |
| `aft`              | yes      | yes        | yes            | no        |
| `multinomial`      | yes      | yes        | yes            | no        |
| `ordinal`          | yes      | yes        | yes            | no        |
| `nbinom`           | yes      | yes        | yes            | no        |
| `hazard`           | yes      | yes        | no             | no        |
| `hurdle.lognormal` | yes      | no         | no             | no        |

A family the entry point does not fit is refused by name. See
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
`family` argument for what each family models and what it supports.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)

## Examples

``` r
f <- dbartsFamilies$student(df = 4)
f
#> dbarts response family: student(df = 4)

set.seed(99)
n <- 100L
x <- matrix(rnorm(n * 2L), n, 2L)
y <- x[, 1L] + rt(n, 4)

fit <- bart(
  x, y,
  family = student(df = 4),
  n.samples = 50L, n.burn = 50L, n.chains = 1L, verbose = FALSE
)
```
