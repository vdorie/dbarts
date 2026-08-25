# Discrete Bayesian Additive Regression Trees Sampler Specification

Resolves the `control`/`model`/`data` triple and family token that
describe a sampler, without constructing one. Intended for packages that
embed dbarts - in particular `LinkingTo: dbarts` packages that create
and hold a sampler through the C API in `dbarts.h` - and that supply
their own design matrix.

## Usage

``` r
dbartsSpec(
    data, control = dbarts::dbartsControl(),
    tree.prior = cgm, node.prior = normal,
    resid.prior = chisq, resid.dist = gaussian,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL, interactions = NULL, blocks = NULL,
    variance = NULL,
    forests = NULL,
    sigma = NA_real_, seed = NA_integer_,
    family = c("auto", "gaussian", "probit", "logistic", "aft",
               "multinomial", "ordinal", "nbinom"),
    dispersion = NA_real_, survival = NULL,
    parentEnv = parent.frame())
```

## Arguments

- data:

  A `dbartsData` object, as built by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md).
  Unlike
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md), this
  function does no response ingestion of its own: the response must
  already be materialized. A data object carrying an \\n \times K\\
  `counts` matrix resolves to `family = "multinomial"` even from
  `"auto"` - the slot IS the declaration - and is refused under any
  other family.

- control, tree.prior, node.prior, resid.prior, resid.dist,
  proposal.probs, monotone, interactions, blocks, variance, forests,
  sigma, seed, family, dispersion:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  The prior arguments are evaluated in dbarts's own prior vocabulary, so
  bare expressions such as `normal(k = chi(1.25, Inf))` resolve
  regardless of what the caller has attached. A
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  `basis` given as a one-sided formula is evaluated in `parentEnv`, this
  surface performing no data ingestion of its own; the data object's
  rows are already whatever its own `subset` kept. A `basis` declared
  here always reaches the sampler: unlike
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md), which
  refuses a basis declaration once `formula` is already a `dbartsData`,
  this function's first argument is always one, so the declaration is
  installed and REPLACES whatever bases the data object carried.

- survival:

  For `family = "aft"` only: a vector of per-observation status codes, 1
  for an observed event and 0 for a right-censored observation, with the
  response holding log event/censoring times.
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) reads
  this off a `Surv` or two-column response; a caller supplying log-times
  directly supplies it here. `NULL` for every other family.

- parentEnv:

  The environment in which the prior arguments are evaluated; defaults
  to the caller's frame, which is almost always what is wanted. Pass
  explicitly when forwarding arguments from another function.

## Details

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) does two
separable things: it ingests a formula or a pair of matrices into a
`dbartsData`, and it then resolves that data into a full sampler
specification - choosing the response family, coding the response,
applying the family's weight policy, estimating a starting `sigma`,
parsing the priors, and attaching the attributes that carry
monotonicity, interaction and block constraints, a variance forest,
residual degrees of freedom, and survival status.
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
exposes the first half. This function exposes the second, and returns
the pieces rather than a sampler.

Both entry points call the same resolution on the `dbartsData` they end
up with, so a sampler built from this function's output draws the same
samples as the equivalent
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) call
given the SAME data object. The two can still resolve a given family
differently, because
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) performs
its own response ingestion first - coding a factor response into a count
matrix, for instance - while this function does none: a call this
function refuses for the raw response
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) can
still accept, by ingesting it into the shape this function requires.

### Use with the C API

The `family` element of the returned list is the token to pass as the
fourth argument of `dbarts_sampler_create`. Passing the empty string
there instead asks the library to dispatch on the shape of the response,
which is correct for `"gaussian"`, `"probit"`, `"ordinal"`, and
`"nbinom"` - each is inferable from the response coding or from an
attribute the library reads - but is *wrong* for `"aft"` and
`"logistic"`, which are indistinguishable by shape from `"gaussian"` and
`"probit"` respectively. Pass the returned token and the question does
not arise.

### Differences from dbarts()

`sigma` defaults to leaving the data object's own value alone, rather
than overwriting it: a caller that computed a starting estimate keeps
it, and an unset (`NA`) value is still estimated during resolution. Cut
points are taken from `control` only when the data does not already
carry resolved per-column counts.

Families that require ingestion this function does not perform are
unavailable: `"hazard"` needs person-period expansion, and
`"hurdle.lognormal"` describes more than one sampler. `"multinomial"` is
in this function's own vocabulary, but only for a `data` object whose
response is already an \\n \times K\\ `counts` matrix (see `data` above,
and
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)'s
`counts` argument); given a factor-coded or plain numeric response
instead, it is refused with a message naming the counts route, even
though [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
accepts the same factor by coding it into that shape first.

## Value

A list with components

- control:

  a `dbartsControl`, carrying any family and variance-forest attributes.

- model:

  a `dbartsModel`, carrying the parsed priors and any monotone,
  interaction, block, and residual-degrees-of-freedom attributes.

- data:

  a `dbartsData`, with the response coded and `sigma` resolved.

- family:

  the resolved family token, never `"auto"`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md),
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md),
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)

## Examples

``` r
f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
      10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

data <- dbartsData(x, y)
control <- dbartsControl(n.chains = 1L, n.samples = 20L, n.burn = 10L,
                         n.trees = 25L, updateState = FALSE)

## an ordinary specification, then the sampler it describes
spec <- dbartsSpec(data, control = control)
spec$family
#> [1] "gaussian"

sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
samples <- sampler$run()

## additive-only fit, monotone in the fourth predictor
constrained <- dbartsSpec(data, control = control,
                          interactions = interactions(max.order = 1L),
                          monotone = c(rep(0, 3), 1, rep(0, 6)))
attr(constrained$model, "interaction.max.order")
#> [1] 1
```
