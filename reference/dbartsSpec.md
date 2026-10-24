# Discrete Bayesian Additive Regression Trees Sampler Specification

Resolves the `control`/`model`/`data` triple and family token that
describe a sampler, without constructing one. Intended for packages that
embed dbarts - in particular `LinkingTo: dbarts` packages that build a
sampler in R from the returned pieces and hold its handle for the C API
in `dbarts.h` - and that supply their own design matrix.

## Usage

``` r
dbartsSpec(
    data, control = dbarts::dbartsControl(),
    tree.prior = cgm, leaf.prior = normal,
    proposal.probs = c(
        birth_death = 0.6, swap = 0, change = 0.4, perturb = 0,
        rule_gibbs = 0, birth = 0.5),
    monotone = NULL, interactions = NULL, blocks = NULL,
    variance = NULL,
    forests = NULL,
    sigest = NULL, seed = NULL,
    family = c("auto", "gaussian", "student", "probit", "logistic", "aft",
               "multinomial", "ordinal", "nbinom"),
    survival = NULL,
    parentEnv = parent.frame(), ...)
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
  other family. Where `family = "auto"` resolves, the family is reported
  in a one-line message when `control` was built with `verbose = TRUE`;
  this function has no `verbose` argument of its own.

- control, tree.prior, leaf.prior, proposal.probs, monotone,
  interactions, blocks, variance, forests, sigest, seed, family, ...:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  (except that `proposal.probs` is a formal of this function rather than
  of `control`), including the frozen `proposal.probs` mixture - all
  five structural probabilities zero, which proposes no structure and
  draws only the leaf values, `sigma` and the family's latents. The
  prior arguments are evaluated in dbarts's own prior vocabulary, so
  bare expressions such as `normal(k = chi(1.25, Inf))` resolve
  regardless of what the caller has attached; so do the forest
  constructors inside `monotone`, `interactions`, `blocks`, `variance`
  and `forests`, evaluated in `parentEnv` (see
  [`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md)).
  The residual prior rides `family` here as everywhere
  (`family = gaussian(sigma = chisq(df, quant))`); the 0.9-x spellings
  `resid.prior`, `sigma` and `node.prior` were never arguments of this
  function and are refused by name, naming
  `family = gaussian(sigma = )`, `sigest` and `leaf.prior` respectively.
  A [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  `basis` given as a one-sided formula is evaluated in `parentEnv`, this
  surface performing no data ingestion of its own; the data object's
  rows are already whatever its own `subset` kept. A `basis` declared
  here always reaches the sampler: unlike
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md), which
  refuses a basis declaration once `formula` is already a `dbartsData`,
  this function's first argument is always one, so the declaration is
  installed and REPLACES whatever bases the data object carried.
  `family`'s eight-token vocabulary is narrower than
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s by
  design - multinomial's K-forest engine and `hurdle.lognormal`'s
  two-sampler composition build through
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md), never
  through this construct-only bypass; the wider set lives there.

- survival:

  For `family = "aft"` only: a vector of per-observation status codes, 1
  for an observed event and 0 for a right-censored observation, with the
  response holding log event/censoring times.
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) reads
  this off a `Surv` or two-column response; a caller supplying log-times
  directly supplies it here. `NULL` for every other family.

- parentEnv:

  The environment in which the prior and forest constructor arguments
  are evaluated; defaults to the caller's frame, which is almost always
  what is wanted. Pass explicitly when forwarding arguments from another
  function.

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

The `family` element of the returned list names the family the sampler
is built with, and the `dbarts_family` enumerator
`dbarts_sampler_family` reports for it (`"gaussian"` names
`DBARTS_FAMILY_GAUSSIAN`, and so on). `"auto"` asks this function to
dispatch on the shape of the response, which is correct for
`"gaussian"`, `"probit"`, and `"ordinal"` - each is inferable from the
response coding or from an attribute - but is *wrong* for `"aft"` and
`"logistic"`, which are indistinguishable by shape from `"gaussian"` and
`"probit"` respectively, and for `"nbinom"`, which is never inferred (a
count vector resolves to `"gaussian"`). Name the family and the question
does not arise.

### Differences from dbarts()

`sigest` defaults to `NULL`, which like `"auto"` leaves the data
object's own value alone, rather than overwriting it: a caller that
computed a starting estimate keeps it, and an unset (`NA`) value is
still estimated during resolution. A number replaces the data object's
value, and `"dense"` or `"sparse"` estimates over it. An explicit `NA`
for `sigest` is a missing value and is refused. Cut points are taken
from `control` only when the data does not already carry resolved
per-column counts.

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

- active:

  the active-row mask a `"probit"`, `"ordinal"` or `"nbinom"` fit's 0/1
  case weights resolved to, `NULL` when they resolved to none; the
  `data` object's weights slot is cleared either way, since neither
  family carries one, so a caller building a sampler from this list must
  install the mask on it. For a count (`"nbinom"`) list,
  `sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)`
  followed by `sampler$setActiveRows(spec$active, updateScale = TRUE)`
  makes the sampler `dbarts` makes: the center is taken over the rows
  the mask keeps, and the two are
  [`identical()`](https://rdrr.io/r/base/identical.html) in the leaf
  prior, `$model`, `$data`, `$activeRows`, the stored state and every
  array `run()` returns from the same seed. A caller who redraws the
  mask each sweep and wants the center re-derived in burn-in and held
  afterward passes `updateScale = i <= n.burn` to that call. Where the
  mask keeps fewer than two distinct counts the call keeps the center
  and warns, so the by-hand route has warned twice, here and at
  `dbartsSpec`, and has spent the sampler's one warning: a later such
  call is silent on it and warns once on `dbarts`'s sampler. For a
  probit or ordinal list `sampler$setActiveRows(spec$active)` with no
  argument makes it `dbarts`'s, there being nothing to re-derive, and a
  list that carries several forests is refused the argument. A list
  whose leaves take covariates (`linear`, `gp`) takes the mask at
  creation,
  `new("dbartsSampler", spec$control, spec$model, spec$data, active = spec$active)`,
  which installs it as `dbarts` does and reads the covariates' centers
  and scales over the rows it keeps; created without it, the sampler
  reads them over every row. A list of several forests, one of them on a
  numeric column that states no `sd`, takes the mask at creation,
  `new("dbartsSampler", spec$control, spec$model, spec$data, active = spec$active)`,
  which installs it as `dbarts` does and reads the column's scale over
  the rows it keeps; a mask installed after creation moves no forest's
  prior. A Gaussian or Student-t list has no mask, its weights riding
  the data object into creation, and needs no call. Where the rows its
  weights keep hold fewer than two distinct values of the response and a
  forest on a numeric column states no `sd`, that route warns twice,
  here for the response's scale and center and at creation for those and
  the column's scale, where `dbarts` warns once.

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
