# Treatment Forest Specification for Bayesian Causal Forests

Build the specification of the treatment forest of a Bayesian causal
forest. Pass the result as the `treatmentForest` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
alongside a `treatment` vector. One constructor carries every knob the
second forest and its glue own, so the fitting functions grow three
arguments rather than fourteen.

## Usage

``` r
treatmentForest(
    n.trees = 50L, base = 0.25, power = 3,
    sd.control = 2, sd.moderate = 1, b.prior.variance = 0.5,
    update.a = TRUE, update.b = TRUE,
    interactions = NULL, blocks = NULL)
```

## Arguments

- n.trees:

  Number of trees in the treatment forest. A single integer of at
  least 1. Fewer trees than the prognostic forest is the usual choice:
  the treatment surface is normally the smoother of the two.

- base, power:

  The treatment forest's tree-structure prior, as
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `tree.prior` is for the prognostic forest. The defaults are shallower
  than the prognostic forest's, for the same reason.

- sd.control, sd.moderate:

  Prior scales, in standard deviations of the response. `sd.control` is
  the half-Cauchy median of the prognostic scalar \\a\\, placing the
  prognostic total at `sd.control` response standard deviations;
  `sd.moderate` places the treatment effect \\(b_1 - b_0)\tau\\ at
  `sd.moderate` of them.

- b.prior.variance:

  Prior variance of the glue coefficients \\b_0\\ and \\b_1\\, which
  scale the treatment forest for control and treated observations.

- update.a, update.b:

  Whether the corresponding glue block is redrawn each sweep. `FALSE`
  fixes it at its prior center for the sampler's life; both are fixed at
  creation and cannot be toggled afterwards.

- interactions, blocks:

  Optional
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)
  and [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)
  constraints on the treatment forest. The arguments of the same names
  on [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  constrain the prognostic forest, so the two can be held to different
  structures - the calibrated-additivity idiom is an additive or
  low-order treatment forest beside a free prognostic one. A `blocks`
  partition here covers the columns the forest may split on, i.e. the
  `moderators` subset when one is given.

## Details

A Bayesian causal forest fits \\y = a\\\mu(x) + b_z\\\tau(x) +
\epsilon\\: a prognostic forest \\\mu\\ over every predictor, and a
treatment forest \\\tau\\ over the moderators, joined by a
three-parameter glue \\(a, b_0, b_1)\\. This constructor describes the
second forest and the glue; the prognostic forest takes the fitting
function's own `control@n.trees` and `tree.prior`.

Both forests' leaf scales come from the model's own calibration map
rather than from the node prior, which is why a `k` hyperprior, a
non-default `k`, and a linear or Gaussian-process node prior are refused
when a treatment vector is supplied. Every value here is validated at
fit time.

## Value

A `dbartsTreatmentForest` specification object, resolved when a sampler
is built.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md),
[`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)

## Examples

``` r
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
z <- rbinom(n, 1L, 0.5)
y <- 2 * x[, 1] + z * (1 + 2 * x[, 3]) + rnorm(n, 0, 0.2)

sampler <- dbarts(x, y, treatment = z, moderators = c("x1", "x3"),
                  treatmentForest = treatmentForest(n.trees = 25L,
                                                    sd.moderate = 1.5),
                  control = dbartsControl(n.chains = 1L, n.trees = 25L,
                                          n.samples = 20L, n.burn = 20L))
samples <- sampler$run(20L, 20L)
```
