# Forest Specification for Multi-Forest Models

Build the specification of one forest of a multi-forest model. Pass a
list of these as the `forests` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
which fits the mean as a weighted sum of the declared ensembles. Every
knob is per forest and carried on the one constructor, so the fitting
functions grow exactly one argument however many forests a model has.

## Usage

``` r
forest(
    basis = NULL, vars = NULL,
    n.trees = NULL, base = NULL, power = NULL, sd = NULL,
    interactions = NULL, blocks = NULL,
    amplitude.prior.variance = NULL, update.amplitude = NULL)
```

## Arguments

- basis:

  The data this forest's amplitudes multiply: a one-sided formula
  (`~ factor(z)`), evaluated against the fit's own `data` and then in
  its own environment, or an already-evaluated vector. It expands by R's
  model-matrix rule - a factor becomes its level indicators, one
  amplitude per level - so the forest enters the mean as \\(\sum_k a_k
  B_k(x_i)) f(x_i)\\. Today's engine fits exactly one shape: a two-level
  factor on the SECOND forest, whose two indicators carry the amplitudes
  \\(b_0, b_1)\\. A numeric basis, a factor of any other width, and a
  basis on the first forest are each refused at creation. Level order is
  meaningful: the second level is the one \\b_1\\ scales. The first
  forest takes no `basis` - its is the implicit intercept its single
  amplitude \\a\\ scales.

- vars:

  Optional restriction of this forest to a subset of the model matrix,
  by column name or index; `NULL` leaves it reading every predictor.
  Only a forest with a `basis` may be restricted, today's engine reading
  the whole design in the first forest.

- n.trees, base, power:

  This forest's tree count and tree-structure prior. `NULL` takes the
  engine's default, which for a basis forest is 50 trees at
  `base = 0.25`, `power = 3` - shallower and fewer than a prognostic
  forest's, the modulating surface normally being the smoother of the
  two. On the FIRST forest these are the fit's own `control@n.trees` and
  `tree.prior`, which they restate rather than add to. When a value here
  disagrees with an explicitly supplied `control@n.trees` or tree-prior
  `base`/`power`, this one governs the fit, being the more specific of
  the two declarations.

- sd:

  This forest's prior scale, in standard deviations of the response: the
  total \\a\\f(x)\\ or \\(b_1 - b_0)f(x)\\ is placed at `sd` of them.
  The defaults are 2 for the first forest and 1 for a basis forest.

- interactions, blocks:

  Optional
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)
  and [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)
  constraints on this forest. The arguments of the same names on
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) are
  the FIRST forest's, so declaring both there is refused as one
  constraint given twice; a second forest's are only expressible here.
  Holding the two forests to different structures is the
  calibrated-additivity idiom - an additive or low-order modulating
  forest beside a free prognostic one. A `blocks` partition covers the
  columns the forest may split on, i.e. the `vars` subset when one is
  given.

- amplitude.prior.variance:

  Prior variance of the \\N(0, \cdot)\\ amplitudes on this forest's
  basis indicators, default `0.5`. Legal only on a forest given a
  `basis`: the first forest's amplitude carries the engine's own
  half-Cauchy scale-mixture prior, which has no caller slot.

- update.amplitude:

  Whether this forest's amplitudes are redrawn each sweep, default
  `TRUE`. `FALSE` fixes them at their prior center for the sampler's
  life; the choice is made at creation and cannot be toggled afterwards.
  Declaring it needs a model with amplitudes, so it is refused on a
  single-forest `forests`.

## Details

With two forests, the second carrying a two-level factor basis, this is
the Bayesian causal forest \\y = a\\\mu(x) + b_z\\\tau(x) + \epsilon\\:
a prognostic forest \\\mu\\ over every predictor, a modulating forest
\\\tau\\ over the columns `vars` allows, and the amplitudes \\(a, b_0,
b_1)\\ joining them, read back with `$getForestAmplitudes`. Gaussian
responses only.

Both forests' leaf scales come from the model's own calibration map
rather than from the node prior, which is why a `k` hyperprior, a
non-default `k`, and a linear or Gaussian-process node prior are refused
when a second forest is declared. Every value here is validated at fit
time, and anything today's engine cannot honour is refused there by name
rather than dropped.

## Value

A `dbartsForest` specification object, resolved when a sampler is built.

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

sampler <- dbarts(x, y,
                  forests = list(forest(),
                                 forest(basis = ~ factor(z),
                                        vars = c("x1", "x3"),
                                        n.trees = 25L, sd = 1.5)),
                  control = dbartsControl(n.chains = 1L, n.trees = 25L,
                                          n.samples = 20L, n.burn = 20L))
samples <- sampler$run(20L, 20L)
amplitudes <- sampler$getForestAmplitudes()
```
