# Forest Specification for Multi-Forest Models

Build the specification of one forest of a multi-forest model. Pass a
list of these as the `forests` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
which fits the mean as a weighted sum of the declared ensembles. Every
knob is per forest and carried on the one constructor, so the fitting
functions grow exactly one argument however many forests a model has.

The constructor is not exported: it resolves by bare name inside the
arguments that take it, and elsewhere is written
`dbartsForests$forest(...)`; see
[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md).

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
  its own environment, or an already-evaluated vector or matrix. It
  expands by R's model-matrix rule - a factor becomes its level
  indicators, one amplitude per level, with no reference level dropped,
  since the forest carries no intercept of its own - so the forest
  enters the mean as \\(\sum_k a_k B_k(x_i)) f(x_i)\\. Any forest may
  carry one, of any width: a two-level factor gives the pair whose
  amplitudes are \\(b_0, b_1)\\, a wider factor one amplitude per level,
  and a numeric vector or matrix is already those columns. Level order
  is meaningful: the second level is the one \\b_1\\ scales. At
  creation, a column no observation enters - all zeros, or a factor
  level no row takes - is refused, since its amplitude would move under
  its prior alone, as is a basis so large or small that a row's norm is
  not representable; `$setForestBasis` keeps an empty factor level (see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
  A forest that declares none takes the implicit intercept its single
  amplitude \\a\\ scales; every forest past the first needs one, since
  the amplitudes multiplying it are what distinguish it from the first.
  Reaching
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  through an already-built `dbartsData` `formula` is the one route that
  refuses this argument outright, since the declaration would otherwise
  have nowhere to ride and be silently discarded;
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
  always takes a pre-built data object and installs the declaration
  instead, replacing whatever bases it carried. On the formula
  interface, an already-evaluated vector or matrix must have one row per
  observation of the FULL `data` - the same rule a formula's own
  evaluation already follows - and is then restricted to the rows
  `subset` keeps, exactly as a predictor column is; a value already
  restricted to those rows instead, matching the subset's count rather
  than the full data's, is refused by name rather than aligned to them
  by position. With no `subset`, or one that keeps every row, the two
  counts agree and nothing changes.

  The right-hand side of a formula is evaluated as R code, and the
  forest does not centre or scale a basis, so a multiplier that should
  be standardized is written so: `basis = ~ scale(w)`, or
  `basis = ~ cbind(scale(w), scale(v))` for two columns.

  For a `forest()` term of a model formula, `predict` rebuilds the basis
  at the new rows from what
  [`scale()`](https://rdrr.io/r/base/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()` or `bs()`
  computed on the training rows, as
  [`lm`](https://rdrr.io/r/stats/lm.html) does, so a new row is
  standardized by the training centre and scale and not by the new rows'
  own. This holds for such a call written alone, or as an operand of
  `+`, `-`, `*`, `/`, `^`, parentheses or
  [`cbind()`](https://rdrr.io/r/base/cbind.html); any other expression,
  such as `I(w - mean(w))`, `scale(w)[, 1]` or `abs(scale(w))`, is
  evaluated on the rows given to `predict`, as `lm` evaluates it, so its
  value at a row depends on the other rows predicted with it; see
  [`SafePrediction`](https://rdrr.io/r/stats/makepredictcall.html). With
  `subset`, the centre and scale are those of the rows `subset` kept,
  where `lm` uses all of them; rows dropped for a missing response still
  enter them. A factor level the fit never saw is refused at `predict`.

  A basis declared through `forests =` is not rebuilt: the caller gives
  it at the new rows, and with `subset` its formula was evaluated on
  every row of `data`.

- vars:

  Optional restriction of this forest to a subset of the model matrix,
  by column name or index; `NULL` leaves it reading every predictor. Any
  forest may be restricted, a single declared forest included: given the
  same residual scale estimate, that fit is the fit on the named columns
  alone, and its other columns report no splits. The estimate is taken
  from every column unless `sigest` states it, so the two fits agree
  draw for draw when it does. On a single forest a `dart` tree prior
  beside `vars` lays its Dirichlet over the named columns, every other
  column reporting probability 0, and `split.probs` keeps its ratios
  among the named columns and must give one of them a positive
  probability (see
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)).
  On the single forest of a hazard fit the `period` column the fit
  appends is always allowed, named or not: `vars` restricts the columns
  the caller supplied. A hazard fit with several forests takes each
  forest's `vars` as written, so a restricted forest there splits on
  `period` only when `vars` names it. The restriction holds for the
  sampler's life; `$setModel` refuses a model that states another, or
  none.

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

  This forest's prior scale, in units of the response family's own
  latent scale and per unit of basis row norm: the total \\a\\f(x)\\ or
  \\(b_1 - b_0)f(x)\\ is placed at `sd` of them. The unit is
  \\\mathrm{sd}(y)\\ for a gaussian response, `1` for `"probit"` and
  \\\pi/\sqrt{3}\\ for `"logistic"` - the standard deviation of the
  link's own error law, a latent model having no response standard
  deviation to name. A forest whose `basis` rows have median non-zero
  norm \\c\\ contributes the scale named here, the calibration map
  dividing \\c\\ out, so rescaling a basis column does not silently
  rescale the prior. Which of the two channels carries it depends on
  whether the forest has a `basis`: without one it is the half-Cauchy
  median of the forest's scalar amplitude, with one it scales the
  forest's own leaf prior. The two channels default differently, and
  neither default is a bare constant. A forest with NO basis takes `2`
  under a gaussian response, where the unit is the response's own
  \\\mathrm{sd}(y)\\ and a drawn \\\sigma\\ absorbs the difference, and
  `1` under `"probit"` and `"logistic"`, where the unit is the link's
  fixed error scale and nothing does. A forest WITH a basis takes
  \\\sqrt{2/K}\\ in a model of \\K\\ forests, so that declaring more of
  them does not widen the prior on the combined location without bound.
  \\K = 2\\ is the fixed point of both statements, \\\sqrt{2/2} = 1\\. A
  value declared here overrides its default and keeps its per-forest
  reading at every \\K\\, so `sd = 1` on each basis forest recovers the
  pre-\\K\\-aware model exactly. It must be positive and finite.

  It is the one argument a live sampler restates:
  `$setLeafPrior(forests = list(forest(sd = ), ...))` writes it in the
  channel creation gave the forest, and `$getLeafPrior(f)$leaf.prior`
  reads it back as `forest(sd = )` (see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
  The channel is fixed at creation: a forest created without a basis
  keeps its half-Cauchy median after `$setForestBasis` gives it one, so
  its `forest(sd = )` round-trips through `$setLeafPrior`, while a fresh
  [`dbarts()`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  given the same bases would read the same `sd` as a leaf-scale factor
  and build a different prior.

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
  given, on the first forest as on any other.

- amplitude.prior.variance:

  Prior variance of the \\N(0, \cdot)\\ amplitudes on this forest's
  basis columns, default `0.5`. Legal only on a forest given a `basis`:
  a forest without one carries a plain scalar amplitude under the
  engine's half-Cauchy scale-mixture prior, whose median is `sd` rather
  than a variance. It is a free multiplier on the induced prior: the
  prior standard deviation of the combined location at row \\i\\ is
  \\\sqrt{\sum_f s_f^2 v_f \\B_f(i,\cdot)\\^2}\\ over the forests
  carrying a basis, with \\s_f\\ read from `$getLeafPrior(f)$k.scale`
  and \\v_f\\ this argument; a basis-free forest's own term is Cauchy
  and has no standard deviation. The budget that sum sits in is set by
  `sd`, whose default already divides it among the \\K\\ forests, so
  raising this argument raises the total rather than redistributing it.
  Under `"probit"` and `"logistic"` that location IS the latent index
  and \\\sigma\\ is pinned, so nothing in the sampler absorbs a
  mis-scaled basis; under a gaussian response it is in
  \\\mathrm{sd}(y)\\ units and a drawn \\\sigma\\ partly does. Every
  input to that expression is readable off the fitted sampler: \\v_f\\
  is the `amplitude.prior.variance` entry of `$getLeafPrior(f)` and
  \\B_f\\ is `data@bases[[f]]`, so the induced prior can be checked
  against what is in force rather than against what the call asked for.
  See the example below. A lone forest carrying a `basis` is refused,
  the amplitudes being what distinguish a forest from another. For
  varying coefficients declare an intercept forest plus one basis forest
  per covariate,
  `forests = list(forest(), forest(basis = ~ z1), forest(basis = ~ z2))`,
  or use one forest with `linear()` leaves.

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
b_1)\\ joining them, read back with `$getForestAmplitudes`. Gaussian,
`"probit"` and `"logistic"` responses; under a latent family the
combination is the index rather than the mean, on the link's own fixed
scale, and `"aft"`, `"ordinal"` and `"nbinom"` are refused at creation
naming what each is missing.

What the defaults put on the combined location, since under a latent
family that location IS the index and no drawn \\\sigma\\ stands between
it and the fitted probabilities. At two forests of the shipped shape -
one carrying no basis - a probit model's prior puts \\P(p \< 0.01
\mathrm{~or~} p \> 0.99)\\ at 0.238, which is the shipped single-forest
binary default's own 0.239; before the `sd` defaults above it was 0.376.
The \\\sqrt{2/K}\\ factor is what holds that as \\K\\ grows, and it
holds it in two different senses. When EVERY forest carries a basis the
induced prior standard deviation of the index is 1.484 latent units at
every \\K\\ - 0.989 of the classic \\k = 2\\ binary leaf-scale budget -
because the whole location is then a sum of fixed-variance channels.
When one forest carries none, its amplitude is Cauchy and has no
variance to enter that budget with, so the fixed-variance part is
BOUNDED by 1.484 rather than pinned at it, rising from 0.699 of the
budget at \\K = 2\\ toward 0.989 and never reaching it; without the
factor it would instead grow past twice the budget by ten forests. Read
the values in force off `$getLeafPrior(f)`'s `leaf.scale.factor` and
`amplitude.prior.scale` entries.

Both forests' leaf scales come from the model's own calibration map
rather than from the leaf prior, which is why a `k` hyperprior, a
non-default `k`, a named leaf-prior `sd`, and a linear or
Gaussian-process leaf prior are refused when a second forest is
declared. This `sd` is not the leaf prior's: it states this forest's
share of the combined location's prior, per unit of basis row norm,
where `normal(sd = )` states a single forest's whole spread. Every value
here is validated at fit time, and anything today's engine cannot honour
is refused there by name rather than dropped.

A `forest()` call written INSIDE a
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)/[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
formula - bare, or as one operand of `:` - declares a second forest
without a separate `forests =` list: `z:forest(x1 + x2)` and
`factor(z):forest(x1 + x2)` desugar to `forest(x1 + x2, basis = ~ z)`
and `forest(x1 + x2, basis = ~ factor(z))`, a `(a + b):forest(x1)`
compound left operand (every member numeric or logical) desugars to one
forest with the two-column basis `~ cbind(a, b)`, and `z * forest(x1)`
is refused, naming both explicit spellings. In this position
`forest()`'s UNNAMED slot is a symbolic predictor SET, not `basis`:
`forest(x1 + x2)` means `vars = c("x1", "x2")`, read from the call and
never evaluated, so the docs always spell `basis =` by name inside a
term. Every other named argument - including `basis` itself, for the
general form `forest(x1 + x2, basis = ~ z)` - evaluates normally. A
symbolic name must be among the formula's own right-hand-side terms; a
factor name expands to every indicator column derived from it. A term's
basis is evaluated against the fit's model frame, after `subset` and any
row selection; a `basis` declared directly on `forests =` is instead
evaluated against the raw `data` and then restricted to the same rows -
the two routes reach the same rows by different timing, so a `subset`
lines them up either way.

## Value

A `dbartsForest` specification object, resolved when a sampler is built.

## See also

[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md),
[`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s ‘Formula
Terms’ section is the full term grammar - colon sugar, the factor forms,
the symbolic unnamed slot, and the refusal list - of which the ‘Details’
paragraph above is the constructor-side summary.

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

# the induced prior sd of the combined location, read off the sampler rather
# than recomputed from the call: forest 2 is the one carrying a basis
calibration <- sampler$getLeafPrior(2L)
basis <- sampler$data@bases[[2L]]
indexSd <- sqrt(calibration$k.scale^2 *
                calibration$amplitude.prior.variance *
                rowSums(basis^2))

# a standardized multiplier: predict centres and scales the new rows by the
# training rows' mean and sd of w
d <- data.frame(y = y, x1 = x[, 1], x3 = x[, 3], w = 50 + 10 * rnorm(n))
fit <- bart(y ~ x1 + x3 + forest(x1 + x3, basis = ~ scale(w), n.trees = 10L),
            d, n.trees = 10L, n.samples = 10L, n.burn = 10L,
            n.chains = 1L, n.threads = 1L, keepTrees = TRUE)
#> family = "auto": continuous response detected, fitting family = "gaussian"; set 'family' to override
#> 
#> Running BART with numeric y
#> 
#> number of trees: 10
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 1
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.008161
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  level fibre gibbs step: auto
#>  proposal probabilities: birth/death 0.60, swap 0.00, change 0.40, perturb 0.00, rule_gibbs 0.00; birth 0.50
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 2
#>  init sigma: 0.948004, curr sigma: 0.948004
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) 
#> Running mcmc loop:
#> total seconds in loop: 0.000546
#> 
#> Tree sizes, last iteration:
#> [1] 2 2 3 3 3 1 4 2 2 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 9) (2: 5) 
#> DONE BART
#> 
predict(fit, d[1:3, ])
#>              1        2        3
#>  [1,] 2.665598 1.435529 1.428253
#>  [2,] 2.718696 1.279158 1.786205
#>  [3,] 2.825532 1.220816 1.820275
#>  [4,] 2.639906 1.695266 1.912748
#>  [5,] 2.711947 1.647589 1.116915
#>  [6,] 2.533666 1.730959 1.593382
#>  [7,] 2.268279 2.188020 2.531138
#>  [8,] 2.679228 1.709843 1.749938
#>  [9,] 2.466707 1.802121 2.166133
#> [10,] 2.533577 1.228568 1.915270
```
