# Prior Specification Constructors

A list of constructor functions building the prior specifications that
the `tree.prior`, `node.prior`, and `resid.prior` arguments of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and the
fitting functions accept. Bundling them keeps generic names like
`normal` out of the search path, where another package could mask them
or be masked depending on attach order.

## Format

A list of functions:

- `cgm(power = 2, base = 0.95, split.probs = NULL)`:

  The Chipman, George, and McCulloch tree prior. `split.probs` weights
  the choice of split variable: `NULL` is uniform, a named vector
  assigns by column or term name (with an optional `".default"`
  element), an unnamed vector assigns by position. Named or positional
  probabilities are matched to the data when a sampler is built.

- `dart(power = 2, base = 0.95, a = 0.5, b = 1, rho = NULL, alpha = 1, update.alpha = TRUE, update.delay = NULL)`:

  The CGM structure prior with DART (Linero 2018): a Dirichlet prior
  over the split-variable probabilities inducing variable selection.
  `alpha` is the concentration, optionally sampled (`update.alpha`) on a
  grid with a Beta(`a`, `b`) prior on `alpha / (alpha + rho)`; `rho`
  defaults to the number of predictors. Updates hold until
  `update.delay` iterations have passed (default: half the control's
  burn-in), so the forest is likelihood-informed when counts first enter
  the Dirichlet.

- `normal(k = NULL)`:

  Normal prior on the node means. `k` scales the standard deviation and
  can be a positive scalar, a hyperprior built with `chi`, or `NULL` for
  the default: 2 for continuous responses, `chi(1.25, Inf)` for binary
  ones. The continuous default follows Chipman, George, and McCulloch's
  argument that with node standard deviation
  `sigma_mu = 0.5 / (k * sqrt(m))` for `m` trees, `k` prior standard
  deviations of \\f(x)\\ span the whole response range regardless of
  `m`; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  Details for the response-scaling caveat this relies on.

- `linear(columns, k = NULL)`:

  Each leaf fits an intercept plus a linear term in the designated
  continuous predictor columns instead of a constant, so the forest
  models smoothly-varying coefficients. `columns` names model matrix
  columns (character) or indexes them (numeric) and is matched to the
  data when a sampler is built; factor columns cannot be designated. The
  covariates are standardized internally and every coefficient shares
  the `normal(k)` prior. Reported leaf values keep the intercept;
  `getTrees` adds one `beta.<column>` column per covariate.
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) accepts
  the same specification through its own `node.prior` argument.

- `gp(columns, k = NULL, lengthscale = NULL, max.leaf.size = 256L)`:

  Each leaf fits a smooth Gaussian-process function of the designated
  continuous predictor columns, drawn under a squared-exponential kernel
  whose prior scale ties to `k` exactly as the other leaf priors;
  `columns` resolves as for `linear`, and `k` may be fixed or given the
  `chi` hyperprior, sampled from the drawn functions' standardized
  magnitudes. `lengthscale` fixes the kernel lengthscales on the
  standardized covariate scale, one value per column or one recycled
  scalar; `NULL` uses the median pairwise-distance heuristic per column.
  Leaves holding more than `max.leaf.size` observations fall back to
  constant fits, confining the cubic kernel cost. Function-valued fits
  ride prediction only: `getTrees` reports `NA` leaf values, and
  `keepTrees` storage grows with the leaf sizes.
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) accepts
  the same specification through its own `node.prior` argument.

- `chi(degreesOfFreedom = 1.25, scale = Inf)`:

  Chi hyperprior over `k`, sampled along with the rest of the model.

- `chisq(df = 3, quant = 0.9)`:

  Chi-squared prior on the residual variance, scaled so that an estimate
  of the residual standard deviation falls at the given quantile.

- `fixed(value = 1)`:

  Fixed residual variance.

## Details

Inside the prior arguments of the fitting functions the same
constructors are available by bare name, so
`dbarts(..., node.prior = normal(chi(1.25)))` works regardless of what
packages are attached: those arguments are evaluated with this
vocabulary layered over the calling environment, along with `num.vars`,
the number of predictor columns. A prior object built ahead of time with
`dbartsPriors$...` can be passed to the same arguments.

## References

Gramacy, R.B. and Lee, H.K.H. (2008) Bayesian Treed Gaussian Process
Models With an Application to Computer Modeling. *Journal of the
American Statistical Association*, **103**(483), 1119–1130.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md)

## Examples

``` r
prior <- dbartsPriors$normal(dbartsPriors$chi(1.25))

x <- matrix(runif(200), ncol = 2)
y <- x[, 1] + rnorm(100, 0, 0.5)
sampler <- dbarts(y ~ x, node.prior = prior,
                  tree.prior = dbartsPriors$cgm(power = 1.5))

## DART: a Dirichlet prior over split-variable probabilities, useful
## when only a handful of predictors out of many are relevant
set.seed(0)
n <- 60L
x.dart <- matrix(runif(n * 8), n)
y.dart <- 4 * x.dart[, 1] + rnorm(n, 0, 0.2)
fit.dart <- dbarts(y.dart ~ x.dart, tree.prior = dart(),
                   control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                            n.threads = 1L))
samples.dart <- fit.dart$run(20L, 20L)

## linear leaves: each leaf fits an intercept plus a slope in x1,
## useful for a function that is smooth and roughly linear in a subset
## of the predictors
set.seed(1)
x1 <- runif(n)
x2 <- runif(n)
y.lin <- 3 * x1 + rnorm(n, 0, 0.2)
df.lin <- data.frame(x1, x2, y.lin)
fit.lin <- dbarts(y.lin ~ x1 + x2, df.lin, node.prior = linear("x1"),
                  control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                           n.threads = 1L))
samples.lin <- fit.lin$run(20L, 20L)

## gp leaves: each leaf fits a smooth Gaussian-process function of x1,
## useful for a smoothly-varying, nonlinear function
set.seed(2)
y.gp <- sin(2 * pi * x1) + rnorm(n, 0, 0.2)
df.gp <- data.frame(x1, x2, y.gp)
fit.gp <- dbarts(y.gp ~ x1 + x2, df.gp,
                 node.prior = gp("x1", max.leaf.size = 30L),
                 control = dbartsControl(n.trees = 10L, n.chains = 1L,
                                          n.threads = 1L))
samples.gp <- fit.gp$run(20L, 20L)
```
