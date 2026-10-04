# Bayesian Additive Regression Trees with Random Effects

**Deprecated.** `rbart_vi` and its methods run the implementation from
dbarts 0.9-x, kept for one release, and are removed in dbarts 1.1-0. A
warning (class `dbartsDeprecatedWarning`) says so the first time
`rbart_vi` is called in a session. Grouped random effects live in the
stan4bart package (`stan4bart::stan4bart`), whose prior on the group
spread differs from the one used here, so a refit there moves the
results rather than reproducing them.

Fits a varying intercept/random effect BART model. A fallback this
function makes (running single-threaded, disabling verbose output for
several threads, recycling `group.by`, drawing effects for unseen
levels, a failed rejection sample) warns with class
`dbartsFallbackWarning`, and the deprecated `value` argument and
`type = "post-mean"` of `predict` warn with class
`dbartsDeprecatedWarning`.

## Usage

``` r
rbart_vi(
    formula, data, test, subset, weights, offset, offset.test = offset,
    group.by, group.by.test, prior = cauchy,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = 2.0,
    power = 2.0, base = 0.95,
    n.trees = 75L,
    n.samples = 1500L, n.burn = 1500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = FALSE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 5L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE,
    keepTrees = TRUE, keepCall = TRUE,
    seed = NA_integer_,
    keepSampler = keepTrees,
    keepTestFits = TRUE,
    callback = NULL,
    ...)

# S3 method for class 'rbart'
plot(
    x, plquants = c(0.05, 0.95), cols = c('blue', 'black'), ...)

# S3 method for class 'rbart'
fitted(
    object,
    type = c("ev", "ppd", "bart", "ranef"),
    sample = c("train", "test"),
    ...)

# S3 method for class 'rbart'
extract(
    object,
    type = c("ev", "ppd", "bart", "ranef", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE,
    ...)

# S3 method for class 'rbart'
predict(
    object, newdata, group.by, offset,
    type = c("ev", "ppd", "bart", "ranef"),
    combineChains = TRUE,
    ...)

# S3 method for class 'rbart'
residuals(object, ...)

# S3 method for class 'rbart'
print(x, ...)
```

## Arguments

- group.by:

  Grouping factor. Can be an integer vector/factor, or a reference to
  such in `data`. In `predict`, a character or numeric vector names its
  groups by value, as a factor of it does.

- group.by.test:

  Grouping factor for test data, of the same type as `group.by`. Can be
  missing.

- prior:

  A function or symbolic reference to built-in priors. A built-in is
  named by symbol or string (`cauchy`, `"gamma"`); any other name refers
  to the caller's function, and an unknown string is refused. Determines
  the prior over the standard deviation of the random effects. Supplied
  functions take two arguments, `x` - the standard deviation, and
  `rel.scale` - the standard deviation of the response variable before
  random effects are fit. Built in priors are `cauchy` with a scale of
  2.5 times the relative scale and `gamma` with a shape of 2.5 and scale
  of 2.5 times the relative scale.

- n.thin:

  The number of tree jumps taken for every stored sample, but also the
  number of samples from the posterior of the standard deviation of the
  random effects before one is kept.

- keepTestFits:

  Logical where, if false, test fits are obtained while running but not
  returned. Useful with `callback`.

- callback:

  Optional function of `trainFits`, `testFits`, `ranef`, `sigma`, and
  `tau`. Called after every post-burn-in iteration and the results of
  which are collected and stored in the final object.

- formula, data, test, subset, weights, offset, offset.test, sigest,
  sigdf, sigquant, k, power, base, n.trees, n.samples, n.burn, n.chains,
  n.threads, combineChains, n.cuts, useQuantiles, keepTrainingFits,
  printEvery, printCutoffs, verbose, keepTrees, keepCall, seed,
  keepSampler, ...:

  Same as in
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md), except
  `sigdf`, `sigquant`, `power` and `base`, which `bart` no longer takes
  and which are as in
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md). `k`
  is applied only when supplied; otherwise the defaults are those of
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md). A
  supplied `seed` is used for this call only: the caller's random number
  stream is left as it was found. Unlike `bart`, the thread count
  changes the draws: with more than one thread the chains run in worker
  processes, each seeded from `seed` or, without one, from R's random
  number stream, so `set.seed` reproduces a run at the same `n.threads`,
  while one thread runs the chains in turn on R's stream. `weights` are
  precisions, so a row's residual variance is \\\sigma^2 / w_i\\, and a
  binary response takes only weights of 0 and 1; other weights are
  refused. A `dbartsData` object carrying `bases` (a multi-forest model)
  is refused: `rbart_vi` fits a single forest.

- object:

  A fitted `rbart` model.

- newdata:

  Same as `test`, but named to match
  [`predict`](https://rdrr.io/r/stats/predict.html) generic.

- type:

  One of `"ev"`, `"ppd"`, `"bart"`, `"ranef"`, or `"trees"` for the
  posterior of the expected value, posterior predictive distribution,
  non-parametric/BART component, random effect, or saved trees
  respectively. The expected value is the sum of the BART component and
  the random effects, while the posterior predictive distribution is a
  response sampled with that mean. To synergize with
  [`predict.glm`](https://rdrr.io/r/stats/predict.glm.html),
  `"response"` can be used as a synonym for `"ev"` and `"link"` can be
  used as a synonym for `"bart"`. For additional details on tree
  extraction, see the corresponding subsection in
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md).

- sample:

  One of `"train"` or `"test"`, referring to the training or tests
  samples respectively.

- x, plquants, cols:

  Same as in
  [`plot.bart`](https://vdorie.github.io/dbarts/reference/bartBT.md).

## Details

Fits a BART model with additive random intercepts, one for each factor
level of `group.by`. For continuous responses:

- \\y_i \sim N(f(x_i) + \alpha\_{g\[i\]}, \sigma^2)\\

- \\\alpha_j \sim N(0, \tau^2)\\.

For binary outcomes the response model is changed to \\P(Y_i = 1) =
\Phi(f(x_i) + \alpha\_{g\[i\]})\\. \\i\\ indexes observations,
\\g\[i\]\\ is the group index of observation \\i\\, \\f(x)\\ and
\\\sigma_y\\ come from a BART model, and \\\alpha_j\\ are the
independent and identically distributed random intercepts. Draws from
the posterior of \\tau\\ are made using a slice sampler, with a width
dynamically determined by assessing the curvature of the posterior
distribution at its mode.

### Out Of Sample Groups

Predicting random effects for groups not in the training sample is
supported by sampling from their posterior predictive distribution, that
is a draw is taken from \\p(\alpha \mid y) = \int p(\alpha \mid
\tau)p(\tau \mid y)d\alpha\\. For out-of-sample groups in the test data,
these random effect draws can be kept with the saved object. For those
supplied to `predict`, they cannot and may change for subsequent calls.

### Data

Data are handled as in dbarts 0.9-x and not as
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) does: a
factor predictor is expanded to indicator columns, rows with missing
values are dropped (`na.omit`), and the response must be continuous or
binary.

### Generics

See the generics section of
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md).

## Value

An object of class `rbart`. Contains all of the same elements of an
object of class
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md), as well as
the elements:

- ranef:

  Samples from the posterior of the random effects. A array/matrix of
  posterior samples. The \\(k, l, j)\\ value is the \\l\\th draw of the
  posterior of the random effect for group \\j\\ (i.e. \\\alpha^\*\_j\\)
  corresponding to chain \\k\\. When `n.chains` is one or
  `combineChains` is `TRUE`, the result is a collapsed down to a matrix.

- ranef.mean:

  Posterior mean of random effects, derived by taking mean across group
  index of samples.

- tau:

  Matrix of posterior samples of `tau`, the standard deviation of the
  random effects. Dimensions are equal to the number of chains times the
  numbers of samples unless `n.chains` is one or `combineChains` is
  `TRUE`.

- `first.tau`:

  Burn-in draws of `tau`.

- `callback`:

  Optional results of `callback` function.

## Author

Vincent Dorie: <vdorie@gmail.com>

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)

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

n.g <- 10
g <- sample(n.g, length(y), replace = TRUE)
sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

y <- y + b[g]

df <- as.data.frame(x)
colnames(df) <- paste0("x_", seq_len(ncol(x)))
df$y <- y
df$g <- g

## low numbers to reduce run time
rbartFit <- suppressWarnings(
    rbart_vi(y ~ . - g, df, group.by = g,
             n.samples = 40L, n.burn = 10L, n.thin = 2L,
             n.chains = 1L,
             n.trees = 25L, n.threads = 1L))
#> family = "auto": continuous response detected, fitting family = "gaussian"; set 'family' to override
#> 
#> Running BART with numeric y
#> 
#> number of trees: 25
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 2
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.003078
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  level fibre gibbs step: auto
#>  proposal probabilities: birth/death 0.60, swap 0.00, change 0.40, perturb 0.00, rule_gibbs 0.00; birth 0.50
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 10
#>  init sigma: 3.181799, curr sigma: 3.181799
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) (3: 100) (4: 100) (5: 100) 
#> (6: 100) (7: 100) (8: 100) (9: 100) (10: 100) 
#> 
#> 
#> Running BART with numeric y
#> 
#> number of trees: 25
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 2
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.003078
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  level fibre gibbs step: auto
#>  proposal probabilities: birth/death 0.60, swap 0.00, change 0.40, perturb 0.00, rule_gibbs 0.00; birth 0.50
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 10
#>  init sigma: 3.181799, curr sigma: 3.181799
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) (3: 100) (4: 100) (5: 100) 
#> (6: 100) (7: 100) (8: 100) (9: 100) (10: 100) 
#> 
```
