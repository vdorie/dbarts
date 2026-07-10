# Bayesian Additive Regression Trees with Random Effects

Fits a varying intercept/random effect BART model.

## Usage

``` r
rbart_vi(
    formula, data, test, subset, weights, offset, offset.test = offset,
    group.by, group.by.test, prior = cauchy,
    sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
    k = 2.0,
    power = 2.0, base = 0.95,
    split.probs = 1 / num.vars,
    dart = FALSE,
    n.trees = 75L,
    n.samples = 1500L, n.burn = 1500L,
    n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
    combineChains = TRUE,
    n.cuts = 100L, useQuantiles = FALSE,
    n.thin = 5L, keepTrainingFits = TRUE,
    printEvery = 100L, printCutoffs = 0L,
    verbose = TRUE,
    keepTrees = TRUE, keepCall = TRUE,
    seed = NA_integer_,
    keepSampler = keepTrees,
    keepTestFits = TRUE,
    callback = NULL,
    factors = c("categorical", "indicators"),
    missing = c("incorporate", "error"),
    ...)

# S3 method for class 'rbart'
plot(
    x, plquants = c(0.05, 0.95), cols = c('blue', 'black'), ...)

# S3 method for class 'rbart'
fitted(
    object,
    type = c("ev", "ppd", "bart", "ranef"),
    sample = c("train", "test"),
    ci.level = NULL,
    ...)

# S3 method for class 'rbart'
extract(
    object,
    type = c("ev", "ppd", "bart", "loglik", "ranef", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE,
    ...)

# S3 method for class 'rbart'
predict(
    object, newdata, group.by, offset, weights,
    type = c("ev", "ppd", "bart", "ranef"),
    combineChains = TRUE,
    ci.level = NULL,
    ...)

# S3 method for class 'rbart'
residuals(object, type = "ev", ...)
```

## Arguments

- group.by:

  Grouping factor. Can be an integer vector/factor, or a reference to
  such in `data`.

- group.by.test:

  Grouping factor for test data, of the same type as `group.by`. Can be
  missing.

- prior:

  A function or symbolic reference to built-in priors. Determines the
  prior over the standard deviation of the random effects. Supplied
  functions take two arguments, `x` - the standard deviation, and
  `rel.scale` - the standard deviation of the response variable before
  random effects are fit. Built in priors are `cauchy` with a scale of
  2.5 times the relative scale and `gamma` with a shape of 2.5 and scale
  of 2.5 times the relative scale. With a built-in prior and no
  `callback`, the whole Gibbs sampler runs inside the engine on one
  multi-chain sampler; a custom prior function runs the random effect
  updates in R instead.

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
  sigdf, sigquant, power, base, split.probs, dart, n.trees, n.samples,
  n.burn, n.chains, n.threads, combineChains, n.cuts, useQuantiles,
  keepTrainingFits, printEvery, printCutoffs, verbose, keepTrees,
  keepCall, seed, keepSampler, factors, missing, ...:

  Same as in
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md), with one
  default difference: `keepTrees` defaults to `TRUE` here (`bart2`'s
  default is `FALSE`). With `dart`, split probability samples appear as
  `varprobs` on the fit. Unlike `bart2`, `rbart_vi` has no `family`
  argument - the response family is always resolved automatically as
  gaussian or probit, with no logistic option - and does not accept
  sparse (`Matrix::dgCMatrix`) predictors.

- k:

  As in [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md),
  except that `rbart_vi` fixes the default at `2.0` for both continuous
  and binary responses; unlike `bart2`, a `NULL` default (and thus the
  `chi` hyperprior on binary fits) is not available here.

- object:

  A fitted `rbart` model.

- newdata:

  Same as `test`, but named to match
  [`predict`](https://rdrr.io/r/stats/predict.html) generic.

- type:

  One of `"ev"`, `"ppd"`, `"bart"`, `"loglik"`, `"ranef"`, or `"trees"`
  for the posterior of the expected value, posterior predictive
  distribution, non-parametric/BART component, training log-likelihood
  (`extract` only), random effect, or saved trees respectively. The
  expected value is the sum of the BART component and the random
  effects, while the posterior predictive distribution is a response
  sampled with that mean. `"loglik"` evaluates the log-likelihood of
  each training observation at each posterior draw, conditioning on the
  drawn random intercepts as well as \\\sigma\\ (gaussian) or the fitted
  probability (binary); when chains are combined the result is a
  samples-by-observations matrix directly consumable by WAIC/PSIS-LOO
  implementations such as those in the loo package. To synergize with
  [`predict.glm`](https://rdrr.io/r/stats/predict.glm.html),
  `"response"` can be used as a synonym for `"ev"` and `"link"` can be
  used as a synonym for `"bart"`. For additional details on tree
  extraction, see the corresponding subsection in
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md).

- sample:

  One of `"train"` or `"test"`, referring to the training or tests
  samples respectively.

- ci.level:

  For `fitted` and `predict`, an optional single number in \\(0, 1)\\.
  As in [`bart`](https://vdorie.github.io/dbarts/reference/bart.md):
  `NULL` (the default) returns the posterior mean, while a level returns
  a matrix of `est`, `ci.lower`, and `ci.upper`, with the interval kind
  following `type`.

- x, plquants, cols:

  Same as in
  [`plot.bart`](https://vdorie.github.io/dbarts/reference/bart.md).

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
the posterior of \\tau\\ are made using a slice sampler; the in-engine
sampler for built-in priors steps out with a width tied to the prior's
scale, while the R implementation used with custom priors determines its
width from the curvature of the posterior at its mode.

### Out Of Sample Groups

Predicting random effects for groups not in the training sample is
supported by sampling from their posterior predictive distribution, that
is a draw is taken from \\p(\alpha \mid y) = \int p(\alpha \mid
\tau)p(\tau \mid y)d\alpha\\. For out-of-sample groups in the test data,
these random effect draws can be kept with the saved object. For those
supplied to `predict`, they cannot and may change for subsequent calls.

### Generics

See the generics section of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md).

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
rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                     n.samples = 40L, n.burn = 10L, n.thin = 2L,
                     n.chains = 1L,
                     n.trees = 25L, n.threads = 1L)
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
#>  proposal probabilities: birth/death 0.50, swap 0.10, change 0.40; birth 0.50
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
#> Running mcmc loop:
#> total seconds in loop: 0.000320
#> 
#> Tree sizes, last iteration:
#> [1] 2 2 2 2 2 2 2 3 2 2 3 2 2 3 2 2 3 3 
#> 1 2 2 4 2 3 3 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 4) (2: 3) (3: 3) (4: 5) (5: 3) 
#> (6: 3) (7: 4) (8: 4) (9: 0) (10: 4) 
#> 
#> DONE BART
#> 
#> Running mcmc loop:
#> total seconds in loop: 0.001614
#> 
#> Tree sizes, last iteration:
#> [1] 2 1 2 3 3 2 1 2 2 3 5 2 2 2 3 2 2 2 
#> 3 2 3 4 3 3 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 4) (2: 3) (3: 5) (4: 7) (5: 2) 
#> (6: 1) (7: 4) (8: 6) (9: 3) (10: 1) 
#> 
#> DONE BART
#> 

## with dart = TRUE, split-variable probabilities are sampled under a
## Dirichlet prior, inducing variable selection
rbartFit.dart <- rbart_vi(y ~ . - g, df, group.by = g, dart = TRUE,
                          n.samples = 40L, n.burn = 10L, n.thin = 2L,
                          n.chains = 1L,
                          n.trees = 25L, n.threads = 1L)
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
#>  proposal probabilities: birth/death 0.50, swap 0.10, change 0.40; birth 0.50
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
#> Running mcmc loop:
#> total seconds in loop: 0.000399
#> 
#> Tree sizes, last iteration:
#> [1] 1 2 3 2 2 4 3 2 4 2 3 2 2 2 3 3 1 4 
#> 2 2 3 2 4 2 5 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 9) (2: 3) (3: 11) (4: 6) (5: 4) 
#> (6: 1) (7: 0) (8: 1) (9: 4) (10: 1) 
#> 
#> DONE BART
#> 
#> Running mcmc loop:
#> total seconds in loop: 0.002088
#> 
#> Tree sizes, last iteration:
#> [1] 2 2 3 2 3 3 3 2 3 2 2 2 2 1 2 2 2 2 
#> 3 3 3 3 2 2 5 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 6) (2: 4) (3: 8) (4: 9) (5: 4) 
#> (6: 0) (7: 2) (8: 1) (9: 1) (10: 1) 
#> 
#> DONE BART
#> 
```
