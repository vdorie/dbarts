# Discrete Bayesian Additive Regression Trees Sampler Control

Convenience function to create a control object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

## Usage

``` r
dbartsControl(
    verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
    keepTrees = FALSE, storage = c("double", "single"),
    n.samples = NA_integer_,
    n.cuts = 100L, n.burn = 200L, n.trees = 75L, n.chains = 4L,
    n.threads = dbarts::guessNumCores(), n.thin = 1L, printEvery = 100L,
    printCutoffs = 0L, rngSeed = NA_integer_, updateState = TRUE)
```

## Arguments

- verbose:

  Logical controlling sampler output to console.

- keepTrainingFits:

  Logical controlling whether or not training fits are returned when the
  sampler runs. These are always computed as part of the fitting
  procedure, so disabling will not substantially impact running time.

- useQuantiles:

  Logical to determine if the empirical quantiles of the columns of
  predictors should be used to determine the tree decision rules. If
  `FALSE`, the rules are spaced uniformly throughout the range of
  covariate values.

- keepTrees:

  A logical that determines whether or not trees are cached as they are
  sampled. In all cases, the current state of the sampler is stored as a
  single set of `n.trees`. When `keepTrees` is `TRUE`, a set of
  `n.trees * n.samples` trees are set aside and populated as the sampler
  runs. If the sampler is stopped and restarted, samples proceed from
  the previously stored tree, looping over if necessary.

- storage:

  A character string selecting the precision of the internal running
  residual. The default `"double"` stores it in double precision and
  produces bitwise-identical draws to previous versions. `"single"` opts
  the running residual into single (32-bit) precision, halving the
  memory traffic of the dominant per-sweep gather for a speedup that is
  largest where memory bandwidth binds - very large `n` and multi-chain
  runs (leaf draws and reductions remain in double precision). Currently
  supported only for continuous (gaussian) responses with the default
  constant leaves; other models signal an error. The reduced-precision
  path changes the sampled values slightly and so is opt-in.

- n.samples:

  A non-negative integer giving the default number of samples to return
  each time the sampler is run. Generally specified by
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  instead, and can be overridden on a per-use basis whenever the sampler
  is
  [`run`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

- n.cuts:

  A positive integer or integer vector giving the number of decision
  rules to be used for each given predictor. If of length less than the
  number of predictors, earlier values are recycled. If for any
  predictor more values are specified than are coherent, fewer may be
  used. See the ‘Decision Rules’ section of
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) for how
  the rules themselves are placed.

- n.burn:

  A non-negative integer determining how many samples, if any, are
  thrown away at the beginning of a run of the sampler.

- n.trees:

  A positive integer giving the number of trees used in the sum-of-trees
  formulation. Default 75, dbarts's own historical choice; BayesTree's
  and [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  default is 200.

- n.chains:

  A positive integer detailing the number of independent chains for the
  sampler to use.

- n.threads:

  A positive integer controlling how many threads will be used for
  various internal calculations, as well as the number of chains.
  Internal calculations are highly optimized so that single-threaded
  performance tends to be superior unless the number of observations is
  very large (\>10k), so that it is often not necessary to have the
  number of threads exceed the number of chains.

- n.thin:

  A positive integer determining how many iterations the MCMC chain
  should jump on the decision trees alone before recording a sample.
  Serves to “thin” the samples against serial correlation. `n.samples`
  are returned regardless of the value of `n.thin`.

- printEvery:

  If `verbose` is `TRUE`, every `printEvery` potential samples (after
  thinning) will issue a verbal statement. Must be a positive integer.

- printCutoffs:

  A non-negative integer specifying how many of the decision rules for a
  variable are printed in verbose mode.

- rngSeed:

  Random number generator seed. Every chain runs its own generator; the
  seed drives a dedicated generator that in turn hands each chain its
  own seed, leaving R's stream untouched. Seeded results do not depend
  on the thread count, and a single-chain run with a given seed
  reproduces the first chain of a multi-chain run with the same seed. If
  equal to `NA`, chain generators are seeded from R's stream at
  creation, so [`set.seed`](https://rdrr.io/r/base/Random.html)
  beforehand suffices for reproducibility; sampling itself never
  advances R's stream.

- updateState:

  Logical setting the default behavior for many
  [sampler](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  methods with regards to the immediate updating of the cached state of
  the object. A current, cached state is only useful when
  [saving](https://rdrr.io/r/base/save.html)/[loading](https://rdrr.io/r/base/load.html)
  the sampler.

## Value

An object of class `dbartsControl`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)

## Examples

``` r
## a small, single-chain, reproducible control for embedding a sampler in a
## larger Gibbs loop: one sample per run(), no burn-in, state cached on demand
control <- dbartsControl(n.chains = 1L, n.threads = 1L,
                         n.burn = 0L, n.samples = 1L,
                         n.trees = 25L, rngSeed = 7L,
                         updateState = FALSE)
control
#> An object of class "dbartsControl"
#> Slot "binary":
#> [1] FALSE
#> 
#> Slot "verbose":
#> [1] FALSE
#> 
#> Slot "keepTrainingFits":
#> [1] TRUE
#> 
#> Slot "useQuantiles":
#> [1] FALSE
#> 
#> Slot "keepTrees":
#> [1] FALSE
#> 
#> Slot "storage":
#> [1] "double"
#> 
#> Slot "n.samples":
#> [1] 1
#> 
#> Slot "n.cuts":
#> [1] 100
#> 
#> Slot "n.burn":
#> [1] 0
#> 
#> Slot "n.trees":
#> [1] 25
#> 
#> Slot "n.chains":
#> [1] 1
#> 
#> Slot "n.threads":
#> [1] 1
#> 
#> Slot "n.thin":
#> [1] 1
#> 
#> Slot "printEvery":
#> [1] 100
#> 
#> Slot "printCutoffs":
#> [1] 0
#> 
#> Slot "rngSeed":
#> [1] 7
#> 
#> Slot "updateState":
#> [1] FALSE
#> 
#> Slot "call":
#> `NA`()
#> 

n <- 50L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.2)
sampler <- dbarts(y ~ x, control = control)
samples <- sampler$run()
str(samples$train)
#>  num [1:50, 1] -0.722 -0.251 0.633 -0.722 0.621 ...
```
