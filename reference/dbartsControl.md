# Discrete Bayesian Additive Regression Trees Sampler Control

Convenience function to create a control object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

## Usage

``` r
dbartsControl(
    verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
    keepTrees = FALSE, n.samples = NA_integer_,
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

  Logical to determine if the empirical quantiles of a columns of
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
  used. See details for more information.

- n.burn:

  A non-negative integer determining how many samples, if any, are
  thrown away at the beginning of a run of the sampler.

- n.trees:

  A positive integer giving the number of trees used in the sum-of-trees
  formulation.

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
