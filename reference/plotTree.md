# Plot a Single Tree From a Fitted BART Model

Minimalist visualization of the branching and leaf contents of one tree
in a fitted
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)/[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md)
or [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)
model, or in a
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
A fit-level convenience wrapper around the sampler's `plotTree` method,
so the trees can be plotted without reaching into the sampler stored on
the fit.

## Usage

``` r
plotTree(object, ...)

# S3 method for class 'bart'
plotTree(object, treeNum = 1L, chainNum, sampleNum, ...)

# S3 method for class 'rbart'
plotTree(object, treeNum = 1L, chainNum = 1L, sampleNum, ...)

# S3 method for class 'dbartsSampler'
plotTree(object, ...)
```

## Arguments

- object:

  A fitted model of class `bart` (from
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
  [`bart2`](https://vdorie.github.io/dbarts/reference/bart.md)) or
  `rbart` (from
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)), or
  a
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
  Fits must have been made with the trees kept (`keeptrees`/`keepTrees`
  equal to `TRUE`).

- treeNum:

  An integer, the index of the tree to plot.

- chainNum:

  An integer, the index of the chain to plot from. For a `bart` fit with
  a single chain it may be omitted; for `rbart` it selects the chain and
  defaults to the first.

- sampleNum:

  An integer, the index of the saved sample to plot from. When omitted,
  the last sample is used (or the current working trees when the trees
  were not kept as samples).

- ...:

  Additional arguments, including `treePlotPars` (a named numeric vector
  controlling `nodeHeight`, `nodeWidth`, and `nodeGap`) and further
  arguments passed on to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `NULL`; called for the side effect of drawing.

## Details

The rendering handles constant, linear, and Gaussian-process leaves. For
a multiple-chain fit, `chainNum` selects the chain. Equivalent to
calling the `plotTree` method on the underlying sampler, i.e.
`object$fit$plotTree(...)` for a `bart` fit.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md),
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)

## Examples

``` r
f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

fit <- bart2(x, y, n.samples = 40L, n.burn = 10L, n.trees = 25L,
             n.chains = 1L, keepTrees = TRUE, verbose = FALSE)

plotTree(fit, treeNum = 1L, sampleNum = 40L)
```
