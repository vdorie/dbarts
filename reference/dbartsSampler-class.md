# Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler

A reference class object that contains a Bayesian Additive Regression
Trees sampler in such a way that it can be modified, stopped, and
started all while maintaining its own state.

## Usage

``` r
# S4 method for class 'dbartsSampler'
run(
  numBurnIn, numSamples, updateState = NA, n.threads = control@n.threads
)
# S4 method for class 'dbartsSampler'
sampleTreesFromPrior(updateState = NA)
# S4 method for class 'dbartsSampler'
sampleNodeParametersFromPrior(updateState = NA)
# S4 method for class 'dbartsSampler'
growFromRoot(n.sweeps = 2L, updateState = NA)
# S4 method for class 'dbartsSampler'
copy(shallow = FALSE)
# S4 method for class 'dbartsSampler'
show()
# S4 method for class 'dbartsSampler'
predict(x.test, offset.test, n.threads = control@n.threads)
# S4 method for class 'dbartsSampler'
setControl(control)
# S4 method for class 'dbartsSampler'
setModel(model)
# S4 method for class 'dbartsSampler'
setData(data)
# S4 method for class 'dbartsSampler'
setResponse(y, updateScale = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setOffset(offset, updateScale = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setWeights(weights, updateState = NA)
# S4 method for class 'dbartsSampler'
setSigma(sigma, updateState = NA)
# S4 method for class 'dbartsSampler'
setPredictor(x, column, forceUpdate, updateCutPoints = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setCutPoints(cuts, column, updateState = NA)
# S4 method for class 'dbartsSampler'
setTestPredictor(x.test, column)
# S4 method for class 'dbartsSampler'
setTestPredictorAndOffset(x.test, offset.test)
# S4 method for class 'dbartsSampler'
setTestOffset(offset.test)
# S4 method for class 'dbartsSampler'
printTrees(treeNums, chainNums, sampleNums)
# S4 method for class 'dbartsSampler'
getTrees(
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL
)
# S4 method for class 'dbartsSampler'
getSigmas(result)
# S4 method for class 'dbartsSampler'
getLatents(result)
# S4 method for class 'dbartsSampler'
getSumsOfSquaredResiduals(result)
# S4 method for class 'dbartsSampler'
installTrees(donor, samples = NULL)
# S4 method for class 'dbartsSampler'
setState(newState)
# S4 method for class 'dbartsSampler'
plotTree(
  treeNum, chainNum, sampleNum, treePlotPars = c(
    nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
  ...
)
```

## Note

`dbartsSampler` is a reference class: its methods are not called as free
functions but as `$`-dispatched calls on a sampler instance, e.g.
`sampler$run(100, 100)` or `sampler$setResponse(newY)`. The S4-method
notation in ‘Usage’ below is an artifact of how reference-class methods
are documented and does not reflect the calling syntax; see ‘Examples’.

## Arguments

- numBurnIn:

  A non-negative integer determining how many iterations the sampler
  should skip before storing results. If missing or `NA`, the default is
  filled in from the sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object.

- numSamples:

  A positive integer determining how many posterior samples should be
  returned. If missing or `NA`, the default is also filled in from the
  control object.

- updateState:

  A logical determining if the local cache of the sampler's state should
  be updated after the call completes. Two conventions apply, by method:
  for `run`, `sampleTreesFromPrior`, and
  `sampleNodeParametersFromPrior`, `NA` (the default) fills in the
  sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object's `updateState`, and explicit `TRUE`/`FALSE` override it. For
  the mutators - `setData`, `setResponse`, `setOffset`, `setWeights`,
  `setSigma`, `setPredictor`, and `setCutPoints` - the state is stored
  only on explicit `TRUE`; `NA` (the default) and `FALSE` both store
  nothing, regardless of `control@updateState`. These are typically
  called once per sweep inside a larger Gibbs/MH loop (as
  `dbartsSampler` is designed for), where storing state on every
  mutation would be wasted work whenever the loop only reads `state`
  occasionally (or never); an unforced `state` promise materializes the
  sampler's *current* state on first access regardless, so a
  mutate-then-first-read sequence needs no explicit store. Pass `TRUE`
  explicitly when `state` was already forced (read or saved) earlier and
  a later mutation must be reflected in the next save.

- shallow:

  A logical determining if the copy should retain the underlying data of
  the sampler (`TRUE`) or have its own copies (`FALSE`).

- control:

  An object inheriting from
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).
  When passed to `setControl`, it cannot change `n.trees`, `n.chains`,
  `useQuantiles`, or `rngSeed` from the values the sampler was created
  with, and cannot set `keepTrees = TRUE` without also giving
  `n.samples`; either is an error.

- model:

  An object inheriting from `dbartsModel`. When passed to `setModel`, it
  cannot switch a DART tree prior on or off relative to the sampler's
  creation-time model; recreate the sampler instead.

- data:

  An object inheriting from
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md).

- y:

  A numeric response vector of length equal to that with which the
  sampler was created.

- x:

  A numeric predictor vector of length equal to that with which the
  sampler was created. Can be of a distinct number of rows for
  `setTestPredictor`.

- x.test:

  A new matrix, data frame, or sparse-bearing test set (the column types
  accepted by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
  of the number of columns equal to that in the current model. A
  data-frame/sparse `x.test` is coded against the training levels;
  `setTestPredictor` and `setTestPredictorAndOffset` then install it as
  a resident, whole-object replacement of the test set (see `column`),
  while `predict` materializes it to a numeric matrix before evaluating.

- offset:

  A numeric vector of length equal to that with which the sampler was
  created, or `NULL`. If `offset.test` was set from `offset`, will
  attempt to update that as well.

- updateScale:

  Logical indicating whether BART's internal scale should re-anchor to
  the new offset (`setOffset`) or response (`setResponse`). Defaults to
  `FALSE`, locking the scale set at creation; should only be `TRUE`
  during burn-in, as re-anchoring mid-run makes the fits across
  iterations no longer comparable.

- offset.test:

  A numeric vector of length equal to that of the test matrix, or
  `NULL`. Can be missing for `setTestPredictorAndOffset`.

- n.threads:

  Currently has no effect: `run` and `predict` both execute serially
  regardless of the value passed here. The sampler's own thread count is
  fixed when it is created, from the `n.threads` of its
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object; this argument is reserved for a future per-call override.

- sigma:

  Numeric vector of residual standard deviations, one for each chain.

- weights:

  A numeric vector of non-negative case weights, of length equal to that
  with which the sampler was created. Weights of zero exclude an
  observation from the likelihood while keeping its fitted values.
  `setWeights` only applies to a gaussian-family sampler; calling it on
  a probit- or logistic-family sampler is refused, since a weighted
  probit has no tractable latent-variable form and logistic weights
  (observation counts) are fixed when the sampler is created. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for
  the family-specific weight rules that apply at creation time.

- cuts:

  Vector of cut points for use with `setCutPoints`.

- column:

  An integer or character string vector specifying which column/columns
  of the predictor matrix is to be replaced or for which cuts should be
  used. When updating predictors, it can be missing and the entire
  matrix will be substituted. For `setTestPredictor`, a per-column
  update is refused when the current test set is a data frame or sparse
  container rather than a plain matrix; replace the whole test set
  instead.

- forceUpdate:

  For `setPredictor`, controls how the trees respond to a new predictor
  that would leave a tree with an empty leaf-node. `FALSE` rolls the
  whole update back and leaves the sampler unchanged (rejection
  sampling); `TRUE` forces the change, collapsing any empty nodes; the
  string `"partial"` installs each observation's new value individually
  and rolls back only those observations whose value would empty a leaf,
  returning a per-observation logical of what was installed. `"partial"`
  requires a single `column` and cannot be combined with
  `updateCutPoints`. When missing, defaults to `TRUE` when the whole
  predictor matrix is being replaced (`column` is missing) and `FALSE`
  when a single column is being replaced.

- updateCutPoints:

  For `setPredictor`, a logical; when `TRUE` the cut points (split
  candidate locations) for the replaced column(s) are recomputed from
  the new values, otherwise the existing cut points are kept. Defaults
  to `FALSE` and cannot be combined with `forceUpdate = "partial"`.

- treeNums:

  An integer vector listing the indices of the trees to print or return.

- chainNums:

  An integer vector listing the chains to return from `getTrees`. When
  missing, all chains are returned.

- sampleNums:

  An integer vector listing the saved samples to return from `getTrees`.
  Applies only when `keepTrees` is `TRUE` and `current` is `FALSE`;
  otherwise the live working trees have no sample dimension and it is
  ignored.

- current:

  For `getTrees`, a logical; when `TRUE` the live working trees (the
  current sampler position) are returned even if `keepTrees` is `TRUE`,
  rather than the saved samples.

- newdata:

  For `getTrees`, an optional matrix of predictors (in the form accepted
  by `predict`). When supplied, its observations are routed through each
  tree so the `n` column counts them instead of the training data; the
  tree structure, splits, and leaf values are unchanged.

- result:

  For `getLatents`, an optional pre-allocated numeric vector (or, with
  multiple chains, matrix) that it fills in place and returns instead of
  allocating a new one; its length must equal the number of observations
  times the number of chains. Omit it to let `getLatents` allocate.
  Accepted but not used by `getSigmas` and `getSumsOfSquaredResiduals`.

- treeNum:

  An integer listing the indices of the tree to plot.

- chainNum:

  For `plotTree`, an integer giving the chain to plot from. Required
  when the sampler has more than one chain; defaults to `1` when there
  is only one.

- sampleNum:

  For `plotTree`, an integer giving the saved sample to plot. Applies
  only when `keepTrees` is `TRUE`, in which case it defaults to the most
  recently drawn sample; ignored (with the live working trees plotted
  instead) when `keepTrees` is `FALSE`.

- treePlotPars:

  A named numeric vector containing the quantities `nodeHeight`,
  `nodeWidth`, and `nodeGap`, all of which control aspects of the
  resulting plot.

- donor:

  For `installTrees`, the fit whose forests seed this sampler: another
  `dbartsSampler`, or a `bart` object fit with `keepSampler = TRUE`. The
  donor must share this sampler's predictors, tree count, cut grid, and
  DART setting.

- samples:

  For `installTrees`, an optional integer vector with one entry per
  chain, each a 1-based index into the donor's pool of samples (its
  saved trees when it kept them, else its final trees per chain). `NULL`
  spreads the chains evenly across the pool.

- n.sweeps:

  For `growFromRoot`, a single positive integer giving the number of
  grow-from-root sweeps run in place (default `2L`). Each sweep rebuilds
  every tree in the forest against the current residual, so a small
  handful reaches a good fit.

- newState:

  For `setState`, a state object previously produced by this sampler
  (its `state` field, or the return of `storeState`) over the same
  model. Must inherit from `bartcoreState`.

- ...:

  Extra arguments to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Details

A `dbartsSampler` is a mutable object which contains information
pertaining to fitting a Bayesian additive regression tree model. The
sampler is first created and then, in a separate instruction, run or
modified. In this way, MCMC samplers can be constructed with BART
components filling arbitrary roles.

### Saving

[`save`](https://rdrr.io/r/base/save.html)-ing and
[`load`](https://rdrr.io/r/base/load.html)ing a `dbarts` sampler for
future use requires that R's serialization mechanism be able to access
the state of the sampler which, for memory purposes, is only made
available to R on request. To make it available, call `storeState()` on
the sampler before saving, e.g. for the object `sampler`, execute
`sampler$storeState()`; this captures the sampler's current state into
the serializable `state` field regardless of the `updateState` setting.

The state object is opaque and engine-specific, and carries a format
version and the package version that wrote it. A state loads only within
the format version of the `dbarts` that wrote it; loading it with an
incompatible version refuses cleanly, naming both versions, rather than
risk a silent misread. There is no cross-version migration: re-fit the
model, or restore the state with the `dbarts` release that wrote it.

To restore a saved state into a sampler, call `setState(newState)`: it
validates that `newState` inherits from `bartcoreState`, re-creates the
underlying engine if needed, pushes the state into it, and caches it on
the `state` field. Assigning the field directly
(`sampler$state <- newState`) does *not* restore the sampler - it only
overwrites the R-side cache, leaving the engine untouched, so the next
run continues from the engine's own state rather than the assigned one.
Always route a restore through `setState`.

### Mutation cost

Each accepted predictor mutation does two things: the engine updates its
internal representation, and the R layer collects the accepted values
into the sampler's data object so that replay, saving, and
re-quantization see current data. The collection is by reference for
data-frame-built samplers (only the affected column changes hands) but
requires copying the full predictor matrix when the sampler was built
from a matrix, so tight loops that mutate predictors every iteration are
better served by data-frame input. The remaining per-call overhead is a
few tens of microseconds of R method dispatch and bookkeeping - the
price of R-level consistency. Clients for which that matters can drive
the sampler through the C interface (see the `dbarts.h` header installed
with the package), which invokes the engine directly and performs no
R-side collection; such clients supply their own current predictor
matrix when replaying saved trees.

### Warm starts

`installTrees` seeds the sampler's forests from a `donor` instead of
drawing trees from the prior, for scaling to more chains or embedding a
fit in a larger sampler. Only the donor's trees, `sigma`, and `k`
transfer; each chain keeps its own random-number stream and redraws
everything else, so several chains seeded from one donor stay
overdispersed. A donor with a different tree count, cut grid, or DART
setting is refused rather than silently reshaped. A warm start biases
the early draws toward the donor, so it shortens burn-in rather than
removing it; keep drawing a non-zero number of burn-in samples before
treating the chain as converged.

`growFromRoot` instead builds the sampler's initial forest by
XBART-style root-down stochastic tree construction (He, Yalov and Hahn
2019): each of `n.sweeps` sweeps rebuilds every tree from the root,
sampling each split from the integrated-likelihood weight of every
candidate cut under the tree prior. This reaches a good fit in far fewer
sweeps than the exact sampler, so it is a fast starting point rather
than a posterior sampler - the exact MCMC sweeps own stationarity once
`run` begins, and the posterior is unchanged. It is available for the
constant-leaf model only; calling it on a `linear` or `gp` node prior is
an error, not a silent fall-back, so initialize those forests with
`sampleTreesFromPrior` instead. As with `installTrees`, the grown forest
biases the early draws toward its fit, so shorten burn-in rather than
skipping it. Each chain grows on its own random-number stream, so the
result is independent of `n.threads`. A composable cross-sampler
workflow follows for free: `donor$growFromRoot(k)` then
`target$installTrees(donor)`.

## Value

For `run`, a named-list with contents `sigma`, `train`, `test`, and
`varcount` (plus `k`, `varprobs`, `tau`, and `ranef` when applicable).
`train` is an array of dimension n.obs x n.samples x n.chains, and
likewise `test` (or `NULL` if the sampler has no test data) and
`varcount` are n.predictors x n.samples x n.chains; `sigma` is n.samples
x n.chains. When `n.chains` is `1` the trailing chain dimension is
dropped, so `train` is a plain n.obs x n.samples matrix and `sigma` a
plain vector of length n.samples. A run can be interrupted with
`Ctrl-C`: it stops between iterations - joining any worker threads
first - and signals an error, returning no samples from the interrupted
run. The sampler's chains are left at the iteration they reached, which
is a valid state to run again from.

For `setPredictor`, `TRUE`/`FALSE` depending on whether or not the
operation was successful. The operation can fail if the new predictor
results in a tree with an empty leaf-node. If only single columns were
replaced, the update is rolled back so that the sampler remains in a
valid state. When `forceUpdate` is `"partial"`, instead returns a
logical vector of length equal to the number of observations, `TRUE`
where that observation's new value was installed and `FALSE` where it
was rolled back to its previous value to keep every tree valid.

`predict` keeps the current test matrix in place and uses the current
set of tree splits. This function has two use cases. The first is when
`keepTrees` of
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
is `TRUE`, in which case the sampler should be run to completion and the
function can be used to interrogate the existing fit. When `keepTrees`
is `FALSE`, the function can be used to obtain the likelihood as part of
a proposed new set of covariates in a Metropolis-Hastings step in a
full-Bayes sampler. This would typically be followed by a call to
`setPredictor` if the step is accepted.

For `getTrees`, a `data.frame` with one row per tree node in
depth-first, left-hand-side pre-order, with columns `chain`, `sample`
(present only for saved samples), `tree`, `n` (the number of
observations in the node), `var` (the splitting variable, or -1 at a
leaf), and `value` (the split value, or the leaf prediction). An ordinal
rule's value is its cut point and observations with values less than or
equal to it go left; a categorical rule carries no data value (its
`value` is `NA`) and its split is reported in the `directions` column
instead. When the sampler has any categorical predictors the result
gains a `directions` column decoding each categorical rule into one
`"L"`/`"R"` character per level, in level order (level `k` goes right
when its character is `"R"`); ordinal rules and leaves are `NA`. When
any predictor contains missing values the result gains a `missing`
column giving the branch (`"L"`/`"R"`) each rule sends missing values
down; rules on complete columns and leaves are `NA`. Under a `linear`
node prior each leaf's `value` is its intercept and the result gains one
`beta.<column>` column per designated covariate holding that leaf's
slope on the internal standardized scale; internal nodes are `NA`.

For `getSigmas`, a numeric vector of length equal to the number of
chains, giving each chain's current residual standard deviation on the
original response scale.

For `getSumsOfSquaredResiduals`, a numeric vector of length equal to the
number of chains, giving each chain's residual sum of squares \\\sum
(y - \hat{y})^2\\ on the original response scale; a binary-response
sampler reports on the latent scale instead.

For `getLatents`, `NULL` when the model has no latent-variable
representation (e.g. a gaussian response); otherwise the sampler's
current latent values - a plain vector of length equal to the number of
observations when there is a single chain, or an observations-by-chains
matrix otherwise - written into `result` when one was supplied.

## See also

[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
for applying a shared per-observation predictor update across several
samplers at once.

[`samplePriorPredictive`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
for repeated
`sampleTreesFromPrior`/`sampleNodeParametersFromPrior`/`predict` draws
on a private sampler, for calibrating priors before fitting.

## Examples

``` r
# the embedded-Gibbs pattern: BART as a conditional model inside a larger
# sampler, alternating its own draws with updates to its response between
# calls to run()
n <- 100
x <- matrix(runif(n * 2), n, 2)
y <- x[, 1] - x[, 2] + rnorm(n, 0, 0.1)

sampler <- dbarts(
  y ~ x,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 0L,
    n.samples = 1L,
    updateState = TRUE
  )
)

## first draw: response as given at creation
samples <- sampler$run()
str(samples$train) # a plain n.obs x n.samples matrix (n.chains == 1)
#>  num [1:100, 1] -0.8207 0.8956 0.2326 -0.3679 0.0244 ...

## an outer Gibbs step changes the response (e.g. after updating some
## other part of a joint model); the sampler picks the new target up on
## the next run() without being recreated
newY <- y + rnorm(n, 0, 0.05)
sampler$setResponse(newY)
samples <- sampler$run()
str(samples$train)
#>  num [1:100, 1] -0.6833 0.6015 0.3523 -0.4685 -0.0563 ...
```
