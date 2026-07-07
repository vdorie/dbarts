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
setResponse(y, updateState = NA)
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
setTestPredictor(x.test, column, updateState = NA)
# S4 method for class 'dbartsSampler'
setTestPredictorAndOffset(x.test, offset.test, updateState = NA)
# S4 method for class 'dbartsSampler'
setTestOffset(offset.test, updateState = NA)
# S4 method for class 'dbartsSampler'
printTrees(treeNums)
# S4 method for class 'dbartsSampler'
getTrees(
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL
)
# S4 method for class 'dbartsSampler'
plotTree(
  treeNum, treePlotPars = c(
    nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
  ...
)
# S4 method for class 'dbartsSampler'
startThreads(n.threads = control@n.threads)
# S4 method for class 'dbartsSampler'
stopThreads()
```

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
  be updated after the completion of the run. If `NA`, the default is
  also filled in from the control object.

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

  A new matrix of test predictors, of the number of columns equal to
  that in the current model.

- offset:

  A numeric vector of length equal to that with which the sampler was
  created, or `NULL`. If `offset.test` was set from `offset`, will
  attempt to update that as well.

- updateScale:

  Logical indicating whether BART's internal scale should update with
  the new offset. Should only be `TRUE` during burn-in.

- offset.test:

  A numeric vector of length equal to that of the test matrix, or
  `NULL`. Can be missing for `setTestPredictors`.

- n.threads:

  If greater than one, chain predictions will take place in parallel.

- sigma:

  Numeric vector of residual standard deviations, one for each chain.

- weights:

  A numeric vector of non-negative case weights, of length equal to that
  with which the sampler was created. Weights of zero exclude an
  observation from the likelihood while keeping its fitted values. This
  gaussian interpretation is the general case; a logistic-family sampler
  treats weights as observation counts and requires positive integers,
  and a probit-family sampler does not accept weights at all.

- cuts:

  Vector of cut points for use with `setCutPoints`.

- column:

  An integer or character string vector specifying which column/columns
  of the predictor matrix is to be replaced or for which cuts should be
  used. When updating predictors, it can be missing and the entire
  matrix will be substituted.

- forceUpdate:

  For `setPredictor`, controls how the trees respond to a new predictor
  that would leave a tree with an empty leaf-node. `FALSE` rolls the
  whole update back and leaves the sampler unchanged (rejection
  sampling); `TRUE` forces the change, collapsing any empty nodes; the
  string `"partial"` installs each observation's new value individually
  and rolls back only those observations whose value would empty a leaf,
  returning a per-observation logical of what was installed. `"partial"`
  requires a single `column` and cannot be combined with
  `updateCutPoints`. When missing, defaults to `TRUE` if the replacement
  predictor is a matrix and `FALSE` otherwise.

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

- treeNum:

  An integer listing the indices of the tree to plot.

- treePlotPars:

  A named numeric vector containing the quantities `nodeHeight`,
  `nodeWidth`, and `nodeGap`, all of which control aspects of the
  resulting plot.

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
available to R on request. To do this, one must “touch” the sampler's
state object before saving, e.g. for the object `sampler`, execute
`invisible(sampler$state)`. This is in addition to guaranteeing that the
`state` object is not `NULL`, which can be done by setting the sampler's
control to an object with `updateState` as `TRUE` or passing `TRUE` as
the `updateState` argument to any of the sampler's applicable methods.

The state object is opaque and engine-specific, and carries a format
version and the package version that wrote it. A state loads only within
the format version of the `dbarts` that wrote it; loading it with an
incompatible version refuses cleanly, naming both versions, rather than
risk a silent misread. There is no cross-version migration: re-fit the
model, or restore the state with the `dbarts` release that wrote it.

## Value

For `run`, a named-list with contents `sigma`, `train`, `test`, and
`varcount`. A run can be interrupted with `Ctrl-C`: it stops between
iterations - joining any worker threads first - and signals an error,
returning no samples from the interrupted run. The sampler's chains are
left at the iteration they reached, which is a valid state to run again
from.

For `setPredictor`, `TRUE`/`FALSE` depending on whether or not the
operation was successful. The operation can fail if the new predictor
results in a tree with an empty leaf-node. If only single columns were
replaced, on the update is rolled-back so that the sampler remains in a
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
equal to it go left; a categorical rule's value is its direction mask,
in which bit `k - 1` being set sends the variable's `k`th level right,
except that a rule on a factor of more than 53 levels has no
double-exact mask and its value is `NA` (the `directions` column still
carries its decode). When the sampler has any categorical predictors the
result gains a `directions` column decoding each categorical rule into
one `"L"`/`"R"` character per level, in level order; ordinal rules and
leaves are `NA`. When any predictor contains missing values the result
gains a `missing` column giving the branch (`"L"`/`"R"`) each rule sends
missing values down; rules on complete columns and leaves are `NA`.
Under a `linear` node prior each leaf's `value` is its intercept and the
result gains one `beta.<column>` column per designated covariate holding
that leaf's slope on the internal standardized scale; internal nodes are
`NA`.

## See also

[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
for applying a shared per-observation predictor update across several
samplers at once.
