# Checks shared by the gp- and linear-leaf node.prior test files: they
# differ only in the formula, the node.prior, and a couple of knob values.
# Every source() of this file must pass local = TRUE for the expect_*
# calls below to reach the run's masked expect_* bindings, and for the
# functions defined here to resolve expect_* (and, for
# checkStateRoundTrip(), statesAgree()) through that same environment at
# call time.

# plotTree does not error under a non-constant node prior, and the leaf
# covariate designation is fixed at creation: swapping to a constant prior
# via setModel() is refused.
checkPlotTreeAndFixedPrior <- function(sampler) {
  pdf(NULL)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_silent(sampler$plotTree(1L, sampleNum = 1L))
  dev.off()

  model.const <- sampler$model
  model.const@node.prior <- dbarts:::normal(2)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(
    sampler$setModel(model.const),
    pattern = "fixed when a sampler is created"
  )
}

# the leaf type's node.prior rides the data-handle views: a full-rows view
# matches the raw-data path bitwise, standardizing with the parent's
# constants; a proper fold serves its held-out rows through the gathered
# covariates. Returns the objects the caller's own rm() teardown names, so
# a call site assigns the result with list2env(..., environment()).
checkDataHandleViews <- function(formula, df, node.prior, n.trees, n, mu) {
  control.view <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    updateState = FALSE
  )
  sampler.view <- dbarts(
    formula,
    df,
    node.prior = node.prior,
    control = control.view
  )
  handle <- dbarts:::bartcoreDataHandle(
    sampler.view$control,
    sampler.view$data,
    sampler.view$model@node.prior@columns
  )
  set.seed(7)
  view <- dbarts:::bartcoreSamplerFromHandle(
    handle,
    sampler.view$control,
    sampler.view$model,
    sampler.view$data,
    trainRows = seq_len(n)
  )
  set.seed(7)
  full <- dbarts:::bartcoreSampler(sampler.view)
  samples.view <- bartcoreRun(view, 40L, 20L)
  samples.full <- bartcoreRun(full, 40L, 20L)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(samples.view$sigma, samples.full$sigma)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(samples.view$train, samples.full$train)

  testRows <- seq(1L, n, by = 4L)
  set.seed(11)
  fold <- dbarts:::bartcoreSamplerFromHandle(
    handle,
    sampler.view$control,
    sampler.view$model,
    sampler.view$data,
    setdiff(seq_len(n), testRows),
    testRows
  )
  samples.fold <- bartcoreRun(fold, 150L, 100L)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_true(all(is.finite(samples.fold$test)))
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_true(cor(rowMeans(samples.fold$test), mu[testRows]) > 0.7)

  list(
    control.view = control.view,
    sampler.view = sampler.view,
    handle = handle,
    view = view,
    full = full,
    samples.view = samples.view,
    samples.full = samples.full,
    testRows = testRows,
    fold = fold,
    samples.fold = samples.fold
  )
}

# state serialization carries the leaf-type-specific per-tree arrays: a
# restored sampler reproduces the model. Requires statesAgree()
# (inst/common/stateContinuation.R) already sourced at the call site.
# Returns the objects the caller's own rm() teardown names, so a call site
# assigns the result with list2env(..., environment()).
checkStateRoundTrip <- function(formula, df, node.prior, n.trees) {
  control.state <- dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = n.trees,
    n.samples = 5L,
    updateState = FALSE
  )
  sampler.state <- dbarts(
    formula,
    df,
    node.prior = node.prior,
    control = control.state
  )
  invisible(sampler.state$run(30L, 2L))
  sampler.state$storeState()
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_true(
    "tree.params" %in% names(sampler.state$state[[1L]]$forests[[1L]])
  )
  sampler.restored <- dbarts(
    formula,
    df,
    node.prior = node.prior,
    control = control.state
  )
  sampler.restored$setState(sampler.state$state)
  sampler.restored$storeState()
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  statesAgree(sampler.restored$state, sampler.state$state)

  list(
    control.state = control.state,
    sampler.state = sampler.state,
    sampler.restored = sampler.restored
  )
}
