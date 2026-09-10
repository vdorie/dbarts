# The engine's fixed limits and the control settings that move them: that each
# cap refuses (or routes) where it says it does, that a setting reaches the
# engine as observable behaviour rather than as a stored value, and that every
# default leaves today's draws alone.

set.seed(0)
n <- 100L
x <- matrix(rnorm(2L * n), n)
y <- x[, 1L] + rnorm(n)

## the per-column cut ceiling refuses by name rather than clamping: a caller
## who asks for 100000 cuts must not silently receive 65533
expect_error(
  dbarts(y ~ x, control = dbartsControl(n.cuts = 100000L, n.chains = 1L)),
  "over the cap of 65533"
)
expect_error(
  dbarts(y ~ x, control = dbartsControl(n.cuts = 65534L, n.chains = 1L)),
  "over the cap of 65533"
)

## the last representable request is still accepted, and costs nothing beyond
## the distinct values a column actually has
sampler <- dbarts(
  y ~ x,
  control = dbartsControl(n.cuts = 65533L, n.chains = 1L, updateState = FALSE)
)
expect_true(all(sampler$data@n.cuts == 65533L))
expect_true(all(is.finite(sampler$run(0L, 2L)$train)))

rm(sampler)


## ---------------------------------------------------------------------------
## the four settable limits. Each pair below fixes everything but the setting
## and asks the engine what it did, since three of the four are byte-identical
## either side of their cutoff and so invisible in the fits. That the DEFAULTS
## reproduce today's draws is what the equivalence baselines gate; here each
## setting is also passed at its default and must leave the draws alone.

## -- categoricalExhaustiveCap: proposal law, so the trees themselves move.
## At twelve present levels the default cap emits eleven sorted prefixes and a
## cap of twelve enumerates all 2047 partitions, so grow-from-root draws a
## different forest; at eight present levels both caps enumerate and the
## forests must agree.
growTrees <- function(nLevels, cap) {
  set.seed(20260910L)
  levels <- letters[seq_len(nLevels)]
  g <- factor(sample(levels, 400L, replace = TRUE))
  effect <- seq(-2, 2, length.out = nLevels)
  names(effect) <- levels
  frame <- data.frame(g = g, z = runif(400L))
  response <- unname(effect[g]) + rnorm(400L, 0, 0.3)
  set.seed(7L)
  sampler <- dbarts::dbarts(
    frame,
    response,
    control = dbarts::dbartsControl(
      n.trees = 20L,
      n.chains = 1L,
      updateState = FALSE,
      categoricalExhaustiveCap = cap
    )
  )
  sampler$growFromRoot(2L)
  trees <- sampler$getTrees(current = TRUE)
  trees[, c("var", "value", "directions")]
}

expect_identical(growTrees(12L, 10L), growTrees(12L, 10L))
expect_false(identical(growTrees(12L, 10L), growTrees(12L, 12L)))
# below either cap the two paths are the same code, so nothing moves
expect_identical(growTrees(8L, 10L), growTrees(8L, 12L))

## -- testFitParallelCutoff: the two routing paths are byte-identical, so the
## engine is asked what it resolved. n.threads above n.chains is the documented
## way to feed the pool, and it warns; that warning is not what is under test.
testFitWorkers <- function(cutoff, nTest) {
  set.seed(4L)
  xTrain <- matrix(rnorm(200L * 2L), 200L)
  yTrain <- xTrain[, 1L] + rnorm(200L)
  xTest <- matrix(rnorm(nTest * 2L), nTest)
  control <- dbarts::dbartsControl(
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 4L,
    updateState = FALSE,
    testFitParallelCutoff = cutoff
  )
  sampler <- suppressWarnings(
    dbarts::dbarts(xTrain, yTrain, xTest, control = control)
  )
  fits <- sampler$run(0L, 2L)
  partition <- .Call(dbarts:::C_dbarts_bartcore_lastTestFitPartition)
  list(
    workers = partition[["n.workers"]],
    rows = partition[["n.rows"]],
    test = fits$test
  )
}

atDefault <- testFitWorkers(65536L, 2000L)
lowered <- testFitWorkers(200L, 2000L)
expect_equal(atDefault$rows, 2000L)
# 2000 test rows sit far below the default cutoff and far above the lowered one
expect_equal(atDefault$workers, 1L)
expect_true(lowered$workers > 1L)
# and the routing is a time choice only: the fits are bit-for-bit the same
expect_identical(atDefault$test, lowered$test)

rm(atDefault, lowered)

## -- predictParallelCutoff: same shape, through the replay's own partition
## channel. The default is calibrated near the measured crossover, so a replay
## of a few hundred thousand traversals now fans out where the uncalibrated 1e7
## kept it inline.
set.seed(5L)
n <- 400L
xFit <- matrix(rnorm(n * 2L), n)
yFit <- xFit[, 1L] + rnorm(n)
xNew <- matrix(rnorm(500L * 2L), 500L)

predictWorkers <- function(cutoff, threads) {
  control <- dbarts::dbartsControl(
    n.trees = 25L,
    n.chains = 1L,
    n.threads = threads,
    n.samples = 10L,
    n.burn = 0L,
    keepTrees = TRUE,
    updateState = FALSE,
    predictParallelCutoff = cutoff
  )
  set.seed(6L)
  sampler <- suppressWarnings(
    dbarts::dbarts(xFit, yFit, control = control)
  )
  invisible(sampler$run(0L, 10L))
  values <- sampler$predict(xNew)
  partition <- .Call(dbarts:::C_dbarts_bartcore_lastPredictPartition)
  list(workers = partition[["n.workers"]], values = values)
}

# 500 rows x 25 trees x 10 draws = 125000 traversals: above the calibrated
# default, below the 1e7 estimate it replaces
calibrated <- predictWorkers(50000L, 4L)
uncalibrated <- predictWorkers(10000000L, 4L)
expect_true(calibrated$workers > 1L)
expect_equal(uncalibrated$workers, 1L)

# the partition gives each slab its own output range and reduces nothing across
# workers, so a replay is bitwise identical however it is dealt out; the moved
# default changes which fits are threaded and must not change what they are
expect_identical(calibrated$values, uncalibrated$values)
expect_identical(calibrated$values, predictWorkers(50000L, 1L)$values)

rm(calibrated, uncalibrated)

## -- sparseDensityThreshold: the layout choice a CSC-built column takes at
## build, which moves no proposal and reaches the answers only in the last
## bits, so the store is asked directly.
if (requireNamespace("Matrix", quietly = TRUE)) {
  set.seed(8L)
  nSparse <- 400L
  density <- 0.25
  columns <- lapply(seq_len(3L), function(j) {
    values <- numeric(nSparse)
    values[sample.int(nSparse, round(nSparse * density))] <- runif(
      round(nSparse * density)
    )
    values
  })
  xCsc <- methods::as(
    methods::as(do.call(cbind, columns), "Matrix"),
    "CsparseMatrix"
  )
  ySparse <- as.vector(xCsc[, 1L]) + rnorm(nSparse, 0, 0.2)

  storageAt <- function(threshold) {
    set.seed(9L)
    sampler <- dbarts::dbarts(
      xCsc,
      ySparse,
      control = dbarts::dbartsControl(
        n.trees = 5L,
        n.chains = 1L,
        updateState = FALSE,
        sparseDensityThreshold = threshold
      )
    )
    sparse <- .Call(
      dbarts:::C_dbarts_bartcore_columnStorageIsSparse,
      sampler$getPointer()
    )
    # one sweep: see the tolerance comment below
    samples <- sampler$run(0L, 1L)
    list(sparse = sparse, fits = samples$train, varcount = samples$varcount)
  }

  atDefault <- storageAt(0.2)
  raised <- storageAt(0.3)
  # a quarter-dense column is densified at the default and kept sparse above it
  expect_true(!any(atDefault$sparse))
  expect_true(all(raised$sparse))
  # the layouts propose the same splits: the discrete decisions are bitwise
  # unmoved, which is the claim the threshold is a memory-time trade
  expect_identical(atDefault$varcount, raised$varcount)
  # the fits are not bitwise, and a tolerance is the honest contract. A dense
  # root partition rewrites indices to the identity before it splits
  # (misc_partitionRange) where the rank-bitmap one permutes in place, so a
  # leaf receives the same members in a different order and its sufficient
  # statistic reassociates. One sweep holds that to a single root partition
  # from the identity both index arrays start at; over more sweeps the gap
  # compounds through the residual and no fixed tolerance holds. The tolerance
  # covers reassociation alone - a moved proposal fails varcount above, which
  # no tolerance hides.
  expect_equal(atDefault$fits, raised$fits, tolerance = 1e-14)

  rm(atDefault, raised)
}

## the defaults are the values the engine used before they were settable, and
## passing them explicitly leaves the draws untouched
defaults <- dbarts::dbartsControl()
expect_equal(defaults@categoricalExhaustiveCap, 10L)
expect_equal(defaults@testFitParallelCutoff, 65536L)
expect_equal(defaults@predictParallelCutoff, 50000L)
expect_equal(defaults@sparseDensityThreshold, 0.2)

fitAt <- function(control) {
  set.seed(11L)
  sampler <- dbarts::dbarts(y ~ x, control = control)
  sampler$run(0L, 3L)$train
}
expect_identical(
  fitAt(dbarts::dbartsControl(
    n.trees = 5L,
    n.chains = 1L,
    updateState = FALSE
  )),
  fitAt(dbarts::dbartsControl(
    n.trees = 5L,
    n.chains = 1L,
    updateState = FALSE,
    categoricalExhaustiveCap = 10L,
    testFitParallelCutoff = 65536L,
    predictParallelCutoff = 50000L,
    sparseDensityThreshold = 0.2
  ))
)

## each refuses what it cannot mean
expect_error(
  dbarts::dbartsControl(categoricalExhaustiveCap = 1L),
  "integer >= 2"
)
expect_error(
  dbarts::dbartsControl(categoricalExhaustiveCap = 31L),
  "at most 30"
)
expect_error(
  dbarts::dbartsControl(testFitParallelCutoff = 0L),
  "positive integer"
)
expect_error(
  dbarts::dbartsControl(predictParallelCutoff = 0L),
  "positive integer"
)
expect_error(
  dbarts::dbartsControl(sparseDensityThreshold = 1.5),
  "in \\[0, 1\\]"
)

## and each is fixed once the sampler exists: the engine reads all four at
## creation, so setControl refuses rather than storing a value the engine will
## never see (and that a re-creation from the stored control would then act on)
fixedSampler <- dbarts::dbarts(
  y ~ x,
  control = dbartsControl(n.trees = 5L, n.chains = 1L, updateState = FALSE)
)
for (change in list(
  c("categoricalExhaustiveCap", 12),
  c("testFitParallelCutoff", 1024),
  c("predictParallelCutoff", 1024),
  c("sparseDensityThreshold", 0.5)
)) {
  moved <- fixedSampler$control
  methods::slot(moved, change[1L]) <- if (
    change[1L] == "sparseDensityThreshold"
  ) {
    as.numeric(change[2L])
  } else {
    as.integer(change[2L])
  }
  expect_error(
    fixedSampler$setControl(moved),
    paste0("changing '", change[1L], "' is not available")
  )
}
rm(fixedSampler, moved, change)

rm(defaults, fitAt, growTrees, testFitWorkers, predictWorkers)


## ---------------------------------------------------------------------------
## the Gaussian-process leaf size cap: not a limit a caller can hit by accident
## but a default that silently substitutes a constant leaf, so the substitution
## is counted and a fit that is mostly not a Gaussian process says so.

set.seed(3L)
nGP <- 300L
gpFrame <- data.frame(x1 = runif(nGP), x2 = runif(nGP))
gpY <- sin(3 * gpFrame$x1) + gpFrame$x2 + rnorm(nGP, 0, 0.2)
gpControl <- dbarts::dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE
)

## a cap of 32 on 300 rows over 10 trees leaves most leaves over it: the fit
## warns, naming the share, and carries the counts behind it
degenerate <- dbarts::dbarts(
  gpY ~ x1 + x2,
  gpFrame,
  control = gpControl,
  node.prior = gp("x1", max.leaf.size = 32L)
)
warned <- NULL
degenerateSamples <- withCallingHandlers(
  degenerate$run(0L, 5L),
  dbartsGPFallbackWarning = function(w) {
    warned <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
expect_true(!is.null(warned))
expect_true(grepl("fell back to a constant leaf", warned, fixed = TRUE))
expect_true(grepl("max.leaf.size", warned, fixed = TRUE))
tally <- attr(degenerateSamples, "gp.fallback")
expect_equal(names(tally), c("evaluations", "fallbacks"))
expect_true(tally[["evaluations"]] > 0)
# the share the message reports is the one the counts give
expect_true(
  grepl(
    sprintf("%.1f%%", 100 * tally[["fallbacks"]] / tally[["evaluations"]]),
    warned,
    fixed = TRUE
  )
)
expect_true(tally[["fallbacks"]] / tally[["evaluations"]] > 0.25)

## a cap above every leaf takes the fallback nowhere, so nothing warns
healthy <- dbarts::dbarts(
  gpY ~ x1 + x2,
  gpFrame,
  control = gpControl,
  node.prior = gp("x1", max.leaf.size = 4096L)
)
quiet <- TRUE
healthySamples <- withCallingHandlers(
  healthy$run(0L, 5L),
  dbartsGPFallbackWarning = function(w) {
    quiet <<- FALSE
    invokeRestart("muffleWarning")
  }
)
expect_true(quiet)
expect_equal(attr(healthySamples, "gp.fallback")[["fallbacks"]], 0)
expect_true(attr(healthySamples, "gp.fallback")[["evaluations"]] > 0)

## a leaf model with no size cap has nothing to report
plainSamples <- dbarts::dbarts(
  gpY ~ x1 + x2,
  gpFrame,
  control = gpControl
)$run(0L, 5L)
expect_null(attr(plainSamples, "gp.fallback"))

## and the census rides the packaged fit, so any other threshold can be
## checked by hand
packaged <- suppressWarnings(dbarts::bart(
  gpY ~ x1 + x2,
  gpFrame,
  node.prior = gp("x1", max.leaf.size = 32L),
  n.trees = 10L,
  n.chains = 1L,
  n.samples = 5L,
  n.burn = 0L,
  verbose = FALSE
))
expect_equal(names(packaged$gp.fallback), c("evaluations", "fallbacks"))
expect_true(packaged$gp.fallback[["fallbacks"]] > 0)

## and it warns ONCE per fit, not once per sampler run: the standard front
## door runs burn-in and sampling as two run() calls
gpWarnings <- 0L
invisible(withCallingHandlers(
  dbarts::bart(
    gpY ~ x1 + x2,
    gpFrame,
    node.prior = gp("x1", max.leaf.size = 32L),
    n.trees = 10L,
    n.chains = 1L,
    n.samples = 5L,
    n.burn = 5L,
    verbose = FALSE
  ),
  dbartsGPFallbackWarning = function(w) {
    gpWarnings <<- gpWarnings + 1L
    invokeRestart("muffleWarning")
  }
))
expect_equal(gpWarnings, 1L)

rm(gpWarnings)

rm(
  degenerate,
  degenerateSamples,
  healthy,
  healthySamples,
  plainSamples,
  packaged,
  tally,
  warned,
  quiet,
  gpControl,
  gpFrame,
  gpY
)
