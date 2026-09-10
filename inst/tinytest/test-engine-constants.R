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
  list(workers = partition[["n.workers"]], rows = partition[["n.rows"]],
       test = fits$test)
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
## build, invisible in the answers, so the store is asked directly.
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
    list(
      sparse = .Call(
        dbarts:::C_dbarts_bartcore_columnStorageIsSparse,
        sampler$getPointer()
      ),
      fits = sampler$run(0L, 2L)$train
    )
  }

  atDefault <- storageAt(0.2)
  raised <- storageAt(0.3)
  # a quarter-dense column is densified at the default and kept sparse above it
  expect_true(!any(atDefault$sparse))
  expect_true(all(raised$sparse))
  # and the two layouts answer identically
  expect_identical(atDefault$fits, raised$fits)

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
  fitAt(dbarts::dbartsControl(n.trees = 5L, n.chains = 1L, updateState = FALSE)),
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
expect_error(dbarts::dbartsControl(categoricalExhaustiveCap = 1L), "integer >= 2")
expect_error(dbarts::dbartsControl(categoricalExhaustiveCap = 31L), "at most 30")
expect_error(dbarts::dbartsControl(testFitParallelCutoff = 0L), "positive integer")
expect_error(dbarts::dbartsControl(predictParallelCutoff = 0L), "positive integer")
expect_error(dbarts::dbartsControl(sparseDensityThreshold = 1.5), "in \\[0, 1\\]")

rm(defaults, fitAt, growTrees, testFitWorkers, predictWorkers)
