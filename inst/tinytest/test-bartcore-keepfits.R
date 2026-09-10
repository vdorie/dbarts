# The run bridge's storage opt-out: keepFits = FALSE keeps no per-observation
# channel, handing the engine a per-chain one-draw scratch buffer instead and
# returning a null slot in its place. The channels are dropped from the
# RESULT, not from the sweep, so every sampled quantity is bit-for-bit what
# the keeping run drew. Driven through the .Call directly - the R surface for
# the argument is not wired up at this layer.

runEntry <- function(sampler, numSamples, keepFits) {
  handle <- dbarts:::bartcoreSampler(sampler)
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    handle$ptr,
    0L,
    as.integer(numSamples),
    NULL,
    NULL,
    keepFits
  )
}

n <- 120L
nTest <- 15L
numSamples <- 6L
numChains <- 2L

# ---- heteroscedastic with a test set: four per-observation channels ----
makeVarianceSampler <- function() {
  set.seed(4242, sample.kind = "Rejection")
  x <- matrix(runif(n * 2L), n, 2L)
  y <- 2 * x[, 1L] + ifelse(x[, 1L] < 0.5, 0.3, 1.5) * rnorm(n)
  x.test <- matrix(runif(nTest * 2L), nTest, 2L)
  dbarts::dbarts(
    x,
    y,
    test = x.test,
    variance = dbarts::varianceForest(n.trees = 10L),
    control = dbarts::dbartsControl(
      n.chains = numChains,
      n.threads = numChains,
      n.trees = 10L,
      n.samples = numSamples,
      n.burn = 0L,
      updateState = FALSE,
      seed = 4242L
    )
  )
}

kept <- runEntry(makeVarianceSampler(), numSamples, TRUE)
dropped <- runEntry(makeVarianceSampler(), numSamples, FALSE)

# the slot list keeps its shape either way: same names, same order
expect_equal(names(kept), names(dropped))

expect_equal(dim(kept$train), c(n, numSamples, numChains))
expect_equal(dim(kept$test), c(nTest, numSamples, numChains))
expect_equal(dim(kept$variance), c(n, numSamples, numChains))
expect_equal(dim(kept$varianceTest), c(nTest, numSamples, numChains))

expect_null(dropped$train)
expect_null(dropped$test)
expect_null(dropped$variance)
expect_null(dropped$varianceTest)

# the scratch run drew exactly what the keeping run drew: the opt-out moves
# the writes, not the sampler
expect_identical(kept$sigma, dropped$sigma)
expect_identical(kept$varcount, dropped$varcount)

# ---- a multi-forest coupling: the per-forest channel drops, its glue stays ----
makeForestSampler <- function() {
  set.seed(5150, sample.kind = "Rejection")
  x <- matrix(runif(n * 3L), n, 3L)
  z <- rbinom(n, 1L, 0.5)
  y <- 2 * sin(pi * x[, 1L]) + z * (1 + x[, 2L]) + rnorm(n, sd = 0.2)
  dbarts::dbarts(
    x,
    y,
    forests = list(
      dbarts::forest(),
      dbarts::forest(basis = ~ factor(z))
    ),
    control = dbarts::dbartsControl(
      n.chains = numChains,
      n.threads = numChains,
      n.trees = 10L,
      n.samples = numSamples,
      n.burn = 0L,
      updateState = FALSE,
      seed = 5150L
    )
  )
}

keptBCF <- runEntry(makeForestSampler(), numSamples, TRUE)
droppedBCF <- runEntry(makeForestSampler(), numSamples, FALSE)

expect_equal(dim(keptBCF$forestFits), c(n, 2L, numSamples, numChains))
expect_null(droppedBCF$forestFits)
expect_null(droppedBCF$train)
# the glue is three doubles a draw, so it is kept either way
expect_false(is.null(droppedBCF$glue))
expect_identical(keptBCF$glue, droppedBCF$glue)
expect_identical(keptBCF$sigma, droppedBCF$sigma)
expect_identical(keptBCF$varcount, droppedBCF$varcount)

# ---- the flag is a real logical, not whatever coerces ----
expect_error(
  runEntry(makeVarianceSampler(), numSamples, NA),
  "'keepFits' must be TRUE or FALSE"
)
expect_error(
  runEntry(makeVarianceSampler(), numSamples, "yes"),
  "'keepFits' must be TRUE or FALSE"
)
