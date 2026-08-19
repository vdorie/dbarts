source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

architectureIsArm <- R.version$arch %in% c("aarch64", "arm", "arm64")

# test that basic Friedman example gets same result with various SIMD sets
maxSIMDLevel <- .Call(dbarts:::C_dbarts_getMaxSIMDInstructionSet)

if (maxSIMDLevel == 0L) {
  # no non-scalar instruction set detected: neither fit below asserts
  # anything, so say so explicitly rather than exiting quietly at 0 tests
  exit_file("no non-scalar SIMD instruction set detected at runtime")
}

# non-scalar levels this host supports, in dispatch order: arm64 has only
# NEON, x86 walks SSE2, SSE4_1, AVX, AVX2 as each becomes available
simdLevels <- if (architectureIsArm) {
  if (maxSIMDLevel >= 1L) 1L else integer(0)
} else {
  Filter(function(level) maxSIMDLevel >= level, c(2L, 5L, 7L, 8L))
}

# Fits x/y at every non-scalar level this host supports and pins each fit
# against the scalar (level 0) draw, bitwise: the engine quantizes predictors
# to cut codes, so nothing short of a genuine kernel defect can move a draw.
# Compares three independently-computed channels (yhat.train, sigma,
# varcount), so a kernel defect confined to one has nowhere to hide. offset,
# when supplied, is forwarded to bart()'s binaryOffset - which dbarts()
# threads straight through as a per-observation offset regardless of family -
# so chain.hpp's trainingFits offset-add pass routes addVectorsInPlace's
# output into yhat.train too. That add is the ONLY addVectorsInPlace call the
# constant-leaf backfitting loop this file otherwise exercises can reach:
# addTreeFitsToTotal's constant-leaf arm is a scalar mu[leaf[i]] gather, so
# without a nonzero offset that kernel is called by no fit here (and a
# structurally-zero offset is nulled out upstream, R/spec.R, so it must be
# genuinely nonzero, not just present).
fitAcrossLevels <- function(
  x,
  y,
  ntree,
  n.sims = 10L,
  n.burn = 0L,
  offset = 0.0
) {
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 0L)
  set.seed(99L)
  scalarFit <- dbarts::bart(
    x,
    y,
    ndpost = n.sims,
    nskip = n.burn,
    ntree = ntree,
    binaryOffset = offset,
    verbose = FALSE
  )
  for (level in simdLevels) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, level)
    set.seed(99L)
    levelFit <- dbarts::bart(
      x,
      y,
      ndpost = n.sims,
      nskip = n.burn,
      ntree = ntree,
      binaryOffset = offset,
      verbose = FALSE
    )
    # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
    expect_equal(scalarFit$yhat.train, levelFit$yhat.train, tolerance = 0)
    # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
    expect_equal(scalarFit$sigma, levelFit$sigma, tolerance = 0)
    # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
    expect_equal(scalarFit$varcount, levelFit$varcount, tolerance = 0)
  }
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, maxSIMDLevel)
  invisible(NULL)
}

# continuous response: the backfitting loop's partition kernels, which every
# fit takes
fitAcrossLevels(testData$x, testData$y, ntree = 50L)

rm(testData)

source(
  system.file("common", "multithreadData.R", package = "dbarts"),
  local = TRUE
)

# a slice of the source set: odd length still straddles every kernel's vector
# width, at a fraction of the cost
testData$x <- testData$x[1:5001L, , drop = FALSE]
testData$y <- testData$y[1:5001L]

# test that long data gets same result with various SIMD sets
fitAcrossLevels(testData$x, testData$y, ntree = 10L)

# same data and tree count, with a genuinely nonzero per-observation offset
# (every element positive, strictly increasing) so the trainingFits offset-add
# pass actually calls addVectorsInPlace and its output lands in yhat.train -
# see fitAcrossLevels' doc comment
testOffset <- 0.1 + 0.01 * seq_len(nrow(testData$x))
fitAcrossLevels(testData$x, testData$y, ntree = 10L, offset = testOffset)

# binary response, same odd-length data: the probit latent-update path
# (z <- -1 + 2*y, then the offset removed) additionally takes
# setVectorToConstant, addVectorsInPlaceWithMultiplier and
# subtractVectorsInPlace, none of which the continuous fits above reach
yBinary <- as.numeric(testData$y > median(testData$y))
fitAcrossLevels(testData$x, yBinary, ntree = 10L)

rm(
  maxSIMDLevel,
  architectureIsArm,
  simdLevels,
  fitAcrossLevels,
  testData,
  testOffset,
  yBinary
)
