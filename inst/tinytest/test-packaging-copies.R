# packageBartResults used to hold three full-size copies of a prediction
# channel at once: the engine's own array, the permuted one it returns, and a
# third that matrix()/t() or apply()'s aperm allocated on the way. It now
# holds two. This file pins both halves of that: the two expressions that
# changed return what the old ones did, bit for bit, over every shape
# packaging feeds them, and they really do allocate one array rather than two.
#
# The old expressions are re-declared here rather than compared against
# recorded draws, which would hold only on the build and instruction set they
# were recorded on while packaging is value-neutral on any build.

# The packaging pair exactly as it stood before the copies came out.
oldConvert <- function(samples, n.chains, combineChains) {
  d <- dim(samples)
  if (is.null(d)) {
    samples
  } else if (!combineChains) {
    if (length(d) == 2L) t(samples) else aperm(samples, c(3L, 2L, 1L))
  } else if (length(d) == 2L) {
    if (n.chains <= 1L) t(samples) else as.vector(samples)
  } else {
    res <- t(matrix(samples, d[1L], prod(d[-1L])))
    if (!is.null(dimnames(samples))) {
      colnames(res) <- dimnames(samples)[[1L]]
    }
    res
  }
}
oldChannelMeans <- function(samples) {
  apply(samples, length(dim(samples)), mean)
}

# One prediction channel in the engine's own layout, through both settings of
# combineChains: the returned array and the posterior mean taken from it must
# both be what the old pair produced.
checkChannel <- function(raw, n.chains) {
  for (combine in c(TRUE, FALSE)) {
    old <- oldConvert(raw, n.chains, combine)
    new <- dbarts:::convertSamplesFromDbartsToBart(raw, n.chains, combine)
    expect_identical(new, old)
    expect_identical(dbarts:::channelMeans(new), oldChannelMeans(old))
  }
}

# the shapes packaging feeds it: one chain (2-D) and several (3-D), and a
# channel carrying names on its parameter margin, which only the combined
# branch threads through
set.seed(17L)
checkChannel(array(rnorm(7L * 5L), c(7L, 5L)), 1L)
checkChannel(array(rnorm(7L * 5L * 3L), c(7L, 5L, 3L)), 3L)
named <- array(rnorm(7L * 5L * 3L), c(7L, 5L, 3L))
dimnames(named) <- list(paste0("obs", seq_len(7L)), NULL, NULL)
checkChannel(named, 3L)

# and the same two expressions on real run channels, gaussian and binary,
# at the chain counts the front door reaches
set.seed(31L)
n <- 100L
p <- 4L
x <- matrix(rnorm(n * p), n, p)
x.test <- matrix(rnorm(20L * p), 20L, p)
y <- 2 * x[, 1L] - x[, 2L] + x[, 3L] * x[, 4L] + rnorm(n)
z <- as.integer(y > 0)

for (case in list(
  list(response = y, n.chains = 2L),
  list(response = z, n.chains = 3L)
)) {
  control <- dbarts::dbartsControl(
    n.samples = 12L,
    n.burn = 5L,
    n.chains = case$n.chains,
    n.trees = 10L,
    n.threads = 1L,
    verbose = FALSE,
    updateState = FALSE
  )
  sampler <- dbarts::dbarts(
    response ~ x,
    data.frame(response = case$response, x = I(x)),
    test = data.frame(x = I(x.test)),
    control = control
  )
  sampler$sampleTreesFromPrior(updateState = FALSE)
  samples <- sampler$run(5L, 12L)
  checkChannel(samples$train, case$n.chains)
  checkChannel(samples$test, case$n.chains)
}

# The allocation count itself: the old pair's extra copy is visible in the
# heap high-water mark, with half an array as the margin. The reduction's own
# summation order is why it runs after the permutation - mean() adds in the
# order it is handed, and only there is an observation's draws contiguous.
maxUsedMiB <- function() gc()[2L, "max used"] * 8 / 1048576
raw <- array(rnorm(300000L), c(10000L, 10L, 3L))
arrayMiB <- 8 * length(raw) / 1048576
invisible(gc(reset = TRUE))
oldCombined <- oldConvert(raw, 3L, TRUE)
oldPeak <- maxUsedMiB()
rm(oldCombined)
invisible(gc(reset = TRUE))
newCombined <- dbarts:::convertSamplesFromDbartsToBart(raw, 3L, TRUE)
newPeak <- maxUsedMiB()
expect_true(newPeak < oldPeak - 0.5 * arrayMiB)

# and the reduction, which used to permute the whole channel again
invisible(gc(reset = TRUE))
oldMeans <- oldChannelMeans(newCombined)
oldMeanPeak <- maxUsedMiB()
invisible(gc(reset = TRUE))
newMeans <- dbarts:::channelMeans(newCombined)
newMeanPeak <- maxUsedMiB()
expect_identical(newMeans, oldMeans)
expect_true(newMeanPeak < oldMeanPeak - 0.5 * arrayMiB)

# --- the same reduction at the summary generics and at a partial-dependence
# replay ---
# Those sites called apply() over the kept axes, so apply() is the reference
# here too. Two shapes the packaging channels never reach: a category-widened
# channel, whose observation margin is next-to-last and whose category margin
# is last (two kept axes), and a predict() replay, whose observations are the
# FIRST margin and whose draws are the trailing ones.
applyMeans <- function(x, trailing) {
  d <- dim(x)
  apply(x, seq.int(length(d) - trailing + 1L, length(d)), mean)
}
reduced <- list(
  list(x = matrix(rnorm(40L * 7L), 40L, 7L), trailing = 1L),
  list(x = array(rnorm(20L * 3L * 5L), c(20L, 3L, 5L)), trailing = 1L),
  list(x = array(rnorm(20L * 3L * 5L), c(20L, 3L, 5L)), trailing = 2L),
  list(x = array(rnorm(11L * 4L * 6L * 3L), c(11L, 4L, 6L, 3L)), trailing = 2L),
  list(x = array(rnorm(9L * 1L * 4L), c(9L, 1L, 4L)), trailing = 1L)
)
for (case in reduced) {
  expect_identical(
    dbarts:::channelMeans(case$x, case$trailing),
    applyMeans(case$x, case$trailing)
  )
}
# and with names on the kept margins, which apply() carries through
named <- array(rnorm(6L * 4L * 3L), c(6L, 4L, 3L))
dimnames(named) <- list(NULL, paste0("o", 1:4), c("a", "b", "c"))
expect_identical(dbarts:::channelMeans(named, 1L), applyMeans(named, 1L))
expect_identical(dbarts:::channelMeans(named, 2L), applyMeans(named, 2L))

# the partial-dependence helper, against the expression it replaced
oldDrawMeans <- function(pred, n.chains) {
  if (n.chains > 1L) {
    as.vector(apply(pred, c(2L, 3L), mean))
  } else {
    apply(pred, 2L, mean)
  }
}
predOne <- matrix(rnorm(50L * 8L), 50L, 8L)
predMany <- array(rnorm(50L * 8L * 3L), c(50L, 8L, 3L))
expect_identical(
  dbarts:::pdbart.drawMeans(predOne, 1L),
  oldDrawMeans(predOne, 1L)
)
expect_identical(
  dbarts:::pdbart.drawMeans(predMany, 3L),
  oldDrawMeans(predMany, 3L)
)

# the allocation, on the shape apply() permutes hardest: a three-margin
# channel reduced over its last
wide <- array(rnorm(200000L), c(10000L, 5L, 4L))
wideMiB <- 8 * length(wide) / 1048576
invisible(gc(reset = TRUE))
oldWide <- applyMeans(wide, 1L)
oldWidePeak <- maxUsedMiB()
invisible(gc(reset = TRUE))
newWide <- dbarts:::channelMeans(wide, 1L)
newWidePeak <- maxUsedMiB()
expect_identical(newWide, oldWide)
expect_true(newWidePeak < oldWidePeak - 0.5 * wideMiB)
