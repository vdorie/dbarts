# packageBartResults used to hold three full-size copies of a prediction
# channel at once: the engine's own array, the permuted one it returns, and a
# third that matrix()/t() or apply()'s aperm allocated on the way. It now
# holds two. This file pins the two things that removal must not cost.
#
# 1. Every element of the packaged fit is what the old pair produced, bit for
#    bit. The comparison is against the pre-change expressions themselves,
#    re-declared below and injected into the shipped packageBartResults
#    through a shadowing environment, rather than against recorded draws:
#    packaging is value-neutral on ANY build, while pinned draws would hold
#    only on the build and instruction set they were recorded on.
# 2. The permutation really does allocate one array rather than two.
#
# The order a mean is summed in is load-bearing, which is why the reduction
# happens after the permutation and not over the engine's own layout: mean()
# adds in the order it is handed, and the returned layout's observation
# margin is the trailing one.

# The packaging pair exactly as it stood before the copies came out.
oldConvert <- function(
  samples,
  n.chains = dim(samples)[length(dim(samples))],
  combineChains = FALSE
) {
  d <- dim(samples)
  if (!combineChains) {
    if (is.null(d)) {
      samples
    } else if (length(d) == 2L) {
      t(samples)
    } else {
      aperm(samples, c(3L, 2L, 1L))
    }
  } else if (is.null(d)) {
    samples
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

# the shipped packager, run with the old pair in scope of its own body; the
# helpers it reaches through resolve there too, so nothing on the packaging
# path is left on the new expressions
oldEnv <- new.env(parent = asNamespace("dbarts"))
oldEnv$convertSamplesFromDbartsToBart <- oldConvert
oldEnv$channelMeans <- oldChannelMeans
for (helper in c("nameVarcount", "shapeMultinomialChannel")) {
  copied <- get(helper, envir = asNamespace("dbarts"))
  environment(copied) <- oldEnv
  assign(helper, copied, envir = oldEnv)
}
oldPackage <- dbarts:::packageBartResults
environment(oldPackage) <- oldEnv

set.seed(31L)
n <- 100L
p <- 4L
x <- matrix(rnorm(n * p), n, p)
x.test <- matrix(rnorm(20L * p), 20L, p)
y <- 2 * x[, 1L] - x[, 2L] + x[, 3L] * x[, 4L] + rnorm(n)
z <- as.integer(y > 0)

runAndPackage <- function(response, n.chains, combineChains) {
  control <- dbarts::dbartsControl(
    n.samples = 12L,
    n.burn = 5L,
    n.chains = n.chains,
    n.trees = 10L,
    n.threads = 1L,
    verbose = FALSE,
    updateState = FALSE
  )
  sampler <- dbarts::dbarts(
    response ~ x,
    data.frame(response = response, x = I(x)),
    test = data.frame(x = I(x.test)),
    control = control
  )
  sampler$sampleTreesFromPrior(updateState = FALSE)
  samples <- sampler$run(5L, 12L)
  list(
    new = dbarts:::packageBartResults(
      sampler,
      samples,
      NULL,
      NULL,
      combineChains,
      FALSE
    ),
    old = oldPackage(sampler, samples, NULL, NULL, combineChains, FALSE)
  )
}

# gaussian and probit at combineChains = TRUE (the folded chain margin, where
# both channels and the posterior means all changed), and both families at
# combineChains = FALSE, which keeps the chain margin - the binary case there
# is the one that already peaked at two copies, so it must come through
# untouched
for (case in list(
  list(response = y, n.chains = 2L, combineChains = TRUE),
  list(response = y, n.chains = 3L, combineChains = FALSE),
  list(response = z, n.chains = 2L, combineChains = TRUE),
  list(response = z, n.chains = 3L, combineChains = FALSE)
)) {
  packaged <- runAndPackage(
    case$response,
    case$n.chains,
    case$combineChains
  )
  expect_identical(packaged$new, packaged$old)
}

# The allocation count itself. A channel this size dwarfs the collector noise
# a fit leaves behind, so the old pair's extra copy is visible in the
# high-water mark; half of one array is the margin.
raw <- array(rnorm(300000L), c(10000L, 10L, 3L))
arrayMb <- 8 * length(raw) / 1048576
invisible(gc(reset = TRUE))
oldCombined <- oldConvert(raw, 3L, TRUE)
oldPeak <- gc()[2L, 6L]
rm(oldCombined)
invisible(gc(reset = TRUE))
newCombined <- dbarts:::convertSamplesFromDbartsToBart(raw, 3L, TRUE)
newPeak <- gc()[2L, 6L]
expect_true(newPeak < oldPeak - 0.5 * arrayMb)

# and the reduction, which used to permute the whole channel again
invisible(gc(reset = TRUE))
oldMeans <- oldChannelMeans(newCombined)
oldMeanPeak <- gc()[2L, 6L]
invisible(gc(reset = TRUE))
newMeans <- dbarts:::channelMeans(newCombined)
newMeanPeak <- gc()[2L, 6L]
expect_identical(newMeans, oldMeans)
expect_true(newMeanPeak < oldMeanPeak - 0.5 * arrayMb)
