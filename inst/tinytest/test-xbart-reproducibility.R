source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that works with fixed seed; the fold splits are drawn via R's
# sample() even with a chain seed set, so the sampling kind is pinned
# against leakage from other test files
oldSampleKind <- RNGkind()[3L]
suppressWarnings(RNGkind(sample.kind = "Rejection"))

x <- testData$x
y <- testData$y

k <- c(4, 8)

runXval <- function(n.threads) {
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 4L,
    n.samples = 20L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    k = k,
    n.threads = n.threads,
    seed = 0L
  )
}

xval.1 <- runXval(1L)
xval.2 <- runXval(1L)

expect_true(all(!is.na(xval.1)))
expect_equal(dim(xval.1), c(4L, length(k)))
expect_equal(xval.1, xval.2)

xval.3 <- runXval(2L)
xval.4 <- runXval(2L)

expect_true(all(!is.na(xval.3)))
expect_equal(dim(xval.3), c(4L, length(k)))
expect_equal(xval.3, xval.4)

expect_true(any(xval.1 != xval.3))

# first-chunk prefix identity: chunk 1's seed is drawn before the chunk
# count is known, so it does not depend on it - replications landing in
# chunk 1 must therefore match bit for bit across thread counts, even
# though replications in later chunks legitimately diverge. A CRAN check
# run refuses more than 2 simultaneous worker processes, so only the
# 1-vs-2-thread split is exercised here.
expect_equal(apply(xval.1 == xval.3, 1L, all), c(TRUE, TRUE, FALSE, FALSE))

rm(xval.4, xval.3, xval.2, xval.1, runXval, k, y, x)

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind)

rm(testData)
