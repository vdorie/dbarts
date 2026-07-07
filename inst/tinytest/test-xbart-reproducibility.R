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

xval.1 <- dbarts::xbart(
  x,
  y,
  method = "k-fold",
  n.reps = 4L,
  n.samples = 20L,
  n.burn = c(10L, 5L, 1L),
  n.test = 5,
  k = k,
  n.threads = 1L,
  seed = 0L
)

xval.2 <- dbarts::xbart(
  x,
  y,
  method = "k-fold",
  n.reps = 4L,
  n.samples = 20L,
  n.burn = c(10L, 5L, 1L),
  n.test = 5,
  k = k,
  n.threads = 1L,
  seed = 0L
)

expect_true(all(!is.na(xval.1)))
expect_equal(dim(xval.1), c(4L, length(k)))
expect_equal(xval.1, xval.2)

xval.3 <- dbarts::xbart(
  x,
  y,
  method = "k-fold",
  n.reps = 4L,
  n.samples = 20L,
  n.burn = c(10L, 5L, 1L),
  n.test = 5,
  k = k,
  n.threads = 2L,
  seed = 0L
)

xval.4 <- dbarts::xbart(
  x,
  y,
  method = "k-fold",
  n.reps = 4L,
  n.samples = 20L,
  n.burn = c(10L, 5L, 1L),
  n.test = 5,
  k = k,
  n.threads = 2L,
  seed = 0L
)

expect_true(all(!is.na(xval.3)))
expect_equal(dim(xval.3), c(4L, length(k)))
expect_equal(xval.3, xval.4)

expect_true(any(xval.1 != xval.3))

rm(xval.4, xval.3, xval.2, xval.1, k, y, x)

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind)

rm(testData)
