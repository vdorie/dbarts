source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that k-fold and random subsample are reproducible, roughly similar
x <- testData$x
y <- testData$y

k <- c(4, 8)

set.seed(0L)
xval.kf <- dbarts::xbart(
  x, y, method = "k-fold",
  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
  k = k,
  n.threads = 1L
)

xval.rs <- dbarts::xbart(
  x, y, method = "random subsample", 
  n.reps = 20L, n.samples = 20L, n.burn = c(10L, 5L, 1L),
  k = k, 
  n.threads = 1L
)

res.kf <- apply(xval.kf, 2L, mean)
res.rs <- apply(xval.rs, 2L, mean)
expect_equal(unname(res.kf), c(2.30094811299725, 4.54475202324197))
expect_equal(unname(res.rs), c(2.35131034003841, 4.57299444101639))

expect_true(all(abs(res.rs - res.kf) < .1))

rm(res.rs, res.kf, xval.rs, xval.kf, k, y, x)


# test that works with fixed seed
x <- testData$x
y <- testData$y

k <- c(4, 8)

xval.1 <- dbarts::xbart(
  x, y, method = "k-fold",
  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
  k = k,
  n.threads = 1L, seed = 0L
)

xval.2 <- dbarts::xbart(
  x, y, method = "k-fold",
  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
  k = k,
  n.threads = 1L, seed = 0L
)

expect_true(all(!is.na(xval.1)))
expect_equal(dim(xval.1), c(4L, length(k)))
expect_equal(xval.1, xval.2)

xval.3 <- dbarts::xbart(
  x, y, method = "k-fold",
  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
  k = k,
  n.threads = 2L, seed = 0L
)

xval.4 <- dbarts::xbart(
  x, y, method = "k-fold",
  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
  k = k,
  n.threads = 2L, seed = 0L
)

expect_true(all(!is.na(xval.3)))
expect_equal(dim(xval.3), c(4L, length(k)))
expect_equal(xval.3, xval.4)

expect_true(any(xval.1 != xval.3))

rm(xval.4, xval.3, xval.2, xval.1, k, y, x)

rm(testData)

