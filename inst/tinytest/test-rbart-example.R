source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart runs example
rbartFit <- dbarts::rbart_vi(
  y ~ x, testData, group.by = g,
  n.samples = 40L, n.burn = 10L, n.thin = 2L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L,
  verbose = FALSE
)
expect_equal(
  dim(rbartFit$yhat.train),
  c(2L, 40L %/% 2L, length(testData$y))
)
expect_equal(
  length(rbartFit$yhat.train.mean),
  length(testData$y)
)
expect_equal(
  dim(rbartFit$ranef),
  c(2L, 40L %/% 2L, length(unique(testData$g)))
)
expect_equal(dim(rbartFit$first.tau), c(2L, 10L %/% 2L))
expect_equal(dim(rbartFit$first.sigma), c(2L, 10L %/% 2L))
expect_equal(dim(rbartFit$tau), c(2L, 40L %/% 2L))
expect_equal(dim(rbartFit$sigma), c(2L, 40L %/% 2L))

expect_true(length(unique(rbartFit$ranef)) > 1L)

# check for one chain
rbartFit <- dbarts::rbart_vi(
  y ~ x, testData, group.by = g,
  n.samples = 80L, n.burn = 20L, n.thin = 2L, n.chains = 1L,
  n.trees = 25L, n.threads = 1L,
  verbose = FALSE
)
expect_equal(
  dim(rbartFit$yhat.train), c(80L %/% 2L, length(testData$y))
)
expect_equal(length(rbartFit$yhat.train.mean), length(testData$y))
expect_equal(
  dim(rbartFit$ranef), c(80L %/% 2L, length(unique(testData$g)))
)
expect_equal(length(rbartFit$first.tau), 20L %/% 2L)
expect_equal(length(rbartFit$first.sigma), 20L %/% 2L)
expect_equal(length(rbartFit$tau), 80L %/% 2L)
expect_equal(length(rbartFit$sigma), 80L %/% 2L)

expect_true(length(unique(rbartFit$ranef)) > 1L)

rm(rbartFit)

rm(testData)

