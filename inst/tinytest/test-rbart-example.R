source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart runs example
# combineChains = FALSE pinned deliberately: the shape assertions below
# expect the raw uncombined (n.chains x n.samples[ x ...]) stored shapes
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  testData,
  group.by = g,
  n.samples = 40L,
  n.burn = 10L,
  n.thin = 2L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE
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
  y ~ x,
  testData,
  group.by = g,
  n.samples = 80L,
  n.burn = 20L,
  n.thin = 2L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(
  dim(rbartFit$yhat.train),
  c(80L %/% 2L, length(testData$y))
)
expect_equal(length(rbartFit$yhat.train.mean), length(testData$y))
expect_equal(
  dim(rbartFit$ranef),
  c(80L %/% 2L, length(unique(testData$g)))
)
expect_equal(length(rbartFit$first.tau), 20L %/% 2L)
expect_equal(length(rbartFit$first.sigma), 20L %/% 2L)
expect_equal(length(rbartFit$tau), 80L %/% 2L)
expect_equal(length(rbartFit$sigma), 80L %/% 2L)

expect_true(length(unique(rbartFit$ranef)) > 1L)

rm(rbartFit)

rm(testData)
