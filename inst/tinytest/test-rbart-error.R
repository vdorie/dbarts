source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

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

expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = g, n.threads = NA_integer_),
  "'n.threads' must be a positive integer, not NA; leave it out"
)

# test that rbart fails with invalid group.by
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = NA, n.threads = 1L),
  "'group.by' must be coercible to factor type"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = not_a_symbol, n.threads = 1L),
  "'group.by' not found"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = testData$g[-1L], n.threads = 1L),
  "group.by' not of length equal to that of data"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = "not a factor", n.threads = 1L),
  "'group.by' not of length equal to that of data"
)

# a NULL power or base is the default
rbartDefault <- function(...) {
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.thin = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 1L,
    ...
  )
}
expect_identical(
  rbartDefault(power = NULL)$yhat.train,
  rbartDefault()$yhat.train
)
expect_identical(
  rbartDefault(base = NULL)$yhat.train,
  rbartDefault()$yhat.train
)
rm(rbartDefault)

rm(testData)
