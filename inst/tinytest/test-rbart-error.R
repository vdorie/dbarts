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

# test that rbart fails with invalid group.by
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = NA, n.threads = 1L),
  "'group.by' must be coercible to factor type"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = not_a_symbol, n.threads = 1L),
  "'group.by' not of length equal to that of data"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = testData$g[-1L], n.threads = 1L),
  "group.by' not of length equal to that of data"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = "not a factor", n.threads = 1L),
  "'group.by' not of length equal to that of data"
)

rm(testData)

