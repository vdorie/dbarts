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

# test that works with multiple threads
x <- testData$x
y <- testData$y
g <- factor(testData$g)

set.seed(0)
expect_inherits(
  dbarts::rbart_vi(
    y ~ x, group.by = g,
    n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
    n.trees = 25L, n.threads = 2L, verbose = FALSE
  ),
  "rbart"
)

rm(testData)

