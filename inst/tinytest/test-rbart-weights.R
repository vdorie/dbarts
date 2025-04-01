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

# test that rbart issues a warning for weights only in training data
n.train <- 80L
data <- as.data.frame(testData)
trainData <- data[seq_len(n.train),]

trainData$w <- 1 / nrow(trainData)

# check that predict works when we've fit with missing levels
expect_warning(rbartFit <- dbarts::rbart_vi(
  y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10,
  data = trainData,
  test = data[seq.int(n.train + 1L, nrow(data)),],
  group.by = g,
  group.by.test = g,
  weights = w,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L,
  verbose = FALSE
))
expect_inherits(rbartFit, "rbart")

rm(rbartFit, trainData, data, n.train)

rm(testData)

