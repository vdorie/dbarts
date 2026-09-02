source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart issues a warning for weights only in training data
n.train <- 80L
data <- as.data.frame(testData)
trainData <- data[seq_len(n.train), ]

trainData$w <- 1 / nrow(trainData)

# check that predict works when we've fit with missing levels
warnings.trainOnlyWeights <- captureWarnings(
  rbartFit <- dbarts::rbart_vi(
    y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10,
    data = trainData,
    test = data[seq.int(n.train + 1L, nrow(data)), ],
    group.by = g,
    group.by.test = g,
    weights = w,
    n.samples = 7L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 2L,
    n.trees = 25L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_equal(length(warnings.trainOnlyWeights), 1L)
expect_match(
  conditionMessage(warnings.trainOnlyWeights[[1L]]),
  "weights specified but not found in test data - ignoring"
)
expect_inherits(warnings.trainOnlyWeights[[1L]], "dbartsIgnoredArgWarning")
expect_inherits(rbartFit, "rbart")

# the aft arm has no weighted latent form in this version, so it refuses
# 'weights' by name rather than fitting and ignoring them
testData$status <- rep_len(c(1L, 0L), length(testData$y))
testData$time <- exp(abs(testData$y) + 1)
aftData <- as.data.frame(testData)
aftData$w <- rep_len(1, nrow(aftData))
expect_error(
  dbarts::rbart_vi(
    cbind(time, status) ~ x.1 + x.2,
    data = aftData,
    group.by = g,
    weights = w,
    family = "aft",
    n.samples = 5L,
    n.burn = 2L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "do not support 'weights'"
)
rm(aftData)

rm(rbartFit, trainData, data, n.train)

rm(testData)
