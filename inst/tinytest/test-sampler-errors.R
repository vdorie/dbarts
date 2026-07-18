source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler settors raise errors
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test <- data.frame(x = testData$x, z = 1 - testData$z)

control <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 1L,
  n.thin = 5L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  verbose = FALSE
)
sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)

expect_error(
  sampler$setControl("not-a-control"),
  "'control' must inherit from dbartsControl"
)
expect_error(
  sampler$setResponse(numeric(0L)),
  paste0("y must be of length equal to ", nrow(train))
)

expect_error(
  sampler$setOffset(numeric(0L)),
  "length of replacement offset is not equal to number of observations"
)
expect_error(
  sampler$setPredictor(numeric(0L), 1L),
  "length of new x does not match y"
)
expect_error(
  sampler$setPredictor(testData$z, 3L),
  "column '3' is out of range"
)
expect_error(
  sampler$setTestPredictor(numeric(0L), 1L),
  "length of new x does not match old x.test"
)
expect_error(
  sampler$setTestPredictor(numeric(0L)),
  "number of columns in 'test' must be equal to that of 'x'"
)
expect_error(
  sampler$setTestPredictor(testData$z, 3L),
  "column '3' is out of range"
)

n <- length(testData$y)
expect_error(
  sampler$setPredictor(matrix(numeric(n * 3L), n)),
  paste0("dimension of x must be equal to ", ncol(sampler$data@x))
)
expect_error(
  sampler$setPredictor(matrix(numeric((n - 1L) * 2L), n - 1L)),
  paste0("dimension of x must be equal to ", nrow(sampler$data@x))
)
expect_error(
  sampler$setPredictor(matrix(numeric(n * 2L), n), 1L),
  "number of columns of new x does not match length of columns to replace"
)

# setResponse/setWeights/setOffset/setSigma reject the same class of bad
# input dbartsData's construction/validity reject, and do so before any
# sampler state is touched (the rejected value never reaches data@... or
# the engine)

yBad <- testData$y
yBad[1L] <- NA
yBefore <- sampler$data@y
expect_error(sampler$setResponse(yBad), "response contains missing values")
expect_identical(sampler$data@y, yBefore)

weightsBefore <- sampler$data@weights
expect_error(
  sampler$setWeights(rep(1, n - 1L)),
  "'weights' must have length equal to that of 'y'"
)
expect_identical(sampler$data@weights, weightsBefore)

weightsBad <- rep(1, n)
weightsBad[1L] <- NA
expect_error(sampler$setWeights(weightsBad), "'weights' cannot be NA")
expect_identical(sampler$data@weights, weightsBefore)

weightsBad <- rep(1, n)
weightsBad[1L] <- -1
expect_error(
  sampler$setWeights(weightsBad),
  "'weights' must all be non-negative"
)
expect_identical(sampler$data@weights, weightsBefore)

offsetBad <- rep(0, n)
offsetBad[1L] <- NA
offsetBefore <- sampler$data@offset
expect_error(sampler$setOffset(offsetBad), "'offset' contains missing values")
expect_identical(sampler$data@offset, offsetBefore)

sigmasBefore <- sampler$getSigmas()
expect_error(sampler$setSigma(c(1, 2)), "'sigma' must be of length 1")
expect_identical(sampler$getSigmas(), sigmasBefore)
expect_error(sampler$setSigma(NA_real_), "'sigma' must be positive")
expect_identical(sampler$getSigmas(), sigmasBefore)
expect_error(sampler$setSigma(-1), "'sigma' must be positive")
expect_identical(sampler$getSigmas(), sigmasBefore)

rm(
  n,
  sampler,
  control,
  test,
  train,
  yBad,
  yBefore,
  weightsBefore,
  weightsBad,
  offsetBad,
  offsetBefore,
  sigmasBefore
)

rm(testData)
