source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler settors raise errors
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test  <- data.frame(x = testData$x, z = 1 - testData$z)

control <- dbarts::dbartsControl(
  n.burn = 0L, n.samples = 1L, n.thin = 5L,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE
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
  
rm(n, sampler, control, test, train)

rm(testData)

