source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler correctly updates R test offsets only when applicable
set.seed(0L)
n <- nrow(testData$X)
control <- dbarts::dbartsControl(
  n.chains = 1L, n.threads = 2L, updateState = FALSE
)

sampler <- dbarts::dbarts(Z ~ X, testData, testData$X, control = control)

sampler$setOffset(0.2)
expect_equal(sampler$data@offset, rep_len(0.2, n))
expect_null(sampler$data@offset.test)

sampler$setOffset(NULL)
expect_null(sampler$data@offset)
expect_null(sampler$data@offset.test)

sampler$setOffset(runif(n))
expect_null(sampler$data@offset.test)


sampler <- dbarts::dbarts(Z ~ X, testData, offset = 0.2, control = control)

expect_null(sampler$data@offset.test)

sampler$setOffset(-0.1)
expect_equal(sampler$data@offset, rep_len(-0.1, n))
expect_null(sampler$data@offset.test)


sampler <- dbarts::dbarts(Z ~ X, testData, testData$X, offset = 0.2, control = control)

sampler$setOffset(0.1)
expect_equal(sampler$data@offset, rep_len(0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.1, n))

sampler$setTestOffset(0.2)
expect_equal(sampler$data@offset.test, rep_len(0.2, n))

sampler$setOffset(-0.1)
expect_equal(sampler$data@offset, rep_len(-0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.2, n))


sampler <- dbarts::dbarts(Z ~ X, testData, testData$X[-1,], offset = 0.2, control = control)

expect_equal(sampler$data@offset.test, rep_len(0.2, n - 1))

sampler$setOffset(0.1)
expect_equal(sampler$data@offset, rep_len(0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.1, n - 1))

sampler$setOffset(rep_len(-0.1, n))
expect_equal(sampler$data@offset, rep_len(-0.1, n))
expect_null(sampler$data@offset.test)

sampler <- dbarts::dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = -0.1, control = control)
sampler$setOffset(0.3)

expect_equal(sampler$data@offset, rep_len(0.3, n))
expect_equal(sampler$data@offset.test, rep_len(-0.1, n))

rm(sampler, control, n)

rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler updates offsets in C++
control <- dbarts::dbartsControl(
  n.burn = 0L, n.samples = 1L,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE
)
sampler <- dbarts::dbarts(
  Z ~ X, testData, testData$X[1L:200L,], subset = 1L:200L,
  control = control
)

sampler$setOffset(0.5)
set.seed(0L)
samples <- sampler$run(25L, 1L)
expect_equal(as.double(samples$train - samples$test), rep_len(0.5, 200L))

sampler$setTestOffset(0.5)
samples <- sampler$run(0L, 1L)
expect_equal(samples$train, samples$test)

rm(samples, sampler, control, testData)

