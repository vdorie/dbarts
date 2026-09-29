source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler correctly updates R test offsets only when applicable
set.seed(0L)
n <- nrow(testData$X)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 2L,
  updateState = FALSE
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


sampler <- dbarts::dbarts(
  Z ~ X,
  testData,
  testData$X,
  offset = 0.2,
  control = control
)

sampler$setOffset(0.1)
expect_equal(sampler$data@offset, rep_len(0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.1, n))

sampler$setTestOffset(0.2)
expect_equal(sampler$data@offset.test, rep_len(0.2, n))

sampler$setOffset(-0.1)
expect_equal(sampler$data@offset, rep_len(-0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.2, n))


sampler <- dbarts::dbarts(
  Z ~ X,
  testData,
  testData$X[-1, ],
  offset = 0.2,
  control = control
)

expect_equal(sampler$data@offset.test, rep_len(0.2, n - 1))

sampler$setOffset(0.1)
expect_equal(sampler$data@offset, rep_len(0.1, n))
expect_equal(sampler$data@offset.test, rep_len(0.1, n - 1))

sampler$setOffset(rep_len(-0.1, n))
expect_equal(sampler$data@offset, rep_len(-0.1, n))
expect_null(sampler$data@offset.test)

sampler <- dbarts::dbarts(
  Z ~ X,
  testData,
  testData$X,
  offset = 0.2,
  offset.test = -0.1,
  control = control
)
sampler$setOffset(0.3)

expect_equal(sampler$data@offset, rep_len(0.3, n))
expect_equal(sampler$data@offset.test, rep_len(-0.1, n))

rm(sampler, control, n)

rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler updates offsets in C++
control <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 1L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  verbose = FALSE
)
sampler <- dbarts::dbarts(
  Z ~ X,
  testData,
  testData$X[1L:200L, ],
  subset = 1L:200L,
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

# NULL is the one way to say "no test offset"; an NA is refused on every
# sampler path that takes one, and a refused setTestOffset leaves the link to
# the regular offset as it was
source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)
naX <- testData$X
naSampler <- dbarts::dbarts(
  naX,
  testData$Z,
  naX[1:3, , drop = FALSE],
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.samples = 2L,
    n.burn = 0L,
    keepTrees = TRUE,
    verbose = FALSE
  )
)
invisible(naSampler$run())
expect_error(naSampler$predict(naX[1:3, ], NA_real_), "use NULL for no offset")
expect_error(
  naSampler$predict(naX[1L, , drop = FALSE], NA_real_),
  "use NULL for no offset"
)
linkBefore <- naSampler$data@testUsesRegularOffset
expect_error(naSampler$setTestOffset(NA), "use NULL for no offset")
expect_identical(naSampler$data@testUsesRegularOffset, linkBefore)
expect_error(
  naSampler$setTestPredictorAndOffset(naX[1:3, ], c(0, NA, 0)),
  "use NULL for no offset"
)
expect_silent(naSampler$setTestOffset(NULL))
rm(naSampler, naX, linkBefore, testData)
