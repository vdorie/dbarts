# A C++ exception thrown under a .Call entry reaches R as an ordinary error,
# whether the entry calls the engine from its own frame or from inside an
# unwind-protected scope, and the session and the sampler carry on. The
# exception is an allocation failure injected into the next test-store build,
# the one engine call every test install makes; the injected failure is
# std::bad_alloc, whose message names it on every standard library.

failNextTestBuild <- function(armed = TRUE) {
  invisible(.Call(dbarts:::C_dbarts_bartcore_setTestStoreHooks, armed))
}

set.seed(20261110L)
n <- 40L
nTest <- 10L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y <- x[, 1L] + rnorm(n, sd = 0.1)
xTest <- matrix(
  runif(2L * nTest),
  nTest,
  2L,
  dimnames = list(NULL, c("x1", "x2"))
)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 2L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 7L
)

# a failed install leaves no test rows; the sampler runs, and takes test rows
# on the next try
expectCarriesOn <- function(sampler, install) {
  expect_null(sampler$run(0L, 1L)$test)
  install()
  expect_equal(dim(sampler$run(0L, 1L)$test), c(nTest, 1L))
}

# --- the dense test setters call the engine from the entry's own frame
sampler <- dbarts(x, y, control = control, n.trees = 5L)
failNextTestBuild()
expect_error(sampler$setTestPredictor(xTest), "bad_alloc")
expectCarriesOn(sampler, function() sampler$setTestPredictor(xTest))

sampler <- dbarts(x, y, control = control, n.trees = 5L)
failNextTestBuild()
expect_error(
  sampler$setTestPredictorAndOffset(xTest, numeric(nTest)),
  "bad_alloc"
)
expectCarriesOn(
  sampler,
  function() sampler$setTestPredictorAndOffset(xTest, numeric(nTest))
)

# --- creation with test rows installs them inside an unwind-protected scope
withTest <- dbarts(x, y, test = xTest, control = control, n.trees = 5L)
create <- function(source, family) {
  .Call(
    dbarts:::C_dbarts_bartcore_create,
    source$control,
    source$model,
    source$data,
    family,
    NULL,
    NULL
  )
}
failNextTestBuild()
expect_error(create(withTest, "gaussian"), "bad_alloc")
expect_inherits(create(withTest, "gaussian"), "externalptr")

# and so does the multinomial route, through its own builder
counts <- matrix(0L, n, 3L)
counts[cbind(seq_len(n), 1L + seq_len(n) %% 3L)] <- 1L
multinomial <- dbarts(
  dbartsData(x, counts = counts, test = xTest),
  family = "multinomial",
  control = control,
  n.trees = 5L
)
failNextTestBuild()
expect_error(create(multinomial, "multinomial"), "bad_alloc")
expect_inherits(create(multinomial, "multinomial"), "externalptr")

# --- a sparse test set is parsed into buffers an unwind-protected scope owns,
# and installed from inside it
if (requireNamespace("Matrix", quietly = TRUE)) {
  sparseTest <- Matrix::Matrix(xTest * (xTest > 0.5), sparse = TRUE)
  sampler <- dbarts(x, y, control = control, n.trees = 5L)
  failNextTestBuild()
  expect_error(sampler$setTestPredictor(sparseTest), "bad_alloc")
  expectCarriesOn(sampler, function() sampler$setTestPredictor(sparseTest))

  sampler <- dbarts(x, y, control = control, n.trees = 5L)
  failNextTestBuild()
  expect_error(
    sampler$setTestPredictorAndOffset(sparseTest, numeric(nTest)),
    "bad_alloc"
  )
  expectCarriesOn(
    sampler,
    function() sampler$setTestPredictorAndOffset(sparseTest, numeric(nTest))
  )
}

failNextTestBuild(FALSE)
