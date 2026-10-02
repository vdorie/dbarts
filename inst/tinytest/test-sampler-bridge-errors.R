# reachable Rf_error/stop paths in the sampler bridge and its R gate that no
# other test file exercises with expect_error: a test offset with no test
# predictors, and setControl attempts to change a creation-fixed count.

set.seed(1)
n <- 100L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] + rnorm(n, 0, 0.3)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE
)

# a sampler with no test predictors refuses a test offset. The R method gates
# it first with a message about the NULL test matrix ...
sampler <- dbarts(y ~ x, control = control)
expect_true(is.null(sampler$data@x.test))
expect_error(
  sampler$setTestOffset(rep(0.1, n)),
  pattern = "test offset must be as well"
)

# ... and reaching the bridge directly hits its own guard with the same intent
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestOffset,
    sampler$getPointer(),
    rep(0.1, 5L)
  ),
  pattern = "cannot set a test offset without test predictors"
)

# setControl cannot change counts fixed at creation: chains, trees, quantile
# use are all rejected
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 10L,
    updateState = FALSE
  )),
  pattern = "changing 'n.chains'"
)
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE
  )),
  pattern = "changing 'n.trees'"
)
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    useQuantiles = TRUE,
    updateState = FALSE
  )),
  pattern = "changing 'useQuantiles'"
)

# retaining trees needs an explicit sample count
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    keepTrees = TRUE,
    updateState = FALSE
  )),
  pattern = "keepTrees requires 'n.samples'"
)

# an active-row mask is 0/1 membership: the R method refuses a fractional
# value, and the bridge behind it refuses one below 1 as it does one above
expect_error(
  sampler$setActiveRows(c(0.5, rep(1, n - 1L))),
  pattern = "must be all 0 or 1"
)
for (value in c(0.5, -1, 2)) {
  expect_error(
    .Call(
      dbarts:::C_dbarts_bartcore_setActiveRows,
      sampler$getPointer(),
      c(value, rep(1, n - 1L))
    ),
    pattern = "exactly 0 or 1",
    info = format(value)
  )
}
rm(value)

rm(sampler, control, x, y, n)
