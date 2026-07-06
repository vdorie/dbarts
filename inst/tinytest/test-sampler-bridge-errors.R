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
bc <- dbarts:::bartcoreSampler(sampler)
expect_error(
  dbarts:::bartcoreSetTestOffset(bc, rep(0.1, 5L)),
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
  pattern = "cannot change 'n.chains'"
)
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE
  )),
  pattern = "cannot change 'n.trees'"
)
expect_error(
  sampler$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    useQuantiles = TRUE,
    updateState = FALSE
  )),
  pattern = "cannot change 'useQuantiles'"
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

rm(sampler, bc, control, x, y, n)
