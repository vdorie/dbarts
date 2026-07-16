# The empty-leaf veto in the MH move kernels (logLikelihoodForBranch) must
# out-penalize any valid branch score, or a change/birth proposal that empties
# a leaf is accepted into the chain state and the resulting tree fails the
# occupancy check when its state is restored into a fresh sampler
# (storeState -> setState, the createStoredBARTSampler path). A valid branch's
# log-likelihood is unbounded below (a -0.5 * centeredSumOfSquares /
# residualVariance term), so a fixed finite penalty is beaten once the branch
# score drops past it - reached here with a small fixed residual variance that
# drives sigma tiny and the score large-negative while the trees grow deep.
# Regression for the -1e7 -> -HUGE_VAL fix.

set.seed(1L)
n <- 2000L
p <- 5L
x <- matrix(runif(n * p), n, p)
y <- 10 *
  sin(pi * x[, 1] * x[, 2]) +
  20 * (x[, 3] - 0.5)^2 +
  10 * x[, 4] +
  5 * x[, 5] +
  rnorm(n)

control <- dbarts::dbartsControl(
  keepTrees = TRUE,
  n.samples = 60L,
  n.burn = 0L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = TRUE,
  rngSeed = 1L
)

# fixed(1e-6): a tiny residual variance forces the veto-crossing regime at
# small n by inflating centeredSumOfSquares / residualVariance
sampler <- dbarts::dbarts(y ~ x, control = control, resid.prior = fixed(1e-6))
invisible(sampler$run(100L, 60L))
sampler$storeState()
state <- sampler$state

# a live tree grown deep enough to have low-count leaves: the configuration a
# crossed veto turns into an empty leaf
expect_true(max(state[[1L]]$forests[[1L]]$tree.sizes) >= 30L)

# restoring into a fresh sampler over the same data must succeed: any empty
# leaf carried in the state trips "state is not consistent with this sampler"
restored <- dbarts::dbarts(y ~ x, control = control, resid.prior = fixed(1e-6))
expect_silent(restored$setState(state))

# the restored trees round-trip exactly (a re-store reproduces the forests)
restored$storeState()
expect_equal(restored$state[[1L]]$forests, state[[1L]]$forests)

# and predict from the restored saved trees matches the original's fits
xTest <- x[1:50L, , drop = FALSE]
expect_equal(
  suppressWarnings(sampler$predict(xTest)),
  suppressWarnings(restored$predict(xTest))
)

rm(sampler, restored, state, control, x, y, xTest, n, p)
