# S0 pins for the latent-subset-mask arc: fixed points the active-rows
# channel must not move without disturbing. Extended at S1.

set.seed(20260812L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L)
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  rngSeed = 7L
)

# Pin: $setWeights on a probit sampler refuses post-creation, unconditionally
# on the value - even the all-ones vector that dbarts(..., weights = rep(1, n))
# accepts and normalizes away at creation is refused here
# (refuseBinaryWeightChange, R_interface_bartcore.cpp, reached from
# bartcore_setWeights before its length or value checks run).
sampler.probit <- dbarts::dbarts(
  x,
  y.binary,
  family = "probit",
  control = control
)
expect_error(
  sampler.probit$setWeights(rep(1, n)),
  "probit models do not support case weights"
)
expect_error(
  sampler.probit$setWeights(runif(n, 0.5, 1.5)),
  "probit models do not support case weights"
)

rm(n, x, y.binary, control, sampler.probit)
