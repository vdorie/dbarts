# The exactness a zero treatment multiplier buys, measured rather than
# derived. Under the fixed glue (a, b0, b1) = (1, 0, 1) every
# control row's treatment multiplier is exactly zero, and with the treatment
# forest restricted to a single binary moderator every control row sharing a
# moderator value shares every treatment leaf. Their reported treatment fits
# must then be BITWISE identical: a fit is finalized as
# forestResponse - runningResidual + lastLeaf, and an excluded row carries a
# forest response of exactly zero rather than a residual divided by a near-zero
# multiplier, so that cancellation is exact instead of accurate to about seven
# digits.
#
# Scope, stated honestly: treated rows are NOT made exact by this. Their
# multiplier is one, their arithmetic is the ordinary backfit cancellation, and
# they keep a spread at the usual 1e-16 scale.

set.seed(29)
n <- 400L
moderator <- rbinom(n, 1L, 0.5)
x <- cbind(m = as.double(moderator), x2 = runif(n), x3 = runif(n))
z <- rbinom(n, 1L, 0.5)
y <- 2 *
  sin(pi * x[, "x2"]) +
  x[, "x3"] +
  z * (1 + 2 * x[, "m"]) +
  rnorm(n, sd = 0.2)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)
bcSampler <- dbarts(
  x,
  y,
  forests = list(
    forest(update.amplitude = FALSE),
    forest(
      basis = ~ factor(z),
      n.trees = 25L,
      vars = "m",
      update.amplitude = FALSE
    )
  ),
  control = control
)
invisible(bcSampler$run(100L, 25L))

expect_equal(
  as.vector(bcSampler$getForestAmplitudes()),
  c(1, 0, 1)
)

tauFits <- as.vector(bcSampler$getForestFits(2L))
expect_true(all(is.finite(tauFits)))

controlCells <- split(tauFits[z == 0L], moderator[z == 0L])
treatedCells <- split(tauFits[z == 1L], moderator[z == 1L])
expect_equal(length(controlCells), 2L)

# not vacuous: the treatment forest really did split on the moderator, so the
# two control cells hold different values
expect_true(controlCells[[1L]][1L] != controlCells[[2L]][1L])

# the payoff - zero spread, bitwise, inside each control cell
expect_true(all(vapply(
  controlCells,
  function(fits) length(unique(fits)) == 1L,
  logical(1L)
)))

# and the honest limit: treated rows are only accurate, not exact
expect_true(all(vapply(
  treatedCells,
  function(fits) diff(range(fits)) < 1e-12,
  logical(1L)
)))
