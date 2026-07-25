# The R-level slice and rejection samplers behind rbart_vi's custom-prior
# loop, pinned directly against known densities. Draws are seeded, so the
# moment checks are deterministic; tolerances leave several standard errors
# of slack.

set.seed(42L)

# log-scale slice sampling of N(2, 1), the sampler finding the mode itself
# (width = NA)
samples <- dbarts:::sliceSample(
  function(x) dnorm(x, 2, 1, log = TRUE),
  start = 0,
  numSamples = 2000L
)
expect_equal(length(samples), 2000L)
expect_equal(mean(samples), 2, tolerance = 0.1)
expect_equal(sd(samples), 1, tolerance = 0.1)

# natural-scale branch with an explicit width
samples <- dbarts:::sliceSample(
  function(x) dnorm(x, -1, 0.5),
  start = -1,
  numSamples = 2000L,
  width = 2,
  log = FALSE
)
expect_equal(mean(samples), -1, tolerance = 0.1)
expect_equal(sd(samples), 0.5, tolerance = 0.2)

# a positive boundary confines the samples
samples <- dbarts:::sliceSample(
  function(x) dgamma(x, 2, 1, log = TRUE),
  start = 1,
  numSamples = 2000L,
  boundary = c(0, Inf)
)
expect_true(all(samples > 0))
expect_equal(mean(samples), 2, tolerance = 0.15)

# rejection sampling of a truncated normal from a uniform envelope, both
# scales
target <- function(x) dnorm(x, log = TRUE)
lconst <- dnorm(0, log = TRUE) - dunif(0, -2, 2, log = TRUE)
samples <- replicate(
  500L,
  dbarts:::rejectionSample(
    target,
    function(x) dunif(x, -2, 2, log = TRUE),
    function() runif(1L, -2, 2),
    lconst,
    c(-2, 2)
  )
)
expect_true(all(abs(samples) < 2))
expect_equal(mean(samples), 0, tolerance = 0.15)

target <- function(x) dnorm(x)
const <- dnorm(0) / dunif(0, -2, 2)
samples <- replicate(
  500L,
  dbarts:::rejectionSample(
    target,
    function(x) dunif(x, -2, 2),
    function() runif(1L, -2, 2),
    const,
    c(-2, 2),
    log = FALSE
  )
)
expect_true(all(abs(samples) < 2))
expect_equal(mean(samples), 0, tolerance = 0.15)

rm(samples, target, lconst, const)

# a draw accepted on the maxIter'th shrink pass is a valid draw: the exhaustion
# check tests acceptance, not the loop counter. A flat density accepts the
# first proposal always, so maxIter = 1 must still return.
draw <- dbarts:::sliceSample(
  function(x) 1,
  start = 0,
  numSamples = 1L,
  width = 0.5,
  maxIter = 1L,
  boundary = c(-1, 1),
  log = FALSE
)
expect_true(is.finite(draw))
expect_true(draw >= -1 && draw <= 1)

# two finite boundaries clamp the stepping-out interval on both sides
samples <- dbarts:::sliceSample(
  function(x) dbeta(x, 2, 3, log = TRUE),
  start = 0.4,
  numSamples = 1000L,
  boundary = c(0, 1)
)
expect_true(all(samples > 0 & samples < 1))
expect_equal(mean(samples), 0.4, tolerance = 0.05)

# a start with negligible density is replaced by a rejection draw off the
# normal approximation at the mode, so the chain still finds the target
samples <- dbarts:::sliceSample(
  function(x) dnorm(x, 0, 1, log = TRUE),
  start = 12,
  numSamples = 1000L
)
expect_equal(mean(samples), 0, tolerance = 0.15)
expect_equal(sd(samples), 1, tolerance = 0.15)

# that recovery needs the mode, which the natural-scale branch does not compute
expect_error(
  dbarts:::sliceSample(
    function(x) dnorm(x),
    start = 5,
    numSamples = 1L,
    width = 1,
    log = FALSE
  ),
  pattern = "not yet implemented"
)

# a target that is non-finite at the start defeats optim and its hand-rolled
# fallback alike; the failure is named rather than propagating an optim error
expect_error(
  dbarts:::sliceSample(
    function(x) if (abs(x) < 1e-12) -Inf else dnorm(x, 3, 1, log = TRUE),
    start = 0,
    numSamples = 10L
  ),
  pattern = "unable to determine initial curvature"
)

# rejection sampling discards proposals outside the boundary rather than
# returning them
samples <- replicate(
  200L,
  dbarts:::rejectionSample(
    function(x) dnorm(x, log = TRUE),
    function(x) dunif(x, -5, 5, log = TRUE),
    function() runif(1L, -5, 5),
    dnorm(0, log = TRUE) - dunif(0, -5, 5, log = TRUE),
    c(-1, 1)
  )
)
expect_true(all(samples > -1 & samples < 1))

# an envelope that never dominates exhausts maxIter, on either scale
expect_error(
  dbarts:::rejectionSample(
    function(x) dnorm(x, 100, 1, log = TRUE),
    function(x) dnorm(x, log = TRUE),
    function() rnorm(1L),
    0,
    c(-Inf, Inf),
    maxIter = 5L
  ),
  pattern = "unable to obtain rejection sample after 5"
)
expect_error(
  dbarts:::rejectionSample(
    function(x) dnorm(x, 100, 1),
    function(x) dnorm(x),
    function() rnorm(1L),
    1,
    c(-Inf, Inf),
    log = FALSE,
    maxIter = 5L
  ),
  pattern = "unable to obtain rejection sample after 5"
)

rm(draw, samples)
