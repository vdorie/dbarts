# The R-level slice and rejection samplers behind rbart_vi's custom-prior
# loop, pinned directly against known densities. Draws are seeded, so the
# moment checks are deterministic; tolerances leave several standard errors
# of slack.

set.seed(42L)

# log-scale slice sampling of N(2, 1), the sampler finding the mode itself
# (width = NA)
samples <- dbarts:::sliceSample(function(x) dnorm(x, 2, 1, log = TRUE),
                                start = 0, numSamples = 2000L)
expect_equal(length(samples), 2000L)
expect_equal(mean(samples), 2, tolerance = 0.1)
expect_equal(sd(samples), 1, tolerance = 0.1)

# natural-scale branch with an explicit width
samples <- dbarts:::sliceSample(function(x) dnorm(x, -1, 0.5), start = -1,
                                numSamples = 2000L, width = 2, log = FALSE)
expect_equal(mean(samples), -1, tolerance = 0.1)
expect_equal(sd(samples), 0.5, tolerance = 0.2)

# a positive boundary confines the samples
samples <- dbarts:::sliceSample(function(x) dgamma(x, 2, 1, log = TRUE),
                                start = 1, numSamples = 2000L,
                                boundary = c(0, Inf))
expect_true(all(samples > 0))
expect_equal(mean(samples), 2, tolerance = 0.15)

# rejection sampling of a truncated normal from a uniform envelope, both
# scales
target <- function(x) dnorm(x, log = TRUE)
lconst <- dnorm(0, log = TRUE) - dunif(0, -2, 2, log = TRUE)
samples <- replicate(500L, dbarts:::rejectionSample(
  target, function(x) dunif(x, -2, 2, log = TRUE),
  function() runif(1L, -2, 2), lconst, c(-2, 2)))
expect_true(all(abs(samples) < 2))
expect_equal(mean(samples), 0, tolerance = 0.15)

target <- function(x) dnorm(x)
const <- dnorm(0) / dunif(0, -2, 2)
samples <- replicate(500L, dbarts:::rejectionSample(
  target, function(x) dunif(x, -2, 2),
  function() runif(1L, -2, 2), const, c(-2, 2), log = FALSE))
expect_true(all(abs(samples) < 2))
expect_equal(mean(samples), 0, tolerance = 0.15)

rm(samples, target, lconst, const)
