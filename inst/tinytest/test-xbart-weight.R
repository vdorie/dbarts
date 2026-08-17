source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that runs correctly with weighted input
x <- testData$x
y <- testData$y
weights <- rep(1, length(y))

xval <- dbarts::xbart(
  x,
  y,
  weights = weights,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  n.reps = 3,
  n.test = 5,
  k = 2,
  n.threads = 2L
)

expect_equal(length(xval), 3L)

rm(weights, y, x, testData)

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# the binary weight policy (R/spec.R's enforceWeightPolicy): xbart reaches
# the same refusals every other fitting function does. A probit has no
# tractable weighted latent-variable form, except weights identically 1,
# treated as absent; a logistic model requires positive integer count
# weights.
X <- testData$X
Z <- testData$Z

expect_error(
  dbarts::xbart(
    X,
    Z,
    weights = rep(2, length(Z)),
    family = "probit",
    n.reps = 1L,
    n.threads = 1L
  ),
  pattern = "probit models do not support weights"
)
expect_silent(dbarts::xbart(
  X,
  Z,
  weights = rep(1, length(Z)),
  family = "probit",
  n.reps = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = c(3L, 2L, 1L)
))

expect_error(
  dbarts::xbart(
    X,
    Z,
    weights = rep(0.5, length(Z)),
    family = "logistic",
    n.reps = 1L,
    n.threads = 1L
  ),
  pattern = "logistic weights are observation counts"
)
expect_silent(dbarts::xbart(
  X,
  Z,
  weights = sample(1:3, length(Z), replace = TRUE),
  family = "logistic",
  n.reps = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = c(3L, 2L, 1L)
))

rm(X, Z, testData)
