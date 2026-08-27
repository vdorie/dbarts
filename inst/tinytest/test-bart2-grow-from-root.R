source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

## warm.start and n.grow.sweeps both request an initialization: refuse
donor <- dbarts::bart2(
  x,
  y,
  n.trees = 8L,
  n.samples = 4L,
  n.burn = 50L,
  n.chains = 1L,
  keepTrees = TRUE,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 1L
)
expect_error(
  dbarts::bart2(
    x,
    y,
    n.trees = 8L,
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE,
    warm.start = donor,
    n.grow.sweeps = 2L
  ),
  "at most one"
)

## a negative sweep count is refused
expect_error(
  dbarts::bart2(
    x,
    y,
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE,
    n.grow.sweeps = -1L
  ),
  "non-negative"
)

## seeded reproducibility of a grow-initialized fit
fitA <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 3L,
  n.grow.sweeps = 2L
)
fitB <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 3L,
  n.grow.sweeps = 2L
)
expect_equal(fitA$yhat.train, fitB$yhat.train)
expect_true(all(is.finite(fitA$yhat.train)))

## the COUNT is forwarded, not just the request: one sweep and two from the
## same seed initialize differently, so their draws cannot coincide
fitOne <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 3L,
  n.grow.sweeps = 1L
)
expect_false(isTRUE(all.equal(fitA$yhat.train, fitOne$yhat.train)))

## the measurable claim: with no burn-in, a grow-from-root start reaches a much
## lower early-iteration training RMSE than the prior-init cold start
earlyRMSEBart2GrowFromRoot <- function(fit) {
  yh <- fit$yhat.train
  perSample <- sqrt(rowMeans(sweep(yh, 2L, y)^2))
  mean(perSample[seq_len(min(10L, length(perSample)))])
}
coldFit <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L
)
growFit <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L,
  n.grow.sweeps = 2L
)
expect_true(
  earlyRMSEBart2GrowFromRoot(growFit) <
    0.9 * earlyRMSEBart2GrowFromRoot(coldFit)
)

## a longer grow-initialized run converges to a sensible fit
convergedFit <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  n.samples = 200L,
  n.burn = 100L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 5L,
  n.grow.sweeps = 2L
)
postMean <- apply(convergedFit$yhat.train, 2L, mean)
expect_true(sqrt(mean((postMean - y)^2)) < sd(y))
