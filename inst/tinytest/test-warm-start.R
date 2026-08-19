source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

## a warm-started forest reproduces the donor sample it was seeded from:
## same split structure, same partition, and the same training fit
donor <- dbarts::bart2(
  x,
  y,
  n.trees = 8L,
  n.samples = 6L,
  n.burn = 100L,
  n.chains = 1L,
  keepTrees = TRUE,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 1L
)

donorTrees <- donor$fit$getTrees(sampleNums = 3L, chainNums = 1L)

dest <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
dest$installTrees(donor, samples = 3L)
destTrees <- dest$getTrees(current = TRUE, chainNums = 1L)

expect_equal(destTrees$var, donorTrees$var)
expect_equal(destTrees$value, donorTrees$value)
expect_equal(destTrees$n, donorTrees$n)
expect_equal(as.vector(dest$predict(x)), donor$yhat.train[3L, ])

## incompatible donors refuse, touching nothing
expect_error(
  dbarts::bart2(
    x,
    y,
    n.trees = 5L,
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE,
    warm.start = donor
  ),
  "shape-compatible"
)
# a donor on a different cut grid is no longer refused: its splits remap onto
# the destination grid (collapsing any the coarser grid starves) and the fit
# runs to a well-formed, finite result
crossGrid <- dbarts::bart2(
  x,
  y,
  n.trees = 8L,
  n.cuts = 20L,
  n.samples = 4L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 5L,
  warm.start = donor
)
expect_true(all(is.finite(as.vector(crossGrid$yhat.train))))
expect_error(
  dbarts::bart2(
    x,
    y,
    n.trees = 8L,
    dart = TRUE,
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE,
    warm.start = donor
  ),
  "DART"
)
expect_error(dest$installTrees(list(1, 2, 3)), "must be a dbarts sampler")
expect_error(dest$installTrees(donor, samples = c(1L, 2L)), "one per chain")

# a raw state donor is read block by block, and a chain declaring no forests at
# all is refused by name: the donor pool reads each chain's first forest
liveSampler <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(n.trees = 8L, n.chains = 1L, n.threads = 1L)
)
invisible(liveSampler$run(5L, 5L))
rawState <- liveSampler$state
expect_silent(dest$installTrees(rawState))
emptyForests <- rawState
emptyForests[[1L]]$forests <- list()
expect_error(
  dest$installTrees(emptyForests),
  "donor chain holds no forests"
)
# the refused donor left the destination running on the state it had
expect_true(all(is.finite(dest$run(0L, 3L)$sigma)))
rm(liveSampler, rawState, emptyForests)

## DART state transfers between two DART fits
donorDart <- dbarts::bart2(
  x,
  y,
  n.trees = 8L,
  dart = TRUE,
  n.samples = 6L,
  n.burn = 100L,
  n.chains = 1L,
  keepTrees = TRUE,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 1L
)
warmDart <- dbarts::bart2(
  x,
  y,
  n.trees = 8L,
  dart = TRUE,
  n.samples = 6L,
  n.burn = 5L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L,
  warm.start = donorDart
)
expect_true(all(is.finite(as.vector(warmDart$yhat.train))))

## the measurable claim: with a converged deep-tree donor and no burn-in, a
## warm start beats a cold start on early-iteration training RMSE
deepDonor <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 200L,
  n.burn = 1000L,
  n.chains = 1L,
  keepTrees = TRUE,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 1L
)
earlyRMSE <- function(fit) {
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
warmFit <- dbarts::bart2(
  x,
  y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L,
  warm.start = deepDonor
)
expect_true(earlyRMSE(warmFit) < 0.8 * earlyRMSE(coldFit))
