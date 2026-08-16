# The negative-binomial dispersion as a per-draw channel: run()$dispersion,
# the $getDispersion() mid-sweep read, and the state-block read they replace.
#
# THE ORACLE is that state read: at every sweep the recorded slot, the getter
# and $state[[chain]]$dispersion are one scalar, for every chain and under both
# r-modes. storeState() returns invisible(NULL), so the state comes off the
# field afterwards.
#
# Every dispersion assertion tests !is.null and the shape BEFORE any value,
# because the channel is NULL on every non-nbinom sampler and
# expect_equal(NULL, NULL) passes silently - a bare value comparison would be
# vacuous on exactly the samplers that carry nothing.
#
# The ESTIMATED-r arm is mandatory and carries a NON-VACUITY measurement: under
# a fixed r every draw repeats one number, so a slot written one sweep late, or
# filled with a constant, is invisible there and visible only here. The fixed r
# is deliberately 8, the shipped grid's median, so a slot filled with that
# median stays green on the fixed arm and moves only on the estimated one.

set.seed(2718L)
n <- 90L
p <- 2L
x <- matrix(runif(n * p), n, p)
yCount <- rnbinom(n, size = 6, mu = exp(0.5 + 0.8 * x[, 1L]))
# the shipped capped positive-integer grid r is drawn on
dispersionGrid <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)

samplerControl <- function(n.chains = 1L, ...) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    seed = 77L,
    n.chains = n.chains,
    ...
  )
}

nbinomSampler <- function(n.chains = 1L, dispersion = NA_real_) {
  dbarts(
    x,
    yCount,
    family = "nbinom",
    dispersion = dispersion,
    control = samplerControl(n.chains),
    verbose = FALSE
  )
}

# One sweep, then read all three. storeSample is the last act of a sweep, so
# the slot the run recorded, the mid-sweep getter and the freshly stored state
# all see the r that sweep settled on.
dispersionCell <- function(sampler) {
  r <- sampler$run(0L, 1L)
  sampler$storeState()
  list(
    slot = as.vector(r$dispersion),
    getter = sampler$getDispersion(),
    state = vapply(sampler$state, function(s) s$dispersion, numeric(1L))
  )
}

# The sweep-by-sweep series the oracle is stated over: one row per sweep, one
# column per chain, for each of the three readers.
dispersionSeries <- function(sampler, n.sweeps = 15L) {
  n.chains <- sampler$control@n.chains
  out <- list(
    slot = matrix(NA_real_, n.sweeps, n.chains),
    getter = matrix(NA_real_, n.sweeps, n.chains),
    state = matrix(NA_real_, n.sweeps, n.chains)
  )
  for (s in seq_len(n.sweeps)) {
    cell <- dispersionCell(sampler)
    out$slot[s, ] <- cell$slot
    out$getter[s, ] <- cell$getter
    out$state[s, ] <- cell$state
  }
  out
}

# --- estimated r, one chain: the arm the late-write and constant-fill
# mutations are only visible on

seriesEst1 <- dispersionSeries(nbinomSampler(1L))
expect_true(!is.null(seriesEst1$slot))
expect_true(!is.null(seriesEst1$getter))
expect_equal(dim(seriesEst1$slot), c(15L, 1L))
expect_equal(dim(seriesEst1$getter), c(15L, 1L))
expect_true(all(!is.na(seriesEst1$slot)))
expect_equal(seriesEst1$slot, seriesEst1$state)
expect_equal(seriesEst1$getter, seriesEst1$state)
expect_true(all(seriesEst1$slot %in% dispersionGrid))
# NON-VACUITY: r must actually move, or the oracle above is a tautology
expect_true(length(unique(seriesEst1$slot[, 1L])) >= 2L)

# --- estimated r, two chains: the per-chain slab stride

seriesEst2 <- dispersionSeries(nbinomSampler(2L))
expect_true(!is.null(seriesEst2$slot))
expect_true(!is.null(seriesEst2$getter))
expect_equal(dim(seriesEst2$slot), c(15L, 2L))
expect_equal(dim(seriesEst2$getter), c(15L, 2L))
expect_true(all(!is.na(seriesEst2$slot)))
expect_equal(seriesEst2$slot, seriesEst2$state)
expect_equal(seriesEst2$getter, seriesEst2$state)
expect_true(all(seriesEst2$slot %in% dispersionGrid))
expect_true(length(unique(seriesEst2$slot[, 1L])) >= 2L)
expect_true(length(unique(seriesEst2$slot[, 2L])) >= 2L)

# --- fixed r: the value the sampler was created with, repeated, on both
# chain counts. 8 is on the grid, so this cell also pins that a fixed r is
# never quietly re-drawn onto a neighbouring grid point.

seriesFix1 <- dispersionSeries(nbinomSampler(1L, dispersion = 8), 5L)
expect_true(!is.null(seriesFix1$slot))
expect_equal(dim(seriesFix1$slot), c(5L, 1L))
expect_equal(seriesFix1$slot, seriesFix1$state)
expect_equal(seriesFix1$getter, seriesFix1$state)
expect_true(all(seriesFix1$slot == 8))

seriesFix2 <- dispersionSeries(nbinomSampler(2L, dispersion = 8), 5L)
expect_true(!is.null(seriesFix2$slot))
expect_equal(dim(seriesFix2$slot), c(5L, 2L))
expect_equal(seriesFix2$slot, seriesFix2$state)
expect_equal(seriesFix2$getter, seriesFix2$state)
expect_true(all(seriesFix2$slot == 8))

# --- a MULTI-sample run: the within-run sample stride and the chain-major
# slab, neither of which a run(0, 1) driver exercises. The last recorded draw
# is the one the post-run state carries.

samplerMulti <- nbinomSampler(2L)
rMulti <- samplerMulti$run(5L, 8L)
samplerMulti$storeState()
expect_true(!is.null(rMulti$dispersion))
expect_equal(dim(rMulti$dispersion), c(8L, 2L))
expect_true(all(rMulti$dispersion %in% dispersionGrid))
expect_equal(
  rMulti$dispersion[8L, ],
  vapply(samplerMulti$state, function(s) s$dispersion, numeric(1L))
)
expect_true(length(unique(as.vector(rMulti$dispersion))) >= 2L)
# a single chain drops the trailing margin, as sigma does
samplerMulti1 <- nbinomSampler(1L)
rMulti1 <- samplerMulti1$run(5L, 8L)
expect_true(!is.null(rMulti1$dispersion))
expect_null(dim(rMulti1$dispersion))
expect_equal(length(rMulti1$dispersion), 8L)

# --- the channel is nbinom-only: NULL slot and NULL getter everywhere else

yGauss <- 3 + 2 * x[, 1L] + rnorm(n)
samplerGauss <- dbarts(x, yGauss, control = samplerControl(), verbose = FALSE)
rGauss <- samplerGauss$run(0L, 2L)
expect_true("dispersion" %in% names(rMulti1))
expect_false("dispersion" %in% names(rGauss))
expect_null(rGauss$dispersion)
expect_null(samplerGauss$getDispersion())

yBinary <- as.numeric(yGauss > median(yGauss))
samplerProbit <- dbarts(x, yBinary, control = samplerControl(), verbose = FALSE)
expect_null(samplerProbit$run(0L, 2L)$dispersion)
expect_null(samplerProbit$getDispersion())

# ordinal owns the conditional slot the dispersion slot is inserted next to, so
# its own channel must be untouched and its getter still NULL
yOrdinal <- cut(
  yGauss,
  quantile(yGauss, c(0, 1 / 3, 2 / 3, 1)),
  labels = FALSE,
  include.lowest = TRUE
)
samplerOrdinal <- dbarts(
  x,
  yOrdinal,
  family = "ordinal",
  control = samplerControl(),
  verbose = FALSE
)
rOrdinal <- samplerOrdinal$run(0L, 2L)
expect_null(rOrdinal$dispersion)
expect_null(samplerOrdinal$getDispersion())
expect_true(!is.null(rOrdinal$cutpoints))
expect_equal(dim(rOrdinal$cutpoints), c(2L, 2L))

# --- the slot INDICES downstream of the insertion. Neither model carries a
# dispersion, which is the point: an off-by-one in varianceTrainSlot or
# forestFitsSlot after the insertion moves these cells and nothing else.

samplerVariance <- dbarts(
  x,
  yGauss,
  variance = TRUE,
  n.trees.variance = 10L,
  control = samplerControl(),
  verbose = FALSE
)
rVariance <- samplerVariance$run(0L, 2L)
expect_null(rVariance$dispersion)
expect_equal(
  names(rVariance),
  c(
    "sigma",
    "train",
    "test",
    "varcount",
    "k",
    "varprobs",
    "tau",
    "ranef",
    "variance",
    "varianceTest"
  )
)
expect_true(!is.null(rVariance$variance))
expect_equal(dim(rVariance$variance), c(n, 2L))
expect_true(all(rVariance$variance > 0))

z <- rbinom(n, 1L, 0.5)
yBcf <- 10 + 4 * x[, 1L] + z * (2 + 3 * x[, 2L]) + rnorm(n, sd = 0.5)
samplerBcf <- dbarts(
  x,
  yBcf,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = samplerControl(),
  verbose = FALSE
)
rBcf <- samplerBcf$run(0L, 2L)
expect_null(rBcf$dispersion)
expect_null(samplerBcf$getDispersion())
expect_equal(
  names(rBcf),
  c(
    "sigma",
    "train",
    "test",
    "varcount",
    "k",
    "varprobs",
    "tau",
    "ranef",
    "forestFits",
    "glue"
  )
)
expect_true(!is.null(rBcf$forestFits))
expect_equal(dim(rBcf$forestFits), c(n, 2L, 2L))
expect_true(!is.null(rBcf$glue))
expect_equal(dim(rBcf$glue), c(3L, 2L))

# --- the nbinom run list itself: the slot lands after ranef, named

expect_equal(
  names(rMulti),
  c(
    "sigma",
    "train",
    "test",
    "varcount",
    "k",
    "varprobs",
    "tau",
    "ranef",
    "dispersion"
  )
)

# --- the host shell of a fit whose model lives elsewhere answers from the
# placeholder, so the read is refused there as $getFitsWithoutOffset() is

samplerShell <- nbinomSampler(1L)
samplerShell$hostFor <- "bart2(family = \"multinomial\")"
expect_error(samplerShell$getDispersion(), pattern = "host sampler of a bart2")

# --- end to end: bart2(family = "nbinom") reads the channel rather than
# serializing state per sweep, and its reported draws still pair with the
# latent psi through mu = r exp(psi)

fit <- bart2(
  x,
  yCount,
  family = "nbinom",
  n.samples = 12L,
  n.burn = 8L,
  n.trees = 15L,
  n.chains = 1L,
  verbose = FALSE
)
expect_true(!is.null(fit$dispersion))
expect_equal(length(fit$dispersion), 12L)
expect_true(all(fit$dispersion %in% dispersionGrid))
expect_true(length(unique(fit$dispersion)) >= 2L)
expect_equal(fit$yhat.train, fit$dispersion * exp(fit$latent.train))
