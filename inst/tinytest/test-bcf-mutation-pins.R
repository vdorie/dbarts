# The BCF two-forest sampler's mutation and reporting surface, pinned one
# assertion at a time. These pins record the sampler's behavior before its
# creation and mutation surface widen, so that widening cannot silently move
# one without a test noticing; each is a single plain assertion, written so it
# flips cleanly to its opposite the day the behavior it pins changes.

set.seed(3)
n <- 300L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)
host <- dbarts(x, y, control = control)
bc <- dbarts:::bartcoreBCFSampler(host, z, n.trees.treatment = 25L)
dbarts:::bartcoreRun(bc, 20L, 0L)

# --- refuses: a whole-data or whole-model mutation rebuilds forest 0 alone ---
expect_error(dbarts:::bartcoreSetData(bc, host$data), "multi-forest")
expect_error(
  dbarts:::bartcoreSetModel(bc, host$model, host$control, host$data),
  "multi-forest"
)

# --- refuses: the test surface is undefined without a test treatment vector
expect_error(dbarts:::bartcorePredict(bc, x[1:5, , drop = FALSE]), "BCF")
expect_error(
  dbarts:::bartcoreSetTestPredictor(bc, x[1:5, , drop = FALSE]),
  "BCF"
)
expect_error(dbarts:::bartcoreSetTestOffset(bc, rep(0, n)), "BCF")

# --- refuses: a transactional predictor update validates only the primary
# forest, and the per-observation session has no force variant
expect_error(
  dbarts:::bartcoreSetPredictor(bc, x, forceUpdate = FALSE),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreUpdatePredictorPerObservation(bc, x[, 1L], 1L),
  "multi-forest"
)

# --- refuses: updateScale = TRUE would re-anchor the response transform while
# both forests keep leaf calibrations stated against the old one
expect_error(
  dbarts:::bartcoreSetResponse(bc, y, updateScale = TRUE),
  "multi-forest"
)

# --- succeeds: the forced whole-matrix predictor swap refreshes every forest
expect_true(dbarts:::bartcoreSetPredictor(bc, x, forceUpdate = TRUE))
# --- succeeds: the treatment swap is the supported multi-forest data swap
expect_silent(dbarts:::bartcoreSetTreatment(bc, z))
# --- succeeds: the scale-pinned response, offset and weight swaps
expect_silent(dbarts:::bartcoreSetResponse(bc, y, updateScale = FALSE))
expect_silent(dbarts:::bartcoreSetOffset(bc, rep(0, n), updateScale = FALSE))
expect_silent(dbarts:::bartcoreSetWeights(bc, rep(1, n)))

# a run stays sane after the accepted mutations above
result <- dbarts:::bartcoreRun(bc, 0L, 5L)
expect_true(all(is.finite(result$train)))

# --- the driver-loop identity, a pinned fact rather than a refusal or a
# success. Per-forest fits are internal-scale; fit.scale (the stored (min,
# max) of y) carries the affine map back to the reported scale. a*mu + b_z*tau
# under that map reconstructs the recorded train draw. ---
reconstructTrain <- function(bcSampler, zVec, chain = 1L) {
  glue <- dbarts:::bartcoreBCFGlue(bcSampler)
  muFits <- dbarts:::bartcoreForestFits(bcSampler, 0L)[, chain]
  tauFits <- dbarts:::bartcoreForestFits(bcSampler, 1L)[, chain]
  fitScale <- dbarts:::bartcoreStoreState(bcSampler)[[chain]]$fit.scale
  scale <- fitScale[2L] - fitScale[1L]
  shift <- scale * 0.5 + fitScale[1L]
  bz <- ifelse(zVec != 0, glue[3L, chain], glue[2L, chain])
  scale * (glue[1L, chain] * muFits + bz * tauFits) + shift
}

reconResult <- dbarts:::bartcoreRun(bc, 0L, 1L)
reconTrain <- reconstructTrain(bc, z)
expect_equal(reconTrain, reconResult$train[, 1L], tolerance = 1e-10)

# a per-sweep run loop is bitwise identical to one batched run of the same
# length: with control@rngSeed set, each chain's rng is independent of R's
# stream, so the two routes to the same posterior draws agree exactly.
control.loop <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 30L,
  updateState = FALSE,
  rngSeed = 71L
)
numSamples <- 10L

makeLoopSampler <- function() {
  loopHost <- dbarts(x, y, control = control.loop)
  dbarts:::bartcoreBCFSampler(loopHost, z, n.trees.treatment = 15L)
}

bcBatch <- makeLoopSampler()
batched <- dbarts:::bartcoreRun(bcBatch, 0L, numSamples)

# the reconstruction above held for one chain; with n.thin = 1 the live trees
# sit at the last recorded sample, so it holds per chain here too
for (chain in seq_len(control.loop@n.chains)) {
  expect_equal(
    reconstructTrain(bcBatch, z, chain = chain),
    batched$train[, numSamples, chain],
    tolerance = 1e-10
  )
}

bcLoop <- makeLoopSampler()
looped.train <- array(0, dim(batched$train))
looped.sigma <- array(0, dim(batched$sigma))
for (s in seq_len(numSamples)) {
  sweep <- dbarts:::bartcoreRun(bcLoop, 0L, 1L)
  looped.train[, s, ] <- sweep$train[, 1L, ]
  looped.sigma[s, ] <- sweep$sigma[1L, ]
}

expect_identical(looped.train, batched$train)
expect_identical(looped.sigma, batched$sigma)
