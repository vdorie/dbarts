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
bc <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
  control = control
)
bc$run(20L, 0L)

# --- refuses: a whole-data or whole-model mutation on a sampler that carries
# forest amplitudes is refused R-side, before either .Call is reached
expect_error(bc$setData(host$data), "carries forest amplitudes")
expect_error(bc$setModel(host$model), "carries forest amplitudes")

# --- refuses: the test surface is undefined without an off-sample basis
expect_error(
  bc$predict(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)
expect_error(
  bc$setTestPredictor(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)
# the bridge's own "have no off-sample basis" refusal is unreachable through
# this method on a sampler with no test matrix at all: $setTestOffset's own
# precondition (data@x.test is NULL) fires first, R-side, before the .Call
# (test-bcf-r5-surface.R notes the same gap)
expect_error(
  bc$setTestOffset(rep(0, n)),
  "when test matrix is NULL, test offset must be as well"
)

# --- succeeds: a transactional predictor update revalidates every forest and
# installs under the empty-leaf veto, rolling the whole change back and
# reporting FALSE if any tree of either forest would lose a leaf. Replacing
# the design with its own values cannot empty a leaf, so this one installs
expect_true(bc$setPredictor(x, forceUpdate = FALSE))
# --- succeeds: the per-observation session's cell guard caches every forest,
# pruned to the trees the column can move, so a row installs only if it empties
# no leaf of either forest and is declined otherwise. Re-installing the
# column's own values moves nothing, so every row installs
expect_true(all(
  bc$setPredictor(x[, 1L], 1L, forceUpdate = "partial")
))
# ... and a column collapsed onto two values of the existing grid empties
# leaves, so the veto declines the rows that would: the per-row rollback, and a
# run afterwards stays finite, which is what says both forests were re-routed
installed.partial <- bc$setPredictor(
  ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
  1L,
  forceUpdate = "partial"
)
expect_true(any(!installed.partial))
expect_true(all(is.finite(bc$run(0L, 5L)$train)))

# --- refuses: updateScale = TRUE would re-anchor the response transform while
# both forests keep leaf calibrations stated against the old one
expect_error(
  bc$setResponse(y, updateScale = TRUE),
  "carries forest amplitudes"
)

# --- succeeds: the forced whole-matrix predictor swap refreshes every forest.
# forceUpdate = TRUE always installs or throws, so the method suppresses the
# return value rather than reporting the veto outcome a transactional update
# would (bartcoreSamplerSetPredictor's own "if (!forceUpdate) ... else
# invisible(NULL)")
expect_silent(bc$setPredictor(x, forceUpdate = TRUE))
# --- succeeds: the treatment swap is the supported multi-forest data swap
expect_silent(bc$setForestBasis(2L, cbind(1 - z, z)))
# --- succeeds: the scale-pinned response, offset and weight swaps
expect_silent(bc$setResponse(y, updateScale = FALSE))
expect_silent(bc$setOffset(rep(0, n), updateScale = FALSE))
expect_silent(bc$setWeights(rep(1, n)))

# a run stays sane after the accepted mutations above
result <- bc$run(0L, 5L)
expect_true(all(is.finite(result$train)))

# --- the driver-loop identity, a pinned fact rather than a refusal or a
# success. Per-forest fits are internal-scale; fit.scale (the stored (min,
# max) of y) carries the affine map back to the reported scale. a*mu + b_z*tau
# under that map reconstructs the recorded train draw. ---
reconstructTrainMutationPins <- function(bcSampler, zVec, chain = 1L) {
  glue <- bcSampler$getForestAmplitudes()
  muFits <- bcSampler$getForestFits(1L)[, chain]
  tauFits <- bcSampler$getForestFits(2L)[, chain]
  bcSampler$storeState()
  fitScale <- bcSampler$state[[chain]]$fit.scale
  scale <- fitScale[2L] - fitScale[1L]
  shift <- scale * 0.5 + fitScale[1L]
  bz <- ifelse(zVec != 0, glue[3L, chain], glue[2L, chain])
  scale * (glue[1L, chain] * muFits + bz * tauFits) + shift
}

reconResult <- bc$run(0L, 1L)
reconTrain <- reconstructTrainMutationPins(bc, z)
expect_equal(reconTrain, reconResult$train[, 1L], tolerance = 1e-10)

# a per-sweep run loop is bitwise identical to one batched run of the same
# length: with control@seed set, each chain's rng is independent of R's
# stream, so the two routes to the same posterior draws agree exactly.
control.loop <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 30L,
  updateState = FALSE,
  seed = 71L
)
numSamples <- 10L

makeLoopSampler <- function() {
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), n.trees = 15L)),
    control = control.loop
  )
}

bcBatch <- makeLoopSampler()
batched <- bcBatch$run(0L, numSamples)

# the reconstruction above held for one chain; with n.thin = 1 the live trees
# sit at the last recorded sample, so it holds per chain here too
for (chain in seq_len(control.loop@n.chains)) {
  expect_equal(
    reconstructTrainMutationPins(bcBatch, z, chain = chain),
    batched$train[, numSamples, chain],
    tolerance = 1e-10
  )
}

bcLoop <- makeLoopSampler()
looped.train <- array(0, dim(batched$train))
looped.sigma <- array(0, dim(batched$sigma))
for (s in seq_len(numSamples)) {
  sweep <- bcLoop$run(0L, 1L)
  looped.train[, s, ] <- sweep$train[, 1L, ]
  looped.sigma[s, ] <- sweep$sigma[1L, ]
}

expect_identical(looped.train, batched$train)
expect_identical(looped.sigma, batched$sigma)
