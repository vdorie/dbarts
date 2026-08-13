# A BCF run's per-draw per-forest reporting channels (docs/design/bcf.md):
# forestFits carries each forest's own internal-scale function values (mu in
# slot 1, tau in slot 2) and glue the (a, b0, b1) that recombines them into that
# draw's location, so one run reports both surfaces for every sample rather than
# only the live values the accessors read after it. A model whose forests
# compose through no such scalars - a single-forest sampler, a K-forest
# multinomial - carries neither channel.

set.seed(17)
n <- 150L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

reportingControl <- function(...) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 40L,
    updateState = FALSE,
    rngSeed = 47L,
    ...
  )
}

numSamples <- 8L

# --- the layout: the forest axis sits between the observations and the samples,
# and the glue is one (a, b0, b1) column per draw ---
sampler <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = reportingControl(n.chains = 1L)
)
batched <- sampler$run(0L, numSamples)

expect_equal(dim(batched$forestFits), c(n, 2L, numSamples))
expect_equal(dim(batched$glue), c(3L, numSamples))
expect_true(all(is.finite(batched$forestFits)))
expect_true(all(is.finite(batched$glue)))

# fit.scale is the stored (min, max) of y; the reported train draw is the
# internal-scale location under that affine map, which is what makes the
# reconstruction below a statement about the channels rather than the scaling
internalTrain <- function(sampler, train, chain = 1L) {
  fitScale <- sampler$state[[chain]]$fit.scale
  scale <- fitScale[2L] - fitScale[1L]
  (train - (scale * 0.5 + fitScale[1L])) / scale
}

# F8: a * mu + b_z * tau over the recorded channels reproduces the recorded
# internal-scale train draw for EVERY sample, not only the last one the live
# per-forest accessors can still see.
reconstructionError <- function(fits, glue, train, zVec) {
  err <- 0
  for (s in seq_len(ncol(glue))) {
    bz <- ifelse(zVec != 0, glue[3L, s], glue[2L, s])
    recon <- glue[1L, s] * fits[, 1L, s] + bz * fits[, 2L, s]
    err <- max(err, max(abs(recon - train[, s])))
  }
  err
}

sampler$storeState()
expect_true(
  reconstructionError(
    batched$forestFits,
    batched$glue,
    internalTrain(sampler, batched$train),
    z
  ) <
    1e-12
)

# --- the channels report exactly what a per-sweep driver loop reports: one run
# of numSamples sweeps and numSamples runs of one sweep agree to the bit, and
# each recorded slab equals what the live per-forest accessors read at that
# point - the loop a front end would otherwise have to own ---
looped <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = reportingControl(n.chains = 1L)
)
loopFits <- array(0, dim(batched$forestFits))
loopGlue <- array(0, dim(batched$glue))
for (s in seq_len(numSamples)) {
  sweep <- looped$run(0L, 1L)
  loopFits[,, s] <- sweep$forestFits[,, 1L]
  loopGlue[, s] <- sweep$glue[, 1L]
  expect_identical(sweep$forestFits[, 1L, 1L], looped$getForestFits(0L)[, 1L])
  expect_identical(sweep$forestFits[, 2L, 1L], looped$getForestFits(1L)[, 1L])
  expect_identical(sweep$glue[, 1L], looped$getForestAmplitudes()[, 1L])
}
expect_identical(loopFits, batched$forestFits)
expect_identical(loopGlue, batched$glue)

# --- several chains: each chain fills its own slab, and the reconstruction
# holds per chain (a mixed-up stride would show as a cross-chain reconstruction
# failure, since a and b_z are per-chain draws) ---
multi <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = reportingControl(n.chains = 2L)
)
multiResult <- multi$run(0L, numSamples)
multi$storeState()

expect_equal(dim(multiResult$forestFits), c(n, 2L, numSamples, 2L))
expect_equal(dim(multiResult$glue), c(3L, numSamples, 2L))
for (chain in seq_len(2L)) {
  expect_true(
    reconstructionError(
      multiResult$forestFits[,,, chain],
      multiResult$glue[,, chain],
      internalTrain(multi, multiResult$train[,, chain], chain),
      z
    ) <
      1e-12
  )
}

# --- null by default: a single-forest sampler's result list is the exact eight
# slots it always was, so nothing is allocated or computed for either channel ---
plain <- dbarts(x, y, control = reportingControl(n.chains = 1L))
plainResult <- plain$run(0L, 2L)
expect_null(plainResult$forestFits)
expect_null(plainResult$glue)
expect_identical(
  names(plainResult),
  c("sigma", "train", "test", "varcount", "k", "varprobs", "tau", "ranef")
)

# --- and a K-forest sampler allocates nothing either: the channels follow the
# coupling (are the forests composed through reportable scalars?), never the
# forest count, so multinomial's K forests stay out of them ---
set.seed(4703)
labels <- rbinom(n, 2L, 0.5)
mnHost <- dbarts(x, y, control = reportingControl(n.chains = 1L))
mn <- dbarts:::bartcoreMultinomialSampler(mnHost, labels, K = 3L)
mnResult <- dbarts:::bartcoreRun(mn, 0L, 2L)
expect_null(mnResult$forestFits)
expect_null(mnResult$glue)
