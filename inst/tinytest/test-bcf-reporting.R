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
    seed = 47L,
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

# the reconstruction identity: a * mu + b_z * tau over the recorded channels
# reproduces the recorded
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
  expect_identical(sweep$forestFits[, 1L, 1L], looped$getForestFits(1L)[, 1L])
  expect_identical(sweep$forestFits[, 2L, 1L], looped$getForestFits(2L)[, 1L])
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

# --- the per-draw per-forest VARIABLE-COUNT channel: on a multi-forest sampler
# run()$varcount carries a forest axis between the predictors and the samples
# (n.predictors x n.forests x n.samples x n.chains, forest-major within a
# sample, prognostic first), the same axis the multinomial coupling's per-
# category counts ride. $getForestVariableCounts reads the CURRENT state; this
# channel is the draw history it never was, and the two meet at the last kept
# draw ---
counting <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = reportingControl(n.chains = 2L, n.thin = 1L)
)
countingResult <- counting$run(0L, numSamples)

expect_equal(dim(countingResult$varcount), c(p, 2L, numSamples, 2L))
expect_true(all(countingResult$varcount >= 0L))

# the per-draw varcount contract, positive half: with n.thin = 1 and no
# sweep since the last storeSample, the last kept draw's slab IS the live
# read, forest by forest and chain by chain - so the forest axis is really
# keyed on the forest and not transposed
expect_identical(
  countingResult$varcount[, 1L, numSamples, ],
  counting$getForestVariableCounts(1L)
)
expect_identical(
  countingResult$varcount[, 2L, numSamples, ],
  counting$getForestVariableCounts(2L)
)
# and the two forests genuinely differ, so the pin above could fail
expect_false(
  identical(
    counting$getForestVariableCounts(1L),
    counting$getForestVariableCounts(2L)
  )
)

# the per-draw varcount contract, NEGATIVE half: one more sweep and the
# equality must break, which is what separates a draw history from a
# repeated live read
invisible(counting$run(0L, 1L))
expect_false(
  identical(
    cbind(
      countingResult$varcount[, 1L, numSamples, ],
      countingResult$varcount[, 2L, numSamples, ]
    ),
    cbind(
      counting$getForestVariableCounts(1L),
      counting$getForestVariableCounts(2L)
    )
  )
)

# the mask-reaches-the-channel contract: a column the treatment forest is
# masked off is structurally zero in EVERY draw of its slab, while the
# prognostic forest is unrestricted - the mask reaches the per-draw channel,
# not merely the live read
restrictedCounts <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), vars = c(3L, 4L))),
  control = reportingControl(n.chains = 2L)
)
restrictedResult <- restrictedCounts$run(0L, numSamples)
expect_true(all(restrictedResult$varcount[1:2, 2L, , ] == 0L))
expect_true(sum(restrictedResult$varcount[3:4, 2L, , ]) > 0L)
expect_true(sum(restrictedResult$varcount[1:2, 1L, , ]) > 0L)

# the mask-reaches-the-channel contract, NEGATIVE half: at the same seed,
# without the mask, the masked columns are non-zero in the treatment
# forest's slab
unrestrictedCounts <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = reportingControl(n.chains = 2L)
)
unrestrictedResult <- unrestrictedCounts$run(0L, numSamples)
expect_true(sum(unrestrictedResult$varcount[1:2, 2L, , ]) > 0L)

# --- and the channel is UNCHANGED on a single-forest sampler: the forest axis
# appears only when there is more than one forest to report, so every existing
# caller's n.predictors x n.samples x n.chains array is what it always was ---
plainCounts <- dbarts(x, y, control = reportingControl(n.chains = 2L))
plainCountsResult <- plainCounts$run(0L, numSamples)
expect_equal(dim(plainCountsResult$varcount), c(p, numSamples, 2L))
expect_identical(
  plainCountsResult$varcount[, numSamples, ],
  plainCounts$getForestVariableCounts(1L)
)

# --- the multi-forest varcount packaging: the bart2 packaging path. A
# dbartsData carrying bases reaches bart2
# today, and its varcount is the one channel that widens, so it is packaged
# through the same K-margin reshape multinomial's per-category counts take:
# draws-first, predictor names on the lead margin, engine-vocabulary forest
# names on the trailing one. n.forests records the count the packaged rank
# alone cannot recover ---
bcfData <- dbartsData(x, y, bases = list(NULL, model.matrix(~ factor(z) - 1)))
bcfFitArgs <- list(
  formula = bcfData,
  n.samples = numSamples,
  n.burn = 3L,
  n.trees = 20L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 51L
)
combinedFit <- do.call(bart2, c(bcfFitArgs, list(n.chains = 2L)))
uncombinedFit <- do.call(
  bart2,
  c(bcfFitArgs, list(n.chains = 2L, combineChains = FALSE))
)
oneChainFit <- do.call(bart2, c(bcfFitArgs, list(n.chains = 1L)))

expect_equal(dim(combinedFit$varcount), c(2L * numSamples, p, 2L))
expect_equal(dim(uncombinedFit$varcount), c(2L, numSamples, p, 2L))
expect_equal(dim(oneChainFit$varcount), c(numSamples, p, 2L))
expect_identical(
  dimnames(combinedFit$varcount),
  list(NULL, colnames(bcfData@x), c("forest1", "forest2"))
)
expect_identical(
  dimnames(uncombinedFit$varcount),
  list(NULL, NULL, colnames(bcfData@x), c("forest1", "forest2"))
)
expect_identical(combinedFit$n.forests, 2L)
expect_identical(uncombinedFit$n.forests, 2L)
expect_null(
  bart2(
    x,
    y,
    n.samples = 2L,
    n.burn = 2L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 51L
  )$n.forests
)

# the forest-1 slab really is the prognostic forest's, and the combined chain
# margin folds chain-major (chain c's last draw is row c * n.samples), pinned
# against the live per-forest read of the sampler that produced it
keptFit <- do.call(
  bart2,
  c(bcfFitArgs, list(n.chains = 2L, keepSampler = TRUE))
)
for (chain in seq_len(2L)) {
  expect_identical(
    keptFit$varcount[chain * numSamples, , 1L],
    keptFit$fit$getForestVariableCounts(1L)[, chain]
  )
  expect_identical(
    keptFit$varcount[chain * numSamples, , 2L],
    keptFit$fit$getForestVariableCounts(2L)[, chain]
  )
}

# the synopsis reports the TRUE per-chain kept-draw count on a fit with no
# kept sampler, where the count can only come off the varcount dimensions
expect_true(any(grepl(
  paste0("kept draws \\(per chain\\): ", numSamples),
  capture.output(print(combinedFit))
)))
expect_true(any(grepl(
  paste0("kept draws \\(per chain\\): ", numSamples),
  capture.output(print(uncombinedFit))
)))
expect_true(any(grepl(
  paste0("kept draws \\(per chain\\): ", numSamples),
  capture.output(print(oneChainFit))
)))

# the multi-forest varcount packaging, NEGATIVE half: without the n.forests
# arm the synopsis would take the
# single-forest branch on the same dimensions and report the PREDICTOR count
# instead of the draw count, so the assertions above pin the arm rather than
# the arithmetic
unarmed <- function(fit, n.chains) {
  d <- dim(fit[["varcount"]])
  if (length(d) == 3L) d[2L] else d[1L] %/% n.chains
}
expect_identical(unarmed(combinedFit, 2L), p)
expect_false(unarmed(combinedFit, 2L) == numSamples)
expect_false(unarmed(uncombinedFit, 2L) == numSamples)
