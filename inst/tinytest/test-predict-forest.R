# Out-of-sample per-forest replay: predict(type = "forest") on an
# amplitude-coupled fit, and the $predictForests method under it
# (docs/design/bcf.md). The in-sample twin is test-bcf-forest-channel.R;
# the contract here is that the same quantity comes back at NEW rows, raw and
# per forest, with the recombination left to the caller.

set.seed(709)
n <- 60L
p <- 3L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
y <- 2 * sin(pi * x[, 1L]) + z * (1 + 2 * x[, 2L]) + rnorm(n, sd = 0.2)

nNew <- 25L
xNew <- matrix(runif(nNew * p), nNew, p)

n.samples <- 6L
twoForests <- list(forest(), forest(basis = ~ factor(z)))

# keepTrees rides the SAMPLING run only, as bart() turns it on: the store is a
# circular buffer, so recording burn-in sweeps too would rotate the saved slots
# away from the recorded draws the identity below reads.
keptControl <- function(n.chains = 1L, keepTrees = FALSE) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 8L,
    n.chains = n.chains,
    n.samples = n.samples,
    n.burn = 2L,
    keepTrees = keepTrees,
    updateState = FALSE,
    seed = 709L
  )
}

fitFrom <- function(forests, n.chains = 1L, combineChains = TRUE) {
  sampler <- dbarts(x, y, forests = forests, control = keptControl(n.chains))
  burn <- dbarts:::runWithBurnIn(sampler, sampler$control, TRUE)
  dbarts:::packageBartResults(
    sampler,
    burn$samples,
    burn$burnInSigma,
    burn$burnInK,
    combineChains,
    TRUE
  )
}

fit <- fitFrom(twoForests)

# --- shape and naming: the in-sample channel's layout at the new rows ---
predicted <- predict(fit, xNew, type = "forest")
expect_equal(dim(predicted), c(n.samples, nNew, 2L))
expect_identical(dimnames(predicted)[[3L]], c("forest1", "forest2"))
expect_true(all(is.finite(predicted)))

# --- the identity: replayed at the TRAINING rows, the per-forest channel is
# the one the run recorded. Not bitwise - the amplitude ridge rescales each
# forest's running total as c * sum(mu_t) while the saved leaves are rescaled
# themselves, so the replay re-sums sum(c * mu_t) - so the same 1e-12 bar
# test-argument-surface.R's reconstruction identity uses.
inSample <- extract(fit, type = "forest")
atTraining <- predict(fit, x, type = "forest")
expect_equal(dim(atTraining), dim(inSample))
expect_true(max(abs(atTraining - inSample)) < 1e-12)

# --- and the whole point of the channel: with the caller's own bases, the
# three-line recombination reproduces a fit at new rows, exactly as the
# in-sample identity does. response.shift and the glue belong to the
# combination, which is why neither is folded into what comes back.
shift <- fit$fit$getCalibration(1L)[1L, "response.shift"]
glueForest <- attr(fit$glue, "forest")
recombine <- function(perForest, bases, nRows) {
  out <- matrix(shift, nrow(fit$glue), nRows)
  for (k in seq_len(fit$n.forests)) {
    basis <- if (is.null(bases[[k]])) matrix(1, nRows, 1L) else bases[[k]]
    g <- fit$glue[, glueForest == dimnames(perForest)[[3L]][k], drop = FALSE]
    out <- out + (g %*% t(basis)) * perForest[,, k]
  }
  out
}
expect_true(
  max(abs(recombine(atTraining, fit$bases, n) - fit$yhat.train)) < 1e-12
)
# the same arithmetic at the new rows is finite and moves with the basis: a
# control-row and a treated-row recombination differ by the treatment forest
zeroBasis <- list(NULL, unname(cbind(rep(1, nNew), rep(0, nNew))))
oneBasis <- list(NULL, unname(cbind(rep(0, nNew), rep(1, nNew))))
newAtZero <- recombine(predicted, zeroBasis, nNew)
newAtOne <- recombine(predicted, oneBasis, nNew)
expect_true(all(is.finite(newAtZero)) && all(is.finite(newAtOne)))
expect_true(max(abs(newAtOne - newAtZero)) > 1e-6)

# --- forest selection reuses extract's vocabulary, and always keeps the
# trailing margin ---
byIndex <- predict(fit, xNew, type = "forest", forest = 2L)
byName <- predict(fit, xNew, type = "forest", forest = "forest2")
expect_equal(dim(byIndex), c(n.samples, nNew, 1L))
expect_identical(byIndex, byName)
expect_equal(byIndex[,, 1L], predicted[,, 2L])
expect_error(
  predict(fit, xNew, type = "forest", forest = "forest3"),
  pattern = "must name one of"
)
expect_error(
  predict(fit, xNew, type = "forest", forest = 3L),
  pattern = "index must be between"
)

# --- chains: combineChains folds and splits the same way the in-sample
# channel does ---
fit2 <- fitFrom(twoForests, n.chains = 2L)
combined <- predict(fit2, xNew, type = "forest")
split <- predict(fit2, xNew, type = "forest", combineChains = FALSE)
expect_equal(dim(combined), c(2L * n.samples, nNew, 2L))
expect_equal(dim(split), c(2L, n.samples, nNew, 2L))
expect_equal(split[1L, , , ], combined[seq_len(n.samples), , ])

# --- refusals ---
# an offset shifts the recombination, never one forest's own total
expect_error(
  predict(fit, xNew, type = "forest", offset = rep(0.5, nNew)),
  pattern = "takes no offset"
)
# a fit with no per-forest reporting, by the same predicate extract reads
plainFit <- fitFrom(NULL)
expect_error(
  predict(plainFit, xNew, type = "forest"),
  pattern = "per-forest reporting"
)
# the engine gate itself, reached through the sampler surface: a single-forest
# sampler and a multinomial one both report no per-forest fits, the latter
# although its raw per-category totals are perfectly well defined
expect_error(
  plainFit$fit$predictForests(xNew),
  pattern = "no per-forest fits"
)
yMulti <- factor(c("lo", "mid", "hi")[1L + (z + as.integer(x[, 3L] > 0.5))])
multiFit <- bart2(
  x,
  yMulti,
  family = "multinomial",
  n.trees = 8L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = n.samples,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_error(
  multiFit$fit$predictForests(xNew),
  pattern = "no per-forest fits"
)

# --- a stored state round trip replays the same forests: the trees ride the
# state, so a reloaded sampler answers identically without re-running ---
sampler <- dbarts(
  x,
  y,
  forests = twoForests,
  control = keptControl(keepTrees = TRUE)
)
sampler$run(2L, n.samples)
before <- sampler$predictForests(xNew)
expect_equal(dim(before), c(nNew, 2L, n.samples))
sampler$storeState()
samplerFile <- tempfile(fileext = ".rds")
saveRDS(sampler, samplerFile)
reloaded <- readRDS(samplerFile)
unlink(samplerFile)
expect_true(max(abs(reloaded$predictForests(xNew) - before)) < 1e-12)

# --- forestFits' own storage shape can be 4-d: n.chains > 1 with
# combineChains = FALSE stores (chains, samples, obs, forests), so the forest
# margin is the LAST axis, not the third - predict(type = "forest") must read
# it off the same axis extract's reshape-first path already does ---
fit3 <- fitFrom(twoForests, n.chains = 2L, combineChains = FALSE)
expect_equal(length(dim(fit3$forestFits)), 4L)
predicted3 <- predict(fit3, xNew, type = "forest")
expect_equal(dim(predicted3), c(2L * n.samples, nNew, 2L))
expect_identical(dimnames(predicted3)[[3L]], c("forest1", "forest2"))

byIndex3 <- predict(fit3, xNew, type = "forest", forest = 2L)
byName3 <- predict(fit3, xNew, type = "forest", forest = "forest2")
expect_identical(byIndex3, byName3)
expect_equal(byIndex3[,, 1L], predicted3[,, 2L])

inSample3 <- extract(fit3, type = "forest")
atTraining3 <- predict(fit3, x, type = "forest")
expect_equal(dim(atTraining3), dim(inSample3))
expect_true(max(abs(atTraining3 - inSample3)) < 1e-12)

# --- newdata may arrive sparse: a plain continuous fixture declares no leaf
# covariate, so refuseSparseLeafCovariate has nothing to refuse and a bare
# dgCMatrix rides resident through the same funnel the "bart" channel uses ---
if (requireNamespace("Matrix", quietly = TRUE)) {
  xNewSparse <- as(as(xNew, "CsparseMatrix"), "generalMatrix")
  predictedSparse <- predict(fit, xNewSparse, type = "forest")
  expect_equal(dim(predictedSparse), dim(predicted))
  expect_true(max(abs(predictedSparse - predicted)) < 1e-12)
}
