# The bartMultinomial generics that used to fall through to a "bart"/default
# method (print's per-chain arithmetic,
# residuals, plot, summary) or had no type argument (predict).

parseKeptDraws <- function(printed) {
  line <- grep("kept draws", printed, value = TRUE)
  as.integer(sub(".*: ", "", line))
}

set.seed(2201)
n <- 150L
p <- 3L
x <- matrix(runif(n * p), n, p)
eta <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 1.5 * (x[, 3L] - 0.5))
probs <- exp(eta) / rowSums(exp(eta))
labels <- vapply(
  seq_len(n),
  function(i) sample.int(3L, 1L, prob = probs[i, ]) - 1L,
  integer(1L)
)
y <- factor(c("lo", "mid", "hi")[labels + 1L], levels = c("lo", "mid", "hi"))

n.trees <- 10L
n.burn <- 5L
n.samples <- 7L

# ---- print.bartMultinomial's "kept draws (per chain)" arithmetic ----
# combineChains = TRUE (the default) must report the per-chain count, not
# the combined total (n.samples * n.chains).

set.seed(3301)
fitCombined <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 3L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  combineChains = TRUE
)
set.seed(3301)
fitSplit <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 3L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  combineChains = FALSE
)
# same fit (same seed, same settings) reported two ways: both must read the
# per-chain count, never the combined total
expect_equal(parseKeptDraws(capture.output(print(fitCombined))), n.samples)
expect_equal(parseKeptDraws(capture.output(print(fitSplit))), n.samples)
expect_false(
  parseKeptDraws(capture.output(print(fitCombined))) == n.samples * 3L
)

# single chain: division by n.chains == 1 is a no-op, so this path was never
# buggy, but must stay correct
fitSingle <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = 11L,
  verbose = FALSE
)
expect_equal(parseKeptDraws(capture.output(print(fitSingle))), 11L)

# ---- residuals.bartMultinomial must not fall through to residuals.default
# (which would silently return NULL) ----

res <- residuals(fitCombined)
expect_false(is.null(res))
expect_equal(dim(res), c(n, 3L))
expect_equal(colnames(res), c("lo", "mid", "hi"))
# 1[y = k] - phat_k sums to 0 across categories, both indicator and phat rows
# being simplexes
expect_true(all(abs(rowSums(res)) < 1e-8))
# matches a hand-rolled indicator-minus-fitted computation exactly
phatCombined <- fitted(fitCombined)
indicator <- matrix(0, n, 3L, dimnames = dimnames(phatCombined))
indicator[cbind(seq_len(n), match(y, levels(y)))] <- 1
expect_equal(res, indicator - phatCombined)

# grouped-count ingestion: the observed-proportion analog (y / rowSums(y))
set.seed(2202)
n.c <- 120L
x.c <- matrix(runif(n.c * p), n.c, p)
eta.c <- cbind(
  2 * (x.c[, 1L] - 0.5),
  x.c[, 2L] - x.c[, 3L],
  1.5 * (x.c[, 3L] - 0.5)
)
probs.c <- exp(eta.c) / rowSums(exp(eta.c))
trials.c <- sample(2:6, n.c, replace = TRUE)
counts.c <- t(vapply(
  seq_len(n.c),
  function(i) rmultinom(1L, trials.c[i], probs.c[i, ])[, 1L],
  integer(3L)
))
colnames(counts.c) <- c("lo", "mid", "hi")
fitCounts <- bart2(
  x.c,
  counts.c,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
resCounts <- residuals(fitCounts)
expect_equal(dim(resCounts), c(n.c, 3L))
expect_equal(resCounts, counts.c / rowSums(counts.c) - fitted(fitCounts))

# ---- plot.bartMultinomial must not fall through to plot.default (which
# errors on "'x' is a list, but does not have components 'x' and 'y'") ----

pdf(NULL)
expect_silent(plot(fitCombined))
expect_silent(plot(fitSingle))
dev.off()

# ---- summary.bartMultinomial must not fall through to summary.default (a
# raw structure dump) ----

sm <- summary(fitCombined)
expect_equal(class(sm), "summary.bart")
expect_equal(sort(sm$stats$variable), sort(paste0("meanProb[", levels(y), "]")))
expect_true(is.data.frame(sm$stats) || inherits(sm$stats, "tbl_df"))
printedSummary <- capture.output(print(sm))
expect_true(any(grepl("meanProb", printedSummary, fixed = TRUE)))

havePosterior <- requireNamespace("posterior", quietly = TRUE)
if (havePosterior) {
  expect_true(sm$posterior)
  expect_true("rhat" %in% names(sm$stats))
} else {
  expect_false(sm$posterior)
  expect_true(is.data.frame(sm$stats))
}

# degrade path, simulated the same way test-convergence-diagnostics.R does
if (havePosterior) {
  oldPosteriorAvailable <- dbarts:::posteriorAvailable
  tryCatch(
    {
      assignInNamespace("posteriorAvailable", function() FALSE, ns = "dbarts")
      smDegraded <- summary(fitCombined)
      expect_false(smDegraded$posterior)
      expect_true(is.data.frame(smDegraded$stats))
      expect_false("rhat" %in% names(smDegraded$stats))
    },
    finally = assignInNamespace(
      "posteriorAvailable",
      oldPosteriorAvailable,
      ns = "dbarts"
    )
  )
  rm(oldPosteriorAvailable, smDegraded)
}

# ---- predict.bartMultinomial's 'type' argument ----

x.test <- x[seq_len(20L), , drop = FALSE]
set.seed(3302)
fitKeep <- bart2(
  x,
  y,
  family = "multinomial",
  test = x.test,
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)

# the default type = "ev" is bitwise identical to a pre-change-style call
# (no type argument at all) and consumes no RNG state
seedBefore <- .Random.seed
predDefault <- predict(fitKeep, x.test)
predEv <- predict(fitKeep, x.test, type = "ev")
expect_identical(predDefault, predEv)
expect_identical(predEv, fitKeep$yhat.test)
expect_identical(.Random.seed, seedBefore)

# type = "ppd" draws valid categories, shaped like "ev" minus the K margin,
# and DOES consume RNG state (only when explicitly requested)
seedBeforePpd <- .Random.seed
predPpd <- predict(fitKeep, x.test, type = "ppd")
expect_false(identical(.Random.seed, seedBeforePpd))
expect_equal(dim(predPpd), dim(predEv)[-length(dim(predEv))])
expect_true(all(predPpd %in% seq_len(3L)))

# predict's ppd matches extract's ppd construction exactly: replaying the
# fit-time test rows reproduces the stored test channel bit for bit (the
# existing reproduction gate), so with the RNG aligned the same way, the two
# ppd draws - one built from the replayed channel, one from the stored one -
# must agree exactly as well
set.seed(4401)
ppdFromExtract <- extract(fitKeep, type = "ppd", sample = "test")
set.seed(4401)
ppdFromPredict <- predict(fitKeep, x.test, type = "ppd")
expect_identical(ppdFromExtract, ppdFromPredict)

# type is validated like every other type= argument in the package, and the
# latent-scale request is refused by NAME rather than by an enum message
expect_error(
  predict(fitKeep, x.test, type = "bart"),
  "non-identified and unrecorded"
)

# type = "class" is fitted()'s own argmax reduction, reached through predict's
# own replay; a class prediction is a label rather than a quantity with a
# credible band, so pairing it with ci.level is refused by name instead of
# silently returning the type = "ev" band in its place
predClass <- predict(fitKeep, x.test, type = "class")
expect_true(is.factor(predClass))
expect_identical(levels(predClass), levels(y))
expect_equal(length(predClass), nrow(x.test))
expect_identical(
  predict(fitKeep, x, type = "class"),
  fitted(fitKeep, type = "class")
)
expect_error(
  predict(fitKeep, x.test, type = "class", ci.level = 0.9),
  "does not support 'ci.level'"
)
expect_error(
  fitted(fitKeep, type = "class", ci.level = 0.9),
  "does not support 'ci.level'"
)

# ---- type= runs through validateType, so the predict.glm synonyms reach
# these methods, and the two latent-scale values they name are refused with
# a reason ----

expect_identical(predict(fitKeep, x.test, type = "response"), predEv)
expect_identical(
  extract(fitKeep, type = "response"),
  extract(fitKeep, type = "ev")
)
expect_identical(fitted(fitKeep, type = "response"), fitted(fitKeep))

latentReason <- paste0(
  "multinomial fits do not support type = \"bart\": the run records ",
  "only the identified softmax probabilities; the raw per-category ",
  "latent fits are non-identified and unrecorded"
)
forestReason <- paste0(
  "multinomial fits do not support type = \"forest\": a category's ",
  "forest is a latent whose level is reproducibly structured yet not ",
  "identified, so a raw replay reads as signal; the identified content ",
  "is the log-ratio of the probabilities predict() reports"
)
# "link" folds onto "bart", so it lands on that value's own refusal
expect_error(
  predict(fitKeep, x.test, type = "link"),
  latentReason,
  fixed = TRUE
)
expect_error(extract(fitKeep, type = "link"), latentReason, fixed = TRUE)
expect_error(fitted(fitKeep, type = "bart"), latentReason, fixed = TRUE)
expect_error(
  predict(fitKeep, x.test, type = "forest"),
  forestReason,
  fixed = TRUE
)
expect_error(extract(fitKeep, type = "forest"), forestReason, fixed = TRUE)
# a value no multinomial method offers refuses against the set it does
expect_error(extract(fitKeep, type = "nonsense"), "type must be in 'ev'")
# loglik is extract-only
expect_error(predict(fitKeep, x.test, type = "loglik"), "type must be in")
expect_error(fitted(fitKeep, type = "loglik"), "type must be in")

# ---- predict on a fit trained with an n x K category offset requires an
# explicit offset argument: the predicted rows are the caller's, so no
# resident offset describes them ----

set.seed(5501)
categoryOffset <- matrix(rnorm(n * 3L, sd = 0.4), n, 3L)
fitOffset <- bart2(
  x,
  y,
  family = "multinomial",
  offset = categoryOffset,
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_error(
  predict(fitOffset, x.test),
  "'offset' is required on a multinomial fit trained with a",
  fixed = TRUE
)
expect_error(
  predict(
    fitOffset,
    x.test,
    offset = matrix(0, nrow(x.test), 2L)
  ),
  "must be a 20 x 3 matrix"
)
# ORACLE: the training rows under the training offset ARE the training
# channel, replayed draw for draw out of the saved trees (measured 3.9e-15)
expect_equal(
  predict(fitOffset, x, offset = categoryOffset),
  fitOffset$yhat.train,
  tolerance = 1e-12
)
# an all-zero matrix asks for the offset-free surface on purpose
predFree <- predict(
  fitOffset,
  x.test,
  offset = matrix(0, nrow(x.test), 3L)
)
expect_equal(dim(predFree), c(n.samples, nrow(x.test), 3L))
expect_true(max(abs(apply(predFree, c(1L, 2L), sum) - 1)) < 1e-12)
# a fit trained without one defaults to NULL, and still takes an explicit
# matrix: the shift enters the raw fits, so it moves the reported surface
expect_identical(predict(fitKeep, x.test), predEv)
expect_false(identical(
  predict(
    fitKeep,
    x.test,
    offset = matrix(c(1, 0, 0), nrow(x.test), 3L, byrow = TRUE)
  ),
  predEv
))

# ---- extract's combineChains formal must be honoured, not ignored (which
# would silently return the fit's own stored layout) ----

fitMC <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 2L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
combinedEv <- extract(fitMC, type = "ev")
splitEv <- extract(fitMC, type = "ev", combineChains = FALSE)
expect_equal(dim(combinedEv), c(2L * n.samples, n, 3L))
expect_equal(dim(splitEv), c(2L, n.samples, n, 3L))

# ---- extract(type = "loglik") is the multinomial log-pmf (its coefficient
# included), checked against an independently coded oracle (dmultinom) for
# both response ingestions; extract-only, sample = "test" refused ----

llLabel <- extract(fitCombined, type = "loglik")
evLabel <- extract(fitCombined, type = "ev")
oracleLabel <- matrix(0, nrow(llLabel), ncol(llLabel))
for (s in seq_len(nrow(llLabel))) {
  for (i in seq_len(ncol(llLabel))) {
    indicatorRow <- as.integer(y[i] == levels(y))
    oracleLabel[s, i] <- dmultinom(
      indicatorRow,
      prob = evLabel[s, i, ],
      log = TRUE
    )
  }
}
expect_equal(llLabel, oracleLabel, tolerance = 1e-12)
expect_equal(dim(llLabel), dim(evLabel)[-3L])

llCounts <- extract(fitCounts, type = "loglik")
evCounts <- extract(fitCounts, type = "ev")
oracleCounts <- vapply(
  seq_len(n.c),
  function(i) dmultinom(counts.c[i, ], prob = evCounts[1L, i, ], log = TRUE),
  numeric(1L)
)
expect_equal(llCounts[1L, ], oracleCounts, tolerance = 1e-12)

expect_error(
  extract(fitCombined, type = "loglik", sample = "test"),
  "no test response exists"
)
llSplit <- extract(fitMC, type = "loglik", combineChains = FALSE)
expect_equal(dim(llSplit), c(2L, n.samples, n))

# ---- fitted/predict's ci.level returns an (obs x K x 3) array on this
# K-margined family: a plain 3-column matrix has no room for the category
# margin, so pooling it would average incomparable probabilities ----

fitCi <- fitted(fitCombined, ci.level = 0.9)
expect_equal(dim(fitCi), c(n, 3L, 3L))
expect_equal(dimnames(fitCi)[[3L]], c("est", "ci.lower", "ci.upper"))
predCi <- predict(fitKeep, x.test, ci.level = 0.9)
expect_equal(dim(predCi), c(20L, 3L, 3L))

# ---- bart-family-only formals are refused by name instead of vanishing
# into '...' - verified against the unmodified build by hand (identical()
# to the call without the argument, the swallow's signature) ----

expect_error(
  extract(fitCombined, forest = 1L),
  "'forest' is not used by extract on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  extract(fitCombined, contribution = TRUE),
  "'contribution' is not used by extract on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  predict(fitKeep, x.test, forest = 1L),
  "'forest' is not used by predict on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  predict(fitKeep, x.test, contribution = TRUE),
  "'contribution' is not used by predict on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  fitted(fitCombined, forest = 1L),
  "'forest' is not used by fitted on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  residuals(fitCombined, forest = 1L),
  "'forest' is not used by residuals on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  fitted(fitCombined, sample = "test"),
  "'sample' is not used by fitted on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  residuals(fitCombined, type = "bart"),
  "'type' is not used by residuals on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  summary(fitCombined, vars = "sigma"),
  "'vars' is not used by summary on a bartMultinomial fit",
  fixed = TRUE
)

# ---- plotTree/survivalProbabilities are refused by name instead of a raw
# 'no applicable method' ----

expect_error(
  plotTree(fitCombined),
  "plotTree is defined for bart, rbart_vi and dbartsSampler fits",
  fixed = TRUE
)
expect_error(
  survivalProbabilities(fitCombined),
  "survivalProbabilities applies to a discrete-time hazard fit",
  fixed = TRUE
)

# ---- plot(object) draws a second panel for both response ingestions, and
# restores the caller's par afterward - assert restoration against a
# sentinel value distinct from both the plot's own layout and the device
# default ----

pdf(NULL)
par(mfrow = c(3L, 3L))
plot(fitCombined)
combinedMfrow <- par("mfrow")
plot(fitCounts)
countsMfrow <- par("mfrow")
dev.off()
expect_equal(combinedMfrow, c(3L, 3L))
expect_equal(countsMfrow, c(3L, 3L))

# ---- as_draws_array/df expose meanProb[<level>], never yhat.train ----

if (havePosterior) {
  meanProbNames <- paste0("meanProb[", levels(y), "]")
  ad <- posterior::as_draws_array(fitCombined)
  expect_equal(sort(dimnames(unclass(ad))[[3L]]), sort(meanProbNames))
  adf <- posterior::as_draws_df(fitCombined)
  expect_true(all(meanProbNames %in% names(adf)))
  rm(meanProbNames, ad, adf)
}

rm(
  parseKeptDraws,
  n,
  p,
  x,
  eta,
  probs,
  labels,
  y,
  n.trees,
  n.burn,
  n.samples,
  fitCombined,
  fitSplit,
  fitSingle,
  res,
  phatCombined,
  indicator,
  n.c,
  x.c,
  eta.c,
  probs.c,
  trials.c,
  counts.c,
  fitCounts,
  resCounts,
  sm,
  printedSummary,
  havePosterior,
  x.test,
  fitKeep,
  seedBefore,
  predDefault,
  predEv,
  seedBeforePpd,
  predPpd,
  ppdFromExtract,
  ppdFromPredict,
  latentReason,
  forestReason,
  categoryOffset,
  fitOffset,
  predFree,
  fitMC,
  combinedEv,
  splitEv,
  llLabel,
  evLabel,
  oracleLabel,
  s,
  i,
  indicatorRow,
  llCounts,
  evCounts,
  oracleCounts,
  llSplit,
  fitCi,
  predCi,
  combinedMfrow,
  countsMfrow
)
