# The dbartsSampler surface, post host-shell removal. Pointer adoption for
# ordinal/nbinom and direct construction for multinomial mean every
# bart2() alternate-family fit's $fit is now the sampler that actually
# ran: no hostFor field, no refuseHostMutation/refuseHostRead guards, no
# host-shell save/reload defect. This file used to pin those defects; it
# now pins the census (the drift detector: any method added or removed
# from dbartsSampler should be caught here) and the capabilities that
# replaced them.

# --- the census: every own method of dbartsSampler. Read the SOURCE rather
# than probe each mutator by call, since most need contrived arguments this
# file has no stake in constructing.
gen <- getRefClass("dbartsSampler")
ownMethods <- gen$def@refMethods
infrastructure <- c(
  "initialize",
  "adoptPointer",
  "reapplyForestWeights",
  "getPointer",
  "show"
)
# envRefClass's own inherited machinery, not part of dbartsSampler's surface
inherited <- c(
  ".objectPackage",
  ".objectParent",
  "callSuper",
  "copy#envRefClass",
  "export",
  "field",
  "getClass",
  "getRefClass",
  "import",
  "initFields",
  "show#envRefClass",
  "trace",
  "untrace",
  "usingMethods"
)
own <- setdiff(names(ownMethods), inherited)
substantiveMethods <- sort(setdiff(own, infrastructure))

expect_equal(length(own), 46L)
expect_equal(length(infrastructure), 5L)
expect_equal(length(substantiveMethods), 41L)

expect_identical(
  substantiveMethods,
  sort(c(
    "run",
    "sampleTreesFromPrior",
    "sampleNodeParametersFromPrior",
    "growFromRoot",
    "setControl",
    "setModel",
    "setData",
    "setResponse",
    "setOffset",
    "setWeights",
    "setCounts",
    "setCategoryOffset",
    "setCategoryTestOffset",
    "setActiveRows",
    "setForestWeights",
    "setForestBasis",
    "setSigma",
    "setPredictor",
    "setCutPoints",
    "setTestPredictor",
    "setTestPredictorAndOffset",
    "setTestOffset",
    "setCalibration",
    "setState",
    "installTrees",
    "getDispersion",
    "getFitsWithoutOffset",
    "copy",
    "predict",
    "predictForests",
    "getLatents",
    "getSigmas",
    "getSumsOfSquaredResiduals",
    "getForestFits",
    "getForestAmplitudes",
    "getForestVariableCounts",
    "getCalibration",
    "storeState",
    "printTrees",
    "getTrees",
    "plotTree"
  ))
)

# --- spot check: a multinomial fit's $fit is the K-forest engine that ran,
# not a placeholder - mutation succeeds where the softmax gives it meaning
# and refuses (by a model reason, never a host-shell one) where it does not;
# predict reports real, non-degenerate probabilities ---
set.seed(94011)
nSpot <- 60L
xSpot <- matrix(runif(nSpot * 2L), nSpot, 2L)
ySpot <- factor(
  c("a", "b", "c")[1L + (xSpot[, 1L] > 0.5) + (xSpot[, 1L] > 0.8)],
  levels = c("a", "b", "c")
)
fitSpot <- bart2(
  xSpot,
  ySpot,
  family = "multinomial",
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = 2L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(fitSpot$fit$setResponse(rnorm(nSpot)), "\\$setCounts")
expect_null(fitSpot$fit$getDispersion())
predSpot <- fitSpot$fit$predict(xSpot)
expect_equal(dim(predSpot), c(nSpot, 3L, 2L))
expect_equal(apply(predSpot[,, 1L], 1L, sum), rep(1.0, nSpot))

# --- fixtures for the remaining pins: one multinomial K = 3, one K = 2, one
# ordinal, one nbinom, each with $fit retained via keepTrees ---
n.trees <- 15L
n.burn <- 15L
n.samples <- 15L

set.seed(9401)
n3 <- 120L
p3 <- 3L
x3 <- matrix(runif(n3 * p3), n3, p3)
eta3 <- cbind(
  2 * (x3[, 1L] - 0.5),
  x3[, 2L] - x3[, 3L],
  1.2 * (x3[, 1L] - x3[, 3L])
)
probs3 <- exp(eta3) / rowSums(exp(eta3))
labels3 <- vapply(
  seq_len(n3),
  function(i) sample.int(3L, 1L, prob = probs3[i, ]) - 1L,
  integer(1L)
)
y3 <- factor(c("a", "b", "c")[labels3 + 1L], levels = c("a", "b", "c"))
set.seed(9402)
fit3 <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(9403)
n2 <- 100L
p2 <- 2L
x2 <- matrix(runif(n2 * p2), n2, p2)
labels2 <- rbinom(n2, 1L, plogis(2 * (x2[, 1L] - 0.5) + x2[, 2L]))
y2 <- factor(c("no", "yes")[labels2 + 1L], levels = c("no", "yes"))
set.seed(9404)
fit2 <- bart2(
  x2,
  y2,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(99)
nOrd <- 80L
xOrd <- matrix(runif(nOrd * 3L), nOrd, 3L)
zOrd <- 2 * (xOrd[, 1L] - 0.5) + rnorm(nOrd)
codesOrd <- 1L + (zOrd > 0) + (zOrd > 0.8)
yOrd <- ordered(c("lo", "mid", "hi")[codesOrd], levels = c("lo", "mid", "hi"))
fitOrd <- bart2(
  xOrd,
  yOrd,
  family = "ordinal",
  n.trees = 10L,
  n.chains = 1L,
  n.burn = 10L,
  n.samples = 10L,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(199)
nNb <- 80L
xNb <- matrix(runif(nNb * 3L), nNb, 3L)
etaNb <- 0.9 * (xNb[, 1L] - 0.5)
yNb <- rnbinom(nNb, size = 5L, mu = exp(etaNb))
fitNb <- bart2(
  xNb,
  yNb,
  family = "nbinom",
  n.trees = 10L,
  n.chains = 1L,
  n.burn = 10L,
  n.samples = 10L,
  verbose = FALSE,
  keepTrees = TRUE
)

# --- $copy() carries an independent, running K-forest engine: mutating one
# side leaves the other untouched ---
dupe3 <- fit3$fit$copy(shallow = TRUE)
originalCounts3 <- fit3$fit$data@counts
swappedCounts3 <- originalCounts3[c(2:n3, 1L), , drop = FALSE]
dupe3$setCounts(swappedCounts3)
expect_identical(dupe3$data@counts, swappedCounts3)
expect_identical(fit3$fit$data@counts, originalCounts3)
expect_silent(invisible(fit3$fit$run(0L, 1L)))

# --- ordinal's $fit accepts mutation on the original and its copy alike ---
dupeOrd <- fitOrd$fit$copy(shallow = TRUE)
expect_silent(dupeOrd$setData(dbartsData(xOrd, as.double(codesOrd))))

# --- $fit is the K-forest engine that ran: predict reports real,
# non-degenerate softmax probabilities, not a constant placeholder ---
predHost3 <- fit3$fit$predict(x3)
expect_equal(dim(predHost3), c(n3, 3L, n.samples))
expect_true(length(unique(round(predHost3[, 1L, 1L], 6L))) > 1L)

# --- $fit is now an ordinary multinomial K-forest sampler regardless of K:
# no more K-dependent placeholder family (gaussian at K >= 3, probit at
# K = 2) - both resolve to "multinomial" directly ---
expect_equal(fit3$fit$model@family, "multinomial")
expect_false(fit3$fit$control@binary)
expect_equal(fit2$fit$model@family, "multinomial")
expect_false(fit2$fit$control@binary)

# --- family = "auto" on a 2-level factor still never reaches multinomial -
# it announces probit directly and returns class "bart" ---
autoWarnings <- 0L
fitAuto <- withCallingHandlers(
  bart2(
    x2,
    y2,
    family = "auto",
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 2L,
    n.samples = 2L,
    verbose = FALSE
  ),
  message = function(m) {
    autoWarnings <<- autoWarnings + 1L
    invokeRestart("muffleMessage")
  }
)
expect_equal(autoWarnings, 1L)
expect_identical(class(fitAuto), "bart")
expect_false(inherits(fitAuto, "bartMultinomial"))

# --- save/reload/predict round-trips for all three host-shell families now:
# $fit is the sampler whose engine actually ran, so getPointer's
# re-creation branch rebuilds it from stored state and predict replays
# through it ---
checkSaveReloadPredicts <- function(fit, xtest) {
  before <- predict(fit, xtest)
  fit$fit$storeState()
  tempFile <- tempfile()
  saveRDS(fit, tempFile)
  reloaded <- readRDS(tempFile)
  unlink(tempFile)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_equal(predict(reloaded, xtest), before)
}
checkSaveReloadPredicts(fit3, x3[1:5, ])
checkSaveReloadPredicts(fit2, x2[1:5, ])
checkSaveReloadPredicts(fitOrd, xOrd[1:5, ])
checkSaveReloadPredicts(fitNb, xNb[1:5, ])
