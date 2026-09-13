# Heteroscedastic BART: a `variance` selector
# adds a second forest modeling s^2(x). The fit recovers f(x) and a plausible
# s(x); a homoscedastic truth does not manufacture spurious heteroscedasticity;
# the surface is gaussian + constant-leaf only, and predict returns s(x).

set.seed(9, sample.kind = "Rejection")

# ---- a step-heteroscedastic truth: the fit recovers f(x) and s(x) ----
n <- 800L
x <- runif(n)
fTrue <- 2 * x
sTrue <- ifelse(x < 0.5, 0.3, 1.5)
y <- fTrue + sTrue * rnorm(n)

fit <- bart(
  x,
  y,
  variance = varianceForest(n.trees = 25L),
  n.trees = 50L,
  n.samples = 400L,
  n.burn = 400L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# s(x) is reported, one posterior slab per training row
expect_true(!is.null(fit$s.train))
expect_equal(dim(fit$s.train), c(400L, n))

sHat <- apply(fit$s.train, 2L, mean)
fHat <- fit$yhat.train.mean

# the mean surface tracks the ramp
expect_true(cor(fHat, fTrue) > 0.9)

# s(x) is larger where the truth is noisier, and tracks the two levels within a
# factor of ~2 (the estimate attenuates with finite trees)
sLow <- mean(sHat[x < 0.5])
sHigh <- mean(sHat[x >= 0.5])
expect_true(sHigh > 2 * sLow)
expect_true(sLow > 0.15 && sLow < 0.6)
expect_true(sHigh > 0.9 && sHigh < 2.2)

# ---- predict returns s(x) alongside f(x) on new data ----
xNew <- matrix(c(0.25, 0.75), 2L, 1L)
pred <- predict(fit, xNew)
s <- attr(pred, "s")
expect_true(!is.null(s))
expect_equal(dim(s), c(400L, 2L))
sNew <- apply(s, 2L, mean)
expect_true(sNew[2L] > 2 * sNew[1L]) # x = 0.75 noisier than x = 0.25

# ---- the SAVED variance trees survive a state round trip ----
# predict addresses the saved variance buffer, never the live trees, so a
# re-created sampler that restored only the live ones replays the identity fill
# and reports s(x) == 0.
fit$fit$storeState()
expect_false(is.null(fit$fit$state[[1L]][["variance.saved.vars"]]))

stateFile <- tempfile()
saveRDS(fit, stateFile)
fitLoaded <- readRDS(stateFile)
predLoaded <- predict(fitLoaded, xNew)
expect_identical(predLoaded, pred)
expect_identical(attr(predLoaded, "s"), s)
expect_true(all(attr(predLoaded, "s") > 0))
unlink(stateFile)

fitCopy <- fit
fitCopy$fit <- fit$fit$copy()
predCopy <- predict(fitCopy, xNew)
expect_identical(predCopy, pred)
expect_identical(attr(predCopy, "s"), s)

# the state-comparison helper reaches the chain-level variance blocks; every
# other caller is homoscedastic, where both sides read NULL and the comparison
# is vacuous
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
fitCopy$fit$storeState()
statesAgree(fitCopy$fit$state, fit$fit$state)
# and it bites: one perturbed saved variance leaf is detected
mutatedState <- fit$fit$state
mutatedState[[1L]][["variance.saved.values"]][1:8] <- as.raw(0L)
expect_false(statesAgree(mutatedState, fit$fit$state, expect = FALSE))

# an ABSENT saved block against a live capacity can only be a state written
# before the channel existed; restoring it would substitute the identity fill
# for the recorded surface and report a plausible constant s(x)
strippedState <- fit$fit$state
for (block in c("vars", "values", "sizes", "flags")) {
  strippedState[[1L]][[paste0("variance.saved.", block)]] <- NULL
}
expect_error(
  fit$fit$setState(strippedState),
  "state is not consistent with this sampler"
)

# a PRESENT but malformed block is named rather than reported generically
badSavedState <- fit$fit$state
badSavedState[[1L]][["variance.saved.sizes"]] <- c(1L, 2L)
expect_error(
  fit$fit$setState(badSavedState),
  "block 'variance.saved.vars' is malformed"
)

# ---- a keepTrees state stored before any recorded sweep is legal ----
# the unwritten slots hold the MULTIPLICATIVE identity 1.0: a 0-valued leaf
# would both annihilate the product predict forms and fail the positivity law
# validation applies to every saved variance tree.
preData <- data.frame(x = x, y = y)
preControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 4L,
  keepTrees = TRUE,
  updateState = FALSE
)
makePre <- function() {
  dbarts(
    y ~ x,
    preData,
    variance = varianceForest(n.trees = 5L),
    control = preControl
  )
}
preDonor <- makePre()
preDonor$storeState()
preDest <- makePre()
preDest$setState(preDonor$state)
preDest$storeState()
statesAgree(preDest$state, preDonor$state)

# ---- a wide (pooled) categorical predictor under a variance forest ----
# a variance tree splitting on a >63-level column keeps its rule's words in a
# side channel rather than the flat record, so flatten, rebuild, validation and
# replay each need it - and the saved buffer carries its own copy.
set.seed(21, sample.kind = "Rejection")
nWide <- 400L
wideData <- data.frame(
  g = factor(sample(80L, nWide, replace = TRUE)),
  z = runif(nWide)
)
wideSd <- ifelse(wideData$z < 0.5, 0.3, 1.2)
wideData$y <- as.numeric(wideData$g) / 40 + wideData$z + wideSd * rnorm(nWide)
fitWide <- bart(
  y ~ g + z,
  wideData,
  variance = varianceForest(n.trees = 8L),
  n.trees = 20L,
  n.samples = 15L,
  n.burn = 15L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
wideNew <- wideData[1:5, c("g", "z")]
predWide <- predict(fitWide, wideNew)
expect_true(all(attr(predWide, "s") > 0))
fitWide$fit$storeState()
wideFile <- tempfile()
saveRDS(fitWide, wideFile)
fitWideLoaded <- readRDS(wideFile)
predWideLoaded <- predict(fitWideLoaded, wideNew)
expect_identical(predWideLoaded, predWide)
expect_identical(attr(predWideLoaded, "s"), attr(predWide, "s"))
unlink(wideFile)

# ---- a homoscedastic truth does not manufacture heteroscedasticity ----
set.seed(11, sample.kind = "Rejection")
xHom <- runif(n)
yHom <- 2 * xHom + 0.8 * rnorm(n)
fitHom <- bart(
  xHom,
  yHom,
  variance = varianceForest(n.trees = 25L),
  n.trees = 50L,
  n.samples = 400L,
  n.burn = 400L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
sHom <- apply(fitHom$s.train, 2L, mean)
# no region should read as wildly more variable than another: the spread of the
# per-observation s(x) around its mean stays modest (no spurious structure)
expect_true(sd(sHom) / mean(sHom) < 0.35)
# and the recovered level is near the truth (0.8)
expect_true(mean(sHom) > 0.55 && mean(sHom) < 1.15)

# ---- a homoscedastic fit (no variance forest) carries no s channel ----
fitPlain <- bart(
  x,
  y,
  n.trees = 50L,
  n.samples = 100L,
  n.burn = 100L,
  n.chains = 1L,
  verbose = FALSE
)
expect_null(fitPlain$s.train)

# ---- the variance forest is gaussian only ----
expect_error(
  bart(
    x,
    as.integer(y > median(y)),
    variance = TRUE,
    family = "probit",
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "variance forest requires"
)

# ---- Student-t residuals: unadjudicated with 'variance' ----
# resid.dist is NSE (parsed in
# dbarts's own vocabulary), so these stay literal calls rather than do.call.
expect_error(
  bart(
    xHom,
    yHom,
    resid.dist = student(3),
    variance = TRUE,
    n.samples = 2L,
    n.burn = 2L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "variance forest does not support"
)
expect_inherits(
  bart(
    xHom,
    yHom,
    resid.dist = student(3),
    n.samples = 2L,
    n.burn = 2L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "bart"
)
# the variance forest alone still constructs
expect_inherits(
  dbarts(
    xHom,
    yHom,
    control = dbartsControl(
      n.chains = 1L,
      n.samples = 2L,
      n.burn = 2L,
      n.trees = 5L,
      n.threads = 1L,
      updateState = FALSE,
      seed = 7L
    ),
    variance = varianceForest(n.trees = 3L)
  ),
  "dbartsSampler"
)

# ---- the variance forest's own prior draw -----------------------------------
#
# sampleTreesFromPrior and sampleNodeParametersFromPrior are MEAN-forest
# entries by contract, so a heteroscedastic chain had no path to a prior-drawn
# s(x) before sampleVarianceForestFromPrior. What is pinned here is the R
# surface: the draw moves the scale surface, leaves it live (positive, and a
# state that restores), leaves the mean forest and sigma where it found them,
# and does nothing at all - not even a generator call - on a homoscedastic
# sampler.

set.seed(41L)
nPrior <- 200L
xPrior <- matrix(
  runif(nPrior * 2L),
  nPrior,
  2L,
  dimnames = list(NULL, c("a", "b"))
)
sPrior <- ifelse(xPrior[, 1L] < 0.5, 0.3, 1.5)
yPrior <- 2 * xPrior[, 2L] + sPrior * rnorm(nPrior)
priorControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 4L,
  updateState = FALSE,
  seed = 41L
)
priorSampler <- dbarts(
  xPrior,
  yPrior,
  test = xPrior[1:10, , drop = FALSE],
  variance = varianceForest(n.trees = 5L),
  control = priorControl
)
beforeDraw <- priorSampler$run(20L, 4L)
beforeTrees <- priorSampler$getTrees()
priorSampler$sampleVarianceForestFromPrior()
afterDraw <- priorSampler$run(0L, 1L)

# the surface moved, and it is still a legal scale surface
expect_true(
  !identical(
    afterDraw$variance[, 1L],
    beforeDraw$variance[, 4L]
  )
)
expect_true(all(afterDraw$variance > 0))
expect_true(all(afterDraw$varianceTest > 0))
# the mean forest and the fixed sigma are untouched by a variance-forest entry
expect_identical(
  beforeTrees[beforeTrees$sample == 4L, c("var", "value")],
  priorSampler$getTrees()[
    priorSampler$getTrees()$sample == 4L,
    c("var", "value")
  ]
)
expect_equal(afterDraw$sigma[[1L]], beforeDraw$sigma[[4L]])

# and the drawn state is live state: the validation a restore runs demands
# every variance leaf a positive scale and every bottom occupied
priorSampler$sampleVarianceForestFromPrior()
priorSampler$storeState()
drawnState <- priorSampler$state
restored <- dbarts(
  xPrior,
  yPrior,
  test = xPrior[1:10, , drop = FALSE],
  variance = varianceForest(n.trees = 5L),
  control = priorControl
)
expect_silent(restored$setState(drawnState))

# a homoscedastic sampler has nothing to draw: not a refusal, and not a
# generator call either, so the next draws are the ones it would have made
homoControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 4L,
  updateState = FALSE,
  seed = 41L
)
homo <- dbarts(xPrior, yPrior, control = homoControl)
homoBefore <- homo$run(5L, 4L)
homoAgain <- dbarts(xPrior, yPrior, control = homoControl)
expect_silent(homoAgain$sampleVarianceForestFromPrior())
expect_identical(homoAgain$run(5L, 4L)$train, homoBefore$train)

# ---- getVariance: the current variance surface, read without a run ----
# The accessor reports exactly what a run records as `variance` and
# `varianceTest`, at the state the read finds rather than at a kept sample, so
# a host driving the sampler one sweep at a time - or drawing the surface from
# its prior - reads s^2(x) here instead of through a recorded channel.

set.seed(23L)
nAcc <- 200L
xAcc <- matrix(
  runif(nAcc * 2L),
  nAcc,
  2L,
  dimnames = list(NULL, c("a", "b"))
)
sAcc <- ifelse(xAcc[, 1L] < 0.5, 0.3, 1.5)
yAcc <- 2 * xAcc[, 2L] + sAcc * rnorm(nAcc)
xAccTest <- xAcc[1:10, , drop = FALSE]
accControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 4L,
  updateState = FALSE,
  seed = 23L
)
accSampler <- dbarts(
  xAcc,
  yAcc,
  test = xAccTest,
  variance = varianceForest(n.trees = 5L),
  control = accControl
)
accRun <- accSampler$run(20L, 4L)

# the state a recorded sweep left: the accessor and the channel agree bitwise,
# one column per chain
expect_equal(dim(accSampler$getVariance()), c(nAcc, 1L))
expect_equal(dim(accSampler$getVariance(test = TRUE)), c(10L, 1L))
expect_identical(as.vector(accSampler$getVariance()), accRun$variance[, 4L])
expect_identical(
  as.vector(accSampler$getVariance(test = TRUE)),
  accRun$varianceTest[, 4L]
)

# several chains report per chain, the channel's own chain margin
multiControl <- dbartsControl(
  n.chains = 3L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 2L,
  updateState = FALSE,
  seed = 24L
)
multiSampler <- dbarts(
  xAcc,
  yAcc,
  test = xAccTest,
  variance = varianceForest(n.trees = 5L),
  control = multiControl
)
multiRun <- multiSampler$run(10L, 2L)
expect_identical(multiSampler$getVariance(), multiRun$variance[, 2L, ])
expect_identical(
  multiSampler$getVariance(test = TRUE),
  multiRun$varianceTest[, 2L, ]
)

# NULL exactly where the channels report nothing: no variance forest, and a
# test read with no test rows
homoAccessor <- dbarts(xAcc, yAcc, test = xAccTest, control = accControl)
expect_null(homoAccessor$getVariance())
expect_null(homoAccessor$getVariance(test = TRUE))
noTestSampler <- dbarts(
  xAcc,
  yAcc,
  variance = varianceForest(n.trees = 5L),
  control = accControl
)
noTestSampler$run(5L, 1L)
expect_null(noTestSampler$getVariance(test = TRUE))
expect_true(all(noTestSampler$getVariance() > 0))

# the test read REBUILDS: it is maintained only at a recorded sweep, so a
# test-predictor swap has to move it. The new rows are training rows, so the
# two reads must agree entry for entry.
accSampler$setTestPredictor(xAcc[21:30, , drop = FALSE])
expect_identical(
  as.vector(accSampler$getVariance(test = TRUE)),
  accSampler$getVariance()[21:30, 1L]
)

# a prior draw moves the surface, and the accessor - unlike predict(), which
# addresses saved samples - answers at the drawn trees
beforePriorDraw <- accSampler$getVariance()
accSampler$sampleVarianceForestFromPrior()
afterPriorDraw <- accSampler$getVariance()
expect_true(!identical(beforePriorDraw, afterPriorDraw))
expect_true(all(afterPriorDraw > 0))

# and the drawn factor is the calibrated one. One variance tree under a
# structure prior that practically never grows is a bare root, so the whole
# surface IS one leaf factor h, whose reciprocal is exactly
# chisq(nu) / (nu lambda^2) at the nu and lambda^2 the sigma prior is
# calibrated to: mean 1 and variance 2 / nu after scaling, so the band below is
# a closed-form standard error rather than a guess.
flatControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 1L,
  updateState = FALSE,
  seed = 25L
)
flatSampler <- dbarts(
  xAcc,
  yAcc,
  variance = varianceForest(n.trees = 1L, base = 1e-10),
  control = flatControl
)
# the seeded surface before any draw is the variance the calibration is stated
# against, so nothing here assumes a response transform
initialVariance <- flatSampler$getVariance()[1L]
residDf <- flatSampler$model@resid.prior@df
rawScale <- qchisq(1 - flatSampler$model@resid.prior@quantile, residDf) /
  residDf
leafScale <- initialVariance * rawScale
numLeafDraws <- 2000L
leafDraws <- numeric(numLeafDraws)
bareRoot <- TRUE
for (i in seq_len(numLeafDraws)) {
  flatSampler$sampleVarianceForestFromPrior()
  surface <- flatSampler$getVariance()
  # a bare root is one factor: every row carries it
  bareRoot <- bareRoot && all(surface == surface[1L])
  leafDraws[i] <- leafScale / surface[1L]
}
expect_true(bareRoot)
standardError <- sqrt(2 / (residDf * numLeafDraws))
expect_true(abs(mean(leafDraws) - 1) < 5 * standardError)
expect_true(abs(2 / var(leafDraws) - residDf) < 0.6)

# ---- samplePriorPredictive(type = "ppd") on a heteroscedastic sampler ----
# The noise is the drawn s(x) itself, read at the rows being predicted. With
# the leaf prior tightened the mean forest's own prior spread is negligible, so
# the ppd's variance IS the prior mean of s^2(x) - which the accessor reports
# directly, and against which the draws are scored here.
ppdControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 4L,
  updateState = FALSE,
  seed = 26L
)
ppdTest <- xAcc[1:4, , drop = FALSE]
ppdSampler <- dbarts(
  xAcc,
  yAcc,
  test = ppdTest,
  variance = varianceForest(n.trees = 5L),
  node.prior = normal(k = 40),
  control = ppdControl
)
set.seed(27L)
ppdDraws <- samplePriorPredictive(
  ppdSampler,
  x.test = ppdTest,
  n.samples = 500L,
  type = "ppd"
)
evDraws <- samplePriorPredictive(
  ppdSampler,
  x.test = ppdTest,
  n.samples = 500L,
  type = "ev"
)
expect_equal(dim(ppdDraws), c(500L, 4L))
expect_true(all(is.finite(ppdDraws)))

priorVariance <- replicate(500L, {
  ppdSampler$sampleVarianceForestFromPrior()
  ppdSampler$getVariance(test = TRUE)[, 1L]
})
# the mean forest contributes almost nothing at k = 40, and what remains is
# the drawn surface: a heavy-tailed mean, hence the factor-of-two band
expect_true(all(apply(evDraws, 2L, var) < 0.1 * apply(ppdDraws, 2L, var)))
varianceRatio <- apply(ppdDraws, 2L, var) / rowMeans(priorVariance)
expect_true(all(varianceRatio > 0.5 & varianceRatio < 2))

# the "ev" surface carries no noise at all, so it is untouched by the lift
expect_true(all(is.finite(evDraws)))
