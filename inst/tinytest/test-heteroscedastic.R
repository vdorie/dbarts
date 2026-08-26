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

fit <- bart2(
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
fitWide <- bart2(
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
fitHom <- bart2(
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
fitPlain <- bart2(
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
  bart2(
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

# ---- Student-t residuals and a grouped fit: unadjudicated with 'variance' ----
# resid.dist is NSE (parsed in
# dbarts's own vocabulary), so these stay literal calls rather than do.call.
expect_error(
  bart2(
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
  bart2(
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
g <- factor(rep(c("a", "b"), length.out = n))
opts2 <- list(
  n.samples = 1L,
  n.burn = 0L,
  n.thin = 1L,
  n.threads = 1L,
  verbose = FALSE
)
# rbart_vi declares no 'variance' formal, and no dots channel to name it
# through: R's own wall
expect_error(
  do.call(rbart_vi, c(list(yHom ~ xHom, group.by = g, variance = TRUE), opts2)),
  "unused argument"
)
expect_inherits(
  do.call(rbart_vi, c(list(yHom ~ xHom, group.by = g, n.trees = 10L), opts2)),
  "rbart"
)

# ---- grouped random effects + a variance forest: one model reached from
# either spelling, refused as an unadjudicated composition (D6) ----
groupedControl <- dbartsControl(
  n.chains = 1L,
  n.samples = 2L,
  n.burn = 2L,
  n.trees = 5L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 7L
)
attr(groupedControl, "bartcore.groups") <- list(
  indices = as.integer(rep(seq_len(4L), length.out = n)),
  n.groups = 4L,
  prior = "cauchy",
  rel.scale = sd(yHom),
  n.steps = 1L
)

# entrance 1: dbarts(), the group attribute already on the control the
# variance formal resolves against - the R check fires
expect_error(
  dbarts(
    xHom,
    yHom,
    control = groupedControl,
    variance = varianceForest(n.trees = 3L)
  ),
  "does not support grouped random effects"
)
# entrance 2: dbartsSpec(), same control, same R-side check
expect_error(
  dbartsSpec(
    dbartsData(xHom, yHom),
    control = groupedControl,
    variance = varianceForest(n.trees = 3L)
  ),
  "does not support grouped random effects"
)
# entrance 3: new("dbartsSampler", ...), built from a dbartsSpec() result
# with the group attribute added AFTERWARDS - resolveSamplerSpec never sees
# it (the other order: variance resolves first, the group attribute arrives
# after), so this is the createHolder backstop's own message
acquiredSpec <- dbartsSpec(
  dbartsData(xHom, yHom),
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
)
attr(acquiredSpec$control, "bartcore.groups") <-
  attr(groupedControl, "bartcore.groups")
expect_error(
  new(
    "dbartsSampler",
    acquiredSpec$control,
    acquiredSpec$model,
    acquiredSpec$data
  ),
  "not supported with a heteroscedastic variance forest"
)

# each half ALONE still constructs: this is a composition check, not a
# regression on either feature
expect_inherits(
  dbarts(xHom, yHom, control = groupedControl),
  "dbartsSampler"
)
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
