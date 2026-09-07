# The level-fibre Gibbs step: the control surface, the settings it is fixed
# against, and what the step must and must not move. The distributional gate
# on the shift's own law - its mean, its covariance, and both poisons - is a
# tests/cpp one, since the shift is not a quantity any R channel reports.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

# ---- the control slot ----

expect_false(dbarts::dbartsControl()@levelGibbs)
expect_true(dbarts::dbartsControl(levelGibbs = TRUE)@levelGibbs)
expect_error(
  dbarts::dbartsControl(levelGibbs = NA),
  "'levelGibbs' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(levelGibbs = "not-a-logical"),
  "'levelGibbs' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(levelGibbs = c(TRUE, TRUE)),
  "'levelGibbs' must be of length 1"
)

# the value rides the sampler's own control, and bart2 carries the formal
onControl <- dbarts::dbartsControl(
  levelGibbs = TRUE,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L
)
onSampler <- dbarts::dbarts(y ~ x, testData, control = onControl)
expect_true(onSampler$control@levelGibbs)

# ---- fixed when the sampler is created ----

# the engine reads it at creation only, so a changed value must be refused
# rather than recorded R-side and dropped: $getPointer's re-creation branch
# rebuilds from the control, and the flag would come back on across a save
# and load
changed <- onSampler$control
changed@levelGibbs <- FALSE
expect_error(
  onSampler$setControl(changed),
  pattern = "changing 'levelGibbs'"
)
# and the same refusal from the off side
offSampler <- dbarts::dbarts(
  y ~ x,
  testData,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 5L
  )
)
turnedOn <- offSampler$control
turnedOn@levelGibbs <- TRUE
expect_error(
  offSampler$setControl(turnedOn),
  pattern = "changing 'levelGibbs'"
)

# a control that only restates the value is accepted
expect_null(onSampler$setControl(onSampler$control))

# ---- the serialized control carries it through a re-creation ----

invisible(onSampler$run(20L, 5L))
onSampler$storeState()
savedState <- onSampler$state
serialized <- tempfile(fileext = ".rds")
saveRDS(onSampler, serialized)
reloaded <- readRDS(serialized)
expect_true(reloaded$control@levelGibbs)
reloaded$storeState()
statesAgree(reloaded$state, savedState)
expect_silent(invisible(reloaded$run(0L, 1L)))
unlink(serialized)
rm(onSampler, offSampler, reloaded, savedState, changed, turnedOn, onControl)

# ---- a fit with the flag on runs, and its trees still sum to its fits ----

fitAt <- function(levelGibbs, ...) {
  dbarts::bart2(
    testData$x,
    testData$y,
    levelGibbs = levelGibbs,
    n.trees = 25L,
    n.samples = 60L,
    n.burn = 60L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 11L,
    ...
  )
}

fitOn <- fitAt(TRUE)
fitOff <- fitAt(FALSE)
expect_true(all(is.finite(fitOn$yhat.train)))
expect_true(all(is.finite(fitOn$sigma)))

# the precondition the shift has to keep: every recorded draw's fit is the
# sum over trees of the recorded leaf values, which is what predict replays.
# A shift written to the leaf tables while the reported fit came from a stale
# aggregate would part the two.
leafSumError <- function(fit) {
  max(abs(predict(fit, testData$x) - fit$yhat.train))
}
expect_true(leafSumError(fitOn) < 1e-10)
expect_true(leafSumError(fitOff) < 1e-10)

# the step is live: the same seed draws a different path with it on
expect_false(identical(fitOn$yhat.train, fitOff$yhat.train))
# and off, naming it changes nothing at all against the default
expect_identical(fitOff$yhat.train, fitAt(FALSE)$yhat.train)
rm(fitOn, fitOff)

# ---- the answer is the same: held-out fits agree within Monte Carlo error ----

heldOut <- function(levelGibbs, seed) {
  fit <- dbarts::bart2(
    testData$x[1:70, ],
    testData$y[1:70],
    test = testData$x[71:100, ],
    levelGibbs = levelGibbs,
    n.trees = 25L,
    n.samples = 500L,
    n.burn = 500L,
    n.chains = 4L,
    n.threads = 1L,
    verbose = FALSE,
    seed = seed
  )
  fit$yhat.test
}
# read against a sham of the same shape - the flag off at a different sampler
# seed - rather than against a standard error, since what separates two arms
# here is chain-to-chain Monte Carlo variation and not sampling noise the
# draws' own spread measures
worstGap <- function(a, b) {
  max(abs(colMeans(a) - colMeans(b)) / apply(b, 2L, sd))
}
reference <- heldOut(FALSE, 21L)
expect_true(
  worstGap(heldOut(TRUE, 21L), reference) <
    3 * worstGap(heldOut(FALSE, 31L), reference)
)
rm(reference, heldOut, worstGap)

# ---- out of scope, and inert: a linear leaf keeps no leaf table to shift ----

df <- data.frame(testData$x[, 1:3], y = testData$y)
names(df)[1:3] <- c("x1", "x2", "x3")
linearAt <- function(levelGibbs) {
  dbarts::bart2(
    y ~ x1 + x2 + x3,
    df,
    node.prior = linear("x2"),
    levelGibbs = levelGibbs,
    n.trees = 10L,
    n.samples = 40L,
    n.burn = 40L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 5L
  )
}
expect_identical(linearAt(TRUE)$yhat.train, linearAt(FALSE)$yhat.train)
rm(linearAt)

# ---- the monotone leaf is in scope, and the constraint survives ----

monotoneFit <- dbarts::bart2(
  testData$x,
  testData$y,
  monotone = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
  keepTrees = TRUE,
  levelGibbs = TRUE,
  n.trees = 20L,
  n.samples = 40L,
  n.burn = 40L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 13L
)
expect_true(all(is.finite(monotoneFit$yhat.train)))
# a shift is common to a tree's occupied leaves, so it leaves every
# within-tree difference alone and the ensemble stays monotone in column 4
grid <- testData$x[rep(1L, 25L), , drop = FALSE]
grid[, 4L] <- seq(0, 1, length.out = 25L)
along <- colMeans(predict(monotoneFit, grid))
expect_true(all(diff(along) >= -1e-8))
rm(monotoneFit, grid, along, df, fitAt, leafSumError, testData)
