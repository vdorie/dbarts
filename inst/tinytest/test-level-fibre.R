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

# the tests below say TRUE and FALSE where the control now says words
treeShiftWord <- function(x) if (x) "always" else "never"

# ---- the control slot ----

# three values: TRUE takes the step every iteration, FALSE never takes it,
# and NA - the default - takes it for a forest exactly where that forest's
# structural mixture is frozen
# (dbartsControl(treeShift = ) spells them "always", "never" and "auto"; the
# slot the bridge reads keeps the tri-state logical)
expect_true(is.na(dbarts::dbartsControl()@levelGibbs))
expect_true(is.na(dbarts::dbartsControl(treeShift = "auto")@levelGibbs))
expect_true(dbarts::dbartsControl(treeShift = "always")@levelGibbs)
expect_false(dbarts::dbartsControl(treeShift = "never")@levelGibbs)
expect_true(dbarts::dbartsControl(treeShift = "al")@levelGibbs)
# anything else, NA and a logical included, is refused by name
for (bad in list(NA, TRUE, "not-a-shift", c("always", "never"), 1L)) {
  expect_error(
    dbarts::dbartsControl(treeShift = bad),
    "'treeShift' must be one of"
  )
}
expect_error(dbarts::dbartsControl(levelGibbs = TRUE), "unused argument")

# the value rides the sampler's own control, and bart carries the formal
onControl <- dbarts::dbartsControl(
  treeShift = "always",
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
  pattern = "changing 'treeShift'"
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
# and the same refusal from the default, which is NA and not FALSE
turnedOn <- offSampler$control
expect_true(is.na(turnedOn@levelGibbs))
turnedOn@levelGibbs <- TRUE
expect_error(
  offSampler$setControl(turnedOn),
  pattern = "changing 'treeShift'"
)

# a control that only restates the value is accepted, NA against NA among
# them - identical() reads two missing values as the same setting
expect_null(onSampler$setControl(onSampler$control))
expect_null(offSampler$setControl(offSampler$control))

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
  dbarts::bart(
    testData$x,
    testData$y,
    control = dbarts::dbartsControl(treeShift = treeShiftWord(levelGibbs)),
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
# and where structure is being proposed - which is every bart fit that does
# not freeze it - the NA default takes no step, so naming FALSE changes
# nothing at all
fitDefault <- dbarts::bart(
  testData$x,
  testData$y,
  n.trees = 25L,
  n.samples = 60L,
  n.burn = 60L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 11L
)
expect_identical(fitDefault$yhat.train, fitOff$yhat.train)
expect_identical(fitOff$yhat.train, fitAt(FALSE)$yhat.train)
rm(fitOn, fitOff, fitDefault)

# ---- automatic: a frozen mixture switches the step on, and nothing else ----

# the mixture is mutable between samples while the levelGibbs slot is fixed at
# creation, so the decision is taken per sweep: a sampler grown under the
# shipped mixture takes no step, and takes one from the sweep its structures
# are frozen at
freeze <- function(control) {
  control@proposal.probs[
    c("birth_death", "swap", "change", "perturb", "rule_gibbs")
  ] <- 0
  control
}
frozenTail <- function(...) {
  control <- dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.burn = 0L,
    n.samples = 25L,
    updateState = FALSE,
    seed = 29L,
    ...
  )
  sampler <- dbarts::dbarts(y ~ x, testData, control = control)
  invisible(sampler$run(50L, 0L))
  sampler$setControl(freeze(sampler$control))
  sampler$run(0L, 25L)$train
}
# the two arms share every growing sweep - the default takes no step while
# structure is proposed, and neither does FALSE - so they stand at one forest
# when the freeze lands, and part only after it
frozenDefault <- frozenTail()
expect_false(identical(frozenDefault, frozenTail(treeShift = "never")))
expect_true(all(is.finite(frozenDefault)))

# frozen from creation instead, where TRUE and the default share every sweep
# rather than only the ones after the freeze, the default draws exactly what
# the step named on draws. It cannot be read against a grown forest from R:
# TRUE steps through the growing sweeps too, so the two arms would part
# before the freeze rather than at it
frozenThroughout <- function(...) {
  control <- dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.burn = 0L,
    n.samples = 25L,
    updateState = FALSE,
    seed = 31L,
    ...
  )
  control@proposal.probs <- c(
    birth_death = 0,
    swap = 0,
    change = 0,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  )
  sampler <- dbarts::dbarts(y ~ x, testData, control = control)
  sampler$run(25L, 25L)$train
}
alwaysFrozen <- frozenThroughout()
expect_identical(alwaysFrozen, frozenThroughout(treeShift = "always"))
expect_false(identical(alwaysFrozen, frozenThroughout(treeShift = "never")))
rm(frozenDefault, frozenTail, freeze, frozenThroughout, alwaysFrozen)

# ---- the answer is the same: held-out fits agree within Monte Carlo error ----

heldOut <- function(levelGibbs, seed) {
  fit <- dbarts::bart(
    testData$x[1:70, ],
    testData$y[1:70],
    test = testData$x[71:100, ],
    control = dbarts::dbartsControl(treeShift = treeShiftWord(levelGibbs)),
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
  dbarts::bart(
    y ~ x1 + x2 + x3,
    df,
    leaf.prior = linear("x2"),
    control = dbarts::dbartsControl(treeShift = treeShiftWord(levelGibbs)),
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

monotoneFit <- dbarts::bart(
  testData$x,
  testData$y,
  monotone = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
  keepTrees = TRUE,
  control = dbarts::dbartsControl(treeShift = "always"),
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

# ---- the setting lives on the control only ----

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# neither structure prior carries it: it was never per forest
expect_false("levelGibbs" %in% names(formals(dbarts:::cgm)))
expect_false("levelGibbs" %in% names(formals(dbarts:::dart)))
expect_false("treeShift" %in% names(formals(dbarts:::cgm)))
expect_error(dbarts::dbartsPriors$cgm(levelGibbs = TRUE), "unused argument")
expect_error(dbarts::dbartsPriors$dart(levelGibbs = TRUE), "unused argument")

# a value on the control reaches the bridge and takes the extra step: the
# draws move against a fit that declared nothing, and "auto" is the default
fitWithControl <- function(...) {
  dbarts::bart(
    testData$x,
    testData$y,
    control = dbarts::dbartsControl(...),
    n.trees = 15L,
    n.samples = 30L,
    n.burn = 30L,
    n.chains = 1L,
    n.threads = 1L,
    keepSampler = TRUE,
    verbose = FALSE,
    seed = 21L
  )
}
shiftDefault <- fitWithControl()
shiftOn <- fitWithControl(treeShift = "always")
shiftOff <- fitWithControl(treeShift = "never")
expect_true(is.na(shiftDefault$fit$control@levelGibbs))
expect_true(shiftOn$fit$control@levelGibbs)
expect_false(shiftOff$fit$control@levelGibbs)
expect_false(identical(shiftDefault$yhat.train, shiftOn$yhat.train))
expect_identical(shiftDefault$yhat.train, shiftOff$yhat.train)

# and from a DART prior, whose control is the same one
dartOn <- dbarts::bart(
  testData$x,
  testData$y,
  tree.prior = dbarts::dbartsPriors$dart(),
  control = dbarts::dbartsControl(treeShift = "always"),
  n.trees = 15L,
  n.samples = 30L,
  n.burn = 30L,
  n.chains = 1L,
  n.threads = 1L,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 21L
)
expect_true(dartOn$fit$control@levelGibbs)

rm(fitWithControl, shiftDefault, shiftOn, shiftOff, dartOn)
