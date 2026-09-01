# A heteroscedastic sampler's variance forest lives outside the forest vector
# the mutation helpers loop, so a predictor or cut-grid change used to leave
# s^2(x) routing observations by the predictors the forest was built with.
# The forced predictor paths and setCutPoints now re-route it.

set.seed(29, sample.kind = "Rejection")

n <- 200L
x <- cbind(x1 = runif(n), x2 = runif(n), x3 = runif(n))
y <- 2 * x[, 1L] + ifelse(x[, 2L] < 0.5, 0.3, 1.5) * rnorm(n)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 25L,
  n.burn = 0L,
  verbose = FALSE
)
buildVarianceSampler <- function(predictors, response) {
  dbarts::dbarts(
    predictors,
    response,
    control = control,
    variance = dbarts::varianceForest(n.trees = 10L)
  )
}

# a fractional varianceForest(vars=) entry is refused, naming the argument,
# rather than silently truncated (coerceOrError's integer branch); a
# whole-number double entry is accepted the same as an integer index
expect_error(
  dbarts::dbarts(
    x,
    y,
    control = control,
    variance = dbarts::varianceForest(vars = 2.9, n.trees = 10L)
  ),
  "'vars' must be a whole number; got '2.9'",
  fixed = TRUE
)
expect_silent(dbarts::dbarts(
  x,
  y,
  control = control,
  variance = dbarts::varianceForest(vars = 2, n.trees = 10L)
))

# distinct reported values, quantized well below any routing difference: the
# reported channels carry per-observation rounding of order 1e-16, which is not
# what these bounds are counting
numDistinct <- function(values) length(unique(signif(as.vector(values), 8)))

# ---- forced predictor swap: an all-identical design admits one value ----
# every row carries the same predictors, so every tree of either forest routes
# every row to the same leaf and both surfaces are constant. The bound is 1
# because the design has one distinct row, not because of any tolerance.
sampler <- buildVarianceSampler(x, y)
invisible(sampler$run(50L, 0L))
identicalRows <- matrix(rep(c(0.4, 0.55, 0.7), each = n), n, 3L)
sampler$setPredictor(identicalRows, forceUpdate = TRUE)
swapped <- sampler$run(0L, 1L)
expect_equal(numDistinct(swapped$variance), 1L)
# non-vacuity: the mean channel, which has always been re-routed, reaches the
# same bound on the same mutated sampler
expect_equal(numDistinct(swapped$train), 1L)
# and the bound is achievable, not merely small: a fresh fit on the
# post-mutation design reports one value too
fresh <- buildVarianceSampler(identicalRows, y)
expect_equal(numDistinct(fresh$run(50L, 1L)$variance), 1L)

# ---- setCutPoints: a coarsened grid caps the surface combinatorially ----
# 2 cuts on each of 3 columns quantizes the design into at most 3^3 = 27 cells,
# and s^2(x_i) is a function of row i's cell alone. A stale partition reports a
# value per training row.
coarse <- buildVarianceSampler(x, y)
invisible(coarse$run(50L, 0L))
coarse$setCutPoints(rep(list(c(1 / 3, 2 / 3)), 3L), 1:3)
quantized <- coarse$run(0L, 1L)
expect_true(numDistinct(quantized$variance) <= 27L)
expect_true(numDistinct(quantized$train) <= 27L)

# every serialized variance threshold must be a member of the new grid: the
# rebuild matches each flattened cut value against the current cut points
# exactly, so a state carrying a stale threshold cannot restore
coarse$storeState()
expect_silent(coarse$setState(coarse$state))
expect_true(numDistinct(coarse$run(0L, 1L)$variance) <= 27L)

# ---- statistical agreement with a from-scratch fit: a repair that re-routes
# ---- but scatters the wrong factors fails here
set.seed(37, sample.kind = "Rejection")
nAgree <- 300L
xOld <- cbind(runif(nAgree), runif(nAgree))
xNew <- cbind(runif(nAgree), runif(nAgree))
yAgree <- 2 * xNew[, 1L] + ifelse(xNew[, 2L] < 0.5, 0.4, 1.6) * rnorm(nAgree)
mutated <- buildVarianceSampler(xOld, yAgree)
invisible(mutated$run(50L, 0L))
mutated$setPredictor(xNew, forceUpdate = TRUE)
sMutated <- sqrt(apply(mutated$run(100L, 100L)$variance, 1L, mean))
scratch <- buildVarianceSampler(xNew, yAgree)
sScratch <- sqrt(apply(scratch$run(150L, 100L)$variance, 1L, mean))
expect_true(cor(sMutated, sScratch) > 0.9)
expect_true(abs(mean(sMutated) / mean(sScratch) - 1) < 0.15)

# ---- setData: a replacement data set, observation count included ----
# Seven n-sized allocations are pinned at creation n - meanWeights_ and the
# variance forest's indexBuffer, factorByTree, combinedVariance, meanResidual,
# divisor and treeResidual - so a changed count used to overrun them. setData
# itself returned cleanly; the fault landed on the FOLLOWING run, which
# segfaulted in 4 tries out of 5 and reported a non-finite variance in the
# fifth. Assert the run completes and every reported value is finite.
set.seed(41, sample.kind = "Rejection")
nSmall <- 200L
nLarge <- 5000L
xSmall <- cbind(runif(nSmall), runif(nSmall))
ySmall <- 2 *
  xSmall[, 1L] +
  ifelse(xSmall[, 2L] < 0.5, 0.3, 1.5) * rnorm(nSmall)
xLarge <- cbind(runif(nLarge), runif(nLarge))
yLarge <- 2 *
  xLarge[, 1L] +
  ifelse(xLarge[, 2L] < 0.5, 0.3, 1.5) * rnorm(nLarge)
resized <- buildVarianceSampler(xSmall, ySmall)
invisible(resized$run(50L, 0L))
resized$setData(dbarts::dbartsData(xLarge, yLarge))
grownRun <- resized$run(0L, 5L)
expect_equal(nrow(grownRun$variance), nLarge)
expect_true(all(is.finite(grownRun$variance)))
expect_true(all(grownRun$variance > 0))
resized$setData(dbarts::dbartsData(xSmall, ySmall))
shrunkRun <- resized$run(0L, 5L)
expect_equal(nrow(shrunkRun$variance), nSmall)
expect_true(all(is.finite(shrunkRun$variance)))
expect_true(all(shrunkRun$variance > 0))

# the routing half, isolated at a FIXED count so no buffer moves: replacing the
# design with an all-identical one admits exactly one reported value
identicalDesign <- matrix(rep(c(0.35, 0.65), each = nSmall), nSmall, 2L)
replaced <- buildVarianceSampler(xSmall, ySmall)
invisible(replaced$run(50L, 0L))
replaced$setData(dbarts::dbartsData(identicalDesign, ySmall))
replacedRun <- replaced$run(0L, 1L)
expect_equal(numDistinct(replacedRun$variance), 1L)
expect_equal(numDistinct(replacedRun$train), 1L)

# ---- statistical agreement with a from-scratch fit on the replacement data
set.seed(43, sample.kind = "Rejection")
nDataOld <- 200L
nDataNew <- 400L
xDataOld <- cbind(runif(nDataOld), runif(nDataOld))
yDataOld <- 2 *
  xDataOld[, 1L] +
  ifelse(xDataOld[, 2L] < 0.5, 0.4, 1.6) * rnorm(nDataOld)
xDataNew <- cbind(runif(nDataNew), runif(nDataNew))
yDataNew <- 2 *
  xDataNew[, 1L] +
  ifelse(xDataNew[, 2L] < 0.5, 0.4, 1.6) * rnorm(nDataNew)
swappedData <- buildVarianceSampler(xDataOld, yDataOld)
invisible(swappedData$run(50L, 0L))
swappedData$setData(dbarts::dbartsData(xDataNew, yDataNew))
sSwapped <- sqrt(apply(swappedData$run(100L, 100L)$variance, 1L, mean))
scratchData <- buildVarianceSampler(xDataNew, yDataNew)
sScratchData <- sqrt(apply(scratchData$run(150L, 100L)$variance, 1L, mean))
expect_true(cor(sSwapped, sScratchData) > 0.9)
expect_true(abs(mean(sSwapped) / mean(sScratchData) - 1) < 0.15)

# ---- the transactional predictor surface, one entry at a time: the two-phase
# revalidation and the per-observation session's cell guard now reach the
# variance forest, so a transactional change either re-routes every variance
# tree with the mean ones or is refused - the whole transaction rolling back,
# or the row declining - and never installs a partition s^2(x) would misroute.
# The forced entries and setCutPoints reach it as before (the sections above
# measure that they do). Its own sampler, at the end of the file, so no
# assertion above sees a different rng stream.
set.seed(47, sample.kind = "Rejection")
nPin <- 200L
xPin <- cbind(x1 = runif(nPin), x2 = runif(nPin), x3 = runif(nPin))
yPin <- 2 * xPin[, 1L] + ifelse(xPin[, 2L] < 0.5, 0.3, 1.5) * rnorm(nPin)
pinned <- buildVarianceSampler(xPin, yPin)
invisible(pinned$run(50L, 0L))
# a jitter small enough that no leaf of any tree - mean or variance - empties,
# so the ACCEPT arm is what runs
xPinJitter <- pmin(
  pmax(xPin + matrix(rnorm(nPin * 3L, 0, 0.005), nPin), 0),
  1
)
# the joint entry resolves its column by NAME; runs before the whole-matrix
# replacement below, which now keeps data@x's dimnames, so both orderings
# would resolve - this one is simplest to seed
installedJointly <- dbarts::updatePredictorPerObservationJointly(
  pinned,
  xPinJitter[, 3L],
  "x3"
)
expect_true(is.logical(installedJointly) && length(installedJointly) == nPin)
expect_true(any(installedJointly))
installed <- pinned$setPredictor(
  xPinJitter[, 2L],
  column = 2L,
  forceUpdate = "partial"
)
expect_true(is.logical(installed) && length(installed) == nPin)
expect_true(any(installed))
expect_true(pinned$setPredictor(
  xPinJitter[, 1L],
  column = 1L,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
))
expect_true(pinned$setPredictor(xPinJitter, forceUpdate = FALSE))
# the DECLINE arm, beside each: a two-level replacement column empties leaves,
# so the whole transaction rolls back and the per-row session declines the rows
# that would empty one
xPinTwoLevel <- ifelse(seq_len(nPin) %% 2L == 0L, 0.25, 0.75)
expect_false(pinned$setPredictor(
  xPinTwoLevel,
  column = 1L,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
))
expect_true(any(
  !pinned$setPredictor(
    xPinTwoLevel,
    column = 1L,
    forceUpdate = "partial"
  )
))
# ... while the forced entries and the cut-grid change stay supported
xPinNew <- xPin * 0.9 + 0.05
expect_silent(pinned$setPredictor(xPinNew, forceUpdate = TRUE))
expect_silent(pinned$setPredictor(xPin[, 1L], column = 1L, forceUpdate = TRUE))
expect_silent(pinned$setCutPoints(rep(list(c(1 / 3, 2 / 3)), 3L), 1:3))
pinnedRun <- pinned$run(0L, 5L)
expect_true(all(is.finite(pinnedRun$variance)))
expect_true(all(pinnedRun$variance > 0))

# ---- setState now validates variance tree occupancy ----
# The transactional entries above install a row only if it empties no leaf of
# any tree, the variance forest included. setState imposed that criterion on
# every mean tree and on no variance tree, so the route the old refusal message
# recommended admitted states the veto refuses: a variance bottom no row reaches
# carries no drawn scale. It is checked now, which is a behavior CHANGE - such
# a state used to install.
set.seed(53, sample.kind = "Rejection")
stateSampler <- buildVarianceSampler(xPin, yPin)
# an explicit grid, so the hand-built splits below sit on known cut values
stateSampler$setCutPoints(list(c(0.25, 0.5, 0.75)), 1L)
invisible(stateSampler$run(20L, 0L))
stateSampler$storeState()
# replace variance tree 1 with a five-node tree: an ordinal split, then a
# NESTED split on the same column. The pre-order layout is (root, root's left
# child, that child's two leaves, root's right leaf), so the second split
# lands inside the first's left branch.
replaceFirstVarianceTree <- function(state, values) {
  chainState <- state[[1L]]
  sizes <- chainState$variance.sizes
  keptNodes <- -seq_len(sizes[1L])
  keptBytes <- -seq_len(sizes[1L] * 8L)
  chainState$variance.vars <- c(
    c(1L, 1L, -1L, -1L, -1L),
    chainState$variance.vars[keptNodes]
  )
  chainState$variance.values <- c(
    writeBin(values, raw()),
    chainState$variance.values[keptBytes]
  )
  # bits 1-2 of the flag byte hold the node kind; 2 is an ordinal split
  chainState$variance.flags <- c(
    as.raw(c(2L, 2L, 0L, 0L, 0L)),
    chainState$variance.flags[keptNodes]
  )
  sizes[1L] <- 5L
  chainState$variance.sizes <- sizes
  state[[1L]] <- chainState
  state
}
# x1 <= 0.25 and then x1 > 0.75: a region no row can reach
strandedState <- replaceFirstVarianceTree(
  stateSampler$state,
  c(0.25, 0.75, 1.3, 0.7, 1.1)
)
expect_error(stateSampler$setState(strandedState), "not consistent")
# non-vacuity: the same hand-built shape with the nesting the other way round
# leaves every bottom occupied and installs
occupiedState <- replaceFirstVarianceTree(
  stateSampler$state,
  c(0.5, 0.25, 1.3, 0.7, 1.1)
)
expect_silent(stateSampler$setState(occupiedState))
stateRun <- stateSampler$run(0L, 3L)
expect_true(all(is.finite(stateRun$variance)))
expect_true(all(stateRun$variance > 0))

rm(
  n,
  x,
  y,
  control,
  buildVarianceSampler,
  numDistinct,
  sampler,
  identicalRows,
  swapped,
  fresh,
  coarse,
  quantized,
  nAgree,
  xOld,
  xNew,
  yAgree,
  mutated,
  sMutated,
  scratch,
  sScratch,
  nSmall,
  nLarge,
  xSmall,
  ySmall,
  xLarge,
  yLarge,
  resized,
  grownRun,
  shrunkRun,
  identicalDesign,
  replaced,
  replacedRun,
  nDataOld,
  nDataNew,
  xDataOld,
  yDataOld,
  xDataNew,
  yDataNew,
  swappedData,
  sSwapped,
  scratchData,
  sScratchData,
  nPin,
  xPin,
  yPin,
  pinned,
  xPinJitter,
  installed,
  installedJointly,
  xPinTwoLevel,
  xPinNew,
  pinnedRun,
  stateSampler,
  replaceFirstVarianceTree,
  strandedState,
  occupiedState,
  stateRun
)
