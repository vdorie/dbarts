# The saved state and the case weights. A state carries the latents and never
# the weights or the censoring status, and setState installs the stored
# latents whatever weights the destination holds. Like setWeights, it draws
# nothing: the next sweep draws the latents first, against the weights in
# force. So a restore under other weights is the restore under the stored ones
# followed by setWeights, and a weight change is undone by setWeights and
# setState in either order, bit for bit.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(9021L, sample.kind = "Rejection")
n <- 81L
x <- matrix(runif(n * 2L), n)
f <- 2.5 * x[, 1L] - 1.2
yBinary <- as.double(rbinom(n, 1L, plogis(f)))
yContinuous <- as.double(f + rnorm(n))

# equal in mean and in sum of squares, different in bytes
wA <- rep(c(1, 4, 4), length.out = n)
wB <- rep(c(2, 2, 5), length.out = n)

stateControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 902L
)
logisticAt <- function(weights) {
  dbarts::dbarts(
    x,
    yBinary,
    weights = weights,
    family = "logistic",
    control = stateControl,
    n.trees = 20L
  )
}
gaussianAt <- function(weights) {
  dbarts::dbarts(
    x,
    yContinuous,
    weights = weights,
    control = stateControl,
    n.trees = 20L
  )
}
# zeros go in through $setWeights, creation warning about them
studentAt <- function(weights) {
  sampler <- dbarts::dbarts(
    x,
    yContinuous,
    weights = wA,
    control = stateControl,
    n.trees = 20L,
    family = dbarts:::student(5)
  )
  sampler$setWeights(weights)
  sampler
}
storedFrom <- function(sampler) {
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler$state
}
storedLatents <- function(state) state[[1L]][["latents"]]

# --- what a stored state carries -------------------------------------------
# nothing about the weights or the censoring status
state.wA <- storedFrom(logisticAt(wA))
expect_null(attr(state.wA, "weights.digest"))
expect_null(attr(state.wA, "weights.zero"))
expect_null(attr(state.wA, "survival.digest"))
expect_identical(attr(state.wA, "formatVersion"), 1L)

# --- a restore installs the stored latents, silently -----------------------
for (weights in list(wA, wB)) {
  restored <- logisticAt(weights)
  expect_identical(length(captureWarnings(restored$setState(state.wA))), 0L)
  expect_identical(restored$getLatents(), storedLatents(state.wA))
}

# the continuation contract: a restored chain agrees with the uninterrupted
# one to within the documented reconstructions, and two restores of one state
# are bitwise each other
for (make in list(logisticAt, gaussianAt)) {
  source.wA <- make(wA)
  donor <- storedFrom(source.wA)
  restored <- make(wA)
  restored$setState(donor)
  continued <- source.wA$run(0L, 5L)$train
  ahead <- restored$run(0L, 5L)$train
  expect_true(max(abs(continued - ahead)) <= 1e-13)
  again <- make(wA)
  again$setState(donor)
  expect_identical(again$run(0L, 5L)$train, ahead)
}

# --- a restore under other weights is the restore then setWeights ----------
# on the families whose latents read the weights, with rows at weight zero on
# Student-t entering and leaving the likelihood
for (arm in list(
  list(make = logisticAt, stored = wA, other = wB),
  list(make = studentAt, stored = replace(wA, 1:20, 0), other = wB),
  list(make = studentAt, stored = wA, other = replace(wB, 62:81, 0))
)) {
  donor <- storedFrom(arm$make(arm$stored))
  liveArm <- arm$make(arm$stored)
  liveArm$setState(donor)
  liveArm$setWeights(arm$other)
  stateArm <- arm$make(arm$other)
  # other weights leave the chains the stored ones
  expect_true(stateArm$setState(donor))
  expect_identical(stateArm$getLatents(), storedLatents(donor))
  expect_identical(stateArm$getLatents(), liveArm$getLatents())
  expect_identical(stateArm$run(0L, 3L)$train, liveArm$run(0L, 3L)$train)
}

# --- run, swap, save, load --------------------------------------------------
# $setWeights does not refresh $state, so a save and load revives through the
# state stored under the previous weights; as the swap drew nothing, that is
# the live sampler
headline <- logisticAt(wA)
invisible(headline$run(10L, 5L))
headline$storeState()
preSwap <- headline$getLatents()
headline$setWeights(wB)
expect_identical(headline$getLatents(), preSwap)
headlineFile <- tempfile()
saveRDS(headline, file = headlineFile)
revived <- readRDS(headlineFile)
unlink(headlineFile)
expect_identical(revived$getLatents(), preSwap)
ahead <- headline$run(0L, 3L)$train
expect_true(max(abs(revived$run(0L, 3L)$train - ahead)) <= 1e-13)

# --- undoing a weight change ------------------------------------------------
# Store, propose other weights, sweep, go back by the old weights and the
# stored state in either order. A sampler that installed the state and never
# proposed is the reference.
for (arm in list(
  list(make = logisticAt, stored = wA, proposed = wB),
  list(make = studentAt, stored = replace(wA, 1:20, 0), proposed = wB),
  list(make = studentAt, stored = wA, proposed = replace(wB, 62:81, 0))
)) {
  reference <- arm$make(arm$stored)
  reference$setState(storedFrom(reference))
  reference <- reference$run(0L, 5L)
  for (weightsFirst in c(TRUE, FALSE)) {
    sampler <- arm$make(arm$stored)
    kept <- storedFrom(sampler)
    sampler$setWeights(arm$proposed)
    invisible(sampler$run(0L, 3L))
    if (weightsFirst) {
      sampler$setWeights(arm$stored)
      expect_true(sampler$setState(kept))
    } else {
      expect_true(sampler$setState(kept))
      sampler$setWeights(arm$stored)
    }
    expect_identical(sampler$getLatents(), storedLatents(kept))
    expect_identical(sampler$run(0L, 5L), reference)
  }
}

# --- another censoring status -----------------------------------------------
# A state stored with censored rows goes into a sampler that scores some of
# them events and some events censored. The install draws nothing: an event
# row holds its observed log-time, a censored row the stored latent (a row
# that was an event sits at its bound), and the chain runs as one that
# installed the state under its own status and was then handed the new one.
time <- exp(yContinuous / 3)
bound <- as.double(quantile(time, 0.7))
logTime <- log(pmin(time, bound))
statusStored <- as.double(time <= bound)
statusOther <- statusStored
statusOther[which(statusStored == 0)[1:5]] <- 1
statusOther[which(statusStored == 1)[1:5]] <- 0
aftAt <- function(status) {
  dbarts::dbarts(
    x,
    cbind(pmin(time, bound), status),
    family = "aft",
    control = stateControl,
    n.trees = 20L
  )
}
aftState <- storedFrom(aftAt(statusStored))
installed <- aftAt(statusOther)
expect_true(installed$setState(aftState))
latents <- installed$getLatents()
expect_equal(latents[statusOther == 1], logTime[statusOther == 1])
expect_identical(
  latents[statusOther == 0],
  storedLatents(aftState)[statusOther == 0]
)
handed <- aftAt(statusStored)
handed$setState(aftState)
handed$setResponse(handed$data@y, status = statusOther)
expect_identical(installed$run(0L, 3L)$train, handed$run(0L, 3L)$train)

# --- the warm start ---------------------------------------------------------
# installTrees reads the forests alone and leaves the destination's latents
warmDonor <- logisticAt(wA)
invisible(warmDonor$run(10L, 5L))
warmDonor$storeState()
warmDest <- logisticAt(wB)
invisible(warmDest$run(5L, 2L))
beforeWarm <- warmDest$getLatents()
warmDest$installTrees(warmDonor)
expect_identical(warmDest$getLatents(), beforeWarm)

rm(
  state.wA,
  weights,
  restored,
  make,
  source.wA,
  donor,
  continued,
  ahead,
  again,
  arm,
  liveArm,
  stateArm,
  headline,
  preSwap,
  headlineFile,
  revived,
  reference,
  weightsFirst,
  sampler,
  kept,
  time,
  bound,
  logTime,
  statusStored,
  statusOther,
  aftState,
  installed,
  latents,
  handed,
  warmDonor,
  warmDest,
  beforeWarm
)
