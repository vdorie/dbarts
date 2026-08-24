# The saved state and the case weights. The weights themselves do not ride a
# state; a digest of the ones in force when it was stored does, and setState
# re-derives the weight-dependent latents against the DESTINATION's weights
# when the two disagree - so a state install lands where the same setWeights
# call would, rather than pairing one vector's latents with another's counts.
# When the two agree, the restore is the identity it has always been.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(9021L, sample.kind = "Rejection")
n <- 81L # divisible by 3, so the colliding pair below repeats whole
x <- matrix(runif(n * 2L), n)
f <- 2.5 * x[, 1L] - 1.2
yBinary <- as.double(rbinom(n, 1L, plogis(f)))
yContinuous <- as.double(f + rnorm(n))
yCategory <- factor(letters[1L + (seq_len(n) %% 3L)])

# 1,4,4 against 2,2,5: equal in mean and in sum of squares, different in
# bytes. Every pin that separates them below also fails a moment digest.
wA <- rep(c(1, 4, 4), length.out = n)
wB <- rep(c(2, 2, 5), length.out = n)
expect_identical(c(sum(wA), sum(wA^2)), c(sum(wB), sum(wB^2)))
expect_false(identical(wA, wB))

stateControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  seed = 902L
)
logisticSampler <- function(weights) {
  dbarts::dbarts(
    x,
    yBinary,
    weights = weights,
    family = "logistic",
    control = stateControl
  )
}
gaussianSampler <- function(weights) {
  dbarts::dbarts(x, yContinuous, weights = weights, control = stateControl)
}
studentSampler <- function(weights) {
  dbarts::dbarts(
    x,
    yContinuous,
    weights = weights,
    control = stateControl,
    resid.dist = dbarts:::student(5)
  )
}
varianceSampler <- function(weights) {
  dbarts::dbarts(
    x,
    yContinuous,
    weights = weights,
    control = stateControl,
    variance = dbarts::varianceForest(n.trees = 10L)
  )
}
multinomialSampler <- function() {
  dbarts::dbarts(x, yCategory, family = "multinomial", control = stateControl)
}
storedFrom <- function(sampler) {
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler$state
}

# --- what a stored state carries -------------------------------------------
# an ADDITIVE top-level attribute: 8 raw bytes beside cutPoints, and the
# encoding version does not move for it
source.wA <- logisticSampler(wA)
state.wA <- storedFrom(source.wA)
digest <- attr(state.wA, "weights.digest")
expect_true(is.raw(digest))
expect_identical(length(digest), 8L)
expect_identical(attr(state.wA, "formatVersion"), 3L)

# --- matched round trip: the identity --------------------------------------
# the destination's weights ARE the ones the stored latents were shaped by, so
# nothing is re-derived and the stored latents install unchanged
matched <- logisticSampler(wA)
matched$setState(state.wA)
expect_identical(matched$getLatents(), state.wA[[1L]][["latents"]])

# and the restore is silent: no warning fires on either the matched or the
# mismatched path
quietMatched <- captureWarnings(logisticSampler(wA)$setState(state.wA))
quietMismatched <- captureWarnings(logisticSampler(wB)$setState(state.wA))
expect_identical(length(quietMatched), 0L)
expect_identical(length(quietMismatched), 0L)

# the continuation contract is untouched: a restored chain agrees with the
# uninterrupted one to within the two documented reconstructions, and two
# restores of one state are bitwise each other
continued <- source.wA$run(0L, 5L)$train
restored <- matched$run(0L, 5L)$train
expect_true(max(abs(continued - restored)) <= 1e-13)
restoredAgain <- logisticSampler(wA)
restoredAgain$setState(state.wA)
expect_identical(restoredAgain$run(0L, 5L)$train, restored)

# gaussian, the same three
gsource.wA <- gaussianSampler(wA)
gstate.wA <- storedFrom(gsource.wA)
gmatched <- gaussianSampler(wA)
gmatched$setState(gstate.wA)
gcontinued <- gsource.wA$run(0L, 5L)$train
grestored <- gmatched$run(0L, 5L)$train
expect_true(max(abs(gcontinued - grestored)) <= 1e-13)
gagain <- gaussianSampler(wA)
gagain$setState(gstate.wA)
expect_identical(gagain$run(0L, 5L)$train, grestored)

# --- the repair: a mismatched restore lands where setWeights lands ----------
# arm 1 restores into the weights the state was stored under and then swaps
# through the LIVE conduit; arm 2 restores into the swapped weights and lets
# the seam reconcile. The two are the same sampler, bitwise.
liveArm <- logisticSampler(wA)
liveArm$setState(state.wA)
liveArm$setWeights(wB)
stateArm <- logisticSampler(wB)
stateArm$setState(state.wA)
expect_identical(stateArm$getLatents(), liveArm$getLatents())
expect_identical(stateArm$run(0L, 3L)$train, liveArm$run(0L, 3L)$train)

# not vacuous: what the repair installs is NOT what the state carries, and the
# stored latents are what a destination at the matching weights gets
expect_false(identical(stateArm$getLatents(), state.wA[[1L]][["latents"]]))

# --- the headline: run, swap, save, load ------------------------------------
# $setWeights does not refresh $state, so an ordinary save/load revives
# through a state stored under the PREVIOUS weights. The revived latents are
# the post-swap ones, not the pre-swap ones rewound.
headline <- logisticSampler(wA)
invisible(headline$run(10L, 5L))
headline$storeState()
preSwap <- headline$getLatents()
headline$setWeights(wB)
liveSwap <- headline$getLatents()
headlineFile <- tempfile()
saveRDS(headline, file = headlineFile)
rm(headline)
revived <- readRDS(headlineFile)
unlink(headlineFile)
revivedLatents <- revived$getLatents()
expect_false(identical(revivedLatents, preSwap))
expect_true(max(abs(revivedLatents - liveSwap)) <= 1e-13)
expect_true(max(abs(revivedLatents - preSwap)) > 1e-3)

# --- neutrality: the repair is a measured no-op elsewhere -------------------
# Every family below takes a MISMATCHED transplant twice - once with the
# digest present (the repair fires) and once with it stripped (it cannot) -
# and runs byte-identically either way. Student-t is pinned on the NON-UNIFORM
# pair: a uniform w = 1 against w = 8 is invariant by construction (sigma^2
# rescales) and would pin nothing.
withoutDigest <- function(state) {
  attr(state, "weights.digest") <- NULL
  state
}
for (make in list(gaussianSampler, studentSampler, varianceSampler)) {
  donor <- storedFrom(make(wA))
  repaired <- make(wB)
  repaired$setState(donor)
  inert <- make(wB)
  inert$setState(withoutDigest(donor))
  expect_identical(repaired$getLatents(), inert$getLatents())
  expect_identical(repaired$run(0L, 3L)$train, inert$run(0L, 3L)$train)
}

# a multinomial sampler refuses weights at creation, so its digest can only be
# made to differ by hand; the repair still fires and is still inert, because
# the family reads no weights to re-derive anything from
mstate <- storedFrom(multinomialSampler())
mforced <- mstate
attr(mforced, "weights.digest") <- as.raw(rep(255L, 8L))
mrepaired <- multinomialSampler()
mrepaired$setState(mforced)
minert <- multinomialSampler()
minert$setState(mstate)
expect_identical(mrepaired$run(0L, 3L)$train, minert$run(0L, 3L)$train)

# no weights and rep(1, n) are the same sampler - to a logistic count,
# lround(1) either way - so they must digest the same and NOT repair each
# other. A null token of its own would fire here.
nullSource <- dbarts::dbarts(
  x,
  yBinary,
  family = "logistic",
  control = stateControl
)
nullState <- storedFrom(nullSource)
onesDest <- logisticSampler(rep(1, n))
onesDest$setState(nullState)
expect_identical(onesDest$getLatents(), nullState[[1L]][["latents"]])
nullDest <- dbarts::dbarts(
  x,
  yBinary,
  family = "logistic",
  control = stateControl
)
nullDest$setState(withoutDigest(nullState))
expect_identical(onesDest$run(0L, 3L)$train, nullDest$run(0L, 3L)$train)

# --- additivity and refusal -------------------------------------------------
# a state with NO digest - one written before the attribute existed - loads on
# every family and restores exactly what it did before it existed
for (make in list(logisticSampler, gaussianSampler, studentSampler)) {
  donor <- storedFrom(make(wA))
  stripped <- make(wA)
  stripped$setState(withoutDigest(donor))
  expect_identical(stripped$getLatents(), donor[[1L]][["latents"]])
  expect_identical(attr(donor, "formatVersion"), 3L)
}

# present but not 8 raw bytes is named, in the malformed-block voice
malformed <- state.wA
attr(malformed, "weights.digest") <- as.raw(1:4)
expect_error(
  logisticSampler(wA)$setState(malformed),
  pattern = "malformed weights digest in bartcore state"
)
attr(malformed, "weights.digest") <- 1:8
expect_error(
  logisticSampler(wA)$setState(malformed),
  pattern = "malformed weights digest in bartcore state"
)

# --- the warm start is exempt ----------------------------------------------
# installTrees reads version and forests, never latents: a warm start
# transplants across data by design, so a weight difference is not a mismatch
# to reconcile and the destination's own latents are left where they were
warmDonor <- logisticSampler(wA)
invisible(warmDonor$run(10L, 5L))
warmDonor$storeState()
warmDest <- logisticSampler(wB)
invisible(warmDest$run(5L, 2L))
beforeWarm <- warmDest$getLatents()
warmDest$installTrees(warmDonor)
expect_identical(warmDest$getLatents(), beforeWarm)

rm(
  source.wA,
  state.wA,
  digest,
  matched,
  quietMatched,
  quietMismatched,
  continued,
  restored,
  restoredAgain,
  gsource.wA,
  gstate.wA,
  gmatched,
  gcontinued,
  grestored,
  gagain,
  liveArm,
  stateArm,
  preSwap,
  liveSwap,
  headlineFile,
  revived,
  revivedLatents,
  make,
  donor,
  repaired,
  inert,
  mstate,
  mforced,
  mrepaired,
  minert,
  nullSource,
  nullState,
  onesDest,
  nullDest,
  stripped,
  malformed,
  warmDonor,
  warmDest,
  beforeWarm
)
