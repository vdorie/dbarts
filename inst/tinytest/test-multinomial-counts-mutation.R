# The multinomial counts mutation channel: a
# softmax sampler's response is the n x K count matrix the combiner borrows,
# not the chain's y, so the response swap every other family reaches through
# setResponse has its own entry here. n and K are fixed at creation.
#
# The oracles are create-vs-swap parities, because there is no state comparator
# that can gate this channel: the counts are not in the serialized state, and
# the combiner caches nothing derived from them. What the parities cannot cover
# is the BURNED-IN swap - a sampler mid-Gibbs has grown trees, a live omega
# column set and a running category mix, and NO freshly built sampler can be
# put in that state - so the self-swap arm below is the only discriminating
# check on it, and its control is the SAME SPLIT run without the swap rather
# than one long run (whether a split run equals a single one is a separate
# question this file deliberately does not rest on).

set.seed(4211)
n <- 120L
p <- 3L
K <- 3L
nTest <- 15L
x <- matrix(runif(n * p), n, p)
x.test <- matrix(runif(nTest * p), nTest, p)
eta <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 0.5 * (x[, 1L] - 0.5))
probs <- exp(eta) / rowSums(exp(eta))
labels <- vapply(
  seq_len(n),
  function(i) sample.int(K, 1L, prob = probs[i, ]) - 1L,
  integer(1L)
)

# A: the one-hot single-trial matrix of the labels. B: grouped counts over the
# same rows, with different row totals - n_i drives the Polya-Gamma draw count,
# so a swap that carried only the successes would still be wrong.
countsA <- matrix(0L, n, K)
countsA[cbind(seq_len(n), labels + 1L)] <- 1L
countsB <- matrix(rpois(n * K, 1.5), n, K)
countsB[rowSums(countsB) == 0L, 1L] <- 1L
storage.mode(countsB) <- "integer"

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE
)

# every channel a multinomial run reports, the set the equivalence fixture
# records: the K softmax train and test probabilities, each category forest's
# raw fits, each category forest's cumulative split counts, and the per-sample
# per-category run varcount. $getForestFits/$getForestVariableCounts index
# forests from 1
recordChannelsCountsMutation <- function(sampler, result) {
  list(
    train = result$train,
    test = result$test,
    forestFits = lapply(seq_len(K), function(k) {
      sampler$getForestFits(k)
    }),
    varcount = lapply(seq_len(K), function(k) {
      sampler$getForestVariableCounts(k)
    }),
    runVarcount = result$varcount
  )
}

buildSamplerCountsMutation <- function(counts, n.chains = 1L) {
  ctrl <- dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    updateState = FALSE
  )
  dbarts(
    dbartsData(x, counts = counts, test = x.test),
    family = "multinomial",
    control = ctrl
  )
}

# --- Create-vs-swap parity. Building over B and building over A then
# swapping in B must be the same sampler, bitwise, on every recorded channel.
# There is no response transform to pin (the multinomial leaf scale is the
# data-independent pi*sqrt(3)/sqrt(2) anchor and sigma is fixed), so the parity
# is exact and unconditional rather than conditional on a fixed residual prior the
# way BCF's weight swap is. ---
parityArmCountsMutation <- function(build, swap, n.chains = 1L) {
  set.seed(707)
  sampler <- buildSamplerCountsMutation(build, n.chains)
  if (!is.null(swap)) {
    sampler$setCounts(swap, updateState = FALSE)
  }
  recordChannelsCountsMutation(sampler, sampler$run(20L, 8L))
}

arm.build <- parityArmCountsMutation(countsB, NULL)
arm.swap <- parityArmCountsMutation(countsA, countsB)
arm.keep <- parityArmCountsMutation(countsA, NULL)

expect_identical(arm.swap$train, arm.build$train)
expect_identical(arm.swap$test, arm.build$test)
expect_identical(arm.swap$forestFits, arm.build$forestFits)
expect_identical(arm.swap$varcount, arm.build$varcount)
expect_identical(arm.swap$runVarcount, arm.build$runVarcount)
# non-vacuity: the two count matrices are not the same posterior, so the arms
# above do not agree because the swap is inert
expect_false(isTRUE(all.equal(arm.keep$train, arm.build$train)))

# the same parity across the two creation entries, which is the enabling case:
# a sampler built from single-trial LABELS is a sampler built from the one-hot
# count matrix, so swapping grouped counts into it lands where building over
# them would have. Both entries share the build path, so only the response
# differs.
set.seed(707)
sampler.labels <- dbarts(
  x,
  factor(labels, levels = seq.int(0L, K - 1L)),
  test = x.test,
  family = "multinomial",
  control = control
)
sampler.labels$setCounts(countsB, updateState = FALSE)
arm.labels <- recordChannelsCountsMutation(
  sampler.labels,
  sampler.labels$run(20L, 8L)
)
expect_identical(arm.labels$train, arm.build$train)
expect_identical(arm.labels$forestFits, arm.build$forestFits)

# every chain sees the swap: a two-chain sampler swapped mid-life is bitwise
# the two-chain sampler built over the new counts, so no chain kept the old
arm.build.chains <- parityArmCountsMutation(countsB, NULL, n.chains = 2L)
arm.swap.chains <- parityArmCountsMutation(countsA, countsB, n.chains = 2L)
expect_identical(arm.swap.chains$train, arm.build.chains$train)
expect_identical(arm.swap.chains$forestFits, arm.build.chains$forestFits)

# --- The burned-in self-swap. Re-installing the counts a sampler is
# already running against, after it has burned in, must change nothing at all -
# the swap must not reseed the omega scratch, redraw anything, or disturb the
# trees. The control is the same split without the swap. Built over the GROUPED
# counts, so the trials the swap re-derives are not all 1 and a recompute that
# lost them would show here. ---
splitArm <- function(swap) {
  set.seed(911)
  sampler <- buildSamplerCountsMutation(countsB)
  sampler$run(25L, 6L)
  if (!is.null(swap)) {
    sampler$setCounts(swap, updateState = FALSE)
  }
  recordChannelsCountsMutation(sampler, sampler$run(0L, 6L))
}

arm.self <- splitArm(countsB)
arm.control <- splitArm(NULL)
expect_identical(arm.self$train, arm.control$train)
expect_identical(arm.self$test, arm.control$test)
expect_identical(arm.self$forestFits, arm.control$forestFits)
expect_identical(arm.self$varcount, arm.control$varcount)
expect_identical(arm.self$runVarcount, arm.control$runVarcount)

# and the burned-in channel is live, not inert: swapping in DIFFERENT counts at
# the same point moves the draws and leaves the sampler sane
arm.burned <- splitArm(countsA)
expect_false(isTRUE(all.equal(arm.burned$train, arm.control$train)))
expect_true(all(is.finite(arm.burned$train)))

# --- The reported probabilities and the reported forest fits are the same
# vintage. storeSample blends the K forests AFTER the level-centering move, and
# the per-forest fits query reads the post-run totalFits, so for the last
# recorded sample of a single-chain run the softmax of the K forest fits must
# be the recorded train channel. A tolerance, not a bitwise check: an R-side
# softmax does not reproduce the engine's reduction order. Pinned here at the
# null offset as the pre-existing invariant it is. ---
set.seed(313)
sampler.vintage <- buildSamplerCountsMutation(countsB)
res.vintage <- sampler.vintage$run(20L, 5L)
fits.vintage <- vapply(
  seq_len(K),
  function(k) sampler.vintage$getForestFits(k)[, 1L],
  numeric(n)
)
softmax.vintage <- exp(fits.vintage - apply(fits.vintage, 1L, max))
softmax.vintage <- softmax.vintage / rowSums(softmax.vintage)
expect_equal(
  softmax.vintage,
  matrix(res.vintage$train[,, 5L], n, K),
  tolerance = 1e-12
)

# --- setState after a setCounts: a restore reinstalls trees against WHATEVER
# counts the sampler holds then, exactly as a single-forest restore does
# against the current y. The counts are data and ride no wire block, so the
# state carries none. ---
set.seed(515)
sampler.state <- buildSamplerCountsMutation(countsA)
sampler.state$run(20L, 4L)
sampler.state$storeState()
state.A <- sampler.state$state
sampler.state$setCounts(countsB, updateState = FALSE)
expect_silent(status <- sampler.state$setState(state.A))
expect_true(status)
res.restored <- sampler.state$run(0L, 4L)
expect_true(all(is.finite(res.restored$train)))
# the restored trees run against B, not against the A they were fitted to: the
# same restore under A draws a different chain
set.seed(515)
sampler.stateA <- buildSamplerCountsMutation(countsA)
sampler.stateA$run(20L, 4L)
sampler.stateA$storeState()
sampler.stateA$setState(sampler.stateA$state)
expect_false(isTRUE(all.equal(
  res.restored$train,
  sampler.stateA$run(0L, 4L)$train
)))

# --- Refusals, the counts half: what the channel refuses, and what the response-side
# conduits now say. ---
sampler.mn <- buildSamplerCountsMutation(countsA)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no counts, and both must name the family
# situation rather than the forest count
sampler.gaussian <- dbarts(x, rnorm(n), control = control)
expect_error(
  sampler.gaussian$setCounts(countsA),
  "no count response"
)
set.seed(17)
z <- rbinom(n, 1L, 0.5)
sampler.bcf <- dbarts(
  x,
  rnorm(n),
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 10L)),
  control = control
)
expect_error(
  sampler.bcf$setCounts(countsA),
  "no count response"
)

# the bridge's own memory-safety backstop on a multi-forest sampler's
# setData/setModel (refuseMultiForestMutation) has no route left through
# either R5 method: a multinomial sampler's own $setData/$setModel refuse
# R-side first (below, "not available on a multinomial sampler"), and a BCF
# sampler's refuse R-side too (refuseAmplitudeMutation, test-bcf-mutation-pins.R).
# Pinned here directly on the raw pointer, the one route still able to reach it.
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setData,
    sampler.mn$getPointer(),
    sampler.mn$data
  ),
  "multi-forest"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setModel,
    sampler.mn$getPointer(),
    sampler.mn$model,
    sampler.mn$data,
    sampler.mn$control
  ),
  "multi-forest"
)

# n and K are out of scope, and the refusal names both. A transposed matrix is
# the case a length test alone would install into the wrong cells.
expect_error(
  sampler.mn$setCounts(countsA[seq_len(n - 1L), ]),
  "same number of rows"
)
expect_error(
  sampler.mn$setCounts(t(countsA)),
  "3 categories"
)

# the count invariants, restated at the entrance rather than inherited from
# creation. $setCounts refuses R-side first, so the bridge's nonnegativity
# check, a memory-safety backstop against a negative trial count, is pinned on
# the raw pointer as well.
counts.negative <- countsA
counts.negative[1L, 1L] <- -1L
expect_error(
  sampler.mn$setCounts(counts.negative),
  "non-negative"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCounts,
    sampler.mn$getPointer(),
    counts.negative
  ),
  "non-negative"
)
# NA_INTEGER is INT_MIN, so the same bridge check catches NA; the R layer
# names it directly
counts.na <- countsA
counts.na[2L, 1L] <- NA_integer_
expect_error(sampler.mn$setCounts(counts.na), "NA")
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCounts,
    sampler.mn$getPointer(),
    counts.na
  ),
  "non-negative"
)
# an empty row is accepted, entering no likelihood, and warns once per session
counts.empty <- countsA
counts.empty[3L, ] <- 0L
emptyKeyEnv <- dbarts:::onceWarnState
emptyKeyEnv[["multinomialZeroTrials"]] <- NULL
expect_warning(
  sampler.mn$setCounts(counts.empty),
  "zero trials"
)
expect_identical(sampler.mn$data@counts, counts.empty)
expect_silent(sampler.mn$setCounts(countsA))
# a row sum past the trial cap, here far enough past it to overflow the int
# the trials are counted in: the accumulation is checked, not wrapped
counts.overflow <- countsA
counts.overflow[4L, 1L] <- 2000000000L
counts.overflow[4L, 2L] <- 2000000000L
expect_error(
  sampler.mn$setCounts(counts.overflow),
  "requires 'counts' row totals no larger than 1000000"
)

# a refusal leaves the sampler byte-identical, including one that fires PART
# WAY THROUGH the validation. The overflow matrix is the case that gets
# furthest: it passes every R-side check and every per-cell check up to the row
# it overflows on, so a build that validated in place would already have
# written new counts into the buffer the combiner borrows.
refusalArmCountsMutation <- function(attempt) {
  set.seed(808)
  sampler <- buildSamplerCountsMutation(countsA)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        sampler$setCounts(attempt, updateState = FALSE)
        NA_character_
      },
      error = conditionMessage
    )
  }
  c(
    list(refused = refused),
    recordChannelsCountsMutation(sampler, sampler$run(15L, 5L))
  )
}
arm.refused <- refusalArmCountsMutation(counts.overflow)
arm.untouched <- refusalArmCountsMutation(NULL)
expect_true(grepl("row totals no larger than 1000000", arm.refused$refused))
expect_identical(arm.refused$train, arm.untouched$train)
expect_identical(arm.refused$forestFits, arm.untouched$forestFits)
expect_identical(arm.refused$runVarcount, arm.untouched$runVarcount)

# the response-side conduits stay refused - a flat y cannot express an n x K
# count matrix, a flat offset points exactly along the softmax's null
# direction, and the case weights an integer weight would express are already
# row-wise count replication - but each refusal now names the channel that
# works instead of reporting a response fixed at creation, which it no longer
# is.
expect_error(
  sampler.mn$setResponse(as.double(labels)),
  "n x K count matrix"
)
expect_error(
  sampler.mn$setOffset(rep(0.5, n)),
  "n x K matrix"
)
expect_error(sampler.mn$setWeights(runif(n, 0.5, 1.5)), "row-wise")
# and a BCF sampler, which DOES opt into the response conduit, refuses an
# updateScale that is neither TRUE nor FALSE ahead of its amplitude guard
# (test-bcf-mutation-pins.R pins that guard at updateScale = TRUE)
expect_error(
  sampler.bcf$setResponse(rnorm(n), updateScale = NA),
  "'updateScale' must be TRUE or FALSE"
)

# the whole-data and whole-model mutations the multinomial battery had never
# pinned, and the pinned-sigma refusal: none of them is opened by the counts
# channel, which replaces the response and nothing else
expect_error(sampler.mn$setData(sampler.mn$data), "not available")
expect_error(
  sampler.mn$setModel(sampler.mn$model),
  "not available"
)
expect_error(
  sampler.mn$setSigma(5),
  "the softmax carries no residual scale to set"
)
