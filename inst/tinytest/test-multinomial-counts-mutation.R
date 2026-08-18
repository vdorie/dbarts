# The multinomial counts mutation channel (docs/design/multinomial.md): a
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
# per-category run varcount
recordChannels <- function(bc, result) {
  list(
    train = result$train,
    test = result$test,
    forestFits = lapply(seq_len(K) - 1L, function(k) {
      dbarts:::bartcoreForestFits(bc, k)
    }),
    varcount = lapply(seq_len(K) - 1L, function(k) {
      dbarts:::bartcoreForestVariableCounts(bc, k)
    }),
    runVarcount = result$varcount
  )
}

buildSampler <- function(counts, n.chains = 1L) {
  ctrl <- dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    updateState = FALSE
  )
  host <- dbarts(x, as.double(labels), test = x.test, control = ctrl)
  dbarts:::bartcoreMultinomialCountSampler(host, counts, K = K)
}

# --- Create-vs-swap parity. Building over B and building over A then
# swapping in B must be the same sampler, bitwise, on every recorded channel.
# There is no response transform to pin (the multinomial leaf scale is the
# data-independent pi*sqrt(3)/sqrt(2) anchor and sigma is fixed), so the parity
# is exact and unconditional rather than conditional on a fixed resid.prior the
# way BCF's weight swap is. ---
parityArm <- function(build, swap, n.chains = 1L) {
  set.seed(707)
  bc <- buildSampler(build, n.chains)
  if (!is.null(swap)) {
    dbarts:::bartcoreSetCounts(bc, swap)
  }
  recordChannels(bc, dbarts:::bartcoreRun(bc, 20L, 8L))
}

arm.build <- parityArm(countsB, NULL)
arm.swap <- parityArm(countsA, countsB)
arm.keep <- parityArm(countsA, NULL)

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
bc.labels <- dbarts:::bartcoreMultinomialSampler(
  dbarts(x, as.double(labels), test = x.test, control = control),
  labels,
  K = K
)
dbarts:::bartcoreSetCounts(bc.labels, countsB)
arm.labels <- recordChannels(
  bc.labels,
  dbarts:::bartcoreRun(bc.labels, 20L, 8L)
)
expect_identical(arm.labels$train, arm.build$train)
expect_identical(arm.labels$forestFits, arm.build$forestFits)

# every chain sees the swap: a two-chain sampler swapped mid-life is bitwise
# the two-chain sampler built over the new counts, so no chain kept the old
arm.build.chains <- parityArm(countsB, NULL, n.chains = 2L)
arm.swap.chains <- parityArm(countsA, countsB, n.chains = 2L)
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
  bc <- buildSampler(countsB)
  dbarts:::bartcoreRun(bc, 25L, 6L)
  if (!is.null(swap)) {
    dbarts:::bartcoreSetCounts(bc, swap)
  }
  recordChannels(bc, dbarts:::bartcoreRun(bc, 0L, 6L))
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
bc.vintage <- buildSampler(countsB)
res.vintage <- dbarts:::bartcoreRun(bc.vintage, 20L, 5L)
fits.vintage <- vapply(
  seq_len(K) - 1L,
  function(k) dbarts:::bartcoreForestFits(bc.vintage, k)[, 1L],
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
bc.state <- buildSampler(countsA)
dbarts:::bartcoreRun(bc.state, 20L, 4L)
state.A <- dbarts:::bartcoreStoreState(bc.state)
dbarts:::bartcoreSetCounts(bc.state, countsB)
expect_silent(dbarts:::bartcoreSetState(bc.state, state.A))
res.restored <- dbarts:::bartcoreRun(bc.state, 0L, 4L)
expect_true(all(is.finite(res.restored$train)))
# the restored trees run against B, not against the A they were fitted to: the
# same restore under A draws a different chain
set.seed(515)
bc.stateA <- buildSampler(countsA)
dbarts:::bartcoreRun(bc.stateA, 20L, 4L)
dbarts:::bartcoreSetState(bc.stateA, dbarts:::bartcoreStoreState(bc.stateA))
expect_false(isTRUE(all.equal(
  res.restored$train,
  dbarts:::bartcoreRun(bc.stateA, 0L, 4L)$train
)))

# --- Refusals, the counts half: what the channel refuses, and what the response-side
# conduits now say. ---
bc.mn <- buildSampler(countsA)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no counts, and both must name the family
# situation rather than the forest count
bc.gaussian <- dbarts:::bartcoreSampler(
  dbarts(x, rnorm(n), control = control)
)
expect_error(
  dbarts:::bartcoreSetCounts(bc.gaussian, countsA),
  "requires a multinomial"
)
set.seed(17)
z <- rbinom(n, 1L, 0.5)
bc.bcf <- dbarts:::bartcoreBCFSampler(
  dbarts(x, rnorm(n), control = control),
  z,
  n.trees.treatment = 10L
)
expect_error(
  dbarts:::bartcoreSetCounts(bc.bcf, countsA),
  "requires a multinomial"
)

# n and K are out of scope, and the refusal names both. A transposed matrix is
# the case a length test alone would install into the wrong cells: it carries
# exactly n * K entries.
expect_error(
  dbarts:::bartcoreSetCounts(bc.mn, countsA[seq_len(n - 1L), ]),
  "120 observations x 3 categories"
)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCounts, bc.mn$ptr, t(countsA)),
  "120 observations x 3 categories"
)
# a real matrix is refused on TYPE rather than rounded: it is the wrong buffer
# to borrow, whatever its values
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCounts,
    bc.mn$ptr,
    matrix(as.double(countsA), n, K)
  ),
  "120 observations x 3 categories"
)

# the count invariants, restated at the entrance rather than inherited from
# creation. The C entry is called directly for each, since the R wrapper
# refuses first; both layers are pinned.
counts.negative <- countsA
counts.negative[1L, 1L] <- -1L
expect_error(
  dbarts:::bartcoreSetCounts(bc.mn, counts.negative),
  "non-negative"
)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCounts, bc.mn$ptr, counts.negative),
  "non-negative"
)
# NA needs no test of its own C-side: NA_INTEGER is INT_MIN, so the
# nonnegativity check catches it. Pinned anyway.
counts.na <- countsA
counts.na[2L, 1L] <- NA_integer_
expect_error(dbarts:::bartcoreSetCounts(bc.mn, counts.na), "missing values")
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCounts, bc.mn$ptr, counts.na),
  "non-negative"
)
# an empty row: PG(0, .) is a point mass at zero and the working response
# divides by omega, so a zero row sum is refused rather than fit
counts.empty <- countsA
counts.empty[3L, ] <- 0L
expect_error(dbarts:::bartcoreSetCounts(bc.mn, counts.empty), "at least one")
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCounts, bc.mn$ptr, counts.empty),
  "at least one"
)
# a row sum that overflows the int the trials are counted in: the accumulation
# is checked, not wrapped
counts.overflow <- countsA
counts.overflow[4L, 1L] <- 2000000000L
counts.overflow[4L, 2L] <- 2000000000L
expect_error(
  dbarts:::bartcoreSetCounts(bc.mn, counts.overflow),
  "fit in an integer"
)

# a refusal leaves the sampler byte-identical, including one that fires PART
# WAY THROUGH the validation. The overflow matrix is the case that gets
# furthest: it passes every R-side check and every per-cell check up to the row
# it overflows on, so a build that validated in place would already have
# written new counts into the buffer the combiner borrows.
refusalArm <- function(attempt) {
  set.seed(808)
  bc <- buildSampler(countsA)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        dbarts:::bartcoreSetCounts(bc, attempt)
        NA_character_
      },
      error = conditionMessage
    )
  }
  c(
    list(refused = refused),
    recordChannels(bc, dbarts:::bartcoreRun(bc, 15L, 5L))
  )
}
arm.refused <- refusalArm(counts.overflow)
arm.untouched <- refusalArm(NULL)
expect_true(grepl("fit in an integer", arm.refused$refused))
expect_identical(arm.refused$train, arm.untouched$train)
expect_identical(arm.refused$forestFits, arm.untouched$forestFits)
expect_identical(arm.refused$runVarcount, arm.untouched$runVarcount)

# the response-side conduits stay refused - a flat y cannot express an n x K
# count matrix, a flat offset points exactly along the softmax's null
# direction, and the case weights an integer weight would express are already
# row-wise count replication - but each refusal now names the channel that
# works instead of reporting a response fixed at creation, which it no longer
# is. The guard is shared with the flat C API, so both surfaces say this.
expect_error(
  dbarts:::bartcoreSetResponse(bc.mn, as.double(labels)),
  "n x K count matrix"
)
expect_error(
  dbarts:::bartcoreSetOffset(bc.mn, rep(0.5, n)),
  "n x K category matrix"
)
expect_error(dbarts:::bartcoreSetWeights(bc.mn, runif(n, 0.5, 1.5)), "n x K")
# and a BCF sampler, which DOES opt into the response conduit, keeps the
# generic wording: the counts hint is conditioned on the capability, not on the
# forest count
expect_error(
  dbarts:::bartcoreSetResponse(bc.bcf, rnorm(n), updateScale = TRUE),
  "multi-forest"
)

# the whole-data and whole-model mutations the multinomial battery had never
# pinned, and the pinned-sigma refusal: none of them is opened by the counts
# channel, which replaces the response and nothing else
host.mn <- dbarts(x, as.double(labels), test = x.test, control = control)
expect_error(dbarts:::bartcoreSetData(bc.mn, host.mn$data), "multi-forest")
expect_error(
  dbarts:::bartcoreSetModel(bc.mn, host.mn$model, host.mn$data),
  "multi-forest"
)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setSigma, bc.mn$ptr, 5),
  "response family fixes the residual standard deviation"
)
