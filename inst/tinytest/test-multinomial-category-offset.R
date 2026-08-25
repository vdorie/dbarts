# The multinomial category offset (docs/design/multinomial.md): an n x K matrix
# entering the latent as f_ik + o_ik, so it shifts the log-sum-exp margins, is
# subtracted back out of each category's working response, and rides the
# reported softmax - never a leaf value, and never the response model's offset,
# which is added to every reported channel AFTER the K forests are blended.
#
# The oracles are create-vs-swap parity and a null-path neutrality gate, plus
# the one check the parities structurally cannot make: a parity cannot tell two
# fits of DIFFERENT VINTAGES apart, since both arms would report the same stale
# column. The same-vintage cross-check below is the oracle for that, and its
# second arm drives a per-observation predictor session, which rewrites the
# category forests' fits between sweeps without ever entering the combiner.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(9021)
n <- 90L
p <- 3L
K <- 3L
x <- matrix(runif(n * p), n, p)
eta <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 0.5 * (x[, 1L] - 0.5))
probs <- exp(eta) / rowSums(exp(eta))
labels <- vapply(
  seq_len(n),
  function(i) sample.int(K, 1L, prob = probs[i, ]) - 1L,
  integer(1L)
)
counts <- matrix(rpois(n * K, 1.2), n, K)
counts[rowSums(counts) == 0L, 1L] <- 1L
storage.mode(counts) <- "integer"

# a genuinely per-category, per-observation offset: no column is a copy of
# another and no row is constant, so nothing about it lies along the softmax's
# null direction
offset <- cbind(
  0.8 * (x[, 1L] - 0.5),
  -0.6 + 0.4 * x[, 2L],
  0.3 * sin(6 * x[, 3L])
)
zeroOffset <- matrix(0, n, K)

control <- function(n.chains = 1L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    updateState = FALSE
  )
}

# the recorded channels of an offset handle built without test data: the train
# blend, the per-category forest fits and the two variable-count channels (the
# test channel and its own offset are test-multinomial-test-offset.R's)
recordChannels <- function(bc, result) {
  list(
    train = result$train,
    forestFits = lapply(seq_len(K) - 1L, function(k) {
      bartcoreForestFits(bc, k)
    }),
    varcount = lapply(seq_len(K) - 1L, function(k) {
      bartcoreForestVariableCounts(bc, k)
    }),
    runVarcount = result$varcount
  )
}

buildSampler <- function(offset = NULL, n.chains = 1L) {
  host <- dbarts(x, as.double(labels), control = control(n.chains))
  dbarts:::bartcoreMultinomialCountSampler(host, counts, K = K, offset = offset)
}

# --- Create-vs-swap parity. Creating with the offset and installing it on an
# offset-free sampler must be the same sampler, bitwise: nothing derived from
# the offset survives a sweep, so there is no state for the two routes to
# disagree about. There is no response transform to pin (the leaf scale is the
# data-independent pi*sqrt(3)/sqrt(2) anchor and sigma is fixed), so the parity
# is exact and unconditional. ---
parityArm <- function(build, swap, n.chains = 1L) {
  set.seed(4242)
  bc <- buildSampler(build, n.chains)
  if (!is.null(swap)) {
    bartcoreSetCategoryOffset(bc, swap)
  }
  recordChannels(bc, bartcoreRun(bc, 20L, 8L))
}

arm.build <- parityArm(offset, NULL)
arm.swap <- parityArm(NULL, offset)
arm.none <- parityArm(NULL, NULL)

expect_identical(arm.swap$train, arm.build$train)
expect_identical(arm.swap$forestFits, arm.build$forestFits)
expect_identical(arm.swap$varcount, arm.build$varcount)
expect_identical(arm.swap$runVarcount, arm.build$runVarcount)
# non-vacuity: the offset moves the answer, so the two arms above do not agree
# because it is inert
expect_false(isTRUE(all.equal(arm.none$train, arm.build$train)))

# the label creation entry takes the offset too, and lands where the count
# entry does: a one-hot count matrix is the single-trial reduction of the same
# engine, offset included
labelCounts <- matrix(0L, n, K)
labelCounts[cbind(seq_len(n), labels + 1L)] <- 1L
set.seed(4242)
bc.labels <- dbarts:::bartcoreMultinomialSampler(
  dbarts(x, as.double(labels), control = control()),
  labels,
  K = K,
  offset = offset
)
set.seed(4242)
bc.onehot <- dbarts:::bartcoreMultinomialCountSampler(
  dbarts(x, as.double(labels), control = control()),
  labelCounts,
  K = K,
  offset = offset
)
expect_identical(
  bartcoreRun(bc.labels, 20L, 8L)$train,
  bartcoreRun(bc.onehot, 20L, 8L)$train
)

# every chain sees the install: a two-chain sampler offset mid-life is bitwise
# the two-chain sampler built with it, so no chain kept the offset-free latent
arm.build.chains <- parityArm(offset, NULL, n.chains = 2L)
arm.swap.chains <- parityArm(NULL, offset, n.chains = 2L)
expect_identical(arm.swap.chains$train, arm.build.chains$train)
expect_identical(arm.swap.chains$forestFits, arm.build.chains$forestFits)

# --- Null-path neutrality, a hard gate rather than a flagged expectation. An
# all-zero matrix must be BITWISE the null path: the margin's log-sum-exp
# compares its arguments before combining them, exp(-0.0) is 1, and the
# Polya-Gamma sampler opens on |psi|, so -0.0 is absorbed at every consumer and
# x + 0.0 is x. Clearing an installed offset must return to the same path, so
# the sampler carries nothing of an offset it no longer has. ---
arm.zero <- parityArm(zeroOffset, NULL)
expect_identical(arm.zero$train, arm.none$train)
expect_identical(arm.zero$forestFits, arm.none$forestFits)
expect_identical(arm.zero$runVarcount, arm.none$runVarcount)

set.seed(4242)
bc.cleared <- buildSampler(offset)
bartcoreSetCategoryOffset(bc.cleared, NULL)
arm.cleared <- recordChannels(
  bc.cleared,
  bartcoreRun(bc.cleared, 20L, 8L)
)
expect_identical(arm.cleared$train, arm.none$train)
expect_identical(arm.cleared$forestFits, arm.none$forestFits)

# and installing the zero matrix on a live sampler is likewise the null path
set.seed(4242)
bc.zeroswap <- buildSampler(NULL)
bartcoreSetCategoryOffset(bc.zeroswap, zeroOffset)
expect_identical(
  bartcoreRun(bc.zeroswap, 20L, 8L)$train,
  arm.none$train
)

# --- Only the row-centred part of the offset is identified: adding a constant
# to every entry of a ROW leaves the margin, the working response and every
# reported probability unchanged in exact arithmetic, which is why the entrance
# does not silently re-centre the input. The two runs diverge only through
# rounding inside the log-sum-exp, so a short run agrees far below any
# tolerance a sign or placement error could hide in. ---
rowShift <- offset + matrix(rep(0.7 * x[, 1L], K), n, K)
columnShift <- offset
columnShift[, 2L] <- columnShift[, 2L] + 0.7
set.seed(1717)
train.plain <- bartcoreRun(buildSampler(offset), 0L, 3L)$train
set.seed(1717)
train.rowshift <- bartcoreRun(buildSampler(rowShift), 0L, 3L)$train
set.seed(1717)
train.colshift <- bartcoreRun(buildSampler(columnShift), 0L, 3L)$train
expect_equal(train.rowshift, train.plain, tolerance = 1e-8)
# non-vacuity: a shift of ONE column is not a null direction and does move it
expect_false(isTRUE(all.equal(train.colshift, train.plain)))

# --- The reported probabilities and the reported forest fits are the same
# vintage, WITH the offset. storeSample blends the K forests after the
# level-centering move and the per-forest query reads the post-run fits, so for
# the last recorded sample of a single-chain run the softmax of the offset fits
# must be the recorded train channel. A tolerance, not a bitwise check: an
# R-side softmax does not reproduce the engine's reduction order.
#
# This is the check a parity cannot make. A blend that refreshed only the
# category it had just drawn would leave the last category holding the previous
# sweep's fit, and BOTH arms of any parity would report that same stale value. ---
softmaxRows <- function(raw) {
  e <- exp(raw - apply(raw, 1L, max))
  e / rowSums(e)
}
sameVintage <- function(bc, result, sampleNum, offsetMatrix) {
  fits <- vapply(
    seq_len(K) - 1L,
    function(k) bartcoreForestFits(bc, k)[, 1L],
    numeric(n)
  )
  max(abs(
    softmaxRows(fits + offsetMatrix) -
      matrix(result$train[,, sampleNum], n, K)
  ))
}

set.seed(313)
bc.vintage <- buildSampler(offset)
res.vintage <- bartcoreRun(bc.vintage, 20L, 5L)
expect_true(sameVintage(bc.vintage, res.vintage, 5L, offset) < 1e-12)

# second arm, against a fits rewrite the combiner never sees: a per-observation
# predictor session rebuilds every category forest's fits from its trees at
# MUTATION time, outside any sweep. A blend that refreshed its offset fits only
# from inside the sweep would report the pre-session fits for one more sample.
set.seed(515)
bc.session <- buildSampler(NULL)
bartcoreRun(bc.session, 20L, 4L)
bartcoreSetCategoryOffset(bc.session, offset)
installed <- bartcoreUpdatePredictorPerObservation(
  bc.session,
  pmin(pmax(x[, 2L] + rnorm(n, 0, 0.02), 0), 1),
  2L
)
expect_true(sum(installed) > 0L)
res.session <- bartcoreRun(bc.session, 0L, 2L)
expect_true(sameVintage(bc.session, res.session, 2L, offset) < 1e-12)

# --- The simplex invariant holds with an offset installed: the offset enters
# BEFORE the softmax, so the reported values are still K probabilities per row.
# (Added after the blend they would not be - the failure a test offset once
# produced on this same coupling.) ---
expect_true(all(is.finite(res.vintage$train)))
expect_true(all(res.vintage$train >= 0 & res.vintage$train <= 1))
expect_equal(
  apply(res.vintage$train, c(1L, 3L), sum),
  matrix(1, n, 5L),
  tolerance = 1e-12
)

# --- Refusals, the offset half. ---
bc.mn <- buildSampler(NULL)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no category offset, and both must name the
# family situation rather than the forest count
bc.gaussian <- dbarts:::bartcoreSampler(
  dbarts(x, rnorm(n), control = control())
)
expect_error(
  bartcoreSetCategoryOffset(bc.gaussian, offset),
  "requires a multinomial"
)
set.seed(17)
bc.bcf <- dbarts:::bartcoreBCFSampler(
  dbarts(x, rnorm(n), control = control()),
  rbinom(n, 1L, 0.5),
  n.trees.treatment = 10L
)
expect_error(
  bartcoreSetCategoryOffset(bc.bcf, offset),
  "requires a multinomial"
)

# n and K are out of scope, and the refusal names both. A transposed matrix is
# the case a length test alone would install cell by cell into the wrong rows:
# it carries exactly n * K entries.
expect_error(
  bartcoreSetCategoryOffset(bc.mn, offset[seq_len(n - 1L), ]),
  "90 x 3"
)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCategoryOffset, bc.mn$ptr, t(offset)),
  "90 observations x 3 categories"
)
# a flat vector is refused rather than recycled across the categories
expect_error(
  bartcoreSetCategoryOffset(bc.mn, rep(0.5, n)),
  "90 x 3"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCategoryOffset,
    bc.mn$ptr,
    rep(0.5, n * K)
  ),
  "90 observations x 3 categories"
)
# an integer matrix is refused on TYPE C-side: it is the wrong buffer to borrow
# whatever its values (the R wrapper coerces first, so the C entry is called
# directly here)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCategoryOffset,
    bc.mn$ptr,
    matrix(0L, n, K)
  ),
  "90 observations x 3 categories"
)

# every non-finite entry is refused: an infinity propagates through the
# log-sum-exp margin into a NaN for every category of its row, and NA is not a
# shift. Both layers are pinned, the C entry called directly for each.
for (bad in c(NA_real_, NaN, Inf, -Inf)) {
  spoiled <- offset
  spoiled[3L, 2L] <- bad
  expect_error(bartcoreSetCategoryOffset(bc.mn, spoiled), "finite")
  expect_error(
    .Call(dbarts:::C_dbarts_bartcore_setCategoryOffset, bc.mn$ptr, spoiled),
    "finite"
  )
}

# a refusal leaves the sampler byte-identical: the entrance validates a whole
# scratch copy and swaps it in only once it holds, because the combiner borrows
# the installed buffer and an in-place write would BE the mutation
refusalArm <- function(attempt) {
  set.seed(808)
  bc <- buildSampler(NULL)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        bartcoreSetCategoryOffset(bc, attempt)
        NA_character_
      },
      error = conditionMessage
    )
  }
  c(
    list(refused = refused),
    recordChannels(bc, bartcoreRun(bc, 15L, 5L))
  )
}
spoiled <- offset
spoiled[n, K] <- Inf
arm.refused <- refusalArm(spoiled)
arm.untouched <- refusalArm(NULL)
expect_true(grepl("finite", arm.refused$refused))
expect_identical(arm.refused$train, arm.untouched$train)
expect_identical(arm.refused$forestFits, arm.untouched$forestFits)

# --- The test surface under a TRAIN offset. Each test-side channel now carries
# its own per-category offset (test-multinomial-test-offset.R), so a sampler
# holding a train offset no longer has to refuse one: test data installs, the
# offset installs against test data, and the recorded test channel is the
# category forests' test blend, which is what the test rows' own offset shifts.
# What none of them does is REUSE the train offset - its rows are the train
# rows - so the test channel here is the surface at a zero test offset, and
# predict, whose rows are the caller's, refuses to guess. ---
nTest <- 12L
x.test <- x[seq_len(nTest), , drop = FALSE]
# built exactly where buildSampler builds its host, so the two draw streams
# line up and the train comparison below is a statement about test data rather
# than about creation order
set.seed(4242)
hostTest <- dbarts(x, as.double(labels), test = x.test, control = control())
bc.createTest <- dbarts:::bartcoreMultinomialCountSampler(
  hostTest,
  counts,
  K = K,
  offset = offset
)
res.createTest <- bartcoreRun(bc.createTest, 20L, 8L)
# the train channels are the no-test-data sampler's, bitwise: test rows consume
# no rng and enter no likelihood
expect_identical(res.createTest$train, arm.build$train)
expect_true(all(is.finite(res.createTest$test)))
expect_equal(
  apply(res.createTest$test, c(1L, 3L), sum),
  matrix(1, nTest, 8L),
  tolerance = 1e-12
)
# the label entry likewise
expect_silent(
  dbarts:::bartcoreMultinomialSampler(hostTest, labels, K = K, offset = offset)
)
# installing the train offset on a sampler that already holds test data
bc.test <- dbarts:::bartcoreMultinomialCountSampler(hostTest, counts, K = K)
expect_silent(bartcoreSetCategoryOffset(bc.test, offset))
# and the reverse direction: test data installs on a sampler already holding a
# train offset, whose test rows the train offset says nothing about
bc.offset <- buildSampler(offset)
expect_silent(bartcoreSetTestPredictor(bc.offset, x.test))
expect_true(all(is.finite(bartcoreRun(bc.offset, 0L, 3L)$test)))
# the combined entry takes the test rows too; its FLAT offset argument stays
# refused, and truthfully - after the blend it leaves the simplex, before it a
# common per-observation shift is inert - and now names the matrix entry
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestPredictorAndOffset,
    bc.offset$ptr,
    x.test,
    rep(0, nTest)
  ),
  "category test offset channel"
)
expect_silent(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestPredictorAndOffset,
    bc.offset$ptr,
    x.test,
    NULL
  )
)
# predict is the one channel that still refuses under a train offset, and for a
# different reason: its rows are the caller's, so no resident offset describes
# them and none is substituted. An explicit matrix is the way through.
expect_error(
  bartcorePredict(bc.offset, x.test),
  "cannot be inferred"
)

# the response-side conduit's offset half now names the entry that works. The
# flat vector stays refused, and truthfully: a common per-observation shift is
# exactly the softmax's null direction, so it could only ever be inert.
expect_error(
  bartcoreSetOffset(bc.mn, rep(0.5, n)),
  "n x K category matrix"
)
expect_error(
  bartcoreSetOffset(bc.mn, rep(0.5, n), updateScale = TRUE),
  "n x K category matrix"
)
# predict's own flat-offset refusal is narrow rather than false, and says which
# form is not: the R wrapper shapes the argument against the handle's K, and
# the C entry refuses a dimensionless vector on its own
expect_error(
  bartcorePredict(bc.mn, x.test, rep(0.5, nTest)),
  "12 x 3 matrix"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_predict,
    bc.mn$ptr,
    x.test,
    rep(0.5, nTest),
    1L
  ),
  "per-category matrix"
)

# --- The counts entrance's range check, which precedes the integer coercion.
# A count past .Machine$integer.max coerces to NA, and the engine would then
# refuse it as negative - a true refusal naming the wrong reason. ---
counts.huge <- matrix(as.double(counts), n, K)
counts.huge[1L, 1L] <- 2^31
expect_error(
  bartcoreSetCounts(bc.mn, counts.huge),
  "representable as integers"
)
expect_error(
  dbarts:::bartcoreMultinomialCountSampler(
    dbarts(x, as.double(labels), control = control()),
    counts.huge,
    K = K
  ),
  "representable as integers"
)
