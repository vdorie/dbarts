# The multinomial category offset: an n x K matrix
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

controlCategoryOffset <- function(n.chains = 1L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    updateState = FALSE
  )
}

# the recorded channels of a sampler built without test data: the train
# blend, the per-category forest fits and the two variable-count channels (the
# test channel and its own offset are test-multinomial-test-offset.R's).
# $getForestFits/$getForestVariableCounts index forests from 1
recordChannelsCategoryOffset <- function(sampler, result) {
  list(
    train = result$train,
    forestFits = lapply(seq_len(K), function(k) {
      sampler$getForestFits(k)
    }),
    varcount = lapply(seq_len(K), function(k) {
      sampler$getForestVariableCounts(k)
    }),
    runVarcount = result$varcount
  )
}

buildSamplerCategoryOffset <- function(offset = NULL, n.chains = 1L) {
  sampler <- dbarts(
    dbartsData(x, counts = counts),
    family = "multinomial",
    control = controlCategoryOffset(n.chains)
  )
  if (!is.null(offset)) {
    sampler$setCategoryOffset(offset, updateState = FALSE)
  }
  sampler
}

# --- Create-vs-swap parity. Creating with the offset and installing it on an
# offset-free sampler must be the same sampler, bitwise: nothing derived from
# the offset survives a sweep, so there is no state for the two routes to
# disagree about. There is no response transform to pin (the leaf scale is the
# data-independent pi*sqrt(3)/sqrt(2) anchor and sigma is fixed), so the parity
# is exact and unconditional. ---
parityArmCategoryOffset <- function(build, swap, n.chains = 1L) {
  set.seed(4242)
  sampler <- buildSamplerCategoryOffset(build, n.chains)
  if (!is.null(swap)) {
    sampler$setCategoryOffset(swap, updateState = FALSE)
  }
  recordChannelsCategoryOffset(sampler, sampler$run(20L, 8L))
}

arm.build <- parityArmCategoryOffset(offset, NULL)
arm.swap <- parityArmCategoryOffset(NULL, offset)
arm.none <- parityArmCategoryOffset(NULL, NULL)

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
sampler.labels <- dbarts(
  x,
  factor(labels, levels = seq.int(0L, K - 1L)),
  family = "multinomial",
  control = controlCategoryOffset()
)
sampler.labels$setCategoryOffset(offset, updateState = FALSE)
set.seed(4242)
sampler.onehot <- dbarts(
  dbartsData(x, counts = labelCounts),
  family = "multinomial",
  control = controlCategoryOffset()
)
sampler.onehot$setCategoryOffset(offset, updateState = FALSE)
expect_identical(
  sampler.labels$run(20L, 8L)$train,
  sampler.onehot$run(20L, 8L)$train
)

# every chain sees the install: a two-chain sampler offset mid-life is bitwise
# the two-chain sampler built with it, so no chain kept the offset-free latent
arm.build.chains <- parityArmCategoryOffset(offset, NULL, n.chains = 2L)
arm.swap.chains <- parityArmCategoryOffset(NULL, offset, n.chains = 2L)
expect_identical(arm.swap.chains$train, arm.build.chains$train)
expect_identical(arm.swap.chains$forestFits, arm.build.chains$forestFits)

# --- Null-path neutrality, a hard gate rather than a flagged expectation. An
# all-zero matrix must be BITWISE the null path: the margin's log-sum-exp
# compares its arguments before combining them, exp(-0.0) is 1, and the
# Polya-Gamma sampler opens on |psi|, so -0.0 is absorbed at every consumer and
# x + 0.0 is x. Clearing an installed offset must return to the same path, so
# the sampler carries nothing of an offset it no longer has. ---
arm.zero <- parityArmCategoryOffset(zeroOffset, NULL)
expect_identical(arm.zero$train, arm.none$train)
expect_identical(arm.zero$forestFits, arm.none$forestFits)
expect_identical(arm.zero$runVarcount, arm.none$runVarcount)

set.seed(4242)
sampler.cleared <- buildSamplerCategoryOffset(offset)
sampler.cleared$setCategoryOffset(NULL, updateState = FALSE)
arm.cleared <- recordChannelsCategoryOffset(
  sampler.cleared,
  sampler.cleared$run(20L, 8L)
)
expect_identical(arm.cleared$train, arm.none$train)
expect_identical(arm.cleared$forestFits, arm.none$forestFits)

# and installing the zero matrix on a live sampler is likewise the null path
set.seed(4242)
sampler.zeroswap <- buildSamplerCategoryOffset(NULL)
sampler.zeroswap$setCategoryOffset(zeroOffset, updateState = FALSE)
expect_identical(
  sampler.zeroswap$run(20L, 8L)$train,
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
train.plain <- buildSamplerCategoryOffset(offset)$run(0L, 3L)$train
set.seed(1717)
train.rowshift <- buildSamplerCategoryOffset(rowShift)$run(0L, 3L)$train
set.seed(1717)
train.colshift <- buildSamplerCategoryOffset(columnShift)$run(0L, 3L)$train
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
sameVintage <- function(sampler, result, sampleNum, offsetMatrix) {
  fits <- vapply(
    seq_len(K),
    function(k) sampler$getForestFits(k)[, 1L],
    numeric(n)
  )
  max(abs(
    softmaxRows(fits + offsetMatrix) -
      matrix(result$train[,, sampleNum], n, K)
  ))
}

set.seed(313)
sampler.vintage <- buildSamplerCategoryOffset(offset)
res.vintage <- sampler.vintage$run(20L, 5L)
expect_true(sameVintage(sampler.vintage, res.vintage, 5L, offset) < 1e-12)

# second arm, against a fits rewrite the combiner never sees: a per-observation
# predictor session rebuilds every category forest's fits from its trees at
# MUTATION time, outside any sweep. A blend that refreshed its offset fits only
# from inside the sweep would report the pre-session fits for one more sample.
set.seed(515)
sampler.session <- buildSamplerCategoryOffset(NULL)
sampler.session$run(20L, 4L)
sampler.session$setCategoryOffset(offset, updateState = FALSE)
installed <- sampler.session$setPredictor(
  pmin(pmax(x[, 2L] + rnorm(n, 0, 0.02), 0), 1),
  2L,
  forceUpdate = "partial"
)
expect_true(sum(installed) > 0L)
res.session <- sampler.session$run(0L, 2L)
expect_true(sameVintage(sampler.session, res.session, 2L, offset) < 1e-12)

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
sampler.mn <- buildSamplerCategoryOffset(NULL)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no category offset, and both must name the
# family situation rather than the forest count
sampler.gaussian <- dbarts(x, rnorm(n), control = controlCategoryOffset())
expect_error(
  sampler.gaussian$setCategoryOffset(offset),
  "no count response"
)
set.seed(17)
z <- rbinom(n, 1L, 0.5)
sampler.bcf <- dbarts(
  x,
  rnorm(n),
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 10L)),
  control = controlCategoryOffset()
)
expect_error(
  sampler.bcf$setCategoryOffset(offset),
  "no count response"
)

# n and K are out of scope, and the refusal names both. A transposed matrix is
# the case a length test alone would install cell by cell into the wrong rows:
# the row-count check catches it, since the row count is what a transpose
# scrambles first
expect_error(
  sampler.mn$setCategoryOffset(offset[seq_len(n - 1L), ]),
  "same number of rows"
)
expect_error(
  sampler.mn$setCategoryOffset(t(offset)),
  "same number of rows"
)
# a flat vector is refused rather than recycled across the categories
expect_error(
  sampler.mn$setCategoryOffset(rep(0.5, n)),
  "numeric matrix"
)

# every non-finite entry is refused: an infinity propagates through the
# log-sum-exp margin into a NaN for every category of its row, and NA is not a
# shift.
for (bad in c(NA_real_, NaN, Inf, -Inf)) {
  spoiled <- offset
  spoiled[3L, 2L] <- bad
  expect_error(sampler.mn$setCategoryOffset(spoiled), "finite")
}

# a refusal leaves the sampler byte-identical: the entrance validates a whole
# scratch copy and swaps it in only once it holds, because the combiner borrows
# the installed buffer and an in-place write would BE the mutation
refusalArmCategoryOffset <- function(attempt) {
  set.seed(808)
  sampler <- buildSamplerCategoryOffset(NULL)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        sampler$setCategoryOffset(attempt, updateState = FALSE)
        NA_character_
      },
      error = conditionMessage
    )
  }
  c(
    list(refused = refused),
    recordChannelsCategoryOffset(sampler, sampler$run(15L, 5L))
  )
}
spoiled <- offset
spoiled[n, K] <- Inf
arm.refused <- refusalArmCategoryOffset(spoiled)
arm.untouched <- refusalArmCategoryOffset(NULL)
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
# built against the same seed placement as buildSamplerCategoryOffset's own
# creation, so the two draw streams line up and the train comparison below is
# a statement about test data rather than about creation order
set.seed(4242)
sampler.createTest <- dbarts(
  dbartsData(x, counts = counts, test = x.test),
  family = "multinomial",
  control = controlCategoryOffset()
)
sampler.createTest$setCategoryOffset(offset, updateState = FALSE)
res.createTest <- sampler.createTest$run(20L, 8L)
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
sampler.labelsTest <- dbarts(
  x,
  factor(labels, levels = seq.int(0L, K - 1L)),
  test = x.test,
  family = "multinomial",
  control = controlCategoryOffset()
)
expect_silent(sampler.labelsTest$setCategoryOffset(offset, updateState = FALSE))
# installing the train offset on a sampler that already holds test data
sampler.test <- dbarts(
  dbartsData(x, counts = counts, test = x.test),
  family = "multinomial",
  control = controlCategoryOffset()
)
expect_silent(sampler.test$setCategoryOffset(offset, updateState = FALSE))
# and the reverse direction: test data installs on a sampler already holding a
# train offset, whose test rows the train offset says nothing about
sampler.offset <- buildSamplerCategoryOffset(offset)
expect_silent(sampler.offset$setTestPredictor(x.test))
expect_true(all(is.finite(sampler.offset$run(0L, 3L)$test)))
# the combined entry takes the test rows too; its FLAT offset argument stays
# refused, and truthfully - after the blend it leaves the simplex, before it a
# common per-observation shift is inert
expect_error(
  sampler.offset$setTestPredictorAndOffset(x.test, rep(0, nTest)),
  "matrix"
)
expect_silent(
  sampler.offset$setTestPredictorAndOffset(x.test, NULL)
)
# predict is the one channel that still refuses under a train offset, and for a
# different reason: its rows are the caller's, so no resident offset describes
# them and none is substituted. An explicit matrix is the way through.
expect_error(
  sampler.offset$predict(x.test),
  "cannot be inferred"
)

# the response-side conduit's offset half now names the entry that works. The
# flat vector stays refused, and truthfully: a common per-observation shift is
# exactly the softmax's null direction, so it could only ever be inert.
expect_error(
  sampler.mn$setOffset(rep(0.5, n)),
  "n x K matrix"
)
expect_error(
  sampler.mn$setOffset(rep(0.5, n), updateScale = TRUE),
  "n x K matrix"
)
# predict's own flat-offset refusal is narrow rather than false, and says which
# form is not
expect_error(
  sampler.mn$predict(x.test, rep(0.5, nTest)),
  "matrix"
)

# --- The counts entrance's range check, which precedes the integer coercion.
# A count past .Machine$integer.max coerces to NA, and the engine would then
# refuse it as negative - a true refusal naming the wrong reason. ---
counts.huge <- matrix(as.double(counts), n, K)
counts.huge[1L, 1L] <- 2^31
expect_error(
  sampler.mn$setCounts(counts.huge),
  "representable as integers"
)
expect_error(
  dbarts(dbartsData(x, counts = counts.huge), family = "multinomial"),
  "representable as integers"
)
