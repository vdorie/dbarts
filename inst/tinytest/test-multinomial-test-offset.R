# The multinomial category offset on the OUT-OF-SAMPLE side
# (docs/design/multinomial.md): an nTest x K matrix entering the reported test
# blend where the train offset enters the reported train blend, and a per-call
# nNew x K matrix entering each predict replay's raw fits before the softmax.
#
# The one fact these oracles are built around is that the two sides are
# SEPARATE objects. The test rows are other rows, so no resident offset
# describes them: the test channel reports the surface at whatever test offset
# is installed (none meaning zero), and predict, whose rows are neither the
# train rows nor the test rows, refuses to guess rather than substitute one.
# The parity arms below therefore also check what must NOT move - the train
# channels are bitwise those of a sampler with no test offset at all, since the
# test fits enter no likelihood.

set.seed(9021)
n <- 90L
p <- 3L
K <- 3L
nTest <- 14L
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
x.test <- x[seq_len(nTest), , drop = FALSE]

# per-category, per-observation and with no common row component, so nothing
# about either offset lies along the softmax's null direction
offset <- cbind(
  0.8 * (x[, 1L] - 0.5),
  -0.6 + 0.4 * x[, 2L],
  0.3 * sin(6 * x[, 3L])
)
testOffset <- cbind(
  -0.5 + 0.9 * x.test[, 2L],
  0.4 * cos(5 * x.test[, 1L]),
  0.7 * (x.test[, 3L] - 0.5)
)
zeroTestOffset <- matrix(0, nTest, K)

control <- function(n.chains = 1L, keepTrees = FALSE, n.samples = NA_integer_) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    keepTrees = keepTrees,
    n.samples = n.samples,
    updateState = FALSE
  )
}

buildSampler <- function(
  offset = NULL,
  offset.test = NULL,
  n.chains = 1L,
  keepTrees = FALSE,
  n.samples = NA_integer_
) {
  host <- dbarts(
    x,
    as.double(labels),
    test = x.test,
    control = control(n.chains, keepTrees, n.samples)
  )
  dbarts:::bartcoreMultinomialCountSampler(
    host,
    counts,
    K = K,
    offset = offset,
    offset.test = offset.test
  )
}

recordChannels <- function(bc, result) {
  list(
    train = result$train,
    test = result$test,
    forestFits = lapply(seq_len(K) - 1L, function(k) {
      dbarts:::bartcoreForestFits(bc, k)
    }),
    runVarcount = result$varcount
  )
}

# --- Create-vs-swap parity on the test channel. Installing the test offset
# after creation must be the sampler built with it, bitwise: nothing derived
# from it is cached, and the blend rematerializes its offset fits at every
# report. ---
parityArm <- function(build, swap, n.chains = 1L) {
  set.seed(4242)
  bc <- buildSampler(offset, build, n.chains)
  if (!is.null(swap)) {
    dbarts:::bartcoreSetCategoryTestOffset(bc, swap)
  }
  recordChannels(bc, dbarts:::bartcoreRun(bc, 20L, 8L))
}

arm.build <- parityArm(testOffset, NULL)
arm.swap <- parityArm(NULL, testOffset)
arm.none <- parityArm(NULL, NULL)

expect_identical(arm.swap$test, arm.build$test)
expect_identical(arm.swap$train, arm.build$train)
expect_identical(arm.swap$forestFits, arm.build$forestFits)
expect_identical(arm.swap$runVarcount, arm.build$runVarcount)
# non-vacuity: the test offset moves the reported test probabilities, so the
# two arms above do not agree because it is inert
expect_false(isTRUE(all.equal(arm.none$test, arm.build$test)))
# and it moves NOTHING else: the test fits enter no likelihood, so every train
# channel is bitwise the no-test-offset sampler's
expect_identical(arm.build$train, arm.none$train)
expect_identical(arm.build$forestFits, arm.none$forestFits)
expect_identical(arm.build$runVarcount, arm.none$runVarcount)

# every chain sees the install
arm.build.chains <- parityArm(testOffset, NULL, n.chains = 2L)
arm.swap.chains <- parityArm(NULL, testOffset, n.chains = 2L)
expect_identical(arm.swap.chains$test, arm.build.chains$test)
expect_identical(arm.swap.chains$train, arm.build.chains$train)

# --- Null-path neutrality, the same hard gate the train offset carries: an
# all-zero matrix is BITWISE the null path (x + 0.0 is x, and the softmax's
# log-sum-exp compares before it combines), and clearing an installed offset
# returns to it. ---
arm.zero <- parityArm(zeroTestOffset, NULL)
expect_identical(arm.zero$test, arm.none$test)
expect_identical(arm.zero$train, arm.none$train)

set.seed(4242)
bc.cleared <- buildSampler(offset, testOffset)
dbarts:::bartcoreSetCategoryTestOffset(bc.cleared, NULL)
expect_identical(
  dbarts:::bartcoreRun(bc.cleared, 20L, 8L)$test,
  arm.none$test
)

set.seed(4242)
bc.zeroswap <- buildSampler(offset, NULL)
dbarts:::bartcoreSetCategoryTestOffset(bc.zeroswap, zeroTestOffset)
expect_identical(
  dbarts:::bartcoreRun(bc.zeroswap, 20L, 8L)$test,
  arm.none$test
)

# --- The simplex invariant with a test offset installed: the offset enters
# BEFORE the softmax, so the reported test values are still K probabilities per
# row. Added after the blend they would not be - the failure a test offset once
# produced on this same coupling. ---
expect_true(all(is.finite(arm.build$test)))
expect_true(all(arm.build$test >= 0 & arm.build$test <= 1))
expect_equal(
  apply(arm.build$test, c(1L, 3L), sum),
  matrix(1, nTest, 8L),
  tolerance = 1e-12
)

# a row-constant shift of the test offset is the softmax's null direction, and
# so leaves every reported test probability where it was; a one-column shift is
# not, and moves them
rowShift <- testOffset + matrix(rep(0.7 * x.test[, 1L], K), nTest, K)
columnShift <- testOffset
columnShift[, 2L] <- columnShift[, 2L] + 0.7
arm.rowshift <- parityArm(rowShift, NULL)
arm.colshift <- parityArm(columnShift, NULL)
expect_equal(arm.rowshift$test, arm.build$test, tolerance = 1e-8)
expect_false(isTRUE(all.equal(arm.colshift$test, arm.build$test)))

# --- Predict-vs-run agreement, the oracle for the second half of the channel.
# The two replays never touch the combiner: they sum the saved trees into their
# own raw slab and softmax it directly. So a predict on the resident test rows,
# handed the same matrix the resident test offset holds, must reproduce the
# run's recorded test channel exactly - which it can only do if the offset was
# threaded into the replays and not only into the blend. ---
set.seed(313)
bc.keep <- buildSampler(offset, testOffset, keepTrees = TRUE, n.samples = 6L)
res.keep <- dbarts:::bartcoreRun(bc.keep, 20L, 6L)
pred.keep <- dbarts:::bartcorePredict(bc.keep, x.test, testOffset)
expect_identical(dim(pred.keep), dim(res.keep$test))
expect_identical(pred.keep, res.keep$test)

# and the replay reads its ARGUMENT, never the resident test offset: predicting
# the same rows with an all-zero matrix gives the offset-free surface, which is
# a different answer, and one that agrees with the same replay off a sampler
# that holds no test offset at all
pred.zero <- dbarts:::bartcorePredict(bc.keep, x.test, zeroTestOffset)
expect_false(isTRUE(all.equal(pred.zero, pred.keep)))
set.seed(313)
bc.keep.none <- buildSampler(offset, NULL, keepTrees = TRUE, n.samples = 6L)
res.keep.none <- dbarts:::bartcoreRun(bc.keep.none, 20L, 6L)
expect_identical(
  dbarts:::bartcorePredict(bc.keep.none, x.test, zeroTestOffset),
  pred.zero
)
expect_identical(res.keep.none$test, pred.zero)

# the predict offset is per CALL and per ROW: a subset of the rows takes the
# matching subset of the offset, and the answer is the same rows' answer
half <- seq_len(5L)
expect_identical(
  dbarts:::bartcorePredict(bc.keep, x.test[half, ], testOffset[half, ]),
  res.keep$test[half, , , drop = FALSE]
)

# --- Refusals. ---
bc.plain <- buildSampler(NULL, NULL)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no category test offset, and both name the
# family situation
bc.gaussian <- dbarts:::bartcoreSampler(
  dbarts(x, rnorm(n), test = x.test, control = control())
)
expect_error(
  dbarts:::bartcoreSetCategoryTestOffset(bc.gaussian, testOffset),
  "requires a multinomial"
)
set.seed(17)
bc.bcf <- dbarts:::bartcoreBCFSampler(
  dbarts(x, rnorm(n), test = x.test, control = control()),
  rbinom(n, 1L, 0.5),
  n.trees.treatment = 10L
)
expect_error(
  dbarts:::bartcoreSetCategoryTestOffset(bc.bcf, testOffset),
  "requires a multinomial"
)

# nTest and K are the current test store's, and the refusal names both. A
# transposed matrix carries exactly nTest * K entries, so a length test alone
# would install it cell by cell into the wrong rows.
expect_error(
  dbarts:::bartcoreSetCategoryTestOffset(bc.plain, testOffset[, seq_len(2L)]),
  "category test offset"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCategoryTestOffset,
    bc.plain$ptr,
    t(testOffset)
  ),
  "14 observations x 3 categories"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCategoryTestOffset,
    bc.plain$ptr,
    rep(0.5, nTest * K)
  ),
  "14 observations x 3 categories"
)
# the TRAIN offset is not a test offset even though both are n x K matrices
# here only by the accident of the test rows being a subset
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setCategoryTestOffset,
    bc.plain$ptr,
    offset
  ),
  "14 observations x 3 categories"
)
# every non-finite entry is refused, at both layers
for (bad in c(NA_real_, NaN, Inf, -Inf)) {
  spoiled <- testOffset
  spoiled[3L, 2L] <- bad
  expect_error(
    dbarts:::bartcoreSetCategoryTestOffset(bc.plain, spoiled),
    "finite"
  )
  expect_error(
    .Call(
      dbarts:::C_dbarts_bartcore_setCategoryTestOffset,
      bc.plain$ptr,
      spoiled
    ),
    "finite"
  )
}

# without test rows there is nothing for a per-test-row offset to describe, and
# accepting one would leave it silently unread
host.notest <- dbarts(x, as.double(labels), control = control())
bc.notest <- dbarts:::bartcoreMultinomialCountSampler(
  host.notest,
  counts,
  K = K
)
expect_error(
  dbarts:::bartcoreSetCategoryTestOffset(bc.notest, testOffset),
  "requires test data"
)
expect_error(
  dbarts:::bartcoreMultinomialCountSampler(
    host.notest,
    counts,
    K = K,
    offset.test = testOffset
  ),
  "requires test data"
)
# and clearing on such a sampler is a no-op rather than an error
expect_silent(dbarts:::bartcoreSetCategoryTestOffset(bc.notest, NULL))

# a refusal leaves the sampler byte-identical: the entrance validates a whole
# scratch copy and swaps it in only once it holds, because the combiner borrows
# the installed buffer and an in-place write would BE the mutation
refusalArm <- function(attempt) {
  set.seed(808)
  bc <- buildSampler(offset, testOffset)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        dbarts:::bartcoreSetCategoryTestOffset(bc, attempt)
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
spoiled <- testOffset
spoiled[nTest, K] <- Inf
arm.refused <- refusalArm(spoiled)
arm.untouched <- refusalArm(NULL)
expect_true(grepl("finite", arm.refused$refused))
expect_identical(arm.refused$test, arm.untouched$test)
expect_identical(arm.refused$train, arm.untouched$train)

# --- Replacing the test rows under an installed test offset is refused rather
# than reinterpreted: the offset describes the rows being replaced, and a row
# count that happens to match is not consent. Clearing first is the way
# through, on both test-predictor entries and on the removal form. ---
bc.resident <- buildSampler(offset, testOffset)
expect_error(
  dbarts:::bartcoreSetTestPredictor(bc.resident, x[seq_len(nTest) + nTest, ]),
  "clear it"
)
# including the removal form, which the R wrapper cannot express (it matricizes
# its argument), so the C entry is called directly
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setTestPredictor, bc.resident$ptr, NULL),
  "clear it"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestPredictorAndOffset,
    bc.resident$ptr,
    x[seq_len(nTest) + nTest, ],
    NULL
  ),
  "clear it"
)
dbarts:::bartcoreSetCategoryTestOffset(bc.resident, NULL)
expect_silent(
  dbarts:::bartcoreSetTestPredictor(bc.resident, x[seq_len(nTest) + nTest, ])
)
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.resident, 0L, 3L)$test)))

# --- The FLAT test offset stays refused on a softmax coupling, forever and
# truthfully: after the blend it moves the reported values off the simplex, and
# before it a common per-observation shift is the softmax's own null direction.
# The message now names the matrix entry that does work. ---
expect_error(
  dbarts:::bartcoreSetTestOffset(bc.plain, rep(0.5, nTest)),
  "bartcore_setCategoryTestOffset"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestPredictorAndOffset,
    bc.plain$ptr,
    x.test,
    rep(0.5, nTest)
  ),
  "bartcore_setCategoryTestOffset"
)

# --- Predict takes its offset from its ARGUMENT or refuses. A sampler holding
# a train category offset cannot answer for rows it has never seen without one
# being named; an explicit all-zero matrix is how the offset-free surface is
# asked for. A sampler holding no category offset keeps predicting as before. ---
set.seed(313)
bc.pred <- buildSampler(offset, NULL, keepTrees = TRUE, n.samples = 4L)
dbarts:::bartcoreRun(bc.pred, 15L, 4L)
expect_error(
  dbarts:::bartcorePredict(bc.pred, x.test),
  "cannot be inferred"
)
expect_silent(dbarts:::bartcorePredict(bc.pred, x.test, zeroTestOffset))
# the row count is the PREDICTED rows', not the sampler's
expect_error(
  dbarts:::bartcorePredict(bc.pred, x.test, zeroTestOffset[seq_len(5L), ]),
  "14 x 3 matrix"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_predict,
    bc.pred$ptr,
    x.test,
    zeroTestOffset[seq_len(5L), ]
  ),
  "14 observations x 3 categories"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_predict,
    bc.pred$ptr,
    x.test,
    matrix(NA_real_, nTest, K)
  ),
  "finite"
)
set.seed(313)
bc.pred.none <- buildSampler(NULL, NULL, keepTrees = TRUE, n.samples = 4L)
dbarts:::bartcoreRun(bc.pred.none, 15L, 4L)
expect_silent(dbarts:::bartcorePredict(bc.pred.none, x.test))
# and an offset supplied there is honored all the same
expect_false(isTRUE(all.equal(
  dbarts:::bartcorePredict(bc.pred.none, x.test, testOffset),
  dbarts:::bartcorePredict(bc.pred.none, x.test)
)))

# the refusal keys on EITHER resident offset, not the train one alone: a
# sampler carrying only a resident TEST offset (no train offset) must refuse a
# no-offset predict exactly as the train-offset-only sampler above does,
# rather than silently reporting the offset-free surface.
set.seed(313)
bc.pred.testonly <- buildSampler(
  NULL,
  testOffset,
  keepTrees = TRUE,
  n.samples = 4L
)
dbarts:::bartcoreRun(bc.pred.testonly, 15L, 4L)
expect_error(
  dbarts:::bartcorePredict(bc.pred.testonly, x.test),
  "cannot be inferred"
)
expect_silent(dbarts:::bartcorePredict(
  bc.pred.testonly,
  x.test,
  zeroTestOffset
))
