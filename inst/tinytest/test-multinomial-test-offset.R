# The multinomial category offset on the OUT-OF-SAMPLE side: an nTest x K
# matrix entering the reported test blend where the train offset enters the
# reported train blend, and a per-call
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

controlTestOffset <- function(
  n.chains = 1L,
  keepTrees = FALSE,
  n.samples = NULL
) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = n.chains,
    n.trees = 25L,
    keepTrees = keepTrees,
    n.samples = n.samples,
    updateState = FALSE
  )
}

buildSamplerTestOffset <- function(
  offset = NULL,
  offset.test = NULL,
  n.chains = 1L,
  keepTrees = FALSE,
  n.samples = NULL
) {
  sampler <- dbarts(
    dbartsData(x, counts = counts, test = x.test),
    family = "multinomial",
    control = controlTestOffset(n.chains, keepTrees, n.samples)
  )
  if (!is.null(offset)) {
    sampler$setCategoryOffset(offset, updateState = FALSE)
  }
  if (!is.null(offset.test)) {
    sampler$setCategoryTestOffset(offset.test, updateState = FALSE)
  }
  sampler
}

recordChannelsTestOffset <- function(sampler, result) {
  list(
    train = result$train,
    test = result$test,
    forestFits = lapply(seq_len(K), function(k) {
      sampler$getForestFits(k)
    }),
    runVarcount = result$varcount
  )
}

# --- Create-vs-swap parity on the test channel. Installing the test offset
# after creation must be the sampler built with it, bitwise: nothing derived
# from it is cached, and the blend rematerializes its offset fits at every
# report. ---
parityArmTestOffset <- function(build, swap, n.chains = 1L) {
  set.seed(4242)
  sampler <- buildSamplerTestOffset(offset, build, n.chains)
  if (!is.null(swap)) {
    sampler$setCategoryTestOffset(swap, updateState = FALSE)
  }
  recordChannelsTestOffset(sampler, sampler$run(20L, 8L))
}

arm.build <- parityArmTestOffset(testOffset, NULL)
arm.swap <- parityArmTestOffset(NULL, testOffset)
arm.none <- parityArmTestOffset(NULL, NULL)

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
arm.build.chains <- parityArmTestOffset(testOffset, NULL, n.chains = 2L)
arm.swap.chains <- parityArmTestOffset(NULL, testOffset, n.chains = 2L)
expect_identical(arm.swap.chains$test, arm.build.chains$test)
expect_identical(arm.swap.chains$train, arm.build.chains$train)

# --- Null-path neutrality, the same hard gate the train offset carries: an
# all-zero matrix is BITWISE the null path (x + 0.0 is x, and the softmax's
# log-sum-exp compares before it combines), and clearing an installed offset
# returns to it. ---
arm.zero <- parityArmTestOffset(zeroTestOffset, NULL)
expect_identical(arm.zero$test, arm.none$test)
expect_identical(arm.zero$train, arm.none$train)

set.seed(4242)
sampler.cleared <- buildSamplerTestOffset(offset, testOffset)
sampler.cleared$setCategoryTestOffset(NULL, updateState = FALSE)
expect_identical(
  sampler.cleared$run(20L, 8L)$test,
  arm.none$test
)

set.seed(4242)
sampler.zeroswap <- buildSamplerTestOffset(offset, NULL)
sampler.zeroswap$setCategoryTestOffset(zeroTestOffset, updateState = FALSE)
expect_identical(
  sampler.zeroswap$run(20L, 8L)$test,
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
arm.rowshift <- parityArmTestOffset(rowShift, NULL)
arm.colshift <- parityArmTestOffset(columnShift, NULL)
expect_equal(arm.rowshift$test, arm.build$test, tolerance = 1e-8)
expect_false(isTRUE(all.equal(arm.colshift$test, arm.build$test)))

# --- Predict-vs-run agreement, the oracle for the second half of the channel.
# The two replays never touch the combiner: they sum the saved trees into their
# own raw slab and softmax it directly. So a predict on the resident test rows,
# handed the same matrix the resident test offset holds, must reproduce the
# run's recorded test channel exactly - which it can only do if the offset was
# threaded into the replays and not only into the blend. ---
set.seed(313)
sampler.keep <- buildSamplerTestOffset(
  offset,
  testOffset,
  keepTrees = TRUE,
  n.samples = 6L
)
res.keep <- sampler.keep$run(20L, 6L)
pred.keep <- sampler.keep$predict(x.test, testOffset)
expect_identical(dim(pred.keep), dim(res.keep$test))
expect_identical(pred.keep, res.keep$test)

# and the replay reads its ARGUMENT, never the resident test offset: predicting
# the same rows with an all-zero matrix gives the offset-free surface, which is
# a different answer, and one that agrees with the same replay off a sampler
# that holds no test offset at all
pred.zero <- sampler.keep$predict(x.test, zeroTestOffset)
expect_false(isTRUE(all.equal(pred.zero, pred.keep)))
set.seed(313)
sampler.keep.none <- buildSamplerTestOffset(
  offset,
  NULL,
  keepTrees = TRUE,
  n.samples = 6L
)
res.keep.none <- sampler.keep.none$run(20L, 6L)
expect_identical(
  sampler.keep.none$predict(x.test, zeroTestOffset),
  pred.zero
)
expect_identical(res.keep.none$test, pred.zero)

# the predict offset is per CALL and per ROW: a subset of the rows takes the
# matching subset of the offset, and the answer is the same rows' answer
half <- seq_len(5L)
expect_identical(
  sampler.keep$predict(x.test[half, ], testOffset[half, ]),
  res.keep$test[half, , , drop = FALSE]
)

# --- Refusals. ---
sampler.plain <- buildSamplerTestOffset(NULL, NULL)

# the capability probe is not a forest count: a gaussian sampler and a BCF
# sampler (two forests) both own no category test offset, and both name the
# family situation
sampler.gaussian <- dbarts(
  x,
  rnorm(n),
  test = x.test,
  control = controlTestOffset()
)
expect_error(
  sampler.gaussian$setCategoryTestOffset(testOffset),
  "no count response"
)
set.seed(17)
z <- rbinom(n, 1L, 0.5)
sampler.bcf <- dbarts(
  x,
  rnorm(n),
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 10L)),
  control = controlTestOffset()
)
expect_error(
  sampler.bcf$setCategoryTestOffset(testOffset),
  "no count response"
)

# nTest and K are the current test store's, and the refusal names both. A
# column-count mismatch is caught R-side; a matching column count but wrong
# row count (the TRAIN offset shape, only accidentally n x K here) reaches the
# bridge's own row-count guard, which the R layer does not track.
expect_error(
  sampler.plain$setCategoryTestOffset(testOffset[, seq_len(2L)]),
  "3 categories"
)
expect_error(
  sampler.plain$setCategoryTestOffset(t(testOffset)),
  "3 categories"
)
# a flat vector is refused rather than recycled across the categories
expect_error(
  sampler.plain$setCategoryTestOffset(rep(0.5, nTest * K)),
  "numeric matrix"
)
# the TRAIN offset is not a test offset even though both are n x K matrices
# here only by the accident of the test rows being a subset
expect_error(
  sampler.plain$setCategoryTestOffset(offset),
  "14 observations x 3 categories"
)
# every non-finite entry is refused
for (bad in c(NA_real_, NaN, Inf, -Inf)) {
  spoiled <- testOffset
  spoiled[3L, 2L] <- bad
  expect_error(
    sampler.plain$setCategoryTestOffset(spoiled),
    "finite"
  )
}

# without test rows there is nothing for a per-test-row offset to describe, and
# accepting one would leave it silently unread
sampler.notest <- dbarts(
  dbartsData(x, counts = counts),
  family = "multinomial",
  control = controlTestOffset()
)
expect_error(
  sampler.notest$setCategoryTestOffset(testOffset),
  "requires test data"
)
# and clearing on such a sampler is a no-op rather than an error
expect_silent(sampler.notest$setCategoryTestOffset(NULL))

# a refusal leaves the sampler byte-identical: the entrance validates a whole
# scratch copy and swaps it in only once it holds, because the combiner borrows
# the installed buffer and an in-place write would BE the mutation
refusalArmTestOffset <- function(attempt) {
  set.seed(808)
  sampler <- buildSamplerTestOffset(offset, testOffset)
  refused <- if (is.null(attempt)) {
    NA_character_
  } else {
    tryCatch(
      {
        sampler$setCategoryTestOffset(attempt, updateState = FALSE)
        NA_character_
      },
      error = conditionMessage
    )
  }
  c(
    list(refused = refused),
    recordChannelsTestOffset(sampler, sampler$run(15L, 5L))
  )
}
spoiled <- testOffset
spoiled[nTest, K] <- Inf
arm.refused <- refusalArmTestOffset(spoiled)
arm.untouched <- refusalArmTestOffset(NULL)
expect_true(grepl("finite", arm.refused$refused))
expect_identical(arm.refused$test, arm.untouched$test)
expect_identical(arm.refused$train, arm.untouched$train)

# --- Replacing the test rows under an installed test offset is refused rather
# than reinterpreted: the offset describes the rows being replaced, and a row
# count that happens to match is not consent. Clearing first is the way
# through, on both test-predictor entries and on the removal form. ---
sampler.resident <- buildSamplerTestOffset(offset, testOffset)
expect_error(
  sampler.resident$setTestPredictor(x[seq_len(nTest) + nTest, ]),
  "clear it"
)
expect_error(sampler.resident$setTestPredictor(NULL), "clear it")
expect_error(
  sampler.resident$setTestPredictorAndOffset(NULL, NULL),
  "clear it"
)
expect_error(
  sampler.resident$setTestPredictorAndOffset(
    x[seq_len(nTest) + nTest, ],
    NULL
  ),
  "clear it"
)
sampler.resident$setCategoryTestOffset(NULL, updateState = FALSE)
expect_silent(
  sampler.resident$setTestPredictor(x[seq_len(nTest) + nTest, ])
)
expect_true(all(is.finite(sampler.resident$run(0L, 3L)$test)))

# --- The FLAT test offset stays refused on a softmax coupling, forever and
# truthfully: after the blend it moves the reported values off the simplex, and
# before it a common per-observation shift is the softmax's own null direction.
# The message now names the matrix entry that does work. ---
expect_error(
  sampler.plain$setTestOffset(rep(0.5, nTest)),
  "category test offset channel"
)
expect_error(
  sampler.plain$setTestPredictorAndOffset(x.test, rep(0.5, nTest)),
  "category test offset channel"
)

# --- Predict takes its offset from its ARGUMENT or refuses. A sampler holding
# a train category offset cannot answer for rows it has never seen without one
# being named; an explicit all-zero matrix is how the offset-free surface is
# asked for. A sampler holding no category offset keeps predicting as before. ---
set.seed(313)
sampler.pred <- buildSamplerTestOffset(
  offset,
  NULL,
  keepTrees = TRUE,
  n.samples = 4L
)
sampler.pred$run(15L, 4L)
expect_error(
  sampler.pred$predict(x.test),
  "cannot be inferred"
)
expect_silent(sampler.pred$predict(x.test, zeroTestOffset))
# the row count is the PREDICTED rows', not the sampler's
expect_error(
  sampler.pred$predict(x.test, zeroTestOffset[seq_len(5L), ]),
  "per-category matrix"
)
expect_error(
  sampler.pred$predict(x.test, matrix(NA_real_, nTest, K)),
  "missing values"
)
set.seed(313)
sampler.pred.none <- buildSamplerTestOffset(
  NULL,
  NULL,
  keepTrees = TRUE,
  n.samples = 4L
)
sampler.pred.none$run(15L, 4L)
expect_silent(sampler.pred.none$predict(x.test))
# and an offset supplied there is honored all the same
expect_false(isTRUE(all.equal(
  sampler.pred.none$predict(x.test, testOffset),
  sampler.pred.none$predict(x.test)
)))

# the refusal keys on EITHER resident offset, not the train one alone: a
# sampler carrying only a resident TEST offset (no train offset) must refuse a
# no-offset predict exactly as the train-offset-only sampler above does,
# rather than silently reporting the offset-free surface.
set.seed(313)
sampler.pred.testonly <- buildSamplerTestOffset(
  NULL,
  testOffset,
  keepTrees = TRUE,
  n.samples = 4L
)
sampler.pred.testonly$run(15L, 4L)
expect_error(
  sampler.pred.testonly$predict(x.test),
  "cannot be inferred"
)
expect_silent(sampler.pred.testonly$predict(
  x.test,
  zeroTestOffset
))
