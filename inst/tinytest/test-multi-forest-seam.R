# The combined-output location-count seam and the multi-forest mutation guard
# (docs/plans/multinomial.md, C2). The location count is 1 for every model
# today, so the recorded train/test fits keep their single-location shape; this
# file is the shape guard the bitwise anchors do not provide (equivalence reads
# only summaries and never the train shape). It also pins the mutation guard
# that closes a latent BCF hazard: a whole-data mutation rebuilds forest 0 only.

set.seed(2718)
n <- 120L
p <- 4L
x <- matrix(runif(n * p), n, p)
f <- 3 * sin(pi * x[, 1L]) + 2 * x[, 2L]
y <- f + rnorm(n)
nTest <- 15L
x.test <- matrix(runif(nTest * p), nTest, p)

# --- single-chain train/test fits keep the L = 1 shape (n x numSamples), with
# no dimnames attached ---
control.one <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 40L,
  updateState = FALSE
)
sampler.one <- dbarts(x, y, test = x.test, control = control.one)
bc.one <- dbarts:::bartcoreSampler(sampler.one)
numSamples <- 30L
result.one <- dbarts:::bartcoreRun(bc.one, 50L, numSamples)

expect_equal(dim(result.one$train), c(n, numSamples))
expect_equal(dim(result.one$test), c(nTest, numSamples))
expect_null(dimnames(result.one$train))
expect_null(dimnames(result.one$test))

# --- multi-chain fits gain the trailing chain dimension only
# (n x numSamples x numChains); still no location dimension and no dimnames ---
numChains <- 2L
control.many <- dbartsControl(
  n.chains = numChains,
  n.threads = numChains,
  n.trees = 40L,
  updateState = FALSE
)
sampler.many <- dbarts(x, y, test = x.test, control = control.many)
bc.many <- dbarts:::bartcoreSampler(sampler.many)
result.many <- dbarts:::bartcoreRun(bc.many, 50L, numSamples)

expect_equal(dim(result.many$train), c(n, numSamples, numChains))
expect_equal(dim(result.many$test), c(nTest, numSamples, numChains))
expect_null(dimnames(result.many$train))
expect_null(dimnames(result.many$test))

# --- multi-forest mutation guard. setData, setResponse, and setWeights rebuild
# forest 0 only (applyNewData and the response refresh), so they are refused
# outright on a BCF (two-forest) sampler. The TRANSACTIONAL predictor paths
# (setPredictor / updatePredictor without forceUpdate, and the per-observation
# sessions, which have no force variant) validate through revalidateAllChains -
# also forest 0 only - and are refused; the FORCE paths refresh every forest
# and stay supported (the bartCause propensity swap, test-bcf.R). setTreatment,
# the supported multi-forest data swap, stays allowed ---
set.seed(11)
z <- rbinom(n, 1L, 0.5)
y.bcf <- f + z * (1 + 2 * x[, 3L]) + rnorm(n, sd = 0.2)
sampler.bcf.host <- dbarts(x, y.bcf, control = control.one)
bc.bcf <- dbarts:::bartcoreBCFSampler(
  sampler.bcf.host,
  z,
  n.trees.treatment = 20L
)

expect_error(
  dbarts:::bartcoreSetData(bc.bcf, sampler.bcf.host$data),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetResponse(bc.bcf, y.bcf + 1),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetWeights(bc.bcf, runif(n, 0.5, 1.5)),
  "multi-forest"
)

# transactional (non-force) predictor updates refuse ...
expect_error(
  dbarts:::bartcoreSetPredictor(bc.bcf, x + 0.01),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreUpdatePredictor(bc.bcf, x[, 1L] + 0.01, 1L),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreUpdatePredictorPerObservation(bc.bcf, x[, 1L], 1L),
  "multi-forest"
)
# ... while the force path refreshes every forest and stays supported
expect_true(dbarts:::bartcoreSetPredictor(bc.bcf, x + 0, forceUpdate = TRUE))

# setTreatment is still allowed and a subsequent run stays sane
dbarts:::bartcoreSetTreatment(bc.bcf, rbinom(n, 1L, 0.5))
result.bcf <- dbarts:::bartcoreRun(bc.bcf, 0L, 5L)
expect_equal(dim(result.bcf$train), c(n, 5L))
expect_true(all(is.finite(result.bcf$train)))

# the same guard is inert on a single-forest sampler: these mutations still work
expect_silent(dbarts:::bartcoreSetResponse(bc.one, y + 1))
expect_true(dbarts:::bartcoreSetPredictor(bc.one, x + 0))

# --- the per-forest variable-count query (C3). On a single-forest sampler the
# reported forest is forest 0, so the current-state query equals the recorded
# varcount channel of a length-1 run: with n.thin = 1 the run leaves the trees
# at the recorded sample's state (no sweep past the last storeSample). The
# recorded channel is per-sample and the query is live, so they coincide only
# right after such a run. ---
one.sample <- dbarts:::bartcoreRun(bc.one, 0L, 1L)
vc.one <- dbarts:::bartcoreForestVariableCounts(bc.one, 0L)
expect_equal(dim(vc.one), c(p, 1L))
expect_true(is.integer(vc.one))
expect_equal(vc.one[, 1L], one.sample$varcount[, 1L])
# an out-of-range forest index errors, as for the forest-fits query
expect_error(
  dbarts:::bartcoreForestVariableCounts(bc.one, 1L),
  "out of range"
)
