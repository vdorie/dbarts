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

# --- multi-forest mutation guard. setData rebuilds forest 0 only
# (applyNewData), so it is refused outright on a BCF (two-forest) sampler, as
# is setModel. setResponse is a whole-data
# mutation that is opt-in rather than refused: BCF's combiner re-derives every
# per-forest residual from y each sweep and the gaussian response re-maps y
# through the pinned transform, so at updateScale = FALSE the swap leaves
# nothing stale. updateScale = TRUE would move the transform while both leaf
# calibrations stayed stated against the old one, and NA is not FALSE, so both
# are still refused; so is a coupling that does not opt in (multinomial, whose
# setResponse is an empty override). Every transactional predictor path is
# ACCEPTED: the whole-matrix and column ones revalidate every forest, and the
# per-observation sessions guard every forest's trees that the column can move.
# The FORCE paths refresh every forest and stay supported (the bartCause
# propensity swap, test-bcf.R).
# setTreatment, the supported multi-forest data swap, stays allowed ---
set.seed(11)
z <- rbinom(n, 1L, 0.5)
y.bcf <- f + z * (1 + 2 * x[, 3L]) + rnorm(n, sd = 0.2)
sampler.bcf.host <- dbarts(x, y.bcf, control = control.one)
bc.bcf <- dbarts:::bartcoreBCFSampler(
  sampler.bcf.host,
  z,
  n.trees.treatment = 20L
)

# refused for memory safety too, not only staleness: the combiner indexes
# borrowed z over the live n, and setData has no z channel (see
# docs/plans/runsbcbcf-repair.md "setData door survey")
expect_error(
  dbarts:::bartcoreSetData(bc.bcf, sampler.bcf.host$data),
  "multi-forest"
)
expect_silent(dbarts:::bartcoreSetResponse(bc.bcf, y.bcf + 1))
expect_error(
  dbarts:::bartcoreSetResponse(bc.bcf, y.bcf + 1, updateScale = TRUE),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetResponse(bc.bcf, y.bcf + 1, updateScale = NA),
  "multi-forest"
)
# the case weights ride the same opt-in: Chain::setWeights is a pointer swap
# plus a positive-weight recount, BCF re-derives every per-forest response and
# precision from y and w each sweep, and the leaf calibration both forests are
# stated against is unweighted - so there is nothing stale and no scale to pin
expect_silent(dbarts:::bartcoreSetWeights(bc.bcf, runif(n, 0.5, 1.5)))
result.weights <- dbarts:::bartcoreRun(bc.bcf, 0L, 5L)
expect_true(all(is.finite(result.weights$train)))
# zero weights drop rows from the likelihood, so the swap must re-count the
# positive ones rather than carry the build count over
w.zero <- runif(n, 0.5, 1.5)
w.zero[seq_len(10L)] <- 0
expect_silent(dbarts:::bartcoreSetWeights(bc.bcf, w.zero))
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.bcf, 0L, 5L)$train)))
dbarts:::bartcoreSetWeights(bc.bcf, rep(1, n))
# setModel writes forest 0's leaf scale from the prior alone (Chain::
# setModel), which would silently discard BCF's calibrated mu leaf scale;
# refused regardless of the argument's class
expect_error(
  dbarts:::bartcoreSetModel(
    bc.bcf,
    sampler.bcf.host$model,
    sampler.bcf.host$control,
    sampler.bcf.host$data
  ),
  "multi-forest"
)

# transactional (non-force) predictor updates are ACCEPTED: the two-phase
# revalidation loops every forest, so a whole-matrix or column proposal
# installs only if no leaf of any tree of either forest would empty, and rolls
# the whole change back otherwise. Re-installing the current values cannot
# empty a leaf, so both accept; a run afterwards must stay finite, which is
# the assertion that every forest really was re-routed rather than left
# against stale codes
expect_true(dbarts:::bartcoreSetPredictor(bc.bcf, x))
expect_true(dbarts:::bartcoreUpdatePredictor(bc.bcf, x[, 1L], 1L))
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.bcf, 0L, 5L)$train)))
# a proposal that WOULD empty a leaf rolls back and reports FALSE rather than
# erroring - the veto, not a refusal. Collapsing a column onto two values of
# the existing grid empties leaves in every tree that splits on it
expect_false(
  dbarts:::bartcoreUpdatePredictor(
    bc.bcf,
    ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
    1L
  )
)
# the per-observation session is accepted too: its cell guard caches every
# forest, pruned to the trees the column can move, so a row installs only if it
# empties no leaf anywhere. Re-installing the column's own values moves
# nothing, so every row installs; the two-level collapse declines the rows that
# would empty a leaf - the per-row rollback - and the run stays finite
expect_true(all(
  dbarts:::bartcoreUpdatePredictorPerObservation(bc.bcf, x[, 1L], 1L)
))
expect_true(any(
  !dbarts:::bartcoreUpdatePredictorPerObservation(
    bc.bcf,
    ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
    1L
  )
))
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.bcf, 0L, 5L)$train)))
# ... and the force path refreshes every forest and stays supported
expect_true(dbarts:::bartcoreSetPredictor(bc.bcf, x + 0, forceUpdate = TRUE))

# setOffset rides the same conduit as setResponse under a different pointer
# (setOffset(yBuild - yNew, FALSE) re-maps through the pinned transform exactly
# as setResponse(yNew, FALSE) does), so it carries the same two conditions:
# permitted at FALSE, refused at TRUE and at NA. NULL clears the offset, which
# never moves the transform
expect_silent(dbarts:::bartcoreSetOffset(bc.bcf, rep(0.1, n)))
expect_error(
  dbarts:::bartcoreSetOffset(bc.bcf, rep(0.1, n), updateScale = TRUE),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetOffset(bc.bcf, rep(0.1, n), updateScale = NA),
  "multi-forest"
)
expect_silent(dbarts:::bartcoreSetOffset(bc.bcf, NULL))

# setTreatment is still allowed and a subsequent run stays sane
dbarts:::bartcoreSetTreatment(bc.bcf, rbinom(n, 1L, 0.5))
result.bcf <- dbarts:::bartcoreRun(bc.bcf, 0L, 5L)
expect_equal(dim(result.bcf$train), c(n, 5L))
expect_true(all(is.finite(result.bcf$train)))

# the weight swap is not merely permitted, it lands the same sampler: under
# resid.prior = fixed() (which removes the creation-time sigest, the one
# quantity a post-creation swap cannot reproduce - and cannot on a
# single-forest sampler either), building with w2 and building with w1 then
# swapping in w2 agree BITWISE on every channel the coupling exposes
w1 <- runif(n, 0.5, 1.5)
w2 <- runif(n, 0.5, 1.5)

bcfWeightArm <- function(build, swap) {
  set.seed(505)
  host <- dbarts(
    x,
    y.bcf,
    weights = build,
    resid.prior = fixed(1),
    control = control.one
  )
  bc <- dbarts:::bartcoreBCFSampler(host, z, n.trees.treatment = 20L)
  if (!is.null(swap)) {
    dbarts:::bartcoreSetWeights(bc, swap)
  }
  res <- dbarts:::bartcoreRun(bc, 20L, 10L)
  list(
    train = res$train,
    varcount = res$varcount,
    glue = dbarts:::bartcoreBCFGlue(bc),
    mu = dbarts:::bartcoreForestFits(bc, 0L),
    tau = dbarts:::bartcoreForestFits(bc, 1L),
    fit.scale = dbarts:::bartcoreStoreState(bc)[[1L]]$fit.scale
  )
}

arm.build <- bcfWeightArm(w2, NULL)
arm.swap <- bcfWeightArm(w1, w2)
arm.keep <- bcfWeightArm(w1, NULL)

expect_identical(arm.swap$train, arm.build$train)
expect_identical(arm.swap$varcount, arm.build$varcount)
expect_identical(arm.swap$glue, arm.build$glue)
expect_identical(arm.swap$mu, arm.build$mu)
expect_identical(arm.swap$tau, arm.build$tau)
# sanity: the two weightings are not the same posterior, so the arms above are
# not agreeing because the weights are inert
expect_false(isTRUE(all.equal(arm.keep$train, arm.build$train)))
# and the swap never moves the response transform - rescale() is unweighted, so
# there is no scale for setWeights to pin
expect_identical(arm.swap$fit.scale, arm.keep$fit.scale)

# a multi-forest coupling that does not opt in stays refused at every
# updateScale: multinomial's setResponse is an empty override, and its response
# is the borrowed count matrix a flat double vector cannot express
set.seed(23)
labels <- sample(0L:2L, n, replace = TRUE)
# the host carries test data so the test-offset refusal below is reached on a
# sampler that HAS a test channel to corrupt
sampler.mn.host <- dbarts(
  x,
  as.double(labels),
  test = x.test,
  control = control.one
)
bc.mn <- dbarts:::bartcoreMultinomialSampler(sampler.mn.host, labels, K = 3L)
expect_error(
  dbarts:::bartcoreSetResponse(bc.mn, as.double(labels)),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetResponse(bc.mn, as.double(labels), updateScale = TRUE),
  "multi-forest"
)
# a flat offset points exactly along the softmax's null direction (a common
# per-observation shift), so it has no semantics here at any updateScale
expect_error(dbarts:::bartcoreSetOffset(bc.mn, rep(0.5, n)), "multi-forest")
expect_error(
  dbarts:::bartcoreSetOffset(bc.mn, rep(0.5, n), updateScale = TRUE),
  "multi-forest"
)
expect_error(
  dbarts:::bartcoreSetOffset(bc.mn, rep(0.5, n), updateScale = NA),
  "multi-forest"
)
# case weights are refused at multinomial creation and stay refused after it:
# the opt-in that opens them for BCF is the combiner's, and this one does not
# take it
expect_error(
  dbarts:::bartcoreSetWeights(bc.mn, runif(n, 0.5, 1.5)),
  "multi-forest"
)
# z is defined only as the contrast the BCF glue forms b_{z_i} against; the
# capability probe catches a K-forest multinomial that a forest count would not
expect_error(
  dbarts:::bartcoreSetTreatment(bc.mn, rbinom(n, 1L, 0.5)),
  "requires a BCF sampler"
)
# a test offset is added AFTER the K forests are blended, so it would move the
# reported probabilities off the simplex; refused post-creation exactly as the
# multinomial constructor refuses one
expect_error(
  dbarts:::bartcoreSetTestOffset(bc.mn, rep(0.5, nTest)),
  "multi-forest"
)

# --- the multinomial predictor-mutation surface, one entry at a time: the
# whole-matrix/column entries install under revalidateAllChains, which loops
# every category forest, and the per-observation sessions install under a
# cell guard that caches every forest pruned to the trees the column can move.
# The forced entries and setCutPoints refresh every category forest
# throughout. ---
expect_true(dbarts:::bartcoreSetPredictor(bc.mn, x, forceUpdate = FALSE))
expect_true(
  dbarts:::bartcoreUpdatePredictor(bc.mn, x[, 1L], 1L, forceUpdate = FALSE)
)
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.mn, 0L, 5L)$train)))
expect_true(all(
  dbarts:::bartcoreUpdatePredictorPerObservation(bc.mn, x[, 1L], 1L)
))
expect_true(any(
  !dbarts:::bartcoreUpdatePredictorPerObservation(
    bc.mn,
    ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
    1L
  )
))
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.mn, 0L, 5L)$train)))
# the joint session installs in every sampler or none; a single-element list is
# the smallest case that reaches the same guard, and it now installs rather
# than refusing
expect_true(all(
  dbarts:::bartcoreUpdatePredictorPerObservationJointly(
    list(bc.mn),
    x[, 1L],
    1L
  )
))
expect_true(dbarts:::bartcoreSetPredictor(bc.mn, x, forceUpdate = TRUE))
expect_true(
  dbarts:::bartcoreUpdatePredictor(bc.mn, x[, 1L], 1L, forceUpdate = TRUE)
)
# setCutPoints carries no transactional guard at all: it refreshes every forest
# unconditionally, pruning whatever the coarsened grid orphans
expect_silent(dbarts:::bartcoreSetCutPoints(bc.mn, list(c(1 / 3, 2 / 3)), 1L))
expect_true(all(is.finite(dbarts:::bartcoreRun(bc.mn, 0L, 5L)$train)))

# the same guard is inert on a single-forest sampler: these mutations still
# work, including setOffset at updateScale = TRUE (rbart_vi's warmup rescale)
expect_silent(dbarts:::bartcoreSetResponse(bc.one, y + 1))
expect_silent(dbarts:::bartcoreSetWeights(bc.one, runif(n, 0.5, 1.5)))
expect_silent(dbarts:::bartcoreSetOffset(bc.one, rep(0.2, n)))
expect_silent(
  dbarts:::bartcoreSetOffset(bc.one, rep(0.2, n), updateScale = TRUE)
)
expect_silent(
  dbarts:::bartcoreSetModel(
    bc.one,
    sampler.one$model,
    sampler.one$control,
    sampler.one$data
  )
)
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
