# The caller-supplied per-forest, per-observation weight (the semantics live on
# Chain::setForestWeights, src/bartcore/chain.hpp): a multiplicative precision
# factor on one forest's own leaf conditionals, composing with the case weights
# so forest f draws against w_i * m_f^2 * s_i. It excludes a row from that
# forest's leaf conditionals WITHOUT excluding it from the model - the row keeps
# its occupancy, its place in the combination, and its residual degrees of
# freedom.
#
# Three things must hold, in decreasing order of how badly a failure would
# hurt: the channel is inert until a caller uses it (the null gate is the whole
# neutrality argument for shipping it); it is exact for the consumer it was
# built for; and it stays well defined in the limit where every row is
# excluded. The engine-level refusals, the substituted-response invariance of
# the excluded rows' node statistics, and the sigma degrees of freedom are
# pinned in tests/cpp, which can see them.
#
# samplePriorPredictive and sampleTreesFromPrior are deliberately NOT threaded
# through this weight: a prior tree draw reads the response model's own working
# weights and never forms a per-forest response, so it is membership-blind by
# construction, which is the right law for a draw from the prior. Neither is
# reachable from a BCF handle today in any case.

set.seed(1207)
n <- 300L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu.true <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau.true <- 1 + 2 * x[, 3L]
y <- mu.true + z * tau.true + rnorm(n, sd = 0.2)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)
sampler <- dbarts(x, y, control = control)

# One BCF run at a pinned seed, with an optional per-forest weight installed
# before the first sweep. The returned channels are everything the run reports
# plus both forests' fits, so "identical" here means the whole chain agreed,
# not just a summary.
runBCF <- function(
  seed,
  forest = NULL,
  weights = NULL,
  fixed.glue = FALSE,
  n.samples = 8L
) {
  set.seed(seed)
  bc <- dbarts:::bartcoreBCFSampler(
    sampler,
    z,
    n.trees.treatment = 25L,
    update.a = !fixed.glue,
    update.b = !fixed.glue
  )
  if (!is.null(forest)) {
    dbarts:::bartcoreSetForestWeights(bc, forest, weights)
  }
  result <- dbarts:::bartcoreRun(bc, 0L, n.samples)
  list(
    train = result$train,
    sigma = result$sigma,
    varcount = result$varcount,
    mu = dbarts:::bartcoreForestFits(bc, 0L),
    tau = dbarts:::bartcoreForestFits(bc, 1L),
    glue = dbarts:::bartcoreBCFGlue(bc)
  )
}

# --- the null gate: an all-ones weight is the identity, bitwise ---
# Installing it takes the composed path (a scratch buffer, an O(n) multiply)
# where not installing it passes the combiner's own pointer through, so this is
# the application path and the pass-through in one comparison.
base <- runBCF(4001L)
expect_identical(runBCF(4001L, 1L, rep(1, n)), base)
expect_identical(runBCF(4001L, 0L, rep(1, n)), base)

# --- the consumer case, both halves ---
# Under the fixed glue (1, 0, 1) the treatment forest already sees a multiplier
# of exactly zero at every control row, hence exactly zero weight there, so
# installing z as its per-forest weight changes nothing on ANY sweep, sweep 0
# included: 0 * 0 is +0 at the control rows and w * b1^2 * 1 is unchanged at the
# treated ones. It does not wait for a glue draw.
base.fixed <- runBCF(4002L, fixed.glue = TRUE)
expect_identical(
  runBCF(4002L, 1L, as.double(z), fixed.glue = TRUE),
  base.fixed
)

# The NEGATIVE half is mandatory, or this gets generalized into a phantom:
# under the DEFAULT glue the two routes exclude DIFFERENT information and must
# diverge. They still agree on sweep 0 alone, because b0 is exactly 0 at
# creation; from sweep 1 b0 is almost surely nonzero, so a control row carries
# real weight that the z route zeroes.
expect_identical(
  runBCF(4003L, 1L, as.double(z), n.samples = 1L),
  runBCF(4003L, n.samples = 1L)
)
expect_false(identical(runBCF(4003L, 1L, as.double(z)), runBCF(4003L)))

# --- every row excluded ---
# A legitimate transient while a caller resamples membership, so it must RUN
# rather than refuse. Forest 1's leaf conditionals then carry weight exactly
# zero, which is a well-defined draw from the leaf prior (tests/cpp pins that
# the marginal is exactly 0.0 on both sides of a birth, so the split decision is
# exactly the prior's) - not a NaN, not an Inf, and not a forced zero.
all.zero <- runBCF(4004L, 1L, rep(0, n), n.samples = 20L)
expect_true(all(is.finite(all.zero$train)))
expect_true(all(is.finite(all.zero$mu)) && all(is.finite(all.zero$tau)))
expect_true(all(is.finite(all.zero$glue)))
expect_true(all(all.zero$sigma > 0))
expect_true(sd(as.vector(all.zero$tau)) > 0)

# --- the weight toggles mid-chain ---
# The combiner caches nothing per-forest across sweeps and occupancy is
# count-based, so a leaf that was non-empty under one membership set is merely
# weightless under the next: nothing has to happen, and nothing does. This is
# where the count-based empty-leaf veto earns its keep, since no membership
# change can strand an empty leaf.
set.seed(4005)
bc.toggle <- dbarts:::bartcoreBCFSampler(sampler, z, n.trees.treatment = 25L)
invisible(dbarts:::bartcoreRun(bc.toggle, 20L, 5L))
toggled <- lapply(
  list(as.double(z), rep(1, n), rep(0, n), as.double(1 - z)),
  function(s) {
    dbarts:::bartcoreSetForestWeights(bc.toggle, 1L, s)
    dbarts:::bartcoreRun(bc.toggle, 0L, 5L)
  }
)
expect_true(all(vapply(
  toggled,
  function(r) all(is.finite(r$train)) && all(r$sigma > 0),
  logical(1L)
)))

# --- refusals ---
set.seed(4006)
bc <- dbarts:::bartcoreBCFSampler(sampler, z, n.trees.treatment = 25L)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc, 1L, rep(-1, n)),
  "non-negative"
)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc, 1L, c(NaN, rep(1, n - 1L))),
  "finite"
)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc, 1L, c(Inf, rep(1, n - 1L))),
  "finite"
)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc, 2L, rep(1, n)),
  "out of range"
)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc, 1L, rep(1, n - 1L)),
  "number of observations"
)

# the same refusals belong to the BRIDGE, not only to the R wrapper: a caller
# reaching the entry point directly is validated identically (safe over fast in
# R does not license an unguarded C entry)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setForestWeights, bc$ptr, 1L, rep(-1, n)),
  "non-negative"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setForestWeights,
    bc$ptr,
    1L,
    c(NA_real_, rep(1, n - 1L))
  ),
  "finite"
)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setForestWeights, bc$ptr, 2L, rep(1, n)),
  "out of range"
)

# a sampler with no such coupling is refused by the capability probe, which
# runs BEFORE any forest count - a K-forest multinomial carries forest 1 and
# would sail through a count test
bc.single <- dbarts:::bartcoreSampler(sampler)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc.single, 0L, rep(1, n)),
  "BCF"
)
set.seed(4007)
labels <- rbinom(n, 2L, 0.5)
bc.multinomial <- dbarts:::bartcoreMultinomialSampler(sampler, labels, K = 3L)
expect_error(
  dbarts:::bartcoreSetForestWeights(bc.multinomial, 1L, rep(1, n)),
  "BCF"
)

# an installed weight can never go stale against a changed n, because a
# multi-forest sampler refuses whole-data replacement outright; there is no
# length to re-check and no dangling borrow to invalidate
dbarts:::bartcoreSetForestWeights(bc, 1L, rep(1, n))
expect_error(
  dbarts:::bartcoreSetData(bc, sampler$data),
  "fixes its data at creation"
)

# --- saved trees and replay ---
# A weightless leaf's prior-drawn value is stored and replayed like any other,
# so keepTrees needs no branch of its own. Out-of-sample predict() stays refused
# under BCF for want of a test treatment vector (test-bcf.R), so the saved trees
# are the replay surface here.
control.keep <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  keepTrees = TRUE,
  n.samples = 10L
)
sampler.keep <- dbarts(x, y, control = control.keep)
set.seed(4008)
bc.keep <- dbarts:::bartcoreBCFSampler(
  sampler.keep,
  z,
  n.trees.treatment = 25L
)
dbarts:::bartcoreSetForestWeights(bc.keep, 1L, rep(0, n))
result.keep <- dbarts:::bartcoreRun(bc.keep, 10L, 10L)
expect_true(all(is.finite(result.keep$train)))
trees.tau <- dbarts:::bartcoreGetTrees(
  bc.keep,
  chainNums = 1L,
  sampleNums = seq_len(10L),
  treeNums = seq_len(25L),
  forest = 1L
)
expect_true(nrow(trees.tau) > 0L)
expect_true(all(is.finite(trees.tau$value)))
expect_true(all(trees.tau$n >= 0L))
