# Internal BCF two-forest surface (docs/design/bcf.md; src/bartcore/). Sanity
# only - creation, a short run, sane glue and per-forest fits, setTreatment,
# and the step-4 state refusal. The exact-posterior gate lives in benchmarks/.

set.seed(3)
n <- 300L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)
sampler <- dbarts(x, y, control = control)

bcSampler <- dbarts:::bartcoreBCFSampler(sampler, z, n.trees.treatment = 25L)

result <- dbarts:::bartcoreRun(bcSampler, 100L, 100L)
expect_equal(dim(result$train), c(n, 100L))
expect_true(all(is.finite(result$train)))
expect_true(all(result$sigma > 0))

# both forests moved off zero and stay finite
muFits <- dbarts:::bartcoreForestFits(bcSampler, 0L)
tauFits <- dbarts:::bartcoreForestFits(bcSampler, 1L)
expect_equal(dim(muFits), c(n, 1L))
expect_true(all(is.finite(muFits)) && all(is.finite(tauFits)))
expect_true(sum(muFits^2) > 0 && sum(tauFits^2) > 0)

# glue is finite and the treated and control scales separate
glue <- dbarts:::bartcoreBCFGlue(bcSampler)
expect_equal(dim(glue), c(3L, 1L))
expect_true(all(is.finite(glue)))
expect_true(glue[2L, 1L] != glue[3L, 1L])

# setTreatment re-forms both residuals; a subsequent run stays sane
dbarts:::bartcoreSetTreatment(bcSampler, rep(0L, n))
result.control <- dbarts:::bartcoreRun(bcSampler, 0L, 5L)
expect_true(all(is.finite(result.control$train)))

# out-of-range forest index errors
expect_error(dbarts:::bartcoreForestFits(bcSampler, 2L), "out of range")

# the single-forest test-fit and prediction surface is undefined under BCF (no
# test treatment vector): setTestPredictor and predict are refused, pointing at
# the per-forest channels (bcf-testfits-guard)
expect_error(
  dbarts:::bartcoreSetTestPredictor(bcSampler, x[1:5, , drop = FALSE]),
  "BCF"
)
expect_error(
  dbarts:::bartcorePredict(bcSampler, x[1:5, , drop = FALSE]),
  "BCF"
)

# state round-trip: store, restore into a fresh BCF sampler, continue
dbarts:::bartcoreSetTreatment(bcSampler, z)
dbarts:::bartcoreRun(bcSampler, 0L, 5L)
state <- dbarts:::bartcoreStoreState(bcSampler)
expect_equal(length(state), 1L)
expect_equal(length(state[[1L]]$forests), 2L)
expect_false(is.null(state[[1L]]$bcf))

glueBefore <- dbarts:::bartcoreBCFGlue(bcSampler)
muBefore <- dbarts:::bartcoreForestFits(bcSampler, 0L)
tauBefore <- dbarts:::bartcoreForestFits(bcSampler, 1L)

restored <- dbarts:::bartcoreBCFSampler(sampler, z, n.trees.treatment = 25L)
dbarts:::bartcoreSetState(restored, state)

# the glue rides the state exactly; the forests restore to a continuation
# (structural, not bitwise: the dropped accumulation history is not reproduced)
expect_equal(dbarts:::bartcoreBCFGlue(restored), glueBefore)
expect_equal(
  dbarts:::bartcoreForestFits(restored, 0L),
  muBefore,
  tolerance = 1e-5
)
expect_equal(
  dbarts:::bartcoreForestFits(restored, 1L),
  tauBefore,
  tolerance = 1e-5
)

result.restored <- dbarts:::bartcoreRun(restored, 0L, 50L)
expect_equal(dim(result.restored$train), c(n, 50L))
expect_true(all(is.finite(result.restored$train)))
expect_true(all(result.restored$sigma > 0))

# fixed-glue path: update.a = update.b = FALSE holds the glue at (1, 0, 1)
bcFixed <- dbarts:::bartcoreBCFSampler(
  sampler,
  z,
  n.trees.treatment = 25L,
  update.a = FALSE,
  update.b = FALSE
)
dbarts:::bartcoreRun(bcFixed, 50L, 50L)
expect_equal(dbarts:::bartcoreBCFGlue(bcFixed)[, 1L], c(1, 0, 1))
# the treatment forest still moves under the fixed z * tau model
expect_true(sum(dbarts:::bartcoreForestFits(bcFixed, 1L)^2) > 0)

# bartCause-style driver: pihat is a prognostic column; the treatment and the
# propensity column are both swapped between runs through the mutation surface
set.seed(7)
pihat <- plogis(x[, 1L] - 0.5)
x.pi <- cbind(x, pihat)
sampler.pi <- dbarts(x.pi, y, control = control)
bcPi <- dbarts:::bartcoreBCFSampler(sampler.pi, z, n.trees.treatment = 25L)
dbarts:::bartcoreRun(bcPi, 100L, 20L)

z2 <- rbinom(n, 1L, pihat)
x.pi[, ncol(x.pi)] <- plogis(x[, 2L] - 0.5)
dbarts:::bartcoreSetPredictor(bcPi, x.pi, forceUpdate = TRUE)
dbarts:::bartcoreSetTreatment(bcPi, z2)
result.pi <- dbarts:::bartcoreRun(bcPi, 0L, 20L)
expect_equal(dim(result.pi$train), c(n, 20L))
expect_true(all(is.finite(result.pi$train)))
expect_true(all(result.pi$sigma > 0))

# --- interweaving glue-ridge move (docs/plans/bcf-ridge-interweaving.md) ---
# The move rescales (a, mu) -> (a/c, c mu) along the likelihood ridge after
# every glue draw under update.a = TRUE (the default), so the runs above
# already exercise it. It is posterior-preserving; these checks pin its
# behaviour through the R stack. The off path (update.a = FALSE) consumes no
# rng and was verified bitwise identical to the pre-change build cross-build
# (docs/plans/bcf-ridge-interweaving.md Status); update.a = FALSE sanity is
# covered by bcFixed above. The exact invariance and keepTrees saved-slot
# correctness are the C++ gates (tests/cpp).
set.seed(101)
n.m <- 200L
x.m <- matrix(runif(n.m * 3L), n.m, 3L)
z.m <- rbinom(n.m, 1L, 0.5)
y.m <- (2 * sin(pi * x.m[, 1L]) + x.m[, 2L]) +
  z.m * (1 + 2 * x.m[, 3L]) +
  rnorm(n.m, sd = 0.2)
control.m <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 60L,
  updateState = FALSE
)
sampler.m <- dbarts(x.m, y.m, control = control.m)

# move active: a longer run stays sane and the glue stays finite
bcMove <- dbarts:::bartcoreBCFSampler(sampler.m, z.m, n.trees.treatment = 30L)
res.move <- dbarts:::bartcoreRun(bcMove, 200L, 100L)
expect_true(all(is.finite(res.move$train)))
expect_true(all(res.move$sigma > 0))
expect_true(all(is.finite(dbarts:::bartcoreBCFGlue(bcMove))))
muMove <- dbarts:::bartcoreForestFits(bcMove, 0L)
expect_true(all(is.finite(muMove)) && sum(muMove^2) > 0)

# keepTrees with the move: storing the mu forest's saved slots (each rescaled
# by the sweep's c) leaves the run sane through the R stack
control.k <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 60L,
  updateState = FALSE,
  keepTrees = TRUE,
  n.samples = 50L
)
sampler.k <- dbarts(x.m, y.m, control = control.k)
bcKeep <- dbarts:::bartcoreBCFSampler(sampler.k, z.m, n.trees.treatment = 30L)
res.keep <- dbarts:::bartcoreRun(bcKeep, 100L, 50L)
expect_equal(dim(res.keep$train), c(n.m, 50L))
expect_true(all(is.finite(res.keep$train)))
expect_true(all(res.keep$sigma > 0))
