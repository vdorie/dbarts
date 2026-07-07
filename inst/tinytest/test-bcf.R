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
