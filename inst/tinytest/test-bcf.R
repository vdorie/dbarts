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

# the per-forest variable-count query (C3) works on both BCF forests: each
# forest's counts are nonnegative integers whose total is that forest's split
# count - positive here, both forests grew splits over the run
vcMu <- dbarts:::bartcoreForestVariableCounts(bcSampler, 0L)
vcTau <- dbarts:::bartcoreForestVariableCounts(bcSampler, 1L)
expect_equal(dim(vcMu), c(p, 1L))
expect_equal(dim(vcTau), c(p, 1L))
expect_true(is.integer(vcMu) && is.integer(vcTau))
expect_true(all(vcMu >= 0L) && all(vcTau >= 0L))
expect_true(sum(vcMu) > 0L && sum(vcTau) > 0L)
expect_error(
  dbarts:::bartcoreForestVariableCounts(bcSampler, 2L),
  "out of range"
)

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
# (forceUpdate = TRUE refreshes every forest; the transactional non-force paths
# are refused on multi-forest samplers, test-multi-forest-seam.R)
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

# --- moderators restriction on the treatment forest ---
# a named design so the subset can be given by name; the top-of-file sampler
# stays unnamed to exercise the names-without-colnames guard
x.mod <- x
colnames(x.mod) <- paste0("x", seq_len(p))
sampler.mod <- dbarts(x.mod, y, control = control)

# (a) resolution errors, each R-side before the bridge
expect_error(
  dbarts:::bartcoreBCFSampler(sampler.mod, z, moderators = "nope"),
  "not found"
)
expect_error(
  dbarts:::bartcoreBCFSampler(sampler.mod, z, moderators = 0L),
  "out of range"
)
expect_error(
  dbarts:::bartcoreBCFSampler(sampler.mod, z, moderators = p + 1L),
  "out of range"
)
expect_error(
  dbarts:::bartcoreBCFSampler(sampler.mod, z, moderators = integer(0)),
  "empty"
)
expect_error(
  dbarts:::bartcoreBCFSampler(sampler, z, moderators = "x1"),
  "no column names"
)

# (b) a restricted forest carries the run; sanity only, C2 owns the posterior
bcMod <- dbarts:::bartcoreBCFSampler(
  sampler.mod,
  z,
  n.trees.treatment = 25L,
  moderators = c("x1", "x3")
)
result.mod <- dbarts:::bartcoreRun(bcMod, 100L, 50L)
expect_equal(dim(result.mod$train), c(n, 50L))
expect_true(all(is.finite(result.mod$train)))
expect_true(all(result.mod$sigma > 0))
muMod <- dbarts:::bartcoreForestFits(bcMod, 0L)
tauMod <- dbarts:::bartcoreForestFits(bcMod, 1L)
expect_true(all(is.finite(muMod)) && all(is.finite(tauMod)))
expect_true(sum(muMod^2) > 0 && sum(tauMod^2) > 0)

# (c) default neutrality: an explicit moderators = NULL reproduces the omitted
# default bitwise (a fixed-seed R-side echo of the equivalence gate)
set.seed(20)
bcOmit <- dbarts:::bartcoreBCFSampler(sampler.mod, z, n.trees.treatment = 25L)
fit.omit <- dbarts:::bartcoreRun(bcOmit, 20L, 20L)$train
set.seed(20)
bcNull <- dbarts:::bartcoreBCFSampler(
  sampler.mod,
  z,
  n.trees.treatment = 25L,
  moderators = NULL
)
fit.null <- dbarts:::bartcoreRun(bcNull, 20L, 20L)$train
expect_identical(fit.null, fit.omit)

# (d) the forest selector on getTrees makes the restriction observable: every
# split variable the tau forest reports lies in the moderator set {x1, x3}
# (columns 1, 3), while the unrestricted mu forest splits somewhere outside it
# - proof the selector addresses different forests. bcMod runs live trees
# (no keepTrees), so query current = TRUE. var is 1-based; leaves report -1.
tauTrees <- dbarts:::bartcoreGetTrees(
  bcMod,
  chainNums = 1L,
  treeNums = seq_len(25L),
  current = TRUE,
  forest = 1L
)
tauSplits <- tauTrees$var[tauTrees$var > 0L]
expect_true(length(tauSplits) > 0L)
expect_true(all(tauSplits %in% c(1L, 3L)))

muTrees <- dbarts:::bartcoreGetTrees(
  bcMod,
  chainNums = 1L,
  treeNums = seq_len(50L),
  current = TRUE,
  forest = 0L
)
muSplits <- muTrees$var[muTrees$var > 0L]
expect_true(any(!(muSplits %in% c(1L, 3L))))

# an out-of-range forest index errors cleanly (bridge-side, as for forest fits)
expect_error(
  dbarts:::bartcoreGetTrees(
    bcMod,
    chainNums = 1L,
    treeNums = 1L,
    current = TRUE,
    forest = 2L
  ),
  "out of range"
)

# the tau forest's variable-count query (C3) sees the same column restriction
# the getTrees selector does: counts outside the moderator subset {x1, x3}
# (columns 1, 3; R rows 1, 3) are exactly zero, a sharp mask assertion, while
# the unrestricted mu forest is free to split outside it (mu depends on x2)
vcTauMod <- dbarts:::bartcoreForestVariableCounts(bcMod, 1L)
expect_equal(dim(vcTauMod), c(p, 1L))
expect_true(all(vcTauMod[c(2L, 4L), 1L] == 0L))
expect_true(sum(vcTauMod[c(1L, 3L), 1L]) > 0L)
vcMuMod <- dbarts:::bartcoreForestVariableCounts(bcMod, 0L)
expect_true(sum(vcMuMod[c(2L, 4L), 1L]) > 0L)
