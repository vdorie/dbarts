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

# --- the scale-pinned response swap. A BCF sampler is the one multi-forest
# coupling that admits setResponse (updateScale = FALSE only): the gaussian
# response re-maps y through the transform pinned at build, and the combiner
# re-derives both per-forest residuals from y every sweep, so the swap refits
# the same calibrated model against a new target. Two arms from one seed, one
# swapping and one not: the affine reported/internal map (recovered the way
# benchmarks/R/sbc.R does) must be identical across the swap, while the
# posterior must actually retarget. ---
yNew <- mu - z * tau + rnorm(n, sd = 0.2)

bcfMap <- function(bc) {
  reported <- dbarts:::bartcoreRun(bc, 0L, 1L)$train[, 1L]
  glue <- dbarts:::bartcoreBCFGlue(bc)
  internal <- glue[1L, 1L] *
    dbarts:::bartcoreForestFits(bc, 0L)[, 1L] +
    ifelse(z != 0, glue[3L, 1L], glue[2L, 1L]) *
      dbarts:::bartcoreForestFits(bc, 1L)[, 1L]
  fitScale <- stats::cov(reported, internal) / stats::var(internal)
  c(mean(reported) - fitScale * mean(internal), fitScale)
}

bcfSwapArm <- function(swap) {
  set.seed(101)
  host <- dbarts(x, y, control = control)
  bc <- dbarts:::bartcoreBCFSampler(host, z, n.trees.treatment = 25L)
  dbarts:::bartcoreRun(bc, 50L, 1L)
  before <- bcfMap(bc)
  if (swap) {
    dbarts:::bartcoreSetResponse(bc, yNew)
  }
  res <- dbarts:::bartcoreRun(bc, 50L, 20L)
  list(before = before, after = bcfMap(bc), fits = rowMeans(res$train))
}

arm.swap <- bcfSwapArm(TRUE)
arm.keep <- bcfSwapArm(FALSE)

expect_true(all(is.finite(arm.swap$fits)))
expect_equal(arm.swap$after, arm.swap$before, tolerance = 1e-6)
# the swapped chain tracks the new response, not the build one, and its
# posterior is nowhere near the arm that never swapped
expect_true(
  mean((arm.swap$fits - yNew)^2) < mean((arm.swap$fits - y)^2)
)
expect_true(mean(abs(arm.swap$fits - arm.keep$fits)) > 0.5)

# --- the per-forest leaf scale rides the state (docs/plans/multiforest-
# mutation-gaps.md item 3). BCF derives both forests' leaf scales from the
# response's SHAPE (the sd of the range-scaled y), so a destination built on a
# differently shaped response calibrates differently; the state now carries the
# scale like it already carried k, and both restore paths install it. ---
set.seed(11)
n.ls <- 200L
x.ls <- matrix(runif(n.ls * 3L), n.ls, 3L)
z.ls <- rbinom(n.ls, 1L, 0.5)
base.ls <- runif(n.ls)
# identical endpoints, different interior: the transform is the SAME on both
# (so a fit.scale guard would admit the divergent case) and only the shape,
# hence the leaf scale, moves
y.a <- base.ls
y.a[1L] <- 0
y.a[n.ls] <- 1
y.b <- 0.5 + 0.15 * (base.ls - 0.5)
y.b[1L] <- 0
y.b[n.ls] <- 1

control.ls <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 30L,
  updateState = FALSE
)
makeBC <- function(y) {
  host <- dbarts(x.ls, y, control = control.ls, resid.prior = fixed(0.2))
  dbarts:::bartcoreBCFSampler(host, z.ls, n.trees.treatment = 15L)
}

set.seed(101)
donor.ls <- makeBC(y.a)
dbarts:::bartcoreRun(donor.ls, 30L, 10L)
state.ls <- dbarts:::bartcoreStoreState(donor.ls)

# stored for every forest, positive, and in the calibration map's ratio:
# mu's is s / sqrt(m.mu) and tau's sdModerate s / (0.674 sqrt(m.tau)), so with
# sdModerate = 1 the ratio is sqrt(m.mu / m.tau) / 0.674
scale.mu <- state.ls[[1L]]$forests[[1L]]$leaf.scale
scale.tau <- state.ls[[1L]]$forests[[2L]]$leaf.scale
expect_true(is.numeric(scale.mu) && length(scale.mu) == 1L)
expect_true(scale.mu > 0 && scale.tau > 0)
expect_equal(scale.tau / scale.mu, sqrt(30 / 15) / 0.674, tolerance = 1e-8)

set.seed(101)
dest.same <- makeBC(y.a)
set.seed(101)
dest.shape <- makeBC(y.b)
# the arm is not vacuous: the two destinations calibrate differently, while the
# transform - what a fit.scale guard would compare - is identical
expect_true(
  dbarts:::bartcoreStoreState(dest.shape)[[1L]]$forests[[1L]]$leaf.scale !=
    scale.mu
)
expect_identical(
  dbarts:::bartcoreStoreState(dest.shape)[[1L]]$fit.scale,
  state.ls[[1L]]$fit.scale
)

# THE closure: one donor state into both destinations, responses then equalized
# (the conditioning-conduit pattern), and identical sweeps agree BITWISE
restoreArm <- function(bc, state) {
  dbarts:::bartcoreSetState(bc, state)
  dbarts:::bartcoreSetResponse(bc, y.a, FALSE)
  dbarts:::bartcoreRun(bc, 0L, 30L)$train
}
fits.same <- restoreArm(dest.same, state.ls)
fits.shape <- restoreArm(dest.shape, state.ls)
expect_identical(fits.shape, fits.same)

# cross-range: the anchor is invariant to an affine rescaling of y, so a 10*y
# destination was never miscalibrated - it stays admitted and consistent
set.seed(101)
dest.range <- makeBC(10 * y.a)
expect_equal(
  dbarts:::bartcoreStoreState(dest.range)[[1L]]$forests[[1L]]$leaf.scale,
  scale.mu
)
expect_identical(restoreArm(dest.range, state.ls), fits.same)

# an old state (the block stripped) restores with no error and reproduces
# PRE-change behavior exactly: the same-shape arm is untouched, and the
# different-shape arm diverges again, which is the defect this closes
state.old <- state.ls
for (f.ls in seq_along(state.old[[1L]]$forests)) {
  state.old[[1L]]$forests[[f.ls]]$leaf.scale <- NULL
}
set.seed(101)
old.same <- makeBC(y.a)
set.seed(101)
old.shape <- makeBC(y.b)
old.fits.same <- restoreArm(old.same, state.old)
old.fits.shape <- restoreArm(old.shape, state.old)
expect_identical(old.fits.same, fits.same)
expect_false(identical(old.fits.shape, old.fits.same))
expect_true(max(abs(old.fits.shape - old.fits.same)) > 0.1)

rm(
  n.ls,
  x.ls,
  z.ls,
  base.ls,
  y.a,
  y.b,
  control.ls,
  makeBC,
  donor.ls,
  state.ls,
  scale.mu,
  scale.tau,
  dest.same,
  dest.shape,
  restoreArm,
  fits.same,
  fits.shape,
  dest.range,
  state.old,
  f.ls,
  old.same,
  old.shape,
  old.fits.same,
  old.fits.shape
)
