# Internal BCF two-forest surface (src/bartcore/). Sanity
# only - creation, a short run, sane glue and per-forest fits, setForestBasis,
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

bcSampler <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
  control = control
)

result <- bcSampler$run(100L, 100L)
expect_equal(dim(result$train), c(n, 100L))
expect_true(all(is.finite(result$train)))
expect_true(all(result$sigma > 0))

# both forests moved off zero and stay finite
muFits <- bcSampler$getForestFits(1L)
tauFits <- bcSampler$getForestFits(2L)
expect_equal(dim(muFits), c(n, 1L))
expect_true(all(is.finite(muFits)) && all(is.finite(tauFits)))
expect_true(sum(muFits^2) > 0 && sum(tauFits^2) > 0)

# the per-forest variable-count query works on both BCF forests: each
# forest's counts are nonnegative integers whose total is that forest's split
# count - positive here, both forests grew splits over the run
vcMu <- bcSampler$getForestVariableCounts(1L)
vcTau <- bcSampler$getForestVariableCounts(2L)
expect_equal(dim(vcMu), c(p, 1L))
expect_equal(dim(vcTau), c(p, 1L))
expect_true(is.integer(vcMu) && is.integer(vcTau))
expect_true(all(vcMu >= 0L) && all(vcTau >= 0L))
expect_true(sum(vcMu) > 0L && sum(vcTau) > 0L)
expect_error(
  bcSampler$getForestVariableCounts(3L),
  "out of range"
)

# glue is finite and the treated and control scales separate
glue <- bcSampler$getForestAmplitudes()
expect_equal(dim(glue), c(3L, 1L))
expect_true(all(is.finite(glue)))
expect_true(glue[2L, 1L] != glue[3L, 1L])

# an all-zero basis column has no rows for its amplitude to multiply:
# $setForestBasis refuses it before the engine or data@bases is touched
basesBefore <- bcSampler$data@bases
expect_error(
  bcSampler$setForestBasis(2L, cbind(rep(1, n), rep(0, n))),
  "all zeros"
)
expect_identical(bcSampler$data@bases, basesBefore)

# out-of-range forest index errors
expect_error(bcSampler$getForestFits(3L), "out of range")

# the single-forest test-fit and prediction surface is undefined here (the
# amplitudes have no off-sample basis): setTestPredictor and predict are
# refused, pointing at the per-forest channels
expect_error(
  bcSampler$setTestPredictor(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)
expect_error(
  bcSampler$predict(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)

# state round-trip: store, restore into a fresh BCF sampler, continue
bcSampler$setForestBasis(2L, cbind(1 - z, z))
bcSampler$run(0L, 5L)
bcSampler$storeState()
state <- bcSampler$state
expect_equal(length(state), 1L)
expect_equal(length(state[[1L]]$forests), 2L)
expect_false(is.null(state[[1L]]$glue))

glueBefore <- bcSampler$getForestAmplitudes()
muBefore <- bcSampler$getForestFits(1L)
tauBefore <- bcSampler$getForestFits(2L)

restored <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
  control = control
)
restored$setState(state)

# the glue rides the state exactly; the forests restore to a continuation
# whose fits are re-derived from the leaves, so they differ from the run's
# running totals by the additive rounding those totals accumulated only
expect_equal(restored$getForestAmplitudes(), glueBefore)
expect_equal(
  restored$getForestFits(1L),
  muBefore,
  tolerance = 1e-12
)
expect_equal(
  restored$getForestFits(2L),
  tauBefore,
  tolerance = 1e-12
)

result.restored <- restored$run(0L, 50L)
expect_equal(dim(result.restored$train), c(n, 50L))
expect_true(all(is.finite(result.restored$train)))
expect_true(all(result.restored$sigma > 0))

# fixed-glue path: update.a = update.b = FALSE holds the glue at (1, 0, 1)
bcFixed <- dbarts(
  x,
  y,
  forests = list(
    forest(update.amplitude = FALSE),
    forest(
      basis = ~ factor(z),
      n.trees = 25L,
      update.amplitude = FALSE
    )
  ),
  control = control
)
bcFixed$run(50L, 50L)
expect_equal(bcFixed$getForestAmplitudes()[, 1L], c(1, 0, 1))
# the treatment forest still moves under the fixed z * tau model
expect_true(sum(bcFixed$getForestFits(2L)^2) > 0)

# bartCause-style driver: pihat is a prognostic column; the treatment and the
# propensity column are both swapped between runs through the mutation surface
# (forceUpdate = TRUE refreshes every forest; the transactional non-force paths
# are refused on multi-forest samplers, test-multi-forest-seam.R)
set.seed(7)
pihat <- plogis(x[, 1L] - 0.5)
x.pi <- cbind(x, pihat)
bcPi <- dbarts(
  x.pi,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
  control = control
)
bcPi$run(100L, 20L)

z2 <- rbinom(n, 1L, pihat)
x.pi[, ncol(x.pi)] <- plogis(x[, 2L] - 0.5)
bcPi$setPredictor(x.pi, forceUpdate = TRUE)
bcPi$setForestBasis(2L, cbind(1 - z2, z2))
result.pi <- bcPi$run(0L, 20L)
expect_equal(dim(result.pi$train), c(n, 20L))
expect_true(all(is.finite(result.pi$train)))
expect_true(all(result.pi$sigma > 0))

# --- a longer run under the amplitude draw ---
# a and its half-Cauchy auxiliary are redrawn every sweep under update.a = TRUE
# (the default); a longer run stays sane through the R stack, with and without
# keepTrees. The cache and replay identities are the C++ gates (tests/cpp).
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
bcMove <- dbarts(
  x.m,
  y.m,
  forests = list(forest(), forest(basis = ~ factor(z.m), n.trees = 30L)),
  control = control.m
)
res.move <- bcMove$run(200L, 100L)
expect_true(all(is.finite(res.move$train)))
expect_true(all(res.move$sigma > 0))
expect_true(all(is.finite(bcMove$getForestAmplitudes())))
muMove <- bcMove$getForestFits(1L)
expect_true(all(is.finite(muMove)) && sum(muMove^2) > 0)

control.k <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 60L,
  updateState = FALSE,
  keepTrees = TRUE,
  n.samples = 50L
)
bcKeep <- dbarts(
  x.m,
  y.m,
  forests = list(forest(), forest(basis = ~ factor(z.m), n.trees = 30L)),
  control = control.k
)
res.keep <- bcKeep$run(100L, 50L)
expect_equal(dim(res.keep$train), c(n.m, 50L))
expect_true(all(is.finite(res.keep$train)))
expect_true(all(res.keep$sigma > 0))

# --- moderators restriction on the treatment forest ---
# a named design so the subset can be given by name; the plain 'x' stays
# unnamed to exercise the names-without-colnames guard
x.mod <- x
colnames(x.mod) <- paste0("x", seq_len(p))

# (a) resolution errors, each R-side before the bridge
bcModSpec <- function(x, moderators) {
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), vars = moderators)),
    control = control
  )
}
expect_error(bcModSpec(x.mod, "nope"), "not found")
expect_error(bcModSpec(x.mod, 0L), "out of range")
expect_error(bcModSpec(x.mod, p + 1L), "out of range")
expect_error(bcModSpec(x.mod, integer(0)), "empty")
# the plain, top-of-file 'x' is unnamed, to exercise the
# names-without-colnames guard
expect_error(bcModSpec(x, "x1"), "no column names")

# (b) a restricted forest carries the run; sanity only, posterior correctness
# is checked elsewhere
bcMod <- dbarts(
  x.mod,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z), n.trees = 25L, vars = c("x1", "x3"))
  ),
  control = control
)
result.mod <- bcMod$run(100L, 50L)
expect_equal(dim(result.mod$train), c(n, 50L))
expect_true(all(is.finite(result.mod$train)))
expect_true(all(result.mod$sigma > 0))
muMod <- bcMod$getForestFits(1L)
tauMod <- bcMod$getForestFits(2L)
expect_true(all(is.finite(muMod)) && all(is.finite(tauMod)))
expect_true(sum(muMod^2) > 0 && sum(tauMod^2) > 0)

# (c) default neutrality: a restriction naming every column reproduces the
# unrestricted forest bitwise (a fixed-seed R-side echo of the equivalence
# gate); an explicit vars = NULL would build the identical forest() object
set.seed(20)
bcOmit <- dbarts(
  x.mod,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
  control = control
)
fit.omit <- bcOmit$run(20L, 20L)$train
set.seed(20)
bcNull <- dbarts(
  x.mod,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z), n.trees = 25L, vars = colnames(x.mod))
  ),
  control = control
)
fit.null <- bcNull$run(20L, 20L)$train
expect_identical(fit.null, fit.omit)

# (d) the forest selector on getTrees makes the restriction observable: every
# split variable the tau forest reports lies in the moderator set {x1, x3}
# (columns 1, 3), while the unrestricted mu forest splits somewhere outside it
# - proof the selector addresses different forests. bcMod runs live trees
# (no keepTrees), so query current = TRUE. var is 1-based; leaves report -1.
# treeNums is left at its default on both reads: mu carries 50 trees (the
# control's n.trees) and tau carries its own 25, and each forest defaults to
# its OWN count rather than forest 1's.
tauTrees <- bcMod$getTrees(forest = 2L, chainNums = 1L, current = TRUE)
expect_true(all(tauTrees$forest == 2L))
expect_equal(range(tauTrees$tree), c(1L, 25L))
tauSplits <- tauTrees$var[tauTrees$var > 0L]
expect_true(length(tauSplits) > 0L)
expect_true(all(tauSplits %in% c(1L, 3L)))

muTrees <- bcMod$getTrees(forest = 1L, chainNums = 1L, current = TRUE)
expect_equal(range(muTrees$tree), c(1L, 50L))
muSplits <- muTrees$var[muTrees$var > 0L]
expect_true(any(!(muSplits %in% c(1L, 3L))))

# an out-of-range forest index errors cleanly, in R ahead of the .Call, with
# the same wording the sibling per-forest readers raise
expect_error(
  bcMod$getTrees(forest = 3L, chainNums = 1L, current = TRUE),
  "forest index out of range",
  fixed = TRUE
)

# an explicit treeNums is checked against EACH selected forest's own count:
# 1:50 is forest 1's whole range but exceeds tau's 25
expect_error(
  bcMod$getTrees(forest = 2L, treeNums = seq_len(50L), current = TRUE),
  "'treeNums' must be in [1, 25] for forest 2",
  fixed = TRUE
)

# forest as a vector stacks the named forests forest-major, each read at its
# own default tree count, exactly as rbinding the individual reads would
bothTrees <- bcMod$getTrees(forest = c(1L, 2L), chainNums = 1L, current = TRUE)
expect_equal(unique(bothTrees$forest), c(1L, 2L))
manualStack <- rbind(muTrees, tauTrees)
row.names(manualStack) <- row.names(bothTrees)
expect_equal(bothTrees, manualStack)
rm(bothTrees, manualStack)

# the tau forest's variable-count query sees the same column restriction
# the getTrees selector does: counts outside the moderator subset {x1, x3}
# (columns 1, 3; R rows 1, 3) are exactly zero, a sharp mask assertion, while
# the unrestricted mu forest is free to split outside it (mu depends on x2)
vcTauMod <- bcMod$getForestVariableCounts(2L)
expect_equal(dim(vcTauMod), c(p, 1L))
expect_true(all(vcTauMod[c(2L, 4L), 1L] == 0L))
expect_true(sum(vcTauMod[c(1L, 3L), 1L]) > 0L)
vcMuMod <- bcMod$getForestVariableCounts(1L)
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
  reported <- bc$run(0L, 1L)$train[, 1L]
  glue <- bc$getForestAmplitudes()
  internal <- glue[1L, 1L] *
    bc$getForestFits(1L)[, 1L] +
    ifelse(z != 0, glue[3L, 1L], glue[2L, 1L]) *
      bc$getForestFits(2L)[, 1L]
  fitScale <- stats::cov(reported, internal) / stats::var(internal)
  c(mean(reported) - fitScale * mean(internal), fitScale)
}

bcfSwapArm <- function(swap) {
  set.seed(101)
  bc <- dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), n.trees = 25L)),
    control = control
  )
  bc$run(50L, 1L)
  before <- bcfMap(bc)
  if (swap) {
    bc$setResponse(yNew)
  }
  res <- bc$run(50L, 20L)
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

# --- the per-forest leaf scale rides the state. BCF derives both forests'
# leaf scales from the response's SHAPE (the sd of the range-scaled y), so a
# destination built on a
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
  dbarts(
    x.ls,
    y,
    forests = list(forest(), forest(basis = ~ factor(z.ls), n.trees = 15L)),
    control = control.ls,
    family = gaussian(sigma = fixed(0.2))
  )
}

set.seed(101)
donor.ls <- makeBC(y.a)
donor.ls$run(30L, 10L)
donor.ls$storeState()
state.ls <- donor.ls$state

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
dest.shape$storeState()
expect_true(
  dest.shape$state[[1L]]$forests[[1L]]$leaf.scale != scale.mu
)
expect_identical(
  dest.shape$state[[1L]]$fit.scale,
  state.ls[[1L]]$fit.scale
)

# THE closure: one donor state into both destinations, responses then equalized
# (the conditioning-conduit pattern), and identical sweeps agree BITWISE
restoreArm <- function(bc, state) {
  bc$setState(state)
  bc$setResponse(y.a, updateScale = FALSE)
  bc$run(0L, 30L)$train
}
fits.same <- restoreArm(dest.same, state.ls)
fits.shape <- restoreArm(dest.shape, state.ls)
expect_identical(fits.shape, fits.same)

# cross-range: the anchor is invariant to an affine rescaling of y, so a 10*y
# destination was never miscalibrated - it stays admitted and consistent
set.seed(101)
dest.range <- makeBC(10 * y.a)
dest.range$storeState()
expect_equal(
  dest.range$state[[1L]]$forests[[1L]]$leaf.scale,
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
