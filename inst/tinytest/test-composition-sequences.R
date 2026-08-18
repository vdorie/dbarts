# Mutation-SEQUENCE composition: what a host does to a live sampler BETWEEN
# sweeps, on model classes the single-axis mutation tests never compose with.
# Seeded from the gate-blindspot audit's uncovered combinations; only the
# cells still open against the current suite are here, and each one carries a
# draw-level oracle rather than a shape or finiteness check.
#
# Two facts the arms below depend on and pin in passing:
#   * with control@seed set, a sampler's generators are independent of R's
#     stream, so two constructions at the same seed are bitwise twins and a
#     single-chain run reproduces chain 1 of a multi-chain run;
#   * a restore is SEMANTIC, not bitwise (inst/common/stateContinuation.R):
#     the canonical rebuild resums the accumulated fits, so a continued chain
#     reproduces the MOVES exactly and the numbers to within the resum.

set.seed(3103L)
n <- 120L
X <- data.frame(v1 = runif(n), v2 = runif(n), v3 = runif(n))
z <- rbinom(n, 1L, 0.5)
y <- 2 * X$v1 - X$v2 + z + rnorm(n, sd = 0.5)
wts <- runif(n, 0.5, 2)
groups <- factor(rep_len(seq_len(4L), n))
# built once so it can also ride ... into the warm-start donor builder
linearLeaf <- dbarts:::linear("v3")

ctl <- function(n.chains = 1L, n.trees = 8L, ...) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = n.trees,
    n.samples = 1L,
    n.burn = 0L,
    seed = 3L,
    updateState = FALSE,
    ...
  )
}
groupedCtl <- function(...) {
  # rbart_vi's internal control attribute is the only route to a grouped
  # in-core sampler
  cc <- ctl(...)
  attr(cc, "bartcore.groups") <- list(
    indices = as.integer(groups),
    n.groups = nlevels(groups),
    prior = "cauchy",
    rel.scale = sd(y),
    n.steps = 1L
  )
  cc
}
bcfSampler <- function(...) {
  dbarts(
    X,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = ctl(),
    ...
  )
}
# the trailing sweep dimension is last on every reported channel, whatever the
# rank (n x S, p x S, n x K x S, p x F x S)
sweepSlice <- function(a, idx) {
  d <- dim(a)
  ii <- rep(list(bquote()), length(d))
  ii[[length(d)]] <- idx
  do.call(`[`, c(list(a), ii, list(drop = FALSE)))
}
sweepTail <- function(a, keep) {
  s <- dim(a)[length(dim(a))]
  sweepSlice(a, seq.int(s - keep + 1L, s))
}


## --- 1. per-observation setPredictor x probit / logistic -------------------
## The partial session's finalize step re-derives leaf stats from the FAMILY
## WORKING RESPONSE (the probit latent z, the Polya-Gamma omega z), a path no
## test reached under a latent family. Oracle: the session must land the same
## sampler as a whole-column install of the column it actually reached, so the
## continuation is bitwise. The scan order is drawn from the chain's own
## generator (sampler.hpp updatePredictorPerObservation), so the twin spends
## one permutation of its own to stay generator-aligned.

set.seed(4201L)
nb <- 40L
vb <- rnorm(nb)
yb <- rbinom(nb, 1L, pnorm(vb))
binaryCtl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  n.samples = 1L,
  n.burn = 0L,
  seed = 7L,
  updateState = FALSE
)
binarySampler <- function(fam) {
  s <- dbarts(
    yb ~ vb,
    data.frame(vb = vb, yb = yb),
    control = binaryCtl,
    family = fam
  )
  invisible(s$run(40L, 1L))
  s
}

for (fam in c("probit", "logistic")) {
  live <- binarySampler(fam)
  twin <- binarySampler(fam)
  orig <- as.numeric(live$data@x[, 1L])
  # collapsing the column onto one value empties leaves, so the veto declines
  # exactly the last occupant of each: the session mixes installs and rollbacks
  collapsed <- rep(orig[1L], nb)
  mask <- live$setPredictor(collapsed, 1L, forceUpdate = "partial")
  expect_true(sum(mask) > 0L && sum(!mask) > 0L, info = fam)
  reached <- ifelse(mask, collapsed, orig)
  expect_identical(as.numeric(live$data@x[, 1L]), reached, info = fam)
  invisible(twin$setPredictor(orig, 1L, forceUpdate = "partial"))
  expect_true(
    isTRUE(twin$setPredictor(reached, 1L, forceUpdate = FALSE)),
    info = fam
  )
  liveRun <- live$run(0L, 5L)
  twinRun <- twin$run(0L, 5L)
  expect_identical(liveRun$train, twinRun$train, info = fam)
  expect_identical(liveRun$varcount, twinRun$varcount, info = fam)
  expect_true(all(is.finite(live$getLatents())), info = fam)

  # a session in which every moving row is declined must leave the chain
  # exactly where it stood: re-proposing the collapse can only move the leaf
  # survivors, and each of those would empty its leaf
  rolled <- binarySampler(fam)
  uncut <- binarySampler(fam)
  invisible(rolled$setPredictor(collapsed, 1L, forceUpdate = "partial"))
  held <- as.numeric(rolled$data@x[, 1L])
  invisible(uncut$setPredictor(collapsed, 1L, forceUpdate = "partial"))
  maskAgain <- rolled$setPredictor(collapsed, 1L, forceUpdate = "partial")
  expect_true(sum(!maskAgain) > 0L, info = fam)
  expect_identical(as.numeric(rolled$data@x[, 1L]), held, info = fam)
  # the control spends the same second scan-order draw on a no-op proposal
  invisible(uncut$setPredictor(held, 1L, forceUpdate = "partial"))
  expect_identical(
    rolled$run(0L, 5L)$train,
    uncut$run(0L, 5L)$train,
    info = fam
  )
}

# the per-observation entry point ADVANCES the generator even when it installs
# nothing (it draws its scan order); the whole-column entry point does not.
# A Gibbs host that reads the mask as an accept mask inherits that advance.
scanA <- binarySampler("probit")
scanB <- binarySampler("probit")
scanC <- binarySampler("probit")
scanOrig <- as.numeric(scanA$data@x[, 1L])
invisible(scanA$setPredictor(scanOrig, 1L, forceUpdate = "partial"))
invisible(scanC$setPredictor(scanOrig, 1L, forceUpdate = FALSE))
scanBase <- scanB$run(0L, 3L)$train
expect_false(identical(scanA$run(0L, 3L)$train, scanBase))
expect_identical(scanC$run(0L, 3L)$train, scanBase)


## --- 2. installTrees warm start x {multi-chain, linear leaves, missing} ----
## test-warm-start.R is single-chain, constant-leaf and complete-data
## throughout. The three cells here route the donor state through the
## per-chain sample map, the vector leaf-parameter layout, and the NA-carrying
## cut grid respectively.

warmDonor <- function(design, ...) {
  bart2(
    design,
    y,
    ...,
    n.trees = 8L,
    n.samples = 4L,
    n.burn = 20L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    keepSampler = TRUE,
    verbose = FALSE,
    seed = 1L
  )
}
donorSample <- function(d, s) d$fit$getTrees(sampleNums = s, chainNums = 1L)$var

donor <- warmDonor(X)

# multi-chain: one donor sample per chain, and each chain gets ITS sample
warmChains <- dbarts(X, y, control = ctl(n.chains = 2L))
warmChains$installTrees(donor, samples = c(1L, 3L))
donor1 <- donorSample(donor, 1L)
donor3 <- donorSample(donor, 3L)
expect_identical(
  warmChains$getTrees(current = TRUE, chainNums = 1L)$var,
  donor1
)
expect_identical(
  warmChains$getTrees(current = TRUE, chainNums = 2L)$var,
  donor3
)
expect_false(identical(donor1, donor3))
warmChainsAgain <- dbarts(X, y, control = ctl(n.chains = 2L))
warmChainsAgain$installTrees(donor, samples = c(1L, 3L))
expect_identical(
  warmChains$run(0L, 3L)$train,
  warmChainsAgain$run(0L, 3L)$train
)

# linear leaves: the donor's leaf parameters are vectors, and a constant-leaf
# donor is refused rather than silently reinterpreted
linDonor <- warmDonor(X, node.prior = linearLeaf)
warmLinear <- dbarts(X, y, node.prior = linearLeaf, control = ctl())
warmLinear$installTrees(linDonor, samples = 2L)
expect_identical(
  warmLinear$getTrees(current = TRUE)$var,
  donorSample(linDonor, 2L)
)
expect_true(all(is.finite(warmLinear$run(0L, 3L)$train)))
expect_error(
  warmLinear$installTrees(donor, samples = 2L),
  pattern = "malformed leaf parameters"
)

# missing predictors: the donor's splits remap onto the NA-carrying grid
Xna <- X
Xna$v2[c(3L, 17L, 40L)] <- NA_real_
naDonor <- warmDonor(Xna)
warmMissing <- dbarts(Xna, y, control = ctl())
warmMissing$installTrees(naDonor, samples = 2L)
expect_identical(
  warmMissing$getTrees(current = TRUE)$var,
  donorSample(naDonor, 2L)
)
expect_true(all(is.finite(warmMissing$run(0L, 3L)$train)))

# a two-forest destination cannot be seeded from a single-forest donor
expect_error(
  bcfSampler()$installTrees(donor, samples = 2L),
  pattern = "shape-compatible"
)


## --- 3. setCutPoints x {BCF, linear leaf covariate, missing} ---------------
## Existing setCutPoints coverage is single-forest constant-leaf complete-data
## (plus a categorical REFUSAL). The pin is draw-level: after the swap, every
## split drawn on the mutated column over many sweeps must lie on the NEW
## grid, and the column must still be drawable (the count is the non-vacuity
## arm - a column no longer split would satisfy the first claim vacuously).

drawnCuts <- function(sampler, column) {
  trees <- sampler$getTrees(current = TRUE)
  sort(unique(trees$value[trees$var == column]))
}

# a single-forest control first: the OLD grid becomes undrawable
gridBefore <- dbarts(X, y, control = ctl(n.trees = 20L))
invisible(gridBefore$run(50L, 1L))
oldCuts <- drawnCuts(gridBefore, 1L)
expect_true(length(oldCuts) > 1L)
gridBefore$setCutPoints(list(0.5), 1L)
gridRun <- gridBefore$run(0L, 40L)
expect_identical(drawnCuts(gridBefore, 1L), 0.5)
expect_true(sum(gridRun$varcount[1L, ]) > 0L)
# and a cut OUTSIDE the old grid becomes drawable
gridWiden <- dbarts(X, y, control = ctl(n.trees = 20L, n.cuts = 5L))
invisible(gridWiden$run(50L, 1L))
widened <- c(0.011, 0.987)
# the engine's own uniform grid for the column it is replacing: the new cuts
# lie outside it entirely, so drawing one is proof the swap took
oldGrid <- min(X$v2) + seq_len(5L) * ((max(X$v2) - min(X$v2)) / 6)
expect_true(all(widened < min(oldGrid) | widened > max(oldGrid)))
gridWiden$setCutPoints(list(widened), 2L)
widenRun <- gridWiden$run(0L, 60L)
expect_true(all(drawnCuts(gridWiden, 2L) %in% widened))
expect_true(sum(widenRun$varcount[2L, ]) > 0L)

# BCF: the grid is shared by both forests, and revalidation must leave both
cutBcf <- bcfSampler()
invisible(cutBcf$run(20L, 1L))
cutBcf$setCutPoints(list(0.4), 1L)
bcfCutRun <- cutBcf$run(0L, 20L)
expect_identical(drawnCuts(cutBcf, 1L), 0.4)
expect_true(sum(bcfCutRun$varcount[1L, , ]) > 0L)
expect_true(all(is.finite(bcfCutRun$train)))
expect_true(all(is.finite(cutBcf$getForestFits(2L))))

# a leaf-covariate column carries a split grid AND a leaf regressor; coarsening
# it must not disturb the leaf model
cutLinear <- dbarts(X, y, node.prior = linearLeaf, control = ctl())
invisible(cutLinear$run(20L, 1L))
cutLinear$setCutPoints(list(0.4), 3L)
linCutRun <- cutLinear$run(0L, 20L)
expect_identical(drawnCuts(cutLinear, 3L), 0.4)
expect_true(sum(linCutRun$varcount[3L, ]) > 0L)
expect_true(all(is.finite(linCutRun$train)))

# NA-carrying column: the imputed rows ride the new grid like any other
cutMissing <- dbarts(Xna, y, control = ctl())
invisible(cutMissing$run(20L, 1L))
cutMissing$setCutPoints(list(0.4), 2L)
naCutRun <- cutMissing$run(0L, 20L)
expect_identical(drawnCuts(cutMissing, 2L), 0.4)
expect_true(sum(naCutRun$varcount[2L, ]) > 0L)
expect_true(all(is.finite(naCutRun$train)))


## --- 4. multi-chain posterior x {BCF, DART, linear, GP} --------------------
## The multi-chain tests for these three assert shapes (or nothing). The
## engine's contract is that a single-chain run at seed S reproduces chain 1
## of any multi-chain run at seed S, which turns them into draw-level tests.

chainArms <- list(
  bcf = function(k) {
    dbarts(
      X,
      y,
      forests = list(forest(), forest(basis = ~ factor(z))),
      control = ctl(n.chains = k)
    )
  },
  dart = function(k) {
    dbarts(X, y, tree.prior = dart(), control = ctl(n.chains = k))
  },
  linear = function(k) {
    dbarts(X, y, node.prior = linearLeaf, control = ctl(n.chains = k))
  },
  gp = function(k) {
    dbarts(X, y, node.prior = gp("v3", k = 2), control = ctl(n.chains = k))
  }
)
for (arm in names(chainArms)) {
  many <- chainArms[[arm]](2L)$run(10L, 4L)
  one <- chainArms[[arm]](1L)$run(10L, 4L)
  expect_identical(
    as.numeric(many$train[,, 1L]),
    as.numeric(one$train),
    info = arm
  )
  expect_identical(
    as.numeric(many$sigma[, 1L]),
    as.numeric(one$sigma),
    info = arm
  )
  expect_false(identical(many$train[,, 1L], many$train[,, 2L]), info = arm)
}


## --- 5. BCF x mutation axes ------------------------------------------------
## setData is refused at both levels and setWeights has its bitwise
## build-vs-swap oracle in test-multi-forest-seam.R; the gaps are setResponse
## and setOffset at draw level, and $setWeights on the public R5 surface.

# setOffset mid-chain: the reported location is exactly the forest total plus
# the installed offset, on a coupling whose reported scale is glued
offsetBcf <- bcfSampler()
invisible(offsetBcf$run(20L, 1L))
offsets <- runif(n, -0.3, 0.3)
offsetBcf$setOffset(offsets)
offsetRun <- offsetBcf$run(0L, 1L)
expect_identical(
  as.numeric(offsetRun$train[, 1L]),
  as.numeric(offsetBcf$getFitsWithoutOffset()) + offsets
)
expect_true(sd(offsets) > 0.1)

# setResponse at the pinned scale: re-installing the same response is exactly
# a no-op, and a response with the treatment contrast removed shrinks the
# treatment forest rather than merely moving the total
responseBcf <- bcfSampler()
responseControl <- bcfSampler()
invisible(responseBcf$run(20L, 1L))
invisible(responseControl$run(20L, 1L))
responseBcf$setResponse(y)
expect_identical(
  responseBcf$run(0L, 3L)$train,
  responseControl$run(0L, 3L)$train
)

flatBcf <- bcfSampler()
keptBcf <- bcfSampler()
invisible(flatBcf$run(20L, 1L))
invisible(keptBcf$run(20L, 1L))
flatBcf$setResponse(y - z)
invisible(flatBcf$run(200L, 1L))
invisible(keptBcf$run(200L, 1L))
expect_true(
  mean(abs(flatBcf$getForestFits(2L))) < mean(abs(keptBcf$getForestFits(2L)))
)
expect_true(all(is.finite(flatBcf$getForestFits(1L))))

# $setWeights on the public two-forest surface: accepted, and it moves the
# posterior rather than being swallowed by the coupling
weightBcf <- bcfSampler()
weightControl <- bcfSampler()
invisible(weightBcf$run(20L, 1L))
invisible(weightControl$run(20L, 1L))
weightBcf$setWeights(rep(1, n))
weightBase <- weightControl$run(0L, 3L)$train
expect_identical(weightBcf$run(0L, 3L)$train, weightBase)
weightBcf$setWeights(wts)
weightRun <- weightBcf$run(0L, 3L)
expect_false(isTRUE(all.equal(
  weightRun$train,
  weightControl$run(0L, 3L)$train
)))
expect_true(all(is.finite(weightRun$train)))
expect_true(all(is.finite(weightBcf$getForestFits(2L))))


## --- 6. state round-trip DEPTH: posterior CONTINUATION --------------------
## The suite's round trips end at the join (state-object equality, or predict
## equality across serialization). This runs PAST it: N sweeps, storeState, a
## FRESH sampler, setState, M more sweeps, against an uninterrupted N+M run at
## the same seed. Restore is semantic by decision, so the claim is exact on
## the channel the rebuild reproduces exactly - the sequence of tree moves,
## read off varcount - and to within the resum on the numbers. A bitwise
## numeric assertion here would be a false gate.

nBefore <- 3L
nAfter <- 5L
continuationArm <- function(
  make,
  runIt = function(s, k) s$run(0L, k),
  storeIt = function(s) {
    s$storeState()
    s$state
  },
  restoreIt = function(s, st) s$setState(st)
) {
  whole <- runIt(make(), nBefore + nAfter)
  split <- make()
  invisible(runIt(split, nBefore))
  fresh <- make()
  restoreIt(fresh, storeIt(split))
  cont <- runIt(fresh, nAfter)
  wholeTrain <- sweepTail(whole$train, nAfter)
  list(
    moves = identical(sweepTail(whole$varcount, nAfter), cont$varcount),
    numeric = max(abs(wholeTrain - cont$train)),
    scale = max(abs(wholeTrain)),
    mixes = max(abs(
      sweepSlice(cont$train, 1L) - sweepSlice(cont$train, nAfter)
    ))
  )
}

mnLabels <- rbinom(n, 2L, 0.4)
continuationArms <- list(
  gaussian = function() {
    continuationArm(function() dbarts(X, y, control = ctl()))
  },
  weighted = function() {
    continuationArm(function() dbarts(X, y, weights = wts, control = ctl()))
  },
  grouped = function() {
    continuationArm(function() dbarts(X, y, control = groupedCtl()))
  },
  bcf = function() continuationArm(bcfSampler),
  multinomial = function() {
    continuationArm(
      function() {
        dbarts:::bartcoreMultinomialSampler(
          dbarts(X, as.double(mnLabels), control = ctl()),
          mnLabels,
          K = 3L
        )
      },
      runIt = function(s, k) dbarts:::bartcoreRun(s, 0L, k),
      storeIt = dbarts:::bartcoreStoreState,
      restoreIt = dbarts:::bartcoreSetState
    )
  }
)
for (arm in names(continuationArms)) {
  got <- continuationArms[[arm]]()
  expect_true(got$moves, info = arm)
  expect_true(got$numeric < 1e-9 * max(1, got$scale), info = arm)
  # non-vacuity: the continuation is a live chain, not a frozen one
  expect_true(got$mixes > 0, info = arm)
}
