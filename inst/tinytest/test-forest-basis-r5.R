# M4.3 of docs/plans/multiforest-extension-surface.md: the K-length forest
# spec surface. setTreatment has retired as a mutator - basis synthesis is
# construction-only and $setForestBasis is the SOLE mutation route - so what is
# pinned here is what that leaves: both orderings of a widen and a swap at the
# surface (arm 6), the persistence contract across a re-creation with the
# amplitude VALUES asserted rather than the surviving width (arm 7), the ragged
# state block's per-forest widths, the ragged run-result glue slot at K > 2,
# and the capability probe that succeeds data@treatment.
#
# Not re-pinned here: the ordering CONTRACT itself (tests/cpp
# testForestBasisOrdering owns the engine half, including the draw-path value
# predicate and the z divergence), the creation surface (test-bcf-creation.R),
# and the flat C route (inst/tinytest/capi/consumer.c).

set.seed(41)
n <- 180L
p <- 4L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

seededControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    n.samples = 6L,
    updateState = FALSE,
    rngSeed = 41L,
    ...
  )
}

twoForests <- function(...) {
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = seededControl(...)
  )
}

wide <- cbind(1, x[, 4L])

# --- arm 6, both orderings at the surface. Installing on one forest never
# disturbs the other, and the two orderings leave the same per-forest fits and
# the same amplitudes. RED against a route that re-synthesizes on install. ---
widenThenSwap <- twoForests()
widenThenSwap$run(5L, 2L)
before <- widenThenSwap$getForestAmplitudes()
widenThenSwap$setForestBasis(1L, wide)
widened <- widenThenSwap$getForestAmplitudes(1L)
# the layout moved, so the added coordinate entered at the neutral 1
expect_equal(nrow(widened), 2L)
expect_equal(widened[1L, 1L], before[1L, 1L])
expect_equal(widened[2L, 1L], 1)
widenThenSwap$setForestBasis(2L, factor(1L - z))
expect_equal(ncol(widenThenSwap$data@bases[[1L]]), 2L)
# a width-preserving swap is the bitwise identity on every amplitude
expect_identical(widenThenSwap$getForestAmplitudes(1L), widened)

swapThenWiden <- twoForests()
swapThenWiden$run(5L, 2L)
swapThenWiden$setForestBasis(2L, factor(1L - z))
swapThenWiden$setForestBasis(1L, wide)
expect_identical(
  swapThenWiden$getForestAmplitudes(),
  widenThenSwap$getForestAmplitudes()
)
expect_identical(
  swapThenWiden$run(0L, 3L)$forestFits,
  widenThenSwap$run(0L, 3L)$forestFits
)

# --- arm 7, persistence. The bases ride CREATION (data@bases), so a widening
# applied after a state restore preserves and remaps the RESTORED amplitudes
# rather than the constructed ones - the restore-then-widen half of the
# contract, asserted on the amplitude VALUES and not only on the surviving
# width. ---
persisted <- twoForests()
persisted$run(5L, 2L)
persisted$setForestBasis(1L, wide)
persisted$setForestBasis(2L, factor(1L - z))
restoredAmplitudes <- persisted$getForestAmplitudes()
persisted$storeState()

persistedFile <- tempfile(fileext = ".rds")
saveRDS(persisted, persistedFile)
reloaded <- readRDS(persistedFile)
unlink(persistedFile)
# forces getPointer()'s transparent re-creation
reloaded$getPointer()
expect_equal(ncol(reloaded$data@bases[[1L]]), 2L)
expect_identical(reloaded$getForestAmplitudes(), restoredAmplitudes)
# and a widening applied AFTER that restore keeps the restored values, which is
# the half a reapply running after setState would lose
reloaded$setForestBasis(2L, cbind(1 - z, z, x[, 1L]))
widerAmplitudes <- reloaded$getForestAmplitudes()
expect_equal(nrow(widerAmplitudes), 5L)
expect_identical(widerAmplitudes[1:4, , drop = FALSE], restoredAmplitudes)
expect_equal(widerAmplitudes[5L, 1L], 1)

# $copy() carries the mirror too, on the same terms setForestWeights' does
copied <- persisted$copy()
expect_equal(ncol(copied$data@bases[[1L]]), 2L)
expect_identical(copied$getForestAmplitudes(), restoredAmplitudes)

# --- the state block carries PER-FOREST WIDTHS, not just a total. q = (1, 3)
# and q = (2, 2) both serialize four amplitudes, so a reader keyed on the total
# would accept either into either and silently permute the blocks. ---
ragged <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~g)),
  control = seededControl()
)
ragged$run(5L, 2L)
ragged$storeState()
raggedState <- ragged$state

square <- twoForests()
square$run(5L, 2L)
square$setForestBasis(1L, wide)
square$storeState()
squareState <- square$state
# same total (4), different layout: each refuses the other
expect_error(square$setState(raggedState), "not consistent")
expect_error(ragged$setState(squareState), "not consistent")
# its own state restores, and the amplitudes come back bitwise
raggedAmplitudes <- ragged$getForestAmplitudes()
ragged$run(0L, 1L)
ragged$setState(raggedState)
expect_identical(ragged$getForestAmplitudes(), raggedAmplitudes)

# --- the run-result glue slot is the RAGGED amplitude vector, sum q_f rows,
# forest-major within a draw; bcf's q = (1, 2) leaves it the 3 rows it shipped
# with, and the reconstruction identity holds at K = 3. ---
threeForests <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z)), forest(basis = ~g)),
  control = seededControl()
)
threeResult <- threeForests$run(5L, 4L)
expect_equal(dim(threeResult$glue), c(6L, 4L))
expect_equal(dim(threeResult$forestFits), c(n, 3L, 4L))
# the last draw's reported train channel is the combination of the three
# forests through their own amplitude blocks
# the reported train channel is on the RESPONSE scale; map it back through the
# stored (min, max) of y, as test-bcf-reporting.R's own helper does
threeForests$storeState()
fitScale <- threeForests$state[[1L]]$fit.scale
scale <- fitScale[2L] - fitScale[1L]
internalTrain <- (threeResult$train - (scale * 0.5 + fitScale[1L])) / scale
bases <- threeForests$data@bases
reconstructionError <- 0
for (sampleNum in seq_len(ncol(threeResult$glue))) {
  amplitudes <- threeResult$glue[, sampleNum]
  fits <- threeResult$forestFits[,, sampleNum]
  recon <- amplitudes[1L] *
    fits[, 1L] +
    as.vector(bases[[2L]] %*% amplitudes[2:3]) * fits[, 2L] +
    as.vector(bases[[3L]] %*% amplitudes[4:6]) * fits[, 3L]
  reconstructionError <-
    max(reconstructionError, max(abs(recon - internalTrain[, sampleNum])))
}
expect_true(reconstructionError < 1e-12)

# --- $getForestAmplitudes(forest) is per forest, and the whole vector at the
# default; the ragged widths are its own bases' ---
expect_equal(nrow(threeForests$getForestAmplitudes(1L)), 1L)
expect_equal(nrow(threeForests$getForestAmplitudes(2L)), 2L)
expect_equal(nrow(threeForests$getForestAmplitudes(3L)), 3L)
expect_identical(
  threeForests$getForestAmplitudes(),
  rbind(
    threeForests$getForestAmplitudes(1L),
    threeForests$getForestAmplitudes(2L),
    threeForests$getForestAmplitudes(3L)
  )
)
expect_error(threeForests$getForestAmplitudes(4L), "out of range")
expect_error(threeForests$getForestAmplitudes(0L), "single positive integer")

# --- the per-forest ASIS ridge is DERIVED from the amplitude prior's kind at
# the transport's own site (R_interface_bartcore.cpp applyBCFSpec:
# forest.ridge = forest.amplitudePriorScale > 0.0), which is what reproduces
# bcf's a-move on and b-move off. Two halves. First the transported scale
# itself: a forest carrying a basis has a FIXED-variance amplitude, so its
# half-Cauchy scale is 0 and its ridge is off; a forest carrying none has the
# scale mixture, so its scale is its `sd` and its ridge is on. ---
ridgeSpec <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, cbind(1 - z, z))),
  seededControl(),
  forests = list(forest(sd = 2.5), forest(sd = 1.25))
)
ridgeParams <- attr(ridgeSpec$control, "bartcore.bcf")$params
# forest 1: no basis -> variance 1, half-Cauchy median 2.5 -> ridge ON
expect_equal(ridgeParams[[1L]][6:7], c(1, 2.5))
# forest 2: a basis -> fixed variance 0.5, scale 0 -> ridge OFF
expect_equal(ridgeParams[[2L]][6:7], c(0.5, 0))

# and second, the derivation is really READ: flipping the transported scale on
# either forest flips that forest's ridge, and a ridge that travels consumes a
# GIG draw per sweep, so every subsequent draw moves. A bridge that set the
# flag by any other rule would leave these three fits identical.
buildFromSpec <- function(spec) {
  new("dbartsSampler", spec$control, spec$model, spec$data)
}
baseline <- buildFromSpec(ridgeSpec)$run(0L, 4L)$train

ridgeOnBasis <- ridgeSpec
attr(ridgeOnBasis$control, "bartcore.bcf")$params[[2L]][7L] <- 2
expect_false(identical(buildFromSpec(ridgeOnBasis)$run(0L, 4L)$train, baseline))

ridgeOffFirst <- ridgeSpec
attr(ridgeOffFirst$control, "bartcore.bcf")$params[[1L]][7L] <- 0
expect_false(identical(
  buildFromSpec(ridgeOffFirst)$run(0L, 4L)$train,
  baseline
))

# --- the single-forest queries (numTrees, and everything the control's tree
# count is cross-checked against) address FOREST 0 whatever the K-length spec
# says, which is what the sampler's own options carry. Pinned on a spec whose
# forests carry DIFFERENT counts, so a read that answered for any other forest
# would name a different number. ---
unevenTrees <- dbarts(
  x,
  y,
  forests = list(
    forest(n.trees = 13L),
    forest(basis = ~ factor(z), n.trees = 31L),
    forest(basis = ~g, n.trees = 7L)
  ),
  control = seededControl()
)
expect_identical(unevenTrees$control@n.trees, 13L)
# setControl cross-checks the control's count against the sampler's shape, so
# it accepts forest 0's and refuses either of the others
treeCountControl <- function(numTrees) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = numTrees,
    n.samples = 6L,
    updateState = FALSE,
    rngSeed = 41L
  )
}
expect_silent(unevenTrees$setControl(treeCountControl(13L)))
expect_error(
  unevenTrees$setControl(treeCountControl(31L)),
  "cannot change .n.trees. on an existing sampler"
)
expect_error(
  unevenTrees$setControl(treeCountControl(7L)),
  "cannot change .n.trees. on an existing sampler"
)

# --- dbartsData(bases = ) takes a CONTINUOUS basis, the positive route the
# "not a numeric basis column" refusal relaxed into: creation succeeds, the
# column lands in the engine, and its single amplitude really multiplies it -
# asserted through the fitted-value reconstruction rather than by dimensions
# alone, so a basis accepted and dropped would still fail. ---
continuous <- x[, 1L] - 0.5
continuousData <- dbartsData(x, y, bases = list(NULL, cbind(continuous)))
continuousFit <- dbarts(continuousData, control = seededControl())
expect_equal(ncol(continuousFit$data@bases[[2L]]), 1L)
expect_equal(continuousFit$data@bases[[2L]][, 1L], continuous)
continuousResult <- continuousFit$run(5L, 3L)
expect_equal(dim(continuousResult$glue), c(2L, 3L))
expect_equal(nrow(continuousFit$getForestAmplitudes(2L)), 1L)
continuousFit$storeState()
continuousScale <- continuousFit$state[[1L]]$fit.scale
continuousSpan <- continuousScale[2L] - continuousScale[1L]
continuousInternal <- (continuousResult$train -
  (continuousSpan * 0.5 + continuousScale[1L])) /
  continuousSpan
continuousError <- 0
for (sampleNum in seq_len(ncol(continuousResult$glue))) {
  amplitudes <- continuousResult$glue[, sampleNum]
  fits <- continuousResult$forestFits[,, sampleNum]
  recon <- amplitudes[1L] *
    fits[, 1L] +
    amplitudes[2L] * continuous * fits[, 2L]
  continuousError <-
    max(continuousError, max(abs(recon - continuousInternal[, sampleNum])))
}
expect_true(continuousError < 1e-12)

# --- the capability probe that succeeds data@treatment: it reads "carries
# amplitudes", and is deliberately NOT a forest count, which a K-forest
# multinomial (which carries no amplitudes at all) would defeat ---
expect_true(dbarts:::isBCFSampler(threeForests))
plain <- dbarts(x, y, control = seededControl())
expect_false(dbarts:::isBCFSampler(plain))
# the five refusals it gates stay red on a multi-forest sampler
expect_error(threeForests$setResponse(y, updateScale = TRUE), "does not supp")
expect_error(
  threeForests$setOffset(rep(0, n), updateScale = TRUE),
  "does not supp"
)
expect_error(threeForests$setData(dbartsData(x, y)), "does not support")
expect_error(threeForests$setModel(threeForests$model), "does not support")
expect_error(
  threeForests$setCalibration(prior.scale = 1, forest = 1L),
  "does not support"
)
# and the multinomial route is NOT misidentified by it
mnLabels <- sample(0:2, n, replace = TRUE)
mnSampler <- dbarts(x, y, control = seededControl())
mn <- dbarts:::bartcoreMultinomialSampler(mnSampler, mnLabels, 3L)
expect_error(
  dbarts:::bartcoreForestAmplitudes(mn, 0L),
  "forests carry them"
)
expect_error(
  dbarts:::bartcoreSetForestBasis(mn, 1L, cbind(0.5, 0.5)),
  "forests carry amplitudes"
)

# --- M4.4: the leaf scale a swap leaves behind. The calibration map divides
# out the median nonzero row norm of each forest's basis, and that divisor is
# otherwise a construction-time constant - no mutation re-derives a K-forest
# leaf scale - so $setForestBasis owns the staleness. Pinned on a PROBIT
# K-forest, where the anchor is the literal 1 and $getCalibration's
# prior.scale IS the map's node scale, so the assertion is exact. Nothing
# else in the suite or in the equivalence trio calls this mutator at all. ---
yBinary <- as.double(y > median(y))
probitForests <- function() {
  dbarts(
    x,
    yBinary,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = seededControl()
  )
}
priorScale <- function(sampler, forest) {
  unname(sampler$getCalibration(forest)[1L, "prior.scale"])
}

# (i) STALENESS: a basis whose median nonzero row norm is 4x the old one moves
# the divisor with it. RED against a divisor left at its construction value.
stale <- probitForests()
expect_equal(priorScale(stale, 2L), 1 / 0.674, tolerance = 1e-12)
stale$run(5L, 1L)
stale$setForestBasis(2L, 4 * cbind(1 - z, z))
expect_equal(priorScale(stale, 2L), 1 / (0.674 * 4), tolerance = 1e-12)
# and the untouched forest keeps its own
expect_equal(priorScale(stale, 1L), 1, tolerance = 1e-12)

# (ii) BITWISE: a re-install of the SAME values is the identity, on the
# calibration and on every draw after it. That is the claim retaining the
# construction-time anchor rests on - recomputing it would recalibrate this
# forest onto a scale the others are not on - and only this pins it.
swapped <- probitForests()
twin <- probitForests()
swapped$run(5L, 1L)
twin$run(5L, 1L)
swapped$setForestBasis(2L, cbind(1 - z, z))
expect_identical(priorScale(swapped, 2L), priorScale(twin, 2L))
swappedResult <- swapped$run(5L, 2L)
twinResult <- twin$run(5L, 2L)
expect_identical(swapped$getForestFits(2L), twin$getForestFits(2L))
expect_identical(swapped$getForestAmplitudes(), twin$getForestAmplitudes())
expect_identical(swappedResult$sigma, twinResult$sigma)
expect_identical(swappedResult$train, twinResult$train)

# (iii) THE WHOLE PRODUCT, re-derived across the swap. (i) and (ii) both run at
# the DEFAULT node scale factor - a literal 1 - so a swap that lost the stored
# factor would pass them both; only a DECLARED forest(sd = ) can see it. Every
# factor of sd_f * s / (0.674 c_f) discriminates here: sd in {2.5, 0.4}, c in
# {3, 5} then 7, divisor 0.674, and s is probit's 1 (the anchor's own arm is
# test-bcf-family.R's, which runs the same fixture under logistic too, where a
# dropped anchor is a 81 percent miss).
ones <- matrix(1, n, 1L)
zBasis <- cbind(1 - z, z)
declared <- dbarts(
  dbartsData(x, yBinary, bases = list(3 * ones, 5 * zBasis)),
  forests = list(forest(sd = 2.5), forest(sd = 0.4)),
  control = seededControl()
)
expect_equal(priorScale(declared, 1L), 2.5 / (0.674 * 3), tolerance = 1e-12)
expect_equal(priorScale(declared, 2L), 0.4 / (0.674 * 5), tolerance = 1e-12)
declared$run(5L, 1L)
declared$setForestBasis(2L, 7 * zBasis)
# the new row norm enters and the DECLARED factor survives: 0.084781687156
expect_equal(priorScale(declared, 2L), 0.4 / (0.674 * 7), tolerance = 1e-12)
expect_equal(priorScale(declared, 1L), 2.5 / (0.674 * 3), tolerance = 1e-12)

# --- the row norm IN FORCE is reported rather than only divided out, and it
# follows the swap. A stored-and-stale value is what this catches: the
# construction-time norm was 5 and the swap makes it 7. The second install
# carries the convention as well - zero rows excluded, an EVEN nonzero count
# averaging the two central order statistics - which R's own median reproduces
# from data@bases, so the assertion is against the R-side rule rather than
# against a number copied out of the engine. ---
mapColumn <- function(sampler, forest, column) {
  unname(sampler$getCalibration(forest)[1L, column])
}
expect_equal(mapColumn(declared, 2L, "basis.row.norm"), 7)
expect_equal(mapColumn(declared, 1L, "basis.row.norm"), 3)
evenBasis <- matrix(c(rep(0, 4L), seq_len(n - 4L)), n, 1L)
evenNorms <- sqrt(rowSums(evenBasis^2))
declared$setForestBasis(2L, evenBasis)
expect_equal(
  mapColumn(declared, 2L, "basis.row.norm"),
  median(evenNorms[evenNorms > 0])
)
expect_equal(
  priorScale(declared, 2L),
  0.4 / (0.674 * median(evenNorms[evenNorms > 0])),
  tolerance = 1e-12
)

# --- the SPEC ECHO's truthfulness across a state install. $setState and
# $installTrees both adopt the donor's leaf scale by design, so a reader that
# echoed the recipient's own spec would print a decomposition of a number the
# recipient no longer runs under. Three rules, four arms. ---
donorForests <- function(sd, variance) {
  dbarts(
    x,
    yBinary,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), sd = sd, amplitude.prior.variance = variance)
    ),
    control = seededControl()
  )
}
donor <- donorForests(2, 0.125)
recipient <- donorForests(0.5, 0.5)
# the two really do differ in the calibration, which is what makes the install
# foreign rather than a re-statement
expect_true(priorScale(donor, 2L) != priorScale(recipient, 2L))
donor$run(5L, 1L)
donor$storeState()
donorState <- donor$state
donorScale <- priorScale(donor, 2L)

# (a) FOREIGN CALIBRATION. The install is accepted - the widths and tree
# counts agree, and neither the state gate nor the forest gate looks at a leaf
# scale - so the recipient runs under the donor's scale. Its stored factor and
# divisor no longer decompose it, and the reader says so with NaN rather than
# printing a decomposition that would recover a wrong anchor.
recipient$setState(donorState)
expect_equal(priorScale(recipient, 2L), donorScale)
expect_true(is.nan(mapColumn(recipient, 2L, "node.scale.factor")))
expect_true(is.nan(mapColumn(recipient, 2L, "node.scale.divisor")))
# the anchor is therefore NOT computable, which is the point of the NaN
expect_true(is.nan(
  mapColumn(recipient, 2L, "prior.scale") *
    mapColumn(recipient, 2L, "node.scale.divisor") *
    mapColumn(recipient, 2L, "basis.row.norm") /
    mapColumn(recipient, 2L, "node.scale.factor")
))
# the row norm needs no rule - bases are not state, so the recipient's own is
# still the one in force
expect_equal(mapColumn(recipient, 2L, "basis.row.norm"), 1)
# and forest 1's columns SURVIVE: both samplers calibrate it identically, so
# the installed scale is bitwise the one in force. The rule compares before it
# assigns, and a rule that cleared on every install would lose this.
expect_equal(mapColumn(recipient, 1L, "node.scale.factor"), 1)
expect_equal(mapColumn(recipient, 1L, "node.scale.divisor"), 1)

# (d) THE AMPLITUDE PRIOR FOLLOWS THE STATE, which is what the next draw will
# use: the recipient reports the DONOR's variance, not its own 0.5. Forest 1
# carries the scale mixture, whose serialized variance is a live auxiliary
# rather than a prior, so its two amplitude columns keep their exclusivity.
expect_equal(mapColumn(recipient, 2L, "amplitude.prior.variance"), 0.125)
expect_true(is.nan(mapColumn(recipient, 2L, "amplitude.prior.scale")))
expect_true(is.nan(mapColumn(recipient, 1L, "amplitude.prior.variance")))
# neither sampler declares that forest's sd, so both carry the family's own
# default median - 1 under this probit fixture, where the anchor is the link's
# error sd and sigma is pinned, rather than gaussian's 2
expect_equal(mapColumn(recipient, 1L, "amplitude.prior.scale"), 1)

# (c) RE-IMPOSITION. $setForestBasis re-derives the leaf scale from the stored
# factor and divisor, so both columns come back and the identity holds again -
# on the RECIPIENT's own calibration, since the map is what was re-imposed.
recipient$setForestBasis(2L, 4 * zBasis)
expect_equal(mapColumn(recipient, 2L, "node.scale.factor"), 0.5)
expect_equal(mapColumn(recipient, 2L, "node.scale.divisor"), 0.674)
expect_equal(mapColumn(recipient, 2L, "basis.row.norm"), 4)
expect_equal(
  mapColumn(recipient, 2L, "prior.scale") *
    mapColumn(recipient, 2L, "node.scale.divisor") *
    mapColumn(recipient, 2L, "basis.row.norm") /
    mapColumn(recipient, 2L, "node.scale.factor"),
  1,
  tolerance = 1e-12
)

# (b) SELF-RESTORE. A store, a run and a restore of a sampler's OWN state
# installs a bitwise-identical scale, so every column survives non-NaN and the
# identity still holds. This is the arm that keeps the rule from being "clear
# on any install".
selfRestore <- donorForests(2, 0.125)
selfRestore$run(5L, 1L)
selfRestore$storeState()
ownState <- selfRestore$state
ownScale <- priorScale(selfRestore, 2L)
selfRestore$run(0L, 1L)
selfRestore$setState(ownState)
expect_equal(priorScale(selfRestore, 2L), ownScale)
expect_equal(mapColumn(selfRestore, 2L, "node.scale.factor"), 2)
expect_equal(mapColumn(selfRestore, 2L, "node.scale.divisor"), 0.674)
expect_equal(mapColumn(selfRestore, 2L, "amplitude.prior.variance"), 0.125)
expect_equal(
  mapColumn(selfRestore, 2L, "prior.scale") *
    mapColumn(selfRestore, 2L, "node.scale.divisor") *
    mapColumn(selfRestore, 2L, "basis.row.norm") /
    mapColumn(selfRestore, 2L, "node.scale.factor"),
  1,
  tolerance = 1e-12
)
