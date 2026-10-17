# The K-length forest spec surface. setTreatment has retired as a mutator -
# basis synthesis is construction-only and $setForestBasis is the SOLE
# mutation route - so what is pinned here is what that leaves: both
# orderings of a swap on each of two forests at the surface (arm 6), the
# persistence contract across a re-creation with the amplitude VALUES asserted
# (arm 7), the ragged
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

seededControlForestBasisR5 <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    n.samples = 6L,
    updateState = FALSE,
    seed = 41L,
    ...
  )
}

twoForests <- function(...) {
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = seededControlForestBasisR5(...)
  )
}

# A forest created with no basis takes none, so the forest swapped below is
# created on one numeric column. No forest of that model is without a basis,
# so the count is stated on a forest and the control names none.
numericFirst <- function() {
  dbarts(
    x,
    y,
    forests = list(
      forest(basis = x[, 4L], n.trees = 25L),
      forest(basis = ~ factor(z))
    ),
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 6L,
      updateState = FALSE,
      seed = 41L
    )
  )
}
other <- x[, 1L]

# --- arm 6, both orderings at the surface. Installing on one forest never
# disturbs the other, and the two orderings leave the same per-forest fits and
# the same amplitudes. RED against a route that re-synthesizes on install. A
# forest keeps its width for life, so every install preserves it, and an
# install that preserves the width is the bitwise identity on every
# amplitude. ---
firstThenSecond <- numericFirst()
firstThenSecond$run(5L, 2L)
before <- firstThenSecond$getForestAmplitudes()
firstThenSecond$setForestBasis(1L, other)
expect_identical(firstThenSecond$getForestAmplitudes(), before)
firstThenSecond$setForestBasis(2L, factor(1L - z))
expect_identical(firstThenSecond$getForestAmplitudes(), before)
# a width the forest was not created with is refused, nothing being installed
expect_error(
  firstThenSecond$setForestBasis(1L, cbind(1, x[, 4L])),
  "forest 1 ('x[, 4L]') multiplies one numeric column, and 'basis' gives it 2",
  fixed = TRUE
)
expect_identical(unname(firstThenSecond$data@bases[[1L]][, 1L]), other)

secondThenFirst <- numericFirst()
secondThenFirst$run(5L, 2L)
secondThenFirst$setForestBasis(2L, factor(1L - z))
secondThenFirst$setForestBasis(1L, other)
expect_identical(
  secondThenFirst$getForestAmplitudes(),
  firstThenSecond$getForestAmplitudes()
)
expect_identical(
  secondThenFirst$run(0L, 3L)$forestFits,
  firstThenSecond$run(0L, 3L)$forestFits
)

# --- arm 7, persistence. The bases ride CREATION (data@bases), so a swap
# applied after a state restore preserves the RESTORED amplitudes rather than
# the constructed ones - the restore-then-swap half of the contract, asserted
# on the amplitude VALUES. ---
persisted <- numericFirst()
persisted$run(5L, 2L)
persisted$setForestBasis(1L, other)
persisted$setForestBasis(2L, factor(1L - z))
restoredAmplitudes <- persisted$getForestAmplitudes()
persisted$storeState()

persistedFile <- tempfile(fileext = ".rds")
saveRDS(persisted, persistedFile)
reloaded <- readRDS(persistedFile)
unlink(persistedFile)
# forces getPointer()'s transparent re-creation
reloaded$getPointer()
expect_identical(unname(reloaded$data@bases[[1L]][, 1L]), other)
expect_identical(reloaded$getForestAmplitudes(), restoredAmplitudes)
# and a swap applied AFTER that restore keeps the restored values, which is
# the half a reapply running after setState would lose
reloaded$setForestBasis(1L, x[, 2L])
expect_identical(reloaded$getForestAmplitudes(), restoredAmplitudes)

# $copy() carries the mirror too, on the same terms setForestWeights' does
copied <- persisted$copy()
expect_identical(unname(copied$data@bases[[1L]][, 1L]), other)
expect_identical(copied$getForestAmplitudes(), restoredAmplitudes)

# --- the state block carries PER-FOREST WIDTHS, not just a total. q = (1, 2)
# and q = (2, 1) both serialize three amplitudes, so a reader keyed on the
# total would accept either into either and silently permute the blocks. ---
ragged <- twoForests()
ragged$run(5L, 2L)
ragged$storeState()
raggedState <- ragged$state

reversed <- dbarts(
  x,
  y,
  forests = list(
    forest(basis = ~ factor(z), n.trees = 25L),
    forest(basis = x[, 4L])
  ),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 6L,
    updateState = FALSE,
    seed = 41L
  )
)
reversed$run(5L, 2L)
reversed$storeState()
reversedState <- reversed$state
# same total (3), different layout: each refuses the other
expect_error(reversed$setState(raggedState), "not consistent")
expect_error(ragged$setState(reversedState), "not consistent")
# its own state restores, and the amplitudes come back bitwise
raggedAmplitudes <- ragged$getForestAmplitudes()
ragged$run(0L, 1L)
ragged$setState(raggedState)
expect_identical(ragged$getForestAmplitudes(), raggedAmplitudes)

# --- the run-result glue slot is the RAGGED amplitude vector, sum q_f rows,
# forest-major within a draw; bcf's q = (1, 2) leaves it the 3 rows it shipped
# with, and the reconstruction identity holds at K = 3, q = (1, 2, 2). ---
threeForests <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z)),
    forest(basis = ~ g == "b")
  ),
  control = seededControlForestBasisR5()
)
threeResult <- threeForests$run(5L, 4L)
expect_equal(dim(threeResult$glue), c(5L, 4L))
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
    as.vector(bases[[3L]] %*% amplitudes[4:5]) * fits[, 3L]
  reconstructionError <-
    max(reconstructionError, max(abs(recon - internalTrain[, sampleNum])))
}
expect_true(reconstructionError < 1e-12)

# --- $getForestAmplitudes(forest) is per forest, and the whole vector at the
# default; the ragged widths are its own bases' ---
expect_equal(nrow(threeForests$getForestAmplitudes(1L)), 1L)
expect_equal(nrow(threeForests$getForestAmplitudes(2L)), 2L)
expect_equal(nrow(threeForests$getForestAmplitudes(3L)), 2L)
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

# --- the amplitude prior's kind follows what each forest multiplies. Two
# halves. First the transported statement itself: each forest's `sd` goes to
# the engine as written, in the response's units, and what the forest
# multiplies is read from the data object; a forest carrying a basis has a
# FIXED-variance amplitude, so it reports no half-Cauchy scale; a forest
# carrying none has the scale mixture, so its scale is its `sd` (over the
# response's standard deviation here, with no offset, subset or zero
# weight). ---
priorSpec <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, factor(z))),
  seededControlForestBasisR5(),
  forests = list(forest(sd = 2.5), forest(sd = 1.25))
)
priorStated <- attr(priorSpec$control, "bartcore.forests")
expect_identical(priorStated$sd, c(2.5, 1.25))
expect_null(priorStated$variance)
expect_null(priorStated$hasBasis)
priorSampler <- new(
  "dbartsSampler",
  priorSpec$control,
  priorSpec$model,
  priorSpec$data
)
# forest 1: no basis -> no variance, half-Cauchy median 2.5 in the response's
# units, which is 2.5 / sd(y) of the response's standard deviation
priorFree <- priorSampler$getLeafPrior(1L)
expect_identical(priorFree$sd, 2.5)
expect_identical(priorFree$multiplier, "none")
expect_equal(priorFree$amplitude.prior.scale, 2.5 / sd(y), tolerance = 1e-12)
# forest 2: a basis -> fixed variance, no half-Cauchy scale
priorBasis <- priorSampler$getLeafPrior(2L)
expect_identical(priorBasis$sd, 1.25)
expect_identical(priorBasis$multiplier, "factor")
expect_null(priorBasis$amplitude.prior.scale)

# and second, the statement is really READ: the half-Cauchy median moves every
# subsequent draw. A bridge that ignored the transported statement would leave
# these fits identical.
buildFromSpec <- function(spec) {
  new("dbartsSampler", spec$control, spec$model, spec$data)
}
baseline <- buildFromSpec(priorSpec)$run(0L, 4L)$train
movedMedian <- priorSpec
attr(movedMedian$control, "bartcore.forests")$sd[1L] <- 5
expect_false(identical(
  buildFromSpec(movedMedian)$run(0L, 4L)$train,
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
    forest(basis = ~ g == "b", n.trees = 7L)
  ),
  control = seededControlForestBasisR5()
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
    seed = 41L
  )
}
expect_silent(unevenTrees$setControl(treeCountControl(13L)))
expect_error(
  unevenTrees$setControl(treeCountControl(31L)),
  "changing .n.trees. is not available on an existing sampler"
)
expect_error(
  unevenTrees$setControl(treeCountControl(7L)),
  "changing .n.trees. is not available on an existing sampler"
)

# --- dbartsData(bases = ) takes a CONTINUOUS basis, the positive route the
# "not a numeric basis column" refusal relaxed into: creation succeeds, the
# column lands in the engine, and its single amplitude really multiplies it -
# asserted through the fitted-value reconstruction rather than by dimensions
# alone, so a basis accepted and dropped would still fail. ---
continuous <- x[, 1L] - 0.5
continuousData <- dbartsData(x, y, bases = list(NULL, cbind(continuous)))
continuousFit <- dbarts(continuousData, control = seededControlForestBasisR5())
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
expect_true(dbarts:::samplerCarriesAmplitudes(threeForests))
plain <- dbarts(x, y, control = seededControlForestBasisR5())
expect_false(dbarts:::samplerCarriesAmplitudes(plain))
# a plain sampler has no basis to change at all, so the refusal names that
# rather than claiming a forest index is out of range
expect_error(
  plain$setForestBasis(1L, other),
  paste0(
    "$setForestBasis is not available on a sampler that carries no forest ",
    "amplitudes: amplitudes are fixed at creation; make a new sampler instead"
  ),
  fixed = TRUE
)
# the four refusals it gates stay red on a multi-forest sampler
expect_error(threeForests$setResponse(y, updateScale = TRUE), "does not supp")
expect_error(
  threeForests$setOffset(rep(0, n), updateScale = TRUE),
  "does not supp"
)
expect_error(threeForests$setData(dbartsData(x, y)), "does not support")
expect_error(threeForests$setModel(threeForests$model), "does not support")
# and $setLeafPrior takes the no-op creation also accepts
threeLeafPrior <- threeForests$getLeafPrior()
threeForests$setLeafPrior(normal(k = 2))
expect_identical(threeForests$getLeafPrior(), threeLeafPrior)
# and the multinomial route is NOT misidentified by it
mnLabels <- sample(0:2, n, replace = TRUE)
mn <- dbarts(
  x,
  factor(mnLabels, levels = 0:2),
  family = "multinomial",
  control = seededControlForestBasisR5()
)
expect_error(
  mn$getForestAmplitudes(1L),
  "forests carry them"
)
expect_error(
  mn$setForestBasis(2L, cbind(0.5, 0.5)),
  "carry no amplitudes"
)

# --- the leaf scale a swap leaves behind. The calibration map divides
# out the median nonzero row norm of each forest's basis, and that divisor is
# otherwise a construction-time constant - no mutation re-derives a K-forest
# leaf scale - so $setForestBasis owns the staleness. Pinned on a PROBIT
# K-forest, where the anchor is the literal 1 and $getLeafPrior's
# anchor (k pinned at 1) IS the map's leaf scale, so the assertion is exact. Nothing
# else in the suite or in the equivalence trio calls this mutator at all. ---
yBinary <- as.double(y > median(y))
probitForests <- function() {
  dbarts(
    x,
    yBinary,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = seededControlForestBasisR5()
  )
}
priorScale <- function(sampler, forest) {
  sampler$getLeafPrior(forest)$k.scale
}

# (i) STALENESS: a basis whose median nonzero row norm is 4x the old one moves
# the divisor with it. RED against a divisor left at its construction value.
# A forest of levels takes no numbers, so the forest is created on a number.
stale <- dbarts(
  x,
  yBinary,
  forests = list(forest(), forest(basis = z)),
  control = seededControlForestBasisR5()
)
expect_equal(priorScale(stale, 2L), 1 / 0.674, tolerance = 1e-12)
stale$run(5L, 1L)
stale$setForestBasis(2L, 4 * z)
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
swapped$setForestBasis(2L, factor(z))
expect_identical(priorScale(swapped, 2L), priorScale(twin, 2L))
swappedResult <- swapped$run(5L, 2L)
twinResult <- twin$run(5L, 2L)
expect_identical(swapped$getForestFits(2L), twin$getForestFits(2L))
expect_identical(swapped$getForestAmplitudes(), twin$getForestAmplitudes())
expect_identical(swappedResult$sigma, twinResult$sigma)
expect_identical(swappedResult$train, twinResult$train)

# (iii) THE WHOLE PRODUCT, re-derived across the swap. (i) and (ii) both run at
# the DEFAULT leaf scale factor - a literal 1 - so a swap that lost the stored
# factor would pass them both; only a DECLARED forest(sd = ) can see it. Every
# factor of sd_f * s / (0.674 c_f) discriminates here: sd in {2.5, 0.4}, c in
# {3, 5} then 7, divisor 0.674, and s is probit's 1 (the anchor's own arm is
# test-bcf-family.R's, which runs the same fixture under logistic too, where a
# dropped anchor is a 81 percent miss).
ones <- matrix(1, n, 1L)
# No forest of this model is without a basis, so the count and the tree prior
# are stated on a forest and the control names no count.
declared <- dbarts(
  dbartsData(x, yBinary, bases = list(3 * ones, 5 * ones)),
  forests = list(
    forest(sd = 2.5, n.trees = 25L, base = 0.95, power = 2),
    forest(sd = 0.4)
  ),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 6L,
    updateState = FALSE,
    seed = 41L
  )
)
expect_equal(priorScale(declared, 1L), 2.5 / (0.674 * 3), tolerance = 1e-12)
expect_equal(priorScale(declared, 2L), 0.4 / (0.674 * 5), tolerance = 1e-12)
declared$run(5L, 1L)
declared$setForestBasis(2L, 7 * ones)
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
  sampler$getLeafPrior(forest)[[column]]
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

# --- the SPEC ECHO across a state install. A state carries no leaf scale and
# no fixed amplitude variance - both are model - so the recipient's map is
# still the decomposition in force after $setState, whatever the donor ran
# under. Four arms. ---
donorForests <- function(sd, basis = ~ factor(z)) {
  multiplied <- list(basis = basis, sd = sd)
  dbarts(
    x,
    yBinary,
    forests = list(forest(), do.call(forest, multiplied)),
    control = seededControlForestBasisR5()
  )
}
donor <- donorForests(2)
recipient <- donorForests(0.5)
# the two really do differ in the calibration, which is what makes the install
# foreign rather than a re-statement
expect_true(priorScale(donor, 2L) != priorScale(recipient, 2L))
donor$run(5L, 1L)
donor$storeState()
donorState <- donor$state
recipientScale <- priorScale(recipient, 2L)

# (a) THE RECIPIENT'S CALIBRATION STANDS. The install is accepted - the
# widths and tree counts agree - and the recipient runs under its own scale,
# which its factor and divisor still decompose.
recipient$setState(donorState)
expect_identical(priorScale(recipient, 2L), recipientScale)
expect_identical(mapColumn(recipient, 2L, "leaf.scale.factor"), 0.5)
expect_identical(mapColumn(recipient, 2L, "leaf.scale.divisor"), 0.674)
expect_equal(
  mapColumn(recipient, 2L, "k.scale") *
    mapColumn(recipient, 2L, "leaf.scale.divisor") *
    mapColumn(recipient, 2L, "basis.row.norm") /
    mapColumn(recipient, 2L, "leaf.scale.factor"),
  1,
  tolerance = 1e-12
)
# the row norm needs no rule - bases are not state, so the recipient's own is
# still the one in force
expect_equal(mapColumn(recipient, 2L, "basis.row.norm"), 1)
expect_equal(mapColumn(recipient, 1L, "leaf.scale.factor"), 1)
expect_equal(mapColumn(recipient, 1L, "leaf.scale.divisor"), 1)

# (d) THE AMPLITUDE PRIOR IS THE RECIPIENT'S, which is what the next draw will
# use. Forest 1 carries the scale mixture, whose serialized variance is a live
# auxiliary rather than a prior, and forest 2 none.
expect_null(mapColumn(recipient, 2L, "amplitude.prior.scale"))
# neither sampler declares that forest's sd, so both carry the family's own
# default median - 1 under this probit fixture, where the anchor is the link's
# error sd and sigma is pinned, rather than gaussian's 2
expect_equal(mapColumn(recipient, 1L, "amplitude.prior.scale"), 1)

# (c) RE-DERIVATION. $setForestBasis re-derives the leaf scale from the stored
# factor and divisor and the new basis, and the identity holds on it. The
# recipient here multiplies a number, which is what takes a rescaled column.
donor <- donorForests(2, as.double(z))
donor$run(5L, 1L)
donor$storeState()
recipient <- donorForests(0.5, as.double(z))
recipient$setState(donor$state)
recipient$setForestBasis(2L, 4 * z)
expect_equal(mapColumn(recipient, 2L, "leaf.scale.factor"), 0.5)
expect_equal(mapColumn(recipient, 2L, "leaf.scale.divisor"), 0.674)
expect_equal(mapColumn(recipient, 2L, "basis.row.norm"), 4)
expect_equal(
  mapColumn(recipient, 2L, "k.scale") *
    mapColumn(recipient, 2L, "leaf.scale.divisor") *
    mapColumn(recipient, 2L, "basis.row.norm") /
    mapColumn(recipient, 2L, "leaf.scale.factor"),
  1,
  tolerance = 1e-12
)

# (b) SELF-RESTORE. A store, a run and a restore of a sampler's OWN state
# leaves every entry and the identity as they were.
selfRestore <- donorForests(2)
selfRestore$run(5L, 1L)
selfRestore$storeState()
ownState <- selfRestore$state
ownScale <- priorScale(selfRestore, 2L)
selfRestore$run(0L, 1L)
selfRestore$setState(ownState)
expect_equal(priorScale(selfRestore, 2L), ownScale)
expect_equal(mapColumn(selfRestore, 2L, "leaf.scale.factor"), 2)
expect_equal(mapColumn(selfRestore, 2L, "leaf.scale.divisor"), 0.674)
expect_equal(
  mapColumn(selfRestore, 2L, "k.scale") *
    mapColumn(selfRestore, 2L, "leaf.scale.divisor") *
    mapColumn(selfRestore, 2L, "basis.row.norm") /
    mapColumn(selfRestore, 2L, "leaf.scale.factor"),
  1,
  tolerance = 1e-12
)

# The C entry re-checks what the R method already refuses, and it holds an
# owning transposition buffer while it does: a non-finite value is rejected
# part way through, and the sampler it left untouched still runs.
direct <- twoForests()
direct$run(5L, 2L)
badBasis <- cbind(rep(1, n), c(rep(0, n - 1L), Inf))
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setForestBasis,
    direct$getPointer(),
    1L,
    badBasis
  ),
  pattern = "not finite"
)
expect_true(all(is.finite(direct$run(0L, 2L)$sigma)))
expect_equal(ncol(direct$data@bases[[2L]]), 2L)

# --- edge bases: a factor level no kept row takes is no column at creation,
# in every spelling, while setForestBasis keeps a level of the forest that
# the current data leave empty (test-bcf-r5-surface.R); a row
# norm that overflows poisons the calibration that divides by it, at either
emptyLevel <- factor(z, levels = 0:2)
# written as code, a level no kept row has is no column, as in lm()
expect_identical(
  colnames(
    dbarts(
      x,
      y,
      forests = list(forest(), forest(basis = emptyLevel)),
      control = seededControlForestBasisR5()
    )$data@bases[[2L]]
  ),
  c("emptyLevel0", "emptyLevel1")
)
# and handed over as a value it is the same two columns
fromValue <- dbarts(
  x,
  y,
  forests = list(forest(), do.call(forest, list(basis = emptyLevel))),
  control = seededControlForestBasisR5()
)
expect_identical(fromValue$data@basis.levels, list(NULL, c("0", "1")))
expect_identical(fromValue$data@bases[[2L]], unname(cbind(1 - z, z)))
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = rep(1e300, n))),
    control = seededControlForestBasisR5()
  ),
  "row's norm is not representable"
)
edge <- numericFirst()
expect_error(edge$setForestBasis(1L, rep(1e300, n)), "not representable")
expect_true(all(is.finite(edge$run(0L, 2L)$train)))
