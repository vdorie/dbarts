# The mid-chain half of the named calibration
# (docs/design/nameable-calibration.md): $getCalibration is the authoritative
# reader of the leaf prior in force, and $setCalibration rewrites its scale
# half on every chain. The oracles are the two fidelity directions - a read
# followed by a write must be BITWISE inert, and a write followed by a read
# must return what was written - plus the refusal matrix and every mutation
# channel the reported value must not surprise on.

set.seed(41)
n <- 120L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
y <- 12 * (x[, 1L] - x[, 2L]) + rnorm(n)

midControl <- function(n.chains = 2L, n.trees = 20L, ...) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = n.trees,
    n.samples = 10L,
    updateState = FALSE,
    seed = 23L,
    keepTrees = FALSE,
    ...
  )
}
namedSampler <- function(response = y, ...) {
  dbarts(
    x,
    response,
    control = midControl(),
    node.prior = normal(k = 2, scale = 1.5),
    ...
  )
}
priorScaleOf <- function(sampler, forest = 1L) {
  sampler$getCalibration(forest)[, "prior.scale"]
}

# the reported shape: one row per chain, the twelve documented columns, and the
# leaf-model tag on an attribute. The EXACT-SET form is what makes a reordering
# visible; a subset check would not see one.
plain <- dbarts(x, y, control = midControl())
calibration <- plain$getCalibration()
expect_equal(dim(calibration), c(2L, 12L))
expect_identical(
  colnames(calibration),
  c(
    "prior.scale",
    "prior.sd",
    "prior.mean",
    "k",
    "k.has.hyperprior",
    "response.scale",
    "response.shift",
    "amplitude.prior.variance",
    "amplitude.prior.scale",
    "node.scale.factor",
    "node.scale.divisor",
    "basis.row.norm"
  )
)
expect_identical(attr(calibration, "leaf.model"), "constant")
# and the five calibration-map columns are NaN on a single-forest sampler: its
# leaf scale is not map-derived, which the reader says positively rather than
# by reporting a plausible 1 a caller would multiply by
mapColumns <- c(
  "amplitude.prior.variance",
  "amplitude.prior.scale",
  "node.scale.factor",
  "node.scale.divisor",
  "basis.row.norm"
)
expect_true(all(is.nan(calibration[, mapColumns])))
# it reads the ENGINE, so an unnamed model reports the family-keyed default
# converted to response units: node.scale 0.5 times the response range
expect_equal(
  unname(calibration[1L, "prior.scale"]),
  0.5 * (max(y) - min(y))
)
expect_equal(
  unname(calibration[1L, "prior.sd"]),
  unname(calibration[1L, "prior.scale"] / calibration[1L, "k"])
)
expect_equal(unname(calibration[1L, "prior.mean"]), (max(y) + min(y)) / 2)
expect_equal(unname(calibration[1L, "response.scale"]), max(y) - min(y))
expect_true(all(calibration[, "k.has.hyperprior"] == 0))
# and a named model reports what it named, not the range
expect_equal(unname(priorScaleOf(namedSampler())), c(1.5, 1.5))

# --- a get-then-set is BITWISE inert. The setter derives the internal
# leaf scale from the value handed back and SKIPS the write when either that
# scale or the prior.scale the getter would report reproduces what is in
# force, so a round trip cannot perturb the last bit and move a draw. ---
inertA <- namedSampler()
inertB <- namedSampler()
inertB$setCalibration(prior.scale = priorScaleOf(inertB)[1L])
expect_identical(inertA$run(20L, 10L)$train, inertB$run(20L, 10L)$train)
# the same on the UNNAMED default, whose in-force value is an inherited range
# rather than a round number and so is the harder round trip
inertC <- dbarts(x, y, control = midControl())
inertD <- dbarts(x, y, control = midControl())
inertD$setCalibration(prior.scale = priorScaleOf(inertD)[1L])
expect_identical(inertC$run(20L, 10L)$train, inertD$run(20L, 10L)$train)
# the arms above are inert under any implementation whose response-unit round
# trip happens to be exact, so they cannot see a missing skip on their own.
# This response scaling was CHOSEN because its round trip ROUNDS - measured,
# the derived internal scale lands one ulp below 1.11803398874989481926e-01 -
# which makes the arm a deterministic falsifier rather than a lucky one.
yRounding <- 3 * y
leafScales <- function(sampler) {
  sampler$storeState()
  vapply(
    sampler$state,
    function(chain) chain$forests[[1L]]$leaf.scale,
    numeric(1L)
  )
}
inertF <- dbarts(x, yRounding, control = midControl())
scaleBefore <- leafScales(inertF)
inertF$setCalibration(prior.scale = priorScaleOf(inertF)[1L])
expect_identical(leafScales(inertF), scaleBefore)
inertG <- dbarts(x, yRounding, control = midControl())
inertH <- dbarts(x, yRounding, control = midControl())
inertH$setCalibration(prior.scale = priorScaleOf(inertH)[1L])
expect_identical(inertG$run(20L, 10L)$train, inertH$run(20L, 10L)$train)
# non-vacuity: a write of anything else is not inert at all
inertE <- namedSampler()
inertE$setCalibration(prior.scale = 1.5 * (1 + 1e-9))
expect_false(identical(inertA$run(20L, 10L)$train, inertE$run(20L, 10L)$train))

# --- set-then-get fidelity, on EVERY chain. This is the only member with
# power against a setter error that is a fixed function of the named value: a
# sqrt(m)-forgetting setter reports 10.6066 for a requested 1.5 and passes the
# arm-vs-arm oracles of test-calibration-creation.R untouched.
#
# DEVIATION from an expected "bitwise" round trip: the reported value is the internal
# scale times the transform, and the internal scale is the requested value
# DIVIDED by that transform, so exactness is a property of the particular
# (value, transform) pair rather than of the implementation. MEASURED:
# (P / f) * f != P for 194596 of 2e6 random positive pairs (9.7%). This
# fixture happens to sit on the lucky side - at its own transform the round
# trip is bitwise exact for all four requests at n.trees 20, 50 and 200, and
# the one rounding cell is m = 25 at P = 0.25, half an ulp out. Asserting
# bitwise would therefore pin an accident of this response vector and this
# tree count, which any fixture edit or re-anchoring channel could break with
# no defect behind it. The m = 25 arm below carries that rounding cell on
# purpose, so the tolerance is exercised rather than merely permitted.
# Recovering the last bit in general would mean either caching the named
# value in the engine - which the state format and every re-anchoring channel
# would then have to carry - or nudging the leaf scale in force off the exact
# quotient. The assertion is therefore ulp-level, which is 15 orders of
# magnitude below the error it exists to catch, and the BITWISE half that IS
# implementation-determined - every chain reporting the same bits - is
# asserted as such. ---
fidelity <- namedSampler()
for (requested in c(1.5, 0.25, 3.75, 12)) {
  fidelity$setCalibration(prior.scale = requested)
  reported <- priorScaleOf(fidelity)
  expect_true(max(abs(reported / requested - 1)) < 4 * .Machine$double.eps)
  expect_identical(reported[[1L]], reported[[2L]])
}
# the rounding cell: m = 25 at P = 0.25 reports half an ulp below what was
# asked for, so this arm is the one that distinguishes the shipped tolerance
# from an equality and would fail a bitwise assertion outright
rounding <- dbarts(
  x,
  y,
  control = midControl(n.trees = 25L),
  node.prior = normal(k = 2, scale = 1.5)
)
rounding$setCalibration(prior.scale = 0.25)
roundingReported <- priorScaleOf(rounding)
expect_true(max(abs(roundingReported / 0.25 - 1)) < 4 * .Machine$double.eps)
expect_identical(roundingReported[[1L]], roundingReported[[2L]])
expect_false(identical(unname(roundingReported[[1L]]), 0.25))
# the sd spelling is the same statement at the k in force
fidelity$setCalibration(prior.sd = 0.75)
expect_true(max(abs(priorScaleOf(fidelity) / 1.5 - 1)) < 1e-14)
expect_true(
  max(abs(fidelity$getCalibration()[, "prior.sd"] / 0.75 - 1)) < 1e-14
)

# --- the static m falsifier. Two arms at DIFFERENT tree counts are not
# bitwise comparable even under a correct implementation, so this shape - one
# number read twice, and one drawn twice - is what has power over m. The read
# and the write share a single conversion helper by construction, which makes
# a factor error in ONE of them unrepresentable but leaves a factor error in
# BOTH invisible to the round trip; the two absolute pins below are what close
# that, one on each side. A sqrt(m)-forgetting pair reports 10.60660 at m = 50
# and 21.21320 at m = 200, a ratio of exactly 2. ---
staticSampler <- function(numTrees) {
  dbarts(
    x,
    y,
    control = midControl(n.chains = 1L, n.trees = numTrees),
    node.prior = normal(k = 2)
  )
}
# the READ, absolutely: an unnamed model runs the family-keyed node scale of
# 0.5, whose response-unit reading is 0.5 times the range at every tree count
staticRead <- vapply(
  c(50L, 200L),
  function(numTrees) priorScaleOf(staticSampler(numTrees))[[1L]],
  numeric(1L)
)
expect_true(max(abs(staticRead / (0.5 * (max(y) - min(y))) - 1)) < 1e-12)
# the round trip, at both counts
staticM <- vapply(
  c(50L, 200L),
  function(numTrees) {
    sampler <- staticSampler(numTrees)
    sampler$setCalibration(prior.scale = 1.5)
    priorScaleOf(sampler)[[1L]]
  },
  numeric(1L)
)
expect_true(max(abs(staticM / 1.5 - 1)) < 1e-14)
expect_true(abs(staticM[2L] / staticM[1L] - 1) < 1e-14)
# the WRITE, absolutely: prior draws of the forest total after the write have
# sd prior.scale / k at every tree count, which is the quantity the reported
# number claims and the only check that crosses the engine's own draw law.
# The standard error of an sd estimate is about sd / sqrt(2R), so 600 draws
# support a 10% band.
staticDrawn <- vapply(
  c(50L, 200L),
  function(numTrees) {
    sampler <- staticSampler(numTrees)
    sampler$setCalibration(prior.scale = 1.5)
    set.seed(3)
    draws <- vapply(
      seq_len(600L),
      function(draw) {
        sampler$sampleTreesFromPrior(updateState = FALSE)
        sampler$sampleNodeParametersFromPrior(updateState = FALSE)
        sampler$predict(x[1L, , drop = FALSE])[[1L]]
      },
      numeric(1L)
    )
    sd(draws)
  },
  numeric(1L)
)
expect_true(max(abs(staticDrawn / 0.75 - 1)) < 0.1)

# --- chains that have diverged are REFUSED, not flattened. The write
# itself is total - one prior.scale is one statement in response units and
# every chain takes it - but the prior.sd sugar divides by a per-chain k, and
# once the chains disagree about k there is no single answer. ---
diverged <- namedSampler()
diverged$run(10L, 5L)
diverged$storeState()
divergedState <- diverged$state
divergedState[[2L]]$forests[[1L]]$k <- 4
diverged$setState(divergedState)
expect_equal(unname(diverged$getCalibration()[, "k"]), c(2, 4))
expect_error(
  diverged$setCalibration(prior.sd = 0.5),
  "chains' 'k' have diverged"
)
# nothing was written by the refusal, and the k-free spelling still serves
expect_equal(unname(priorScaleOf(diverged)), c(1.5, 1.5))
diverged$setCalibration(prior.scale = 3)
expect_true(max(abs(priorScaleOf(diverged) / 3 - 1)) < 1e-14)
# a divergent LEAF SCALE is likewise reported rather than hidden, and the
# total write flattens it, which is the documented every-chain rule
divergedScale <- namedSampler()
divergedScale$storeState()
scaleState <- divergedScale$state
scaleState[[2L]]$forests[[1L]]$leaf.scale <-
  2 * scaleState[[2L]]$forests[[1L]]$leaf.scale
divergedScale$setState(scaleState)
expect_equal(unname(priorScaleOf(divergedScale)), c(1.5, 3))
divergedScale$setCalibration(prior.scale = 1.5)
expect_true(max(abs(priorScaleOf(divergedScale) / 1.5 - 1)) < 1e-14)

# --- the calibration refusal matrix, mid-chain half. ---

# a value that is not a positive finite number is an ERROR, not a refusal
badValues <- namedSampler()
expect_error(badValues$setCalibration(prior.scale = -1), "must be positive")
expect_error(badValues$setCalibration(prior.scale = 0), "must be positive")
expect_error(badValues$setCalibration(prior.scale = Inf), "must be positive")
expect_error(badValues$setCalibration(prior.scale = NaN), "must be positive")
expect_error(badValues$setCalibration(prior.scale = NA_real_), "positive")
expect_error(badValues$setCalibration(prior.scale = c(1, 2)), "single number")
expect_error(badValues$setCalibration(prior.sd = -1), "must be positive")
# exactly one spelling
expect_error(badValues$setCalibration(), "exactly one")
expect_error(
  badValues$setCalibration(prior.scale = 1, prior.sd = 1),
  "exactly one"
)
# NaN is refused at creation too, where is.na() would otherwise read it as the
# unnamed spelling and let it reach the bridge
expect_error(dbartsPriors$normal(scale = NaN), "'scale' must be positive")
expect_error(
  validObject(new("dbartsModel", prior.scale = NaN)),
  "prior.scale must be NA"
)

# the forest index is 1-based on this surface, and out of range is refused
expect_error(badValues$getCalibration(2L), "forest index out of range")
expect_error(badValues$setCalibration(prior.scale = 1, forest = 2L), "range")
expect_error(badValues$getCalibration(0L), "single positive integer")

# prior.mean is read-only; the message names the offset lever and the recipe
meanRefusal <- tryCatch(
  badValues$setCalibration(prior.mean = 0),
  error = function(e) conditionMessage(e)
)
expect_true(is.character(meanRefusal))
expect_true(grepl("setOffset", meanRefusal, fixed = TRUE))
expect_true(grepl("prior.mean", meanRefusal, fixed = TRUE))
# the lever the message names is the offset channel, and the reported quantity
# is the transform's shift: an offset issued at the default updateScale = FALSE
# shifts the modelled quantity while leaving the reported mean pinned, which is
# exactly the contract test-calibration-creation.R's offset-recipe test
# exercises end to end.
recipe <- namedSampler()
recipeMean <- recipe$getCalibration()[1L, "prior.mean"]
expect_equal(unname(recipeMean), (max(y) + min(y)) / 2)
recipe$setOffset(rep_len(-recipeMean, n))
expect_identical(recipe$getCalibration()[1L, "prior.mean"], recipeMean)

# the sd spelling under a sampled k, with both remedies named. The binary
# default IS the sampled path, so this is an ordinary probit sampler.
yBinary <- rbinom(n, 1L, pnorm(x[, 1L] - x[, 2L]))
sampledK <- dbarts(x, yBinary, control = midControl())
expect_true(all(sampledK$getCalibration()[, "k.has.hyperprior"] == 1))
sdRefusal <- tryCatch(
  sampledK$setCalibration(prior.sd = 1.5),
  error = function(e) conditionMessage(e)
)
expect_true(is.character(sdRefusal))
expect_true(grepl("'prior.scale'", sdRefusal, fixed = TRUE))
expect_true(grepl("pin 'k'", sdRefusal, fixed = TRUE))
# prior.scale is honored EXACTLY under the same hyperprior - the hyperprior
# scales it rather than replacing it - and survives the k draws of a run
sampledK$setCalibration(prior.scale = 1.5)
expect_true(max(abs(priorScaleOf(sampledK) / 1.5 - 1)) < 1e-14)
sampledK$run(20L, 10L)
expect_true(max(abs(priorScaleOf(sampledK) / 1.5 - 1)) < 1e-14)
expect_true(any(sampledK$getCalibration()[, "k"] != 2))

# a two-forest sampler: the getter serves it per forest, the setter refuses it
# by name, because the calibration map owns both halves
zBCF <- rbinom(n, 1L, 0.5)
yBCF <- 12 * (x[, 1L] - x[, 2L]) + 2 * zBCF + rnorm(n)
bcf <- dbarts(
  x,
  yBCF,
  forests = list(forest(), forest(basis = ~ factor(zBCF))),
  control = midControl()
)
bcfCalibration <- bcf$getCalibration(1L)
expect_equal(dim(bcfCalibration), c(2L, 12L))
expect_true(all(bcfCalibration[, "prior.scale"] > 0))
expect_true(all(bcf$getCalibration(2L)[, "prior.scale"] > 0))
# BCF pins k at 1 per its map, which the sampler-wide option does not say
expect_true(all(bcfCalibration[, "k"] == 1))
expect_true(all(bcfCalibration[, "k.has.hyperprior"] == 0))
# here the five map columns are the ones IN FORCE, and the two amplitude
# columns are EXCLUSIVE per forest: forest 1 declares no basis, so it carries
# the half-Cauchy scale mixture and reports no variance, and forest 2 the
# reverse. Which one it is agrees with the transported params, so the reader
# and the creation route cannot disagree about a forest's amplitude channel.
bcfParams <- attr(bcf$control, "bartcore.forests")$params
expect_true(all(is.nan(bcfCalibration[, "amplitude.prior.variance"])))
expect_equal(
  unname(bcfCalibration[, "amplitude.prior.scale"]),
  rep_len(bcfParams[[1L]][7L], 2L)
)
bcfCalibration2 <- bcf$getCalibration(2L)
expect_equal(
  unname(bcfCalibration2[, "amplitude.prior.variance"]),
  rep_len(bcfParams[[2L]][6L], 2L)
)
expect_true(all(is.nan(bcfCalibration2[, "amplitude.prior.scale"])))
# and the anchor s the map states every node scale against is recoverable from
# the reported decomposition, which is the only route to it under gaussian,
# where s is the data-dependent scaled response sd: the two forests recover
# the SAME s, and it is the one the response transform implies
recoveredAnchor <- function(row) {
  unname(
    row[, "prior.scale"] *
      row[, "node.scale.divisor"] *
      row[, "basis.row.norm"] /
      row[, "node.scale.factor"]
  )
}
expect_equal(
  recoveredAnchor(bcfCalibration2),
  recoveredAnchor(bcfCalibration),
  tolerance = 1e-12
)
expect_true(all(recoveredAnchor(bcfCalibration) > 0))
expect_error(
  bcf$setCalibration(prior.scale = 1.5),
  "multi-forest calibration map"
)
# the engine's own refusal, which the low-level route reaches past the R
# guard, names the map by its coupling: the two-forest map at K = 2, and a
# generic one above it, where no two-forest map owns the scale
handleOf <- function(sampler) list(ptr = sampler$getPointer())
expect_error(
  dbarts:::bartcoreSetForestPriorScale(handleOf(bcf), 0L, 1.5),
  "two-forest calibration map"
)
threeForests <- dbarts(
  x,
  yBCF,
  forests = list(
    forest(),
    forest(basis = ~ factor(zBCF)),
    forest(basis = x[, 3L])
  ),
  control = midControl()
)
threeRefusal <- tryCatch(
  dbarts:::bartcoreSetForestPriorScale(handleOf(threeForests), 0L, 1.5),
  error = function(e) conditionMessage(e)
)
expect_true(grepl("multi-forest calibration map", threeRefusal, fixed = TRUE))
expect_false(grepl("two-forest", threeRefusal, fixed = TRUE))

# the multinomial coupling, through the low-level handle its forests live on
labels <- sample(0:2, n, replace = TRUE)
host <- dbarts(x, y, control = midControl(n.chains = 1L))
multinomial <- dbarts:::bartcoreMultinomialSampler(host, labels, K = 3L)
multinomialCalibration <- dbarts:::bartcoreForestCalibration(multinomial, 0L)
expect_equal(dim(multinomialCalibration), c(1L, 12L))
expect_true(multinomialCalibration[1L, "prior.scale"] > 0)
# a K-forest sampler with no calibration map: the five map columns are NaN,
# which a forest-count test would get wrong in the other direction
expect_true(all(is.nan(multinomialCalibration[, mapColumns])))
# the softmax map works on a unit-scale latent and fixes k
expect_equal(unname(multinomialCalibration[1L, "response.scale"]), 1)
expect_true(multinomialCalibration[1L, "k.has.hyperprior"] == 0)
expect_error(
  dbarts:::bartcoreSetForestPriorScale(multinomial, 0L, 1.5),
  "softmax calibration map"
)

# a host shell carries the design and priors of a fit whose model lives
# elsewhere, so a write through it would never reach that fit
set.seed(43)
multinomialFit <- bart2(
  x,
  factor(labels),
  family = "multinomial",
  keepTrees = TRUE,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.trees = 10L,
  verbose = FALSE
)
expect_error(
  multinomialFit$fit$setCalibration(prior.scale = 1.5),
  "host sampler of a bart2"
)

# DART is NOT refused. $setModel refuses a DART sampler outright, so this is
# the concrete gain: the calibration moves and the fit runs clean afterwards.
dartSampler <- dbarts(
  x,
  y,
  control = midControl(),
  tree.prior = dart(),
  node.prior = normal(k = 2, scale = 1.5)
)
expect_error(dartSampler$setModel(dartSampler$model), "DART")
dartSampler$setCalibration(prior.scale = 0.4)
expect_true(max(abs(priorScaleOf(dartSampler) / 0.4 - 1)) < 1e-14)
dartRun <- dartSampler$run(40L, 10L)
expect_true(all(is.finite(dartRun$train)))
expect_true(sum(dartRun$varcount) > 0)

# the write is TOTAL over the four leaf models - each carries the one scale
# field - and the tag says which, so a reader knows whether the reported
# prior.sd is an equality or a bound
leafArms <- list(
  constant = dbarts(x, y, control = midControl()),
  monotone = dbarts(x, y, control = midControl(), monotone = c(x1 = 1)),
  linear = dbarts(x, y, control = midControl(), node.prior = linear("x1")),
  gp = dbarts(x, y, control = midControl(), node.prior = gp("x1"))
)
for (tag in names(leafArms)) {
  sampler <- leafArms[[tag]]
  expect_identical(attr(sampler$getCalibration(), "leaf.model"), tag)
  sampler$setCalibration(prior.scale = 1.5)
  expect_true(
    max(abs(priorScaleOf(sampler) / 1.5 - 1)) < 1e-14,
    info = tag
  )
  expect_true(all(is.finite(sampler$run(10L, 5L)$train)), info = tag)
}

# --- the mid-chain half of the mutation table. Each row is one assertion
# about what the reported prior.scale does. ---

# the AUTHORITY RULE: the model slot records the named INTENT and the engine
# holds what is in force, and only a re-anchoring channel parts them
authority <- namedSampler()
expect_equal(unname(priorScaleOf(authority)), c(1.5, 1.5))
authority$setResponse(y / 8, updateScale = TRUE)
expect_equal(authority$model@prior.scale, 1.5)
expect_true(all(abs(priorScaleOf(authority) - 1.5) > 1e-8))
expect_equal(unname(priorScaleOf(authority)), rep_len(1.5 / 8, 2L))
# and re-issuing the intent is either spelling of the recipe
authority$setModel(authority$model)
expect_true(max(abs(priorScaleOf(authority) / 1.5 - 1)) < 1e-14)

# setResponse / setOffset at the default updateScale = FALSE hold the
# transform, so the calibration in force does not move: the composition path
held <- namedSampler()
heldScale <- priorScaleOf(held)
held$setResponse(y / 8)
expect_identical(priorScaleOf(held), heldScale)
held$setOffset(rep_len(3, n))
expect_identical(priorScaleOf(held), heldScale)
held$setWeights(runif(n, 0.5, 2))
expect_identical(priorScaleOf(held), heldScale)
held$setSigma(2.5)
expect_identical(priorScaleOf(held), heldScale)
# setOffset at updateScale = TRUE re-anchors and the getter shows it
held$setOffset(rep_len(3, n), updateScale = TRUE)
expect_true(all(priorScaleOf(held) != heldScale))
# setData always re-anchors, so it moves
replaced <- namedSampler()
replaced$setData(dbartsData(x, y / 8))
expect_true(all(abs(priorScaleOf(replaced) - 1.5) > 1e-8))

# $setCalibration touches nothing else: not sigma, not the tree prior, not the
# response transform, not k
isolated <- namedSampler()
isolated$run(10L, 5L)
sigmaBefore <- isolated$getSigmas()
isolatedBefore <- isolated$getCalibration()
isolated$setCalibration(prior.scale = 4)
isolatedAfter <- isolated$getCalibration()
expect_identical(isolated$getSigmas(), sigmaBefore)
expect_identical(
  isolatedAfter[, c("prior.mean", "k", "response.scale", "response.shift")],
  isolatedBefore[, c("prior.mean", "k", "response.scale", "response.shift")]
)

# setModel also re-pins sigma for a gaussian sampler with no variance forest,
# which is the channel a caller reaching for setModel to restate a calibration
# would fire as well; setCalibration is the lever that does not
repinned <- dbarts(
  x,
  y,
  control = midControl(),
  node.prior = normal(k = 2, scale = 1.5),
  resid.prior = dbartsPriors$fixed(1)
)
repinned$setSigma(3.5)
expect_equal(unname(repinned$getSigmas()[1L]), 3.5)
repinned$setModel(repinned$model)
expect_equal(unname(repinned$getSigmas()[1L]), 1)
expect_true(max(abs(priorScaleOf(repinned) / 1.5 - 1)) < 1e-14)

# storeState / setState adopt the calibration from the state, which is what
# makes the getter - rather than the stale model slot - the authoritative
# reader after a restore
adopted <- namedSampler()
adopted$setCalibration(prior.scale = 0.6)
adopted$storeState()
donorState <- adopted$state
recipient <- namedSampler()
recipient$setState(donorState)
expect_identical(priorScaleOf(recipient), priorScaleOf(adopted))
expect_equal(recipient$model@prior.scale, 1.5)

# a warm start ADOPTS the donor's calibration - its trees were drawn under it -
# and the documented recipe is to re-issue the write afterwards
warmDonor <- namedSampler()
warmDonor$setCalibration(prior.scale = 0.6)
warmDonor$run(20L, 10L)
warmDonor$storeState()
warmRecipient <- namedSampler()
warmRecipient$installTrees(warmDonor)
expect_true(max(abs(priorScaleOf(warmRecipient) / 0.6 - 1)) < 1e-14)
warmRecipient$setCalibration(prior.scale = 1.5)
expect_true(max(abs(priorScaleOf(warmRecipient) / 1.5 - 1)) < 1e-14)

# --- the save/load gate: updateState = TRUE captures the write, and the
# calibration survives the serialize/re-create round trip. ---
roundTrip <- function(object) {
  tempFile <- tempfile()
  on.exit(unlink(tempFile))
  saveRDS(object, file = tempFile)
  readRDS(tempFile)
}
saved <- namedSampler()
saved$run(10L, 5L)
saved$storeState()
saved$setCalibration(prior.scale = 0.6, updateState = TRUE)
restored <- roundTrip(saved)
expect_true(max(abs(priorScaleOf(restored) / 0.6 - 1)) < 1e-14)
expect_identical(priorScaleOf(restored), priorScaleOf(saved))
# non-vacuity: without the capture the write does not reach the saved state,
# which is the recorded save/load defect the argument exists for
uncaptured <- namedSampler()
uncaptured$run(10L, 5L)
uncaptured$storeState()
uncaptured$setCalibration(prior.scale = 0.6)
expect_true(max(abs(priorScaleOf(roundTrip(uncaptured)) / 1.5 - 1)) < 1e-14)
