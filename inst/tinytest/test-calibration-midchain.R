# The mid-chain half of the leaf prior named by k or sd: $getLeafPrior reads
# the leaf prior in force in the terms it was named in, and $setLeafPrior
# restates its spread, or the hyperprior it is drawn under, on every chain, in
# the creation vocabulary. The oracles are the two fidelity directions - a
# read followed by a write must be BITWISE inert, and a write followed by a
# read must return what was written - plus the refusal matrix and every
# mutation channel the reported value must not surprise on.

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
    leaf.prior = normal(sd = 0.75),
    ...
  )
}
priorSdOf <- function(sampler, forest = 1L) {
  sampler$getLeafPrior(forest)[, "prior.sd"]
}
# the engine's own reading, beneath the R restatement
engineReading <- function(sampler, forest = 1L) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    forest - 1L
  )
}

# the reported shape: one row per chain, the fourteen documented columns, and
# the leaf-model tag and what prior.sd is the sd of on attributes. The
# EXACT-SET form is what makes a reordering visible; a subset check would not
# see one.
plain <- dbarts(x, y, control = midControl())
calibration <- plain$getLeafPrior()
expect_equal(dim(calibration), c(2L, 14L))
expect_identical(
  colnames(calibration),
  c(
    "prior.sd",
    "prior.sd.df",
    "prior.sd.scale",
    "prior.mean",
    "k",
    "k.has.hyperprior",
    "anchor",
    "response.scale",
    "response.shift",
    "amplitude.prior.variance",
    "amplitude.prior.scale",
    "leaf.scale.factor",
    "leaf.scale.divisor",
    "basis.row.norm"
  )
)
expect_identical(attr(calibration, "leaf.model"), "constant")
expect_identical(attr(calibration, "prior.sd.of"), "leaf value")
# and the five calibration-map columns are NaN on a single-forest sampler: its
# leaf scale is not map-derived, which the reader says positively rather than
# by reporting a plausible 1 a caller would multiply by
mapColumns <- c(
  "amplitude.prior.variance",
  "amplitude.prior.scale",
  "leaf.scale.factor",
  "leaf.scale.divisor",
  "basis.row.norm"
)
expect_true(all(is.nan(calibration[, mapColumns])))
# k is fixed, so no sd law is in force
expect_true(all(is.nan(calibration[, c("prior.sd.df", "prior.sd.scale")])))
# it reads the ENGINE, so an unnamed model reports the family-keyed default
# converted to response units: leaf.scale 0.5 times the response range, and
# its k is bitwise the engine's
expect_equal(unname(calibration[1L, "anchor"]), 0.5 * (max(y) - min(y)))
expect_identical(calibration[, "k"], engineReading(plain)[, "k"])
expect_equal(
  unname(calibration[1L, "prior.sd"]),
  unname(calibration[1L, "anchor"] / calibration[1L, "k"])
)
expect_equal(unname(calibration[1L, "prior.mean"]), (max(y) + min(y)) / 2)
expect_equal(unname(calibration[1L, "response.scale"]), max(y) - min(y))
expect_true(all(calibration[, "k.has.hyperprior"] == 0))
# a named model reports what it named, and its k relative to the data's anchor
named <- namedSampler()
namedReading <- named$getLeafPrior()
expect_equal(unname(namedReading[, "prior.sd"]), c(0.75, 0.75))
expect_equal(unname(namedReading[, "anchor"]), rep(0.5 * (max(y) - min(y)), 2L))
expect_equal(
  unname(namedReading[, "k"]),
  unname(namedReading[, "anchor"] / namedReading[, "prior.sd"])
)
# the engine's own k there is the reference 2, which the reader restates
expect_identical(unname(engineReading(named)[, "k"]), c(2, 2))
# the leaf model qualifies what prior.sd is the sd of
expect_identical(
  attr(
    dbarts(
      x,
      y,
      control = midControl(),
      leaf.prior = linear("x1")
    )$getLeafPrior(),
    "prior.sd.of"
  ),
  "coefficient"
)
expect_identical(
  attr(
    dbarts(x, y, control = midControl(), leaf.prior = gp("x1"))$getLeafPrior(),
    "prior.sd.of"
  ),
  "amplitude"
)

# --- a get-then-set is BITWISE inert. The setter writes the anchor the read
# implies and the engine SKIPS a write reproducing what is in force, so a
# round trip cannot perturb the last bit and move a draw. ---
inertA <- namedSampler()
inertB <- namedSampler()
inertB$setLeafPrior(normal(sd = priorSdOf(inertB)[[1L]]))
expect_identical(inertA$run(20L, 10L)$train, inertB$run(20L, 10L)$train)
# the same on the UNNAMED default, whose in-force value is an inherited range
# rather than a round number and so is the harder round trip
inertC <- dbarts(x, y, control = midControl())
inertD <- dbarts(x, y, control = midControl())
inertD$setLeafPrior(normal(sd = priorSdOf(inertD)[[1L]]))
expect_identical(inertC$run(20L, 10L)$train, inertD$run(20L, 10L)$train)
# a scaling whose response-unit round trip ROUNDS, so the arm needs the skip
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
inertF$setLeafPrior(normal(sd = priorSdOf(inertF)[[1L]]))
expect_identical(leafScales(inertF), scaleBefore)
inertG <- dbarts(x, yRounding, control = midControl())
inertH <- dbarts(x, yRounding, control = midControl())
inertH$setLeafPrior(normal(sd = priorSdOf(inertH)[[1L]]))
expect_identical(inertG$run(20L, 10L)$train, inertH$run(20L, 10L)$train)
# a binary default: the read sd law written back is inert, and so is the k
# spelling of the unnamed default
yBinary <- rbinom(n, 1L, pnorm(x[, 1L] - x[, 2L]))
binaryA <- dbarts(x, yBinary, control = midControl())
binaryB <- dbarts(x, yBinary, control = midControl())
binaryRead <- binaryB$getLeafPrior()
expect_identical(unname(binaryRead[, "prior.sd.df"]), c(1.5, 1.5))
expect_identical(unname(binaryRead[, "prior.sd.scale"]), c(1.5, 1.5))
binaryB$setLeafPrior(normal(
  sd = invchi(binaryRead[1L, "prior.sd.df"], binaryRead[1L, "prior.sd.scale"])
))
expect_identical(binaryA$run(20L, 10L)$train, binaryB$run(20L, 10L)$train)
# the improper limit reads a zero scale and writes back as nothing at all
improperA <- dbarts(
  x,
  yBinary,
  control = midControl(),
  leaf.prior = normal(k = chi(1.25, Inf))
)
improperB <- dbarts(
  x,
  yBinary,
  control = midControl(),
  leaf.prior = normal(k = chi(1.25, Inf))
)
expect_identical(unname(improperB$getLeafPrior()[, "prior.sd.scale"]), c(0, 0))
improperB$setLeafPrior(normal(sd = invchi(1.25, 0)))
expect_identical(improperA$run(20L, 10L)$train, improperB$run(20L, 10L)$train)
# non-vacuity: a write of anything else is not inert at all
inertE <- namedSampler()
inertE$setLeafPrior(normal(sd = 0.75 * (1 + 1e-9)))
expect_false(identical(inertA$run(20L, 10L)$train, inertE$run(20L, 10L)$train))

# --- set-then-get fidelity, on EVERY chain. The reported value is the
# internal scale times the transform, and the internal scale is the requested
# value divided by it, so exactness is a property of the (value, transform)
# pair rather than of the implementation; the assertion is ulp-level, and the
# BITWISE half that IS implementation-determined - every chain reporting the
# same bits - is asserted as such. ---
fidelity <- namedSampler()
for (requested in c(0.75, 0.125, 1.875, 6)) {
  fidelity$setLeafPrior(normal(sd = requested))
  reported <- priorSdOf(fidelity)
  expect_true(max(abs(reported / requested - 1)) < 4 * .Machine$double.eps)
  expect_identical(reported[[1L]], reported[[2L]])
}
# a drawn spread reports the law it was written as
fidelity$setLeafPrior(normal(sd = invchi(2, 0.5)))
expect_identical(unname(fidelity$getLeafPrior()[, "prior.sd.df"]), c(2, 2))
expect_true(
  max(abs(fidelity$getLeafPrior()[, "prior.sd.scale"] / 0.5 - 1)) <
    4 * .Machine$double.eps
)
expect_true(all(fidelity$getLeafPrior()[, "k.has.hyperprior"] == 1))
# and a k spelling returns the reader to k terms, relative to the data
fidelity$setLeafPrior(normal(k = 3))
expect_identical(unname(fidelity$getLeafPrior()[, "k"]), c(3, 3))
expect_true(all(is.na(fidelity$model@prior.scale)))

# --- the static m falsifier. Two arms at DIFFERENT tree counts are not
# bitwise comparable even under a correct implementation, so this shape - one
# number read twice, and one drawn twice - is what has power over m. A
# sqrt(m)-forgetting pair reports a ratio of exactly 2 between m = 50 and
# m = 200. ---
staticSampler <- function(numTrees) {
  dbarts(
    x,
    y,
    control = midControl(n.chains = 1L, n.trees = numTrees),
    leaf.prior = normal(k = 2)
  )
}
# the READ, absolutely: an unnamed model runs the family-keyed leaf scale of
# 0.5, whose response-unit reading is 0.5 times the range at every tree count
staticRead <- vapply(
  c(50L, 200L),
  function(numTrees) staticSampler(numTrees)$getLeafPrior()[1L, "anchor"],
  numeric(1L)
)
expect_true(max(abs(staticRead / (0.5 * (max(y) - min(y))) - 1)) < 1e-12)
# the round trip, at both counts
staticM <- vapply(
  c(50L, 200L),
  function(numTrees) {
    sampler <- staticSampler(numTrees)
    sampler$setLeafPrior(normal(sd = 0.75))
    priorSdOf(sampler)[[1L]]
  },
  numeric(1L)
)
expect_true(max(abs(staticM / 0.75 - 1)) < 1e-14)
expect_true(abs(staticM[2L] / staticM[1L] - 1) < 1e-14)
# the WRITE, absolutely: prior draws of the forest total after the write have
# sd prior.sd at every tree count, the only check that crosses the engine's
# own draw law. 600 draws support a 10% band.
staticDrawn <- vapply(
  c(50L, 200L),
  function(numTrees) {
    sampler <- staticSampler(numTrees)
    sampler$setLeafPrior(normal(sd = 0.75))
    set.seed(3)
    draws <- vapply(
      seq_len(600L),
      function(draw) {
        sampler$sampleTreesFromPrior(updateState = FALSE)
        sampler$sampleLeafParametersFromPrior(updateState = FALSE)
        sampler$predict(x[1L, , drop = FALSE])[[1L]]
      },
      numeric(1L)
    )
    sd(draws)
  },
  numeric(1L)
)
expect_true(max(abs(staticDrawn / 0.75 - 1)) < 0.1)

# --- a divergent LEAF SCALE is reported rather than hidden, and the write
# flattens it, which is the documented every-chain rule ---
divergedScale <- namedSampler()
divergedScale$storeState()
scaleState <- divergedScale$state
scaleState[[2L]]$forests[[1L]]$leaf.scale <-
  2 * scaleState[[2L]]$forests[[1L]]$leaf.scale
divergedScale$setState(scaleState)
expect_equal(unname(priorSdOf(divergedScale)), c(0.75, 1.5))
divergedScale$setLeafPrior(normal(sd = 0.75))
expect_true(max(abs(priorSdOf(divergedScale) / 0.75 - 1)) < 1e-14)

# --- the refusal matrix, mid-chain half. ---

# a value that is not a positive finite number is an ERROR, not a refusal
badValues <- namedSampler()
expect_error(badValues$setLeafPrior(normal(sd = -1)), "must be positive")
expect_error(badValues$setLeafPrior(normal(sd = 0)), "must be positive")
expect_error(badValues$setLeafPrior(normal(sd = Inf)), "must be positive")
expect_error(badValues$setLeafPrior(normal(sd = NaN)), "must be positive")
expect_error(badValues$setLeafPrior(normal(sd = NA_real_)), "must not be NA")
expect_error(badValues$setLeafPrior(normal(sd = c(1, 2))), "single number")
expect_error(badValues$setLeafPrior(normal(k = 1, sd = 1)), "not both")
# a specification, and nothing else
expect_error(badValues$setLeafPrior(), "'leaf.prior' must be given")
expect_error(badValues$setLeafPrior(0.5), "leaf prior specification")
expect_error(badValues$setLeafPrior(chisq()), "leaf prior specification")
# the removed number-valued spellings are unused arguments
expect_error(badValues$setLeafPrior(prior.scale = 1), "unused argument")
expect_error(badValues$setLeafPrior(prior.sd = 1), "unused argument")
expect_error(badValues$setLeafPrior(prior.mean = 0), "unused argument")
expect_error(
  badValues$setLeafPrior(normal(sd = 1), forest = 1L),
  "unused argument"
)
# the leaf model is $setModel's to change
expect_error(
  badValues$setLeafPrior(linear("x1", sd = 1)),
  "leaf model is written normal\\(\\), not linear\\(\\).*\\$setModel"
)
linearSampler <- dbarts(
  x,
  y,
  control = midControl(),
  leaf.prior = linear("x1", sd = 0.5)
)
expect_error(
  linearSampler$setLeafPrior(linear("x2", sd = 1)),
  "columns it names differ"
)
linearSampler$setLeafPrior(linear(sd = 0.25))
expect_equal(unname(priorSdOf(linearSampler)), c(0.25, 0.25))
linearSampler$setLeafPrior(linear("x1", sd = 0.5))
expect_equal(unname(priorSdOf(linearSampler)), c(0.5, 0.5))
gpSampler <- dbarts(
  x,
  y,
  control = midControl(),
  leaf.prior = gp("x1", sd = 0.5, lengthscale = 0.3)
)
expect_error(
  gpSampler$setLeafPrior(gp(sd = 1, lengthscale = 0.4)),
  "lengthscale it names differs"
)
expect_error(
  gpSampler$setLeafPrior(gp(sd = 1, max.leaf.size = 10L)),
  "max.leaf.size it names differs"
)
gpSampler$setLeafPrior(gp(sd = 1, lengthscale = 0.3))
expect_equal(unname(priorSdOf(gpSampler)), c(1, 1))
# an sd hyperprior under a monotone constraint, as at creation
monotoneSampler <- dbarts(x, y, control = midControl(), monotone = c(x1 = 1))
expect_error(
  monotoneSampler$setLeafPrior(normal(sd = invchi(1.5, 1))),
  "monotone"
)
# NaN is refused in a hand-built model, where is.na() would otherwise read it
# as the unnamed spelling and let it reach the bridge
expect_error(
  validObject(new("dbartsModel", prior.scale = NaN)),
  "prior.scale must be NA"
)
# the constructors resolve inside the call without the package attached
detachedWrite <- namedSampler()
evalq(
  sampler$setLeafPrior(normal(sd = invchi(1.5, 0.5))),
  list2env(list(sampler = detachedWrite), parent = baseenv())
)
expect_true(all(detachedWrite$getLeafPrior()[, "k.has.hyperprior"] == 1))

# the forest index is 1-based on the reader, and out of range is refused
expect_error(badValues$getLeafPrior(2L), "forest index out of range")
expect_error(badValues$getLeafPrior(0L), "single positive integer")
expect_error(
  badValues$getLeafPrior(1.5),
  "'forest' must be a whole number; got '1.5'",
  fixed = TRUE
)

# prior.mean is read-only: the lever is the offset channel, and the reported
# quantity is the transform's shift, so an offset issued at the default
# updateScale = FALSE shifts the modelled quantity while leaving the reported
# mean pinned
recipe <- namedSampler()
recipeMean <- recipe$getLeafPrior()[1L, "prior.mean"]
expect_equal(unname(recipeMean), (max(y) + min(y)) / 2)
recipe$setOffset(rep_len(-recipeMean, n))
expect_identical(recipe$getLeafPrior()[1L, "prior.mean"], recipeMean)

# a drawn k: the sd law is written and survives the k draws of a run, which
# move the reported spread and not the law
sampledK <- dbarts(x, yBinary, control = midControl())
expect_true(all(sampledK$getLeafPrior()[, "k.has.hyperprior"] == 1))
sampledK$setLeafPrior(normal(sd = invchi(1.5, 0.75)))
expect_true(
  max(abs(sampledK$getLeafPrior()[, "prior.sd.scale"] / 0.75 - 1)) < 1e-14
)
invisible(sampledK$run(20L, 10L))
expect_true(
  max(abs(sampledK$getLeafPrior()[, "prior.sd.scale"] / 0.75 - 1)) < 1e-14
)
expect_true(any(sampledK$getLeafPrior()[, "prior.sd"] != 0.75))
# and a number fixes it, which is a change of the hyperprior itself
sampledK$setLeafPrior(normal(sd = 1.5))
expect_true(all(sampledK$getLeafPrior()[, "k.has.hyperprior"] == 0))
expect_equal(unname(priorSdOf(sampledK)), c(1.5, 1.5))

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
bcfCalibration <- bcf$getLeafPrior(1L)
expect_equal(dim(bcfCalibration), c(2L, 14L))
expect_true(all(bcfCalibration[, "prior.sd"] > 0))
expect_true(all(bcf$getLeafPrior(2L)[, "prior.sd"] > 0))
# BCF pins k at 1 per its map, which the sampler-wide option does not say
expect_true(all(bcfCalibration[, "k"] == 1))
expect_true(all(bcfCalibration[, "k.has.hyperprior"] == 0))
expect_true(all(is.nan(bcfCalibration[, c("prior.sd.df", "prior.sd.scale")])))
# here the five map columns are the ones IN FORCE, and the two amplitude
# columns are EXCLUSIVE per forest: forest 1 declares no basis, so it carries
# the half-Cauchy scale mixture and reports no variance, and forest 2 the
# reverse
bcfParams <- attr(bcf$control, "bartcore.forests")$params
expect_true(all(is.nan(bcfCalibration[, "amplitude.prior.variance"])))
expect_equal(
  unname(bcfCalibration[, "amplitude.prior.scale"]),
  rep_len(bcfParams[[1L]][7L], 2L)
)
bcfCalibration2 <- bcf$getLeafPrior(2L)
expect_equal(
  unname(bcfCalibration2[, "amplitude.prior.variance"]),
  rep_len(bcfParams[[2L]][6L], 2L)
)
expect_true(all(is.nan(bcfCalibration2[, "amplitude.prior.scale"])))

# forest = NULL stacks every forest's reading with the forest margin LAST;
# a single-forest sampler's NULL read is bitwise its forest = 1 read
expect_identical(plain$getLeafPrior(), plain$getLeafPrior(1L))
bcfCalibrationAll <- bcf$getLeafPrior()
expect_equal(dim(bcfCalibrationAll), c(2L, 14L, 2L))
# a slice drops the array-level attributes below, which is ordinary R
# behavior and not part of what this reader promises
stripAttributes <- function(reading) {
  attr(reading, "leaf.model") <- NULL
  attr(reading, "prior.sd.of") <- NULL
  reading
}
expect_identical(bcfCalibrationAll[,, 1L], stripAttributes(bcfCalibration))
expect_identical(bcfCalibrationAll[,, 2L], stripAttributes(bcfCalibration2))
expect_identical(colnames(bcfCalibrationAll), colnames(bcfCalibration))
expect_identical(
  attr(bcfCalibrationAll, "leaf.model"),
  attr(bcfCalibration, "leaf.model")
)
expect_identical(
  attr(bcfCalibrationAll, "prior.sd.of"),
  attr(bcfCalibration, "prior.sd.of")
)

# and the anchor s the map states every leaf scale against is recoverable from
# the reported decomposition (k is pinned at 1, so prior.sd is the map's leaf
# scale): the two forests recover the SAME s
recoveredAnchor <- function(row) {
  unname(
    row[, "prior.sd"] *
      row[, "leaf.scale.divisor"] *
      row[, "basis.row.norm"] /
      row[, "leaf.scale.factor"]
  )
}
expect_equal(
  recoveredAnchor(bcfCalibration2),
  recoveredAnchor(bcfCalibration),
  tolerance = 1e-12
)
expect_true(all(recoveredAnchor(bcfCalibration) > 0))
expect_error(
  bcf$setLeafPrior(normal(sd = 0.75)),
  "multi-forest calibration map.*forest\\(sd = \\)"
)

# the multinomial coupling, through the public sampler its forests live on
labels <- sample(0:2, n, replace = TRUE)
multinomial <- dbarts(
  x,
  factor(labels),
  family = "multinomial",
  control = midControl(n.chains = 1L)
)
multinomialCalibration <- multinomial$getLeafPrior(1L)
expect_equal(dim(multinomialCalibration), c(1L, 14L))
expect_true(multinomialCalibration[1L, "prior.sd"] > 0)
expect_true(all(is.nan(multinomialCalibration[, mapColumns])))
# the softmax map works on a unit-scale latent and fixes k
expect_equal(unname(multinomialCalibration[1L, "response.scale"]), 1)
expect_true(multinomialCalibration[1L, "k.has.hyperprior"] == 0)
expect_error(
  multinomial$setLeafPrior(normal(k = 3)),
  "softmax calibration map.*normal\\(k = \\)"
)

# $fit is the K-forest engine that ran, not a host shell: setLeafPrior
# is refused for the softmax's own reason
set.seed(43)
multinomialFit <- bart(
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
  multinomialFit$fit$setLeafPrior(normal(k = 3)),
  "softmax calibration map"
)

# DART is NOT refused. $setModel refuses a DART sampler outright; the write
# moves the spread, or the law, and the fit runs clean afterwards
dartSampler <- dbarts(
  x,
  y,
  control = midControl(),
  tree.prior = dart(),
  leaf.prior = normal(sd = 0.75)
)
expect_error(dartSampler$setModel(dartSampler$model), "DART")
dartSampler$setLeafPrior(normal(sd = 0.2))
expect_true(max(abs(priorSdOf(dartSampler) / 0.2 - 1)) < 1e-14)
dartSampler$setLeafPrior(normal(k = 3))
expect_identical(unname(dartSampler$getLeafPrior()[, "k"]), c(3, 3))
dartRun <- dartSampler$run(40L, 10L)
expect_true(all(is.finite(dartRun$train)))
expect_true(sum(dartRun$varcount) > 0)

# the write is TOTAL over the four leaf models, each through its own
# constructor, and the tags say which
leafArms <- list(
  constant = list(
    dbarts(x, y, control = midControl()),
    quote(normal(sd = 0.75))
  ),
  monotone = list(
    dbarts(x, y, control = midControl(), monotone = c(x1 = 1)),
    quote(normal(sd = 0.75))
  ),
  linear = list(
    dbarts(x, y, control = midControl(), leaf.prior = linear("x1")),
    quote(linear(sd = 0.75))
  ),
  gp = list(
    dbarts(x, y, control = midControl(), leaf.prior = gp("x1")),
    quote(gp(sd = 0.75))
  )
)
for (tag in names(leafArms)) {
  sampler <- leafArms[[tag]][[1L]]
  expect_identical(attr(sampler$getLeafPrior(), "leaf.model"), tag)
  eval(bquote(sampler$setLeafPrior(.(leafArms[[tag]][[2L]]))))
  expect_true(max(abs(priorSdOf(sampler) / 0.75 - 1)) < 1e-14, info = tag)
  expect_true(all(is.finite(sampler$run(10L, 5L)$train)), info = tag)
}

# --- the mid-chain half of the mutation table. Each row is one assertion
# about what the reported prior.sd does. ---

# a named sd is absolute: a re-anchoring channel moves the transform, and the
# sampler restates the named anchor against it, so the sd holds
absolute <- namedSampler()
absolute$setResponse(y / 8, updateScale = TRUE)
expect_true(max(abs(priorSdOf(absolute) / 0.75 - 1)) < 1e-14)
# the latest write is what is restated, not the creation value
absolute$setLeafPrior(normal(sd = 1))
absolute$setResponse(y, updateScale = TRUE)
expect_true(max(abs(priorSdOf(absolute) / 1 - 1)) < 1e-14)
# the same through setOffset and setData
absolute$setOffset(rep_len(3, n), updateScale = TRUE)
expect_true(max(abs(priorSdOf(absolute) / 1 - 1)) < 1e-14)
absolute$setData(dbartsData(x, y / 8))
expect_true(max(abs(priorSdOf(absolute) / 1 - 1)) < 1e-14)
# non-vacuity: a k-named sampler is relative to the data, so the same channel
# moves its spread with the response
relative <- dbarts(x, y, control = midControl())
relativeSd <- priorSdOf(relative)
relative$setResponse(y / 8, updateScale = TRUE)
expect_equal(priorSdOf(relative), relativeSd / 8)

# setResponse / setOffset at the default updateScale = FALSE hold the
# transform, so the prior in force does not move: the composition path
held <- dbarts(x, y, control = midControl())
heldSd <- priorSdOf(held)
held$setResponse(y / 8)
expect_identical(priorSdOf(held), heldSd)
held$setOffset(rep_len(3, n))
expect_identical(priorSdOf(held), heldSd)
held$setWeights(runif(n, 0.5, 2))
expect_identical(priorSdOf(held), heldSd)
held$setSigma(2.5)
expect_identical(priorSdOf(held), heldSd)

# $setLeafPrior touches nothing else: not sigma, not the response transform,
# not the drawn k in force
isolated <- namedSampler()
invisible(isolated$run(10L, 5L))
sigmaBefore <- isolated$getSigmas()
isolatedBefore <- isolated$getLeafPrior()
isolated$setLeafPrior(normal(sd = 2))
isolatedAfter <- isolated$getLeafPrior()
expect_identical(isolated$getSigmas(), sigmaBefore)
unmoved <- c("prior.mean", "anchor", "response.scale", "response.shift")
expect_identical(isolatedAfter[, unmoved], isolatedBefore[, unmoved])

# $setModel re-pins a fixed sigma; $setLeafPrior does not, whether it writes
# the spread alone or changes the hyperprior
repinned <- dbarts(
  x,
  y,
  control = midControl(),
  leaf.prior = normal(sd = 0.75),
  family = gaussian(sigma = dbartsPriors$fixed(1))
)
repinned$setSigma(3.5)
repinned$setLeafPrior(normal(sd = 0.5))
expect_equal(unname(repinned$getSigmas()), c(3.5, 3.5))
repinned$setLeafPrior(normal(k = 3))
expect_equal(unname(repinned$getSigmas()), c(3.5, 3.5))
repinned$setModel(repinned$model)
expect_equal(unname(repinned$getSigmas()), c(1, 1))

# storeState / setState adopt the calibration from the state, which is what
# makes the getter - rather than the model slot - the authoritative reader
# after a restore
adopted <- namedSampler()
adopted$setLeafPrior(normal(sd = 0.3))
adopted$storeState()
donorState <- adopted$state
recipient <- namedSampler()
recipient$setState(donorState)
expect_identical(priorSdOf(recipient), priorSdOf(adopted))
expect_equal(recipient$model@prior.scale, 1.5)

# a warm start ADOPTS the donor's calibration - its trees were drawn under it -
# and the documented recipe is to restate the prior afterwards
warmDonor <- namedSampler()
warmDonor$setLeafPrior(normal(sd = 0.3))
invisible(warmDonor$run(20L, 10L))
warmDonor$storeState()
warmRecipient <- namedSampler()
warmRecipient$installTrees(warmDonor)
expect_true(max(abs(priorSdOf(warmRecipient) / 0.3 - 1)) < 1e-14)
warmRecipient$setLeafPrior(normal(sd = 0.75))
expect_true(max(abs(priorSdOf(warmRecipient) / 0.75 - 1)) < 1e-14)

# --- the save/load gate: updateState = TRUE captures the write, and the
# prior survives the serialize/re-create round trip. ---
roundTripCalibrationMidchain <- function(object) {
  tempFile <- tempfile()
  on.exit(unlink(tempFile))
  saveRDS(object, file = tempFile)
  readRDS(tempFile)
}
saved <- namedSampler()
invisible(saved$run(10L, 5L))
saved$storeState()
saved$setLeafPrior(normal(sd = 0.3), updateState = TRUE)
restored <- roundTripCalibrationMidchain(saved)
expect_true(max(abs(priorSdOf(restored) / 0.3 - 1)) < 1e-14)
expect_identical(priorSdOf(restored), priorSdOf(saved))
# non-vacuity: without the capture the write does not reach the saved state
uncaptured <- namedSampler()
invisible(uncaptured$run(10L, 5L))
uncaptured$storeState()
uncaptured$setLeafPrior(normal(sd = 0.3))
expect_true(
  max(abs(priorSdOf(roundTripCalibrationMidchain(uncaptured)) / 0.75 - 1)) <
    1e-14
)
