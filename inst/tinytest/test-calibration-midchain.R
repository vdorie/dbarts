# The mid-chain half of the leaf prior named by k or sd: $getLeafPrior reads
# the leaf prior in force in the terms it was named in, $getK each chain's k, and $setLeafPrior
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
# the spread in force on each chain, which the reader states as anchor / k
priorSdOf <- function(sampler, forest = 1L) {
  sampler$getLeafPrior(forest)$k.scale / sampler$getK(forest)
}
# a specification built in the constructors' own vocabulary
priorOf <- function(expr) eval(substitute(expr), dbartsPriors, parent.frame())
# the engine's own reading, beneath the R restatement
engineReading <- function(sampler, forest = 1L) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    forest - 1L
  )
}

# the reported shape: the prior alone, a list of the seven documented
# elements. The EXACT-SET form is what makes a reordering visible; a subset
# check would not see one.
plain <- dbarts(x, y, control = midControl())
calibration <- plain$getLeafPrior()
expect_identical(
  names(calibration),
  c(
    "leaf.prior",
    "leaf.model",
    "prior.sd.of",
    "prior.mean",
    "k.scale",
    "response.scale",
    "response.shift"
  )
)
expect_identical(calibration$leaf.model, "constant")
expect_identical(calibration$prior.sd.of, "leaf value")
# and the five calibration-map entries are absent on a single-forest sampler:
# its leaf scale is not map-derived, which the reader says positively rather
# than by reporting a plausible 1 a caller would multiply by
mapColumns <- c(
  "amplitude.prior.variance",
  "amplitude.prior.scale",
  "leaf.scale.factor",
  "leaf.scale.divisor",
  "basis.row.norm"
)
expect_false(any(mapColumns %in% names(calibration)))
# the unnamed default states the k it resolved to
expect_identical(calibration$leaf.prior, priorOf(normal(k = 2)))
# it reads the ENGINE, so an unnamed model reports the family-keyed default
# converted to response units: leaf.scale 0.5 times the response range, and
# getK is bitwise the engine's k
expect_equal(calibration$k.scale, 0.5 * (max(y) - min(y)))
expect_identical(plain$getK(), engineReading(plain)[, "k"])
expect_identical(plain$getK(), c(2, 2))
expect_equal(calibration$prior.mean, (max(y) + min(y)) / 2)
expect_equal(calibration$response.scale, max(y) - min(y))
# a named model states what it named; getK is the engine's k, the reference 2
# relative to the named anchor
named <- namedSampler()
namedReading <- named$getLeafPrior()
expect_identical(namedReading$leaf.prior, priorOf(normal(sd = 0.75)))
expect_identical(namedReading$k.scale, 1.5)
expect_identical(named$getK(), c(2, 2))
# the leaf model qualifies what the sd is the sd of
expect_identical(
  dbarts(
    x,
    y,
    control = midControl(),
    leaf.prior = linear("x1")
  )$getLeafPrior()$prior.sd.of,
  "coefficient"
)
expect_identical(
  dbarts(
    x,
    y,
    control = midControl(),
    leaf.prior = gp("x1")
  )$getLeafPrior()$prior.sd.of,
  "amplitude"
)

# --- a get-then-set is BITWISE inert. The setter writes the anchor the read
# implies and the engine SKIPS a write reproducing what is in force, so a
# round trip cannot perturb the last bit and move a draw. ---
leafScales <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_getLeafPrior, sampler$getPointer(), 0L)[,
    "prior.scale"
  ]
}
inertA <- namedSampler()
inertB <- namedSampler()
inertB$setLeafPrior(inertB$getLeafPrior()$leaf.prior)
expect_identical(inertA$run(20L, 10L)$train, inertB$run(20L, 10L)$train)
# the named value itself, written again: the writer derives the internal scale
# with creation's arithmetic, so it lands on creation's bits. This response
# scaling was chosen because the two orders of that arithmetic differ by an
# ulp on it (measured), so the arm falsifies a writer that reorders them
inertSdA <- dbarts(
  x,
  y / 5,
  control = midControl(),
  leaf.prior = normal(sd = 1)
)
inertSdB <- dbarts(
  x,
  y / 5,
  control = midControl(),
  leaf.prior = normal(sd = 1)
)
sdScaleBefore <- leafScales(inertSdB)
inertSdB$setLeafPrior(normal(sd = 1))
expect_identical(leafScales(inertSdB), sdScaleBefore)
expect_identical(inertSdA$run(20L, 10L)$train, inertSdB$run(20L, 10L)$train)
# the same on the UNNAMED default, whose in-force value is an inherited range
# rather than a round number and so is the harder round trip
inertC <- dbarts(x, y, control = midControl())
inertD <- dbarts(x, y, control = midControl())
inertD$setLeafPrior(normal(sd = priorSdOf(inertD)[[1L]]))
expect_identical(inertC$run(20L, 10L)$train, inertD$run(20L, 10L)$train)
# a scaling whose response-unit round trip ROUNDS, so the arm needs the skip
yRounding <- 3 * y
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
expect_identical(binaryRead$leaf.prior, priorOf(normal(k = chi(1.5, 2))))
binaryB$setLeafPrior(normal(sd = invchi(1.5, binaryRead$k.scale / 2)))
expect_identical(binaryA$run(20L, 10L)$train, binaryB$run(20L, 10L)$train)
# and the reader's own output written back, for a drawn k and for each
# spelling of an sd law, which reads its scale off the anchor in force
roundTrip <- function(leafPrior) {
  twins <- lapply(1:2, function(i) {
    eval(bquote(dbarts(
      x,
      yBinary,
      control = midControl(),
      leaf.prior = .(leafPrior)
    )))
  })
  twins[[2L]]$setLeafPrior(twins[[2L]]$getLeafPrior()$leaf.prior)
  expect_identical(
    twins[[1L]]$run(20L, 10L)$train,
    twins[[2L]]$run(20L, 10L)$train,
    info = deparse(leafPrior)
  )
}
roundTrip(quote(normal(k = chi(1.25, 3))))
roundTrip(quote(normal(sd = invchi(2, 0.3))))
roundTrip(quote(normal(sd = invchi(2, 0))))
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
expect_identical(
  improperB$getLeafPrior()$leaf.prior,
  priorOf(normal(k = chi(1.25, Inf)))
)
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
fidelityLaw <- fidelity$getLeafPrior()$leaf.prior@prior.sd
expect_identical(fidelityLaw@df, 2)
expect_true(abs(fidelityLaw@scale / 0.5 - 1) < 4 * .Machine$double.eps)
# and a k spelling returns the reader to k terms, relative to the data
fidelity$setLeafPrior(normal(k = 3))
expect_identical(fidelity$getLeafPrior()$leaf.prior, priorOf(normal(k = 3)))
expect_identical(fidelity$getK(), c(3, 3))
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
  function(numTrees) staticSampler(numTrees)$getLeafPrior()$k.scale,
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

# --- a state holds neither the leaf scale nor a fixed k, so one that names
# them - an earlier writer's blocks, here edited on one chain - installs
# without them and every chain stays under the sampler's prior; the reader
# never reports NA, and an NA spread is still refused on write ---
divergedScale <- namedSampler()
priorBefore <- divergedScale$getLeafPrior()
divergedScale$storeState()
scaleState <- divergedScale$state
expect_null(scaleState[[2L]]$forests[[1L]]$leaf.scale)
scaleState[[2L]]$forests[[1L]]$leaf.scale <- 0.123
divergedScale$setState(scaleState)
expect_identical(divergedScale$getLeafPrior(), priorBefore)
expect_identical(divergedScale$getK(), c(2, 2))
naSd <- divergedScale$getLeafPrior()$leaf.prior
naSd@prior.sd <- NA_real_
expect_error(divergedScale$setLeafPrior(naSd), "'sd' is NA")
naModel <- divergedScale$model
naModel@leaf.prior <- naSd
expect_error(divergedScale$setModel(naModel), "'sd' is NA.*name a value")
expect_error(
  dbarts(x, y, control = midControl(), leaf.prior = naSd),
  "'sd' is NA"
)
expect_error(new("dbartsNormalPrior", prior.sd = NA_real_), "'sd' is NA")
divergedK <- dbarts(x, y, control = midControl())
divergedK$storeState()
kState <- divergedK$state
expect_null(kState[[2L]]$forests[[1L]]$k)
kState[[2L]]$forests[[1L]]$k <- 3
divergedK$setState(kState)
expect_identical(divergedK$getK(), c(2, 2))
naK <- divergedK$getLeafPrior()$leaf.prior
expect_identical(naK@k, 2)
naK@k <- NA_real_
expect_error(divergedK$setLeafPrior(naK), "'k' is NA")
expect_identical(divergedK$getK(), c(2, 2))
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
  "takes no 'forest' index"
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
# the constructors keep requiring columns; only the sampler's own write may
# omit them
expect_error(dbartsPriors$linear(sd = 1), "requires 'columns'")
expect_error(dbartsPriors$gp(sd = 1), "requires 'columns'")
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
expect_true(is(
  detachedWrite$getLeafPrior()$leaf.prior@prior.sd,
  "dbartsSdHyperprior"
))

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
recipeMean <- recipe$getLeafPrior()$prior.mean
expect_equal(recipeMean, (max(y) + min(y)) / 2)
recipe$setOffset(rep_len(-recipeMean, n))
expect_identical(recipe$getLeafPrior()$prior.mean, recipeMean)

# a drawn k: the sd law is written and survives the k draws of a run, which
# move the reported spread and not the law
sampledK <- dbarts(x, yBinary, control = midControl())
sdLawScale <- function(sampler) sampler$getLeafPrior()$leaf.prior@prior.sd@scale
expect_true(is(sampledK$getLeafPrior()$leaf.prior@k, "dbartsChiHyperprior"))
sampledK$setLeafPrior(normal(sd = invchi(1.5, 0.75)))
expect_true(abs(sdLawScale(sampledK) / 0.75 - 1) < 1e-14)
sampledRun <- sampledK$run(20L, 10L)
expect_true(abs(sdLawScale(sampledK) / 0.75 - 1) < 1e-14)
expect_true(any(priorSdOf(sampledK) != 0.75))
# getK is the value the run records, bitwise its last draw
expect_identical(sampledK$getK(), sampledRun$k[10L, ])
# and a number fixes it, which is a change of the hyperprior itself
sampledK$setLeafPrior(normal(sd = 1.5))
expect_identical(sampledK$getLeafPrior()$leaf.prior, priorOf(normal(sd = 1.5)))
expect_equal(unname(priorSdOf(sampledK)), c(1.5, 1.5))

# a two-forest sampler: the getter serves it per forest, and the setter takes
# only a forest's spread, forests = list(forest(sd = )), since the calibration
# map owns the leaf scales
zBCF <- rbinom(n, 1L, 0.5)
yBCF <- 12 * (x[, 1L] - x[, 2L]) + 2 * zBCF + rnorm(n)
bcf <- dbarts(
  x,
  yBCF,
  forests = list(forest(), forest(basis = ~ factor(zBCF))),
  control = midControl()
)
bcfCalibration <- bcf$getLeafPrior(1L)
bcfCalibration2 <- bcf$getLeafPrior(2L)
# BCF pins k at 1 per its map, which the sampler-wide option does not say, so
# each forest states its spread as creation does, forest(sd = ): the
# half-Cauchy median on the forest without a basis, the leaf-scale factor on
# the one with
expect_identical(bcf$getK(1L), c(1, 1))
expect_identical(
  bcfCalibration$leaf.prior,
  dbartsForests$forest(sd = bcfCalibration$amplitude.prior.scale)
)
expect_identical(
  bcfCalibration2$leaf.prior,
  dbartsForests$forest(sd = bcfCalibration2$leaf.scale.factor)
)
expect_identical(bcfCalibration$prior.sd.of, "amplitude scale")
expect_identical(bcfCalibration2$prior.sd.of, "forest total")
expect_true(bcfCalibration$k.scale > 0 && bcfCalibration2$k.scale > 0)
# here the map entries are the ones IN FORCE, and the two amplitude entries
# are EXCLUSIVE per forest: forest 1 declares no basis, so it carries the
# half-Cauchy scale mixture and no variance, and forest 2 the reverse
bcfParams <- attr(bcf$control, "bartcore.forests")$params
expect_identical(
  names(bcfCalibration)[-(1:7)],
  c(
    "amplitude.prior.scale",
    "leaf.scale.factor",
    "leaf.scale.divisor",
    "basis.row.norm"
  )
)
expect_null(bcfCalibration$amplitude.prior.variance)
expect_equal(bcfCalibration$amplitude.prior.scale, bcfParams[[1L]][7L])
expect_equal(bcfCalibration2$amplitude.prior.variance, bcfParams[[2L]][6L])
expect_null(bcfCalibration2$amplitude.prior.scale)

# forest = NULL lists every forest's prior; a single-forest sampler's NULL
# read is its forest = 1 read
expect_identical(plain$getLeafPrior(), plain$getLeafPrior(1L))
expect_identical(bcf$getLeafPrior(), list(bcfCalibration, bcfCalibration2))
expect_identical(plain$getK(), plain$getK(1L))
expect_identical(bcf$getK(), rbind(bcf$getK(1L), bcf$getK(2L)))

# and the anchor s the map states every leaf scale against is recoverable from
# the reported decomposition (k is pinned at 1, so the anchor is the map's leaf
# scale): the two forests recover the SAME s
recoveredAnchor <- function(prior) {
  prior$k.scale *
    prior$leaf.scale.divisor *
    prior$basis.row.norm /
    prior$leaf.scale.factor
}
expect_equal(
  recoveredAnchor(bcfCalibration2),
  recoveredAnchor(bcfCalibration),
  tolerance = 1e-12
)
expect_true(recoveredAnchor(bcfCalibration) > 0)
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
expect_true(priorSdOf(multinomial) > 0)
expect_false(any(mapColumns %in% names(multinomialCalibration)))
# the softmax map works on a unit-scale latent and fixes k
expect_equal(multinomialCalibration$response.scale, 1)
expect_true(is.numeric(multinomialCalibration$leaf.prior@k))
expect_identical(nrow(multinomial$getK()), 3L)
multinomial$setLeafPrior(normal(k = 3))
expect_true(all(multinomial$getK() == 3))
expect_error(
  multinomial$setLeafPrior(normal(sd = 1)),
  "normal\\(k = \\).*softmax calibration map"
)

# $fit is the K-forest engine that ran, not a host shell: a named sd is
# refused for the softmax's own reason
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
  multinomialFit$fit$setLeafPrior(normal(sd = 1)),
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
expect_identical(dartSampler$getK(), c(3, 3))
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
  expect_identical(sampler$getLeafPrior()$leaf.model, tag)
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
kBefore <- isolated$getK()
isolated$setLeafPrior(normal(sd = 2))
isolatedAfter <- isolated$getLeafPrior()
expect_identical(isolated$getSigmas(), sigmaBefore)
expect_identical(isolated$getK(), kBefore)
unmoved <- c("prior.mean", "response.scale", "response.shift")
expect_identical(isolatedAfter[unmoved], isolatedBefore[unmoved])

# $setModel re-pins a fixed sigma at the model's value, which $setSigma
# rewrites; $setLeafPrior does not re-pin, whether it writes the spread alone
# or changes the hyperprior
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
expect_equal(repinned$model@resid.prior@value, 3.5^2)
repinned$setModel(repinned$model)
expect_equal(unname(repinned$getSigmas()), c(3.5, 3.5))

# storeState / setState leave the recipient's calibration: the leaf prior is
# the sampler's model, and the state carries none of it
adopted <- namedSampler()
adopted$setLeafPrior(normal(sd = 0.3))
adopted$storeState()
donorState <- adopted$state
recipient <- namedSampler()
recipient$setState(donorState)
expect_true(max(abs(priorSdOf(recipient) / 0.75 - 1)) < 1e-14)
expect_equal(recipient$model@prior.scale, 1.5)

# a warm start leaves it too: the donor's trees seed the recipient, which runs
# under its own prior
warmDonor <- namedSampler()
warmDonor$setLeafPrior(normal(sd = 0.3))
invisible(warmDonor$run(20L, 10L))
warmDonor$storeState()
warmRecipient <- namedSampler()
warmRecipient$installTrees(warmDonor)
expect_true(max(abs(priorSdOf(warmRecipient) / 0.75 - 1)) < 1e-14)

# --- the save/load gate: the write is recorded on the model, so the prior
# survives the serialize/re-create round trip whether or not the state was
# stored after it. ---
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
uncaptured <- namedSampler()
invisible(uncaptured$run(10L, 5L))
uncaptured$storeState()
uncaptured$setLeafPrior(normal(sd = 0.3))
expect_true(
  max(abs(priorSdOf(roundTripCalibrationMidchain(uncaptured)) / 0.3 - 1)) <
    1e-14
)

# --- a write into a drawn prior keeps each chain's spread in force at the
# call, k re-expressed against the new k.scale, whatever the prior before it:
# drawn, a fixed sd or a fixed k, in either spelling. The new prior acts from
# the next draw of k. Each case runs 200 sweeps on two chains first, so the
# chains' spreads differ from each other and from any named value. ---
spreadAfter <- function(from, to, response = y) {
  sampler <- eval(bquote(dbarts(
    x,
    response,
    control = midControl(),
    leaf.prior = .(from)
  )))
  invisible(sampler$run(200L, 1L))
  before <- priorSdOf(sampler)
  eval(bquote(sampler$setLeafPrior(.(to))))
  list(sampler = sampler, before = before, after = priorSdOf(sampler))
}
keepsSpread <- function(from, to, drawn = TRUE) {
  kept <- spreadAfter(from, to)
  info <- paste(deparse(from), "to", deparse(to))
  # a drawn spread has moved apart on the two chains; a fixed one has not
  expect_true(drawn == (kept$before[[1L]] != kept$before[[2L]]), info = info)
  expect_true(
    max(abs(kept$after / kept$before - 1)) < 1e-14,
    info = info
  )
  kept$sampler
}
keepsSpread(quote(normal(k = chi(1.5, 2))), quote(normal(sd = invchi(3, 2))))
keepsSpread(quote(normal(sd = invchi(3, 2))), quote(normal(k = chi(1.5, 2))))
keepsSpread(quote(normal(sd = invchi(3, 2))), quote(normal(sd = invchi(3, 4))))
# from a fixed sd (dec-B392), and between the spellings out of a fixed prior
# (dec-B393): the fixed spread is the one kept
fixedSd <- keepsSpread(
  quote(normal(sd = 0.5)),
  quote(normal(sd = invchi(3, 2))),
  drawn = FALSE
)
expect_true(max(abs(priorSdOf(fixedSd) / 0.5 - 1)) < 1e-14)
fixedK <- spreadAfter(quote(normal(k = 3)), quote(normal(sd = invchi(3, 2))))
expect_true(max(abs(fixedK$after / fixedK$before - 1)) < 1e-14)
fixedSdToK <- spreadAfter(
  quote(normal(sd = 0.5)),
  quote(normal(k = chi(1.5, 2)))
)
expect_true(max(abs(fixedSdToK$after / 0.5 - 1)) < 1e-14)
# and the new prior acts from the next draw: the spreads move on
invisible(fixedSd$run(5L, 1L))
expect_true(all(priorSdOf(fixedSd) != 0.5))
# under the k spelling a changed chi() leaves k.scale, so k is left as it is
chiScale <- dbarts(
  x,
  y,
  control = midControl(),
  leaf.prior = normal(k = chi(1.5, 2))
)
invisible(chiScale$run(200L, 1L))
chiK <- chiScale$getK()
chiScale$setLeafPrior(normal(k = chi(1.5, 4)))
expect_identical(chiScale$getK(), chiK)
# a fixed value stated sets the spread as written, from a drawn k
toFixed <- spreadAfter(quote(normal(k = chi(1.5, 2))), quote(normal(sd = 0.25)))
expect_true(max(abs(toFixed$after / 0.25 - 1)) < 1e-14)
# a re-anchor under a drawn invchi() is not a leaf-prior write: k is left
reanchored <- dbarts(
  x,
  y,
  control = midControl(),
  leaf.prior = normal(sd = invchi(3, 2))
)
invisible(reanchored$run(200L, 1L))
reanchoredK <- reanchored$getK()
reanchored$setResponse(2 * y + 3, updateScale = TRUE)
expect_identical(reanchored$getK(), reanchoredK)
