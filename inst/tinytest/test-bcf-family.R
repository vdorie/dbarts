# The K-forest coupling under a non-gaussian response: probit and logistic
# are built, the wider families are doors that refuse by name. What is
# pinned here is (i) the
# ALL-BASIS shape, in which every forest is fixed-variance and nothing in the
# sampler can absorb a mis-stated prior, (ii) the CALIBRATION ANCHOR itself,
# which is the only number here that a running, mixing sampler cannot
# reveal to be wrong, (iii) the guards that fire on the plainest possible
# binary two-forest call unless they are keyed to the family, and (iv) the
# family-keyed refusals that flip the moment the sampler reports its own
# family. The anchor assertions are exact rather than statistical: under a
# latent family the response transform is the identity and the map's sqrt(m)
# cancels, so $getLeafPrior()'s anchor (k pinned at 1) IS the map's leaf scale.

set.seed(29)
n <- 240L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
eta <- -0.4 + 1.5 * x[, 1L] + 0.8 * z
yBalanced <- rbinom(n, 1L, pnorm(eta))
yRare <- rbinom(n, 1L, 0.1)
yContinuous <- eta + rnorm(n, sd = 0.3)
counts <- rpois(n, 2)

uncountedControlBcfFamily <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 4L,
    updateState = FALSE,
    seed = 29L,
    ...
  )
}
seededControlBcfFamily <- function(...) {
  uncountedControlBcfFamily(n.trees = 25L, ...)
}
# A model in which every forest has a basis takes no tree count from the
# control and no tree prior from the fitting function. Its first forest states
# the count and the tree prior the other fits of this file run under; where it
# multiplies a constant column, which has no standard deviation to take a
# default from, it states one unit of the response's scale too, the size a
# factor takes by default.
everyForestHasBasis <- function(bases) !any(vapply(bases, is.null, NA))
unitOf <- function(y, family) {
  switch(family, probit = 1, logistic = pi / sqrt(3), sd(y))
}
statedOnFirst <- function(forests = NULL, sd = NULL) {
  if (is.null(forests)) {
    forests <- list(dbarts::dbartsForests$forest(sd = sd))
  }
  forests[[1L]][c("n.trees", "base", "power")] <- list(25L, 0.95, 2)
  forests
}

# every entry non-NULL, so both/all forests take the fixed-variance channel:
# the shape in which the anchor is load-bearing, since no half-Cauchy
# amplitude and no drawn sigma stands between the prior and the index
ones <- matrix(1, n, 1L)
zFactor <- factor(z)
gLevel <- g == "b"
hFactor <- factor(x[, 1L] > 0.5)
unitBases2 <- list(ones, zFactor)
unitBases3 <- list(ones, zFactor, gLevel)
unitBases4 <- list(ones, zFactor, gLevel, hFactor)
# and one cell whose rows are four, which is what pins that a numeric column
# multiplies its forest as it is given, never divided
scaledBases <- list(4 * ones, zFactor)
# the product pin's own fixture, whose every factor discriminates; see (d)
pinBases <- list(3 * ones, 5 * ones)
# stated in the response's units, as multiples 2.5 and 0.4 of the latent
# family's own scale s: the forest's leaf-scale factors are those multiples
pinScalesFor <- function(s) {
  list(
    dbarts::dbartsForests$forest(sd = 2.5 * s),
    dbarts::dbartsForests$forest(sd = 0.4 * s)
  )
}

basisSampler <- function(y, bases, family = "auto", forests = NULL, ...) {
  resolved <- if (family == "auto") "gaussian" else family
  dbarts(
    dbartsData(x, y, bases = bases),
    control = uncountedControlBcfFamily(),
    family = family,
    forests = statedOnFirst(forests, unitOf(y, resolved)),
    ...
  )
}
priorScales <- function(sampler) {
  unname(vapply(
    seq_along(sampler$data@bases),
    function(f) sampler$getLeafPrior(f)$k.scale,
    numeric(1L)
  ))
}
# the engine's own reading of forest `f`: the half-Cauchy median a drawn
# forest with no basis carries its sd in, which the reader does not report
engineSlot <- function(sampler, f, name) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    f - 1L
  )[[1L, name]]
}

# --- the all-basis fixture itself, gaussian, K = 2 and K = 3. Both halves of
# what makes it the honest test of an anchor: every forest resolves to the
# fixed-variance channel (its own scale its sd over the half-normal median,
# no half-Cauchy amplitude scale), and it runs. ---
for (bases in list(unitBases2, unitBases3, scaledBases)) {
  fit <- basisSampler(yContinuous, bases)
  expect_equal(length(fit$data@bases), length(bases))
  for (forestIndex in seq_along(bases)) {
    reported <- fit$getLeafPrior(forestIndex)
    expect_equal(
      reported$k.scale * 0.674,
      unname(reported$sd),
      tolerance = 1e-12
    )
    expect_true(is.nan(engineSlot(fit, forestIndex, "amplitude.prior.scale")))
  }
  expect_true(all(is.finite(fit$run(4L, 2L)$train)))
}

# --- THE ANCHOR'S OWN GATE. Three assertions per latent family over that
# fixture. (a) leads because the naive alternative - keeping the sample sd of
# the cold-start working response - equals probit's anchor at p = 0.5 up to a
# Bessel factor, so a balanced fixture passes the defect it is meant to catch.
anchors <- list(probit = 1.0, logistic = pi / sqrt(3.0))
for (family in names(anchors)) {
  s <- anchors[[family]]

  # (a) BASE-RATE INVARIANCE: two samplers differing only in the response's
  # base rate report the same calibration, to the bit
  balanced <- basisSampler(yBalanced, scaledBases, family)
  rare <- basisSampler(yRare, scaledBases, family)
  expect_identical(priorScales(balanced), priorScales(rare))

  # (b) THE ABSOLUTE ANCHOR: the map's own expression, with s written as a
  # literal; the column of 4s stating one unit per unit of it is not divided
  expect_equal(priorScales(balanced), c(s, s) / 0.674, tolerance = 1e-12)

  # (c) THE INDUCED INDEX: sqrt(sum_f prior.sd_f^2 v_f ||B_f(i,.)||^2), a
  # numeric coefficient's variance 1 and a level's 1/2. Each forest with a
  # basis defaults to one unit whatever the number of forests (dec-B291), so
  # the map disperses as 1.04912 sqrt(K + 1) with a numeric forest among them:
  # the index prior grows with every forest the caller adds, which is the
  # ruling's stated cost. Three counts, because the claim is the growth in K
  # rather than a value at one K.
  for (bases in list(unitBases2, unitBases3, unitBases4)) {
    scales <- priorScales(basisSampler(yBalanced, bases, family))
    variances <- c(1, rep(0.5, length(bases) - 1L))
    expect_equal(
      sqrt(sum(scales^2 * variances)),
      1.04912 * sqrt(length(bases) + 1) * s,
      tolerance = 1e-4
    )
  }
  # and a stated sd is per unit of the column as given: the same forest on a
  # column four times as large induces four times the prior on the index
  scaled <- priorScales(basisSampler(yBalanced, scaledBases, family))
  expect_equal(
    sqrt(sum(scaled^2 * c(1, 0.5) * c(4, 1)^2)),
    1.04912 * sqrt(2 * 16 + 1) * s,
    tolerance = 1e-4
  )

  # (d) THE PRODUCT ITSELF, at a DECLARED forest(sd = ). Each factor
  # discriminates: the multiple of s the forest states (sd = 2.5 s and 0.4 s,
  # in the response's units, so the leaf scale factor is the multiple) in
  # {2.5, 0.4}, s in {1, 1.813799} and the divisor 0.674, none 1, on columns
  # of 3s and 5s that a stated sd never reads.
  #
  # BOTH families, because probit ALONE CANNOT SEE A DROPPED ANCHOR - s = 1
  # there, so sd / 0.674 and sd s / 0.674 are the same number. Under logistic
  # they differ by 81 percent. Deleting the logistic arm leaves the anchor
  # unpinned, which is the near-miss this test guards against.
  declaredMap <- basisSampler(
    yBalanced,
    pinBases,
    family,
    forests = pinScalesFor(s)
  )
  expect_equal(
    priorScales(declaredMap),
    c(2.5, 0.4) * s / 0.674,
    tolerance = 1e-12
  )
  # (e) and the reader's sd is the number stated, its own scale that over the
  # half-normal median, with no half-Cauchy median on a forest with a basis
  for (f in seq_along(pinBases)) {
    reported <- declaredMap$getLeafPrior(f)
    expect_equal(unname(reported$sd), c(2.5, 0.4)[f] * s, tolerance = 1e-12)
    expect_equal(
      reported$k.scale * 0.674,
      unname(reported$sd),
      tolerance = 1e-12
    )
    expect_true(is.nan(engineSlot(declaredMap, f, "amplitude.prior.scale")))
  }
}

# --- positive builds on the shipped two-forest shape, where forest 1 carries
# the half-Cauchy amplitude and forest 2 the fixed-variance pair ---
twoForests <- list(
  dbarts::dbartsForests$forest(),
  dbarts::dbartsForests$forest(basis = ~ factor(z))
)
kForest <- function(y, family = "auto", ...) {
  dbarts(
    x,
    y,
    forests = twoForests,
    family = family,
    control = seededControlBcfFamily(),
    ...
  )
}
probit <- kForest(yBalanced)
logistic <- kForest(yBalanced, "logistic")
# a numeric 0/1 response resolves to probit with nothing said, which is the
# single-forest behavior this route now inherits
expect_equal(probit$model@family, "probit")
expect_equal(logistic$model@family, "logistic")
for (fit in list(probit, logistic)) {
  s <- anchors[[fit$model@family]]
  # forest 1 declares no basis, so its leaf scale stays at the anchor itself
  expect_equal(priorScales(fit), c(s, s / 0.674), tolerance = 1e-12)
  # a basis-free forest carries the half-Cauchy scale mixture and its median,
  # a basis forest none, and which a forest gets is the engine's, from what it
  # multiplies: the control carries no sd for either forest, so the numbers
  # the reader reports are the engine's own defaults.
  transported <- attr(fit$control, "bartcore.forests")
  expect_identical(transported$sd, c(NA_real_, NA_real_))
  free <- fit$getLeafPrior(1L)
  carried <- fit$getLeafPrior(2L)
  # and the median that reaches it is the FAMILY's own default, one anchor unit
  # rather than gaussian's two: under a pinned sigma nothing absorbs the
  # difference, and at 1 the induced index prior sits on the shipped
  # single-forest binary default's coverage rather than 1.6x outside it
  expect_equal(engineSlot(fit, 1L, "amplitude.prior.scale"), 1)
  expect_equal(free$sd, s)
  expect_identical(free$multiplier, "none")
  expect_identical(carried$multiplier, "factor")
  expect_true(is.nan(engineSlot(fit, 2L, "amplitude.prior.scale")))
  result <- fit$run(4L, 2L)
  expect_true(all(is.finite(result$train)))
  # sigma is pinned by the family, not merely fixed by a prior
  expect_true(all(result$sigma == 1))
  expect_equal(dim(fit$getForestAmplitudes()), c(3L, 1L))
  expect_true(all(is.finite(fit$getForestFits(2L))))
  # the latent channel opens with the family; it is empty under gaussian
  expect_equal(length(fit$getLatents()), n)
}
expect_equal(length(kForest(yContinuous)$getLatents()), 0L)

# forest = NULL stacks both per-forest readers with the forest margin between
# the rows and the chains, on this two-forest sampler; predictor rownames stay
# on margin 1 in both the single-forest and the stacked shape
fitsAll <- probit$getForestFits()
expect_equal(dim(fitsAll), c(n, 2L, 1L))
expect_equal(c(fitsAll[, 1L, ]), c(probit$getForestFits(1L)))
expect_equal(c(fitsAll[, 2L, ]), c(probit$getForestFits(2L)))

countsAll <- probit$getForestVariableCounts()
expect_equal(dim(countsAll), c(p, 2L, 1L))
expect_equal(
  unname(countsAll[, 1L, ]),
  unname(c(probit$getForestVariableCounts(1L)))
)
expect_equal(
  unname(countsAll[, 2L, ]),
  unname(c(probit$getForestVariableCounts(2L)))
)
expect_identical(rownames(countsAll), colnames(x))

# --- THE TWO DEFAULTS THEMSELVES, as the engine resolves them: the control
# carries no number for either. The GAUSSIAN arm is MANDATORY rather than
# decorative: under a latent family a basis-free forest's median and unit are
# the literal 1, so a latent-only fixture cannot see a confusion between the
# two channels. Only gaussian, where the median is 2, discriminates. ---
mapSampler <- function(y, bases, family, forests = NULL, ...) {
  allBasis <- everyForestHasBasis(bases)
  dbarts(
    dbartsData(x, y, bases = bases),
    control = if (allBasis) {
      uncountedControlBcfFamily()
    } else {
      seededControlBcfFamily()
    },
    family = family,
    forests = if (allBasis) {
      statedOnFirst(forests, unitOf(y, family))
    } else {
      forests
    },
    ...
  )
}
# one entry of the reader for each forest, NA where the forest has none
readerSlot <- function(sampler, name) {
  vapply(
    seq_along(sampler$data@bases),
    function(f) {
      value <- sampler$getLeafPrior(f)[[name]]
      if (is.null(value)) NA_real_ else unname(value)
    },
    numeric(1L)
  )
}
freeBases2 <- list(NULL, zFactor)
freeBases3 <- list(NULL, zFactor, gLevel)
freeBases4 <- list(NULL, zFactor, gLevel, hFactor)

# the basis-free channel's half-Cauchy median is the FAMILY's: 2 where the
# anchor is the response's own sd and sigma is DRAWN, 1 where the anchor is the
# link's error sd and sigma is PINNED
expect_equal(
  engineSlot(
    mapSampler(yContinuous, freeBases2, "gaussian"),
    1L,
    "amplitude.prior.scale"
  ),
  2
)
expect_equal(
  engineSlot(
    mapSampler(yBalanced, freeBases2, "probit"),
    1L,
    "amplitude.prior.scale"
  ),
  1
)
expect_equal(
  engineSlot(
    mapSampler(yBalanced, freeBases2, "logistic"),
    1L,
    "amplitude.prior.scale"
  ),
  1
)

# a factor's sd is one unit at K = 2, 3 and 4: a forest with a basis defaults
# to one unit whatever the number of forests (dec-B291), all-basis or beside a
# forest with none, and the default is a multiple of the response's scale
# reported as that multiple times the unit
for (family in c("gaussian", "probit", "logistic")) {
  response <- if (family == "gaussian") yContinuous else yBalanced
  for (bases in list(unitBases2, unitBases3, unitBases4)) {
    sampler <- mapSampler(response, bases, family)
    unit <- sampler$data@forest.calibration[["unit"]]
    expect_identical(
      readerSlot(sampler, "sd")[-1L],
      rep(unit, length(bases) - 1L)
    )
  }
  for (bases in list(freeBases2, freeBases3, freeBases4)) {
    sampler <- mapSampler(response, bases, family)
    unit <- sampler$data@forest.calibration[["unit"]]
    expect_identical(
      readerSlot(sampler, "sd")[-1L],
      rep(unit, length(bases) - 1L)
    )
  }
}

# a DECLARED sd keeps its per-forest reading at every K: the rule is a DEFAULT
# and not map algebra, so a caller who states the previous values (here the
# multiples of 1, the probit unit) gets the previous model back, bitwise, at any
# count
declaredAll <- mapSampler(
  yBalanced,
  unitBases3,
  "probit",
  forests = list(
    dbarts::dbartsForests$forest(sd = 1),
    dbarts::dbartsForests$forest(sd = 1),
    dbarts::dbartsForests$forest(sd = 1)
  )
)
expect_equal(readerSlot(declaredAll, "sd"), rep(1, 3L))
declaredShape <- mapSampler(
  yBalanced,
  freeBases3,
  "probit",
  forests = list(
    dbarts::dbartsForests$forest(sd = 2),
    dbarts::dbartsForests$forest(sd = 1),
    dbarts::dbartsForests$forest(sd = 1)
  )
)
expect_equal(readerSlot(declaredShape, "sd"), c(2, 1, 1))
expect_equal(engineSlot(declaredShape, 1L, "amplitude.prior.scale"), 2)

# --- the two adjacent guards that fire on the plainest binary two-forest call
# unless they are keyed to the family. Both are silent above, which is what
# makes those builds positive evidence; here is the other half, that an
# explicit non-default is still refused by name. ---
expect_true(is(probit$model@leaf.hyperprior, "dbartsFixedHyperprior"))
expect_equal(probit$model@leaf.hyperprior@k, 2)
expect_equal(probit$model@leaf.scale, 3.0)
expect_equal(logistic$model@leaf.scale, pi * sqrt(3.0))
for (family in c("probit", "logistic")) {
  # the prior vocabulary resolves inside the call, so these are stated here
  # rather than through the helper above
  expect_error(
    dbarts(
      x,
      yBalanced,
      forests = twoForests,
      family = family,
      leaf.prior = normal(chi(1.5, 2)),
      control = seededControlBcfFamily()
    ),
    "a 'k' hyperprior"
  )
  expect_error(
    dbarts(
      x,
      yBalanced,
      forests = twoForests,
      family = family,
      leaf.prior = normal(sd = 0.75),
      control = seededControlBcfFamily()
    ),
    "a model of several forests states each forest's size on the forest.*leave 'sd' out of 'leaf.prior'"
  )
  expect_error(
    dbarts(
      x,
      yBalanced,
      forests = twoForests,
      family = family,
      leaf.prior = normal(sd = invchi(1.5, 1)),
      control = seededControlBcfFamily()
    ),
    "a model of several forests states each forest's size on the forest"
  )
  # leaf.scale is written by the family switch rather than by the caller, so
  # the bridge's own backstop is where a non-default one can be stated at all
  spec <- dbartsSpec(
    dbartsData(x, yBalanced),
    seededControlBcfFamily(),
    forests = twoForests,
    family = family
  )
  spec$model@leaf.scale <- 0.9
  expect_error(
    new("dbartsSampler", spec$control, spec$model, spec$data),
    "a non-default leaf scale"
  )
}

# --- the family-keyed predicates a K-forest sampler newly answers for itself.
# Each was gaussian's answer while family_ was pinned. ---
binaryFits <- list(probit = probit, logistic = logistic)
for (name in names(binaryFits)) {
  fit <- binaryFits[[name]]
  expect_error(fit$setSigma(2), "fixes the residual standard deviation")
  # the weight conduit is family-keyed, not coupling-keyed: probit carries no
  # weights under any forest count, while a logistic sampler's are its
  # observation counts and a swap redraws the Polya-Gamma latents against them
  if (name == "probit") {
    expect_error(
      fit$setWeights(rep(2, n)),
      "probit models do not support case weights"
    )
  } else {
    latentsBefore <- fit$getLatents()
    fit$setWeights(rep(2, n))
    expect_false(identical(fit$getLatents(), latentsBefore))
    expect_error(fit$setWeights(rep(2.5, n)), "must be positive integers")
    rm(latentsBefore)
  }
  # the response-side conduit OPENS: the latents are refreshed against the
  # combined location, which is what used to make this unsafe
  fit$setResponse(1 - yBalanced)
  fit$setOffset(rep(0.1, n), updateScale = FALSE)
  expect_true(all(is.finite(fit$run(2L, 1L)$train)))
  expect_error(fit$setResponse(yContinuous), "must be coded 0 or 1")
  expect_error(
    fit$setOffset(rep(0.1, n), updateScale = TRUE),
    "does not support a sampler that carries forest amplitudes"
  )
}

# --- the whole binary weight policy arrives with the family and costs no new
# code: probit refuses weights outright, logistic reads them as counts ---
expect_error(
  kForest(yBalanced, "probit", weights = rep(2, n)),
  "weights"
)
expect_true(all(is.finite(
  kForest(yBalanced, "logistic", weights = rep(2, n))$run(2L, 1L)$train
)))

# --- the doors, on all three routes. The R surface refuses first, the two
# bridge entries are the backstops a direct consumer meets. ---
expect_error(
  dbarts(
    x,
    counts,
    forests = twoForests,
    family = "nbinom",
    control = seededControlBcfFamily()
  ),
  "does not support family \"nbinom\""
)
expect_error(
  dbarts(
    x,
    cbind(exp(yContinuous), 1),
    forests = twoForests,
    family = "aft",
    control = seededControlBcfFamily()
  ),
  "does not support family \"aft\""
)
# createHolder's backstop: a hand-built model states the family the R surface
# would have refused
doorSpec <- dbartsSpec(
  dbartsData(x, yContinuous),
  seededControlBcfFamily(),
  forests = twoForests
)
doorSpec$model@family <- "aft"
attr(doorSpec$control, "bartcore.survival") <- rep(1, n)
expect_error(
  new("dbartsSampler", doorSpec$control, doorSpec$model, doorSpec$data),
  "does not support an AFT"
)

rm(
  anchors,
  balanced,
  basisSampler,
  bases,
  binaryFits,
  counts,
  declaredAll,
  declaredMap,
  declaredShape,
  doorSpec,
  eta,
  everyForestHasBasis,
  family,
  fit,
  freeBases2,
  freeBases3,
  freeBases4,
  g,
  gLevel,
  hFactor,
  kForest,
  logistic,
  mapSampler,
  n,
  name,
  ones,
  p,
  pinBases,
  pinScalesFor,
  priorScales,
  probit,
  rare,
  engineSlot,
  readerSlot,
  reported,
  response,
  result,
  s,
  scaled,
  scaledBases,
  sampler,
  scales,
  seededControlBcfFamily,
  spec,
  statedOnFirst,
  unitOf,
  twoForests,
  uncountedControlBcfFamily,
  unit,
  unitBases2,
  unitBases3,
  unitBases4,
  variances,
  x,
  yBalanced,
  yContinuous,
  yRare,
  z,
  zFactor
)
