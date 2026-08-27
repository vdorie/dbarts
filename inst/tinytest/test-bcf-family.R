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
# cancels, so $getCalibration()'s prior.scale IS the map's node scale.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

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

seededControlBcfFamily <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    n.samples = 4L,
    updateState = FALSE,
    seed = 29L,
    ...
  )
}

# every entry non-NULL, so both/all forests take the fixed-variance channel:
# the shape in which the anchor is load-bearing, since no half-Cauchy
# amplitude and no drawn sigma stands between the prior and the index
ones <- matrix(1, n, 1L)
zBasis <- cbind(1 - z, z)
gBasis <- unname(model.matrix(~ g - 1L))
hBasis <- unname(model.matrix(~ factor(x[, 1L] > 0.5) - 1L))
unitBases2 <- list(ones, zBasis)
unitBases3 <- list(ones, zBasis, gBasis)
unitBases4 <- list(ones, zBasis, gBasis, hBasis)
# and one cell whose rows are NOT unit norm, which is what pins the map's
# per-unit-of-row-norm divisor rather than assuming it. LOAD-BEARING: delete
# the 4 and the dropped-divisor mutation goes green.
scaledBases <- list(ones, 4 * zBasis)
# the product pin's own fixture, whose every factor discriminates; see (d)
pinBases <- list(3 * ones, 5 * zBasis)
pinScales <- list(forest(sd = 2.5), forest(sd = 0.4))

basisSampler <- function(y, bases, family = "auto", ...) {
  dbarts(
    dbartsData(x, y, bases = bases),
    control = seededControlBcfFamily(),
    family = family,
    ...
  )
}
priorScales <- function(sampler) {
  unname(vapply(
    seq_along(sampler$data@bases),
    function(f) sampler$getCalibration(f)[1L, "prior.scale"],
    numeric(1L)
  ))
}
medianRowNorm <- function(basis) {
  norms <- sqrt(rowSums(basis^2))
  median(norms[norms > 0])
}

# --- the all-basis fixture itself, gaussian, K = 2 and K = 3. Both halves of
# what makes it the honest test of an anchor: every forest resolves to the
# fixed-variance channel (divisor the half-normal median, no half-Cauchy
# amplitude scale), and it runs. ---
for (bases in list(unitBases2, unitBases3, scaledBases)) {
  params <- attr(
    dbartsSpec(
      dbartsData(x, yContinuous, bases = bases),
      seededControlBcfFamily()
    )$control,
    "bartcore.forests"
  )$params
  expect_equal(length(params), length(bases))
  for (forestParams in params) {
    expect_equal(forestParams[5L], 0.674)
    expect_equal(forestParams[7L], 0)
  }
  fit <- basisSampler(yContinuous, bases)
  expect_equal(length(fit$data@bases), length(bases))
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
  # literal and the row-norm divisor recomputed here rather than read back
  expect_equal(
    priorScales(balanced),
    1.0 * s / (0.674 * vapply(scaledBases, medianRowNorm, numeric(1L))),
    tolerance = 1e-12
  )

  # (c) THE INDUCED INDEX: sqrt(sum_f prior.scale_f^2 v_f ||B_f(i,.)||^2) at
  # the default amplitude prior variance and unit row norms. The map itself
  # still disperses as 1.04912 sqrt(K) - it carries no per-K renormalization -
  # but the DEFAULT node scale factor is now sqrt(2/K),
  # whose product with that sqrt(K) is sqrt(2) at every count. So the assertion
  # is a CONSTANT rather than a function of K, which is strictly stronger and is
  # the design statement itself: the all-basis index prior is 1.4837 s whether
  # the caller writes the same model as two forests or as four. Three counts,
  # because the claim is invariance in K rather than a value at one K, and K = 2
  # is the fixed point every shipped configuration sits on.
  for (bases in list(unitBases2, unitBases3, unitBases4)) {
    scales <- priorScales(basisSampler(yBalanced, bases, family))
    expect_equal(
      sqrt(sum(scales^2 * 0.5)),
      1.04912 * sqrt(2) * s,
      tolerance = 1e-4
    )
  }
  # and that index is INVARIANT to a rescaled basis, which is the whole point
  # of dividing the row norm out: the same two forests, one of them at four
  # times the norm, induce the same prior on the index
  scaled <- priorScales(basisSampler(yBalanced, scaledBases, family))
  expect_equal(
    sqrt(sum(scaled^2 * 0.5 * c(1, 4)^2)),
    1.04912 * sqrt(2) * s,
    tolerance = 1e-4
  )

  # (d) THE PRODUCT ITSELF, at a DECLARED forest(sd = ). Every assertion above
  # runs at the default node scale factor, the literal 1, so three of the map's
  # four factors are pinned and the fourth is vacated - a recorded
  # hazard (unit values silently vacate pins) sitting on the expression whose
  # default is about to move. Here each of the four discriminates: sd in
  # {2.5, 0.4}, s in {1, 1.813799}, c in {3, 5}, divisor 0.674, no two equal
  # and none 1. The oracle, re-derived: probit 1.236399604352 and
  # 0.118694362018, logistic 2.242580816313 and 0.215287758366.
  #
  # BOTH families, because probit ALONE CANNOT SEE A DROPPED ANCHOR - s = 1
  # there, so sd / (0.674 c) and sd s / (0.674 c) are the same number. Under
  # logistic they differ by 81 percent. Deleting the logistic arm leaves the
  # anchor unpinned, which is the near-miss this test guards against.
  expect_equal(
    priorScales(basisSampler(yBalanced, pinBases, family, forests = pinScales)),
    c(2.5, 0.4) * s / (0.674 * vapply(pinBases, medianRowNorm, numeric(1L))),
    tolerance = 1e-12
  )

  # (e) THE FACTORS THEMSELVES, which the reader now reports beside the
  # product. Each column against the value this fixture declares, and then the
  # ANCHOR recovered by the identity the Rd states - prior.scale * divisor *
  # row norm / factor - which is the only route to s and is what the
  # decomposition exists for. The recovery is what the logistic arm
  # discriminates: a map that dropped the anchor recovers 1 here rather than
  # 1.813799, while probit cannot tell the two apart.
  declaredMap <- basisSampler(
    yBalanced,
    pinBases,
    family,
    forests = pinScales
  )
  for (f in seq_along(pinBases)) {
    reported <- declaredMap$getCalibration(f)[1L, ]
    expect_equal(unname(reported["node.scale.factor"]), c(2.5, 0.4)[f])
    expect_equal(unname(reported["node.scale.divisor"]), 0.674)
    expect_equal(
      unname(reported["basis.row.norm"]),
      medianRowNorm(pinBases[[f]])
    )
    expect_equal(
      unname(
        reported["prior.scale"] *
          reported["node.scale.divisor"] *
          reported["basis.row.norm"] /
          reported["node.scale.factor"]
      ),
      s,
      tolerance = 1e-12
    )
    # every forest here carries a basis, so every one takes the fixed-variance
    # channel and reports no half-Cauchy median
    expect_equal(unname(reported["amplitude.prior.variance"]), 0.5)
    expect_true(is.nan(reported["amplitude.prior.scale"]))
  }
}

# --- positive builds on the shipped two-forest shape, where forest 1 carries
# the half-Cauchy amplitude and forest 2 the fixed-variance pair ---
twoForests <- list(forest(), forest(basis = ~ factor(z)))
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
  # forest 1 declares no basis, so its node scale stays at the anchor itself
  expect_equal(priorScales(fit), c(s, s / 0.674), tolerance = 1e-12)
  # and the reader SHOWS that rather than leaving it to be inferred: factor,
  # divisor and row norm are each 1 on that forest. The two amplitude columns
  # are EXCLUSIVE - a basis-free forest carries the half-Cauchy scale mixture
  # and reports its median, a basis forest the fixed variance - and which
  # spelling a forest gets agrees with the transported per-forest params, so
  # the reader and the creation route cannot disagree about the channel.
  fitParams <- attr(fit$control, "bartcore.forests")$params
  free <- fit$getCalibration(1L)[1L, ]
  carried <- fit$getCalibration(2L)[1L, ]
  expect_equal(
    unname(free[c(
      "node.scale.factor",
      "node.scale.divisor",
      "basis.row.norm"
    )]),
    c(1, 1, 1)
  )
  expect_equal(unname(free["amplitude.prior.scale"]), fitParams[[1L]][7L])
  # and the median that reaches it is the FAMILY's own default, one anchor unit
  # rather than gaussian's two: under a pinned sigma nothing absorbs the
  # difference, and at 1 the induced index prior sits on the shipped
  # single-forest binary default's coverage rather than 1.6x outside it
  expect_equal(unname(free["amplitude.prior.scale"]), 1)
  expect_true(is.nan(free["amplitude.prior.variance"]))
  expect_equal(unname(carried["amplitude.prior.variance"]), fitParams[[2L]][6L])
  expect_true(is.nan(carried["amplitude.prior.scale"]))
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

# --- THE TWO DEFAULTS THEMSELVES, in their own transport slots. Everything
# above reads them only after the calibration map has folded factor, anchor,
# divisor and row norm into one number, so a confusion between the two channels
# or between the families would have to survive as a coincidence of products.
#
# The GAUSSIAN arm is MANDATORY rather than decorative: under a latent family a
# basis-free forest's vector is now c(50, 0.25, 3, 1, 1, 1, 1, 1) - positions 4
# through 8 all the literal 1 - so a latent-only fixture cannot see a transport
# confusion among those slots at all. Only gaussian, where slot 7 is 2,
# discriminates position. That is the same "unit values silently vacate
# pins" hazard arriving through the new default introduced here. ---
transportParams <- function(y, bases, family, ...) {
  attr(
    dbartsSpec(
      dbartsData(x, y, bases = bases),
      seededControlBcfFamily(),
      family = family,
      ...
    )$control,
    "bartcore.forests"
  )$params
}
paramSlot <- function(params, index) {
  vapply(params, function(vector) vector[index], numeric(1L))
}
freeBases2 <- list(NULL, zBasis)
freeBases3 <- list(NULL, zBasis, gBasis)
freeBases4 <- list(NULL, zBasis, gBasis, hBasis)

# slot 7, the basis-free channel's half-Cauchy median, is the FAMILY's: 2 where
# the anchor is the response's own sd and sigma is DRAWN, 1 where the anchor is
# the link's error sd and sigma is PINNED
expect_equal(
  paramSlot(transportParams(yContinuous, freeBases2, "gaussian"), 7L)[1L],
  2
)
expect_equal(
  paramSlot(transportParams(yBalanced, freeBases2, "probit"), 7L)[1L],
  1
)
expect_equal(
  paramSlot(transportParams(yBalanced, freeBases2, "logistic"), 7L)[1L],
  1
)

# slot 4, the fixed-variance channel's node scale factor, is K-AWARE - 1,
# 0.816497, 0.707107 at K = 2, 3, 4 - on the RESOLVED forest count and NOT on
# the count of basis forests: the shipped shape's K - 1 basis forests take
# sqrt(2/K), so a law normalized on the basis count would read 1 at K = 3 here
# rather than 0.816497, and would move the all-basis K = 2 arm the Option L
# anchor was argued on. The basis-free forest keeps the literal 1 in that slot -
# the law touches the withBasis branch only, a Cauchy channel having no finite
# variance to enter the budget with.
for (family in c("gaussian", "probit", "logistic")) {
  response <- if (family == "gaussian") yContinuous else yBalanced
  for (bases in list(unitBases2, unitBases3, unitBases4)) {
    expect_equal(
      paramSlot(transportParams(response, bases, family), 4L),
      rep(sqrt(2 / length(bases)), length(bases))
    )
  }
  for (bases in list(freeBases2, freeBases3, freeBases4)) {
    factors <- paramSlot(transportParams(response, bases, family), 4L)
    expect_equal(factors[1L], 1)
    expect_equal(
      factors[-1L],
      rep(sqrt(2 / length(bases)), length(bases) - 1L)
    )
  }
}

# a DECLARED sd keeps its per-forest reading at every K: the law is a DEFAULT
# and not map algebra, so a caller who states the previous values gets the
# previous model back, bitwise, at any count
declaredAll <- transportParams(
  yBalanced,
  unitBases3,
  "probit",
  forests = list(forest(sd = 1), forest(sd = 1), forest(sd = 1))
)
expect_equal(paramSlot(declaredAll, 4L), rep(1, 3L))
declaredShape <- transportParams(
  yBalanced,
  freeBases3,
  "probit",
  forests = list(forest(sd = 2), forest(sd = 1), forest(sd = 1))
)
expect_equal(paramSlot(declaredShape, 4L), c(1, 1, 1))
expect_equal(paramSlot(declaredShape, 7L)[1L], 2)

# --- the two adjacent guards that fire on the plainest binary two-forest call
# unless they are keyed to the family. Both are silent above, which is what
# makes those builds positive evidence; here is the other half, that an
# explicit non-default is still refused by name. ---
expect_true(is(probit$model@node.hyperprior, "dbartsFixedHyperprior"))
expect_equal(probit$model@node.hyperprior@k, 2)
expect_equal(probit$model@node.scale, 3.0)
expect_equal(logistic$model@node.scale, pi * sqrt(3.0))
for (family in c("probit", "logistic")) {
  # the prior vocabulary resolves inside the call, so these are stated here
  # rather than through the helper above
  expect_error(
    dbarts(
      x,
      yBalanced,
      forests = twoForests,
      family = family,
      node.prior = normal(chi(1.5, 2)),
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
      node.prior = normal(2, scale = 1.5),
      control = seededControlBcfFamily()
    ),
    "a named 'prior.scale'"
  )
  # node.scale is written by the family switch rather than by the caller, so
  # the bridge's own backstop is where a non-default one can be stated at all
  spec <- dbartsSpec(
    dbartsData(x, yBalanced),
    seededControlBcfFamily(),
    forests = twoForests,
    family = family
  )
  spec$model@node.scale <- 0.9
  expect_error(
    new("dbartsSampler", spec$control, spec$model, spec$data),
    "a non-default node scale"
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

# --- the internal creation route derives its family from the model it is
# handed, so the gate harness can construct one; a supplied family = names it
# at the call site by writing that same slot ---
host <- dbarts(x, yBalanced, control = seededControlBcfFamily())
internalProbit <- dbarts:::bartcoreBCFSampler(host, z)
expect_equal(
  unname(bartcoreForestCalibration(internalProbit, 0L)[
    1L,
    "prior.scale"
  ]),
  1.0,
  tolerance = 1e-12
)
internalLogistic <- dbarts:::bartcoreBCFSampler(host, z, family = "logistic")
expect_equal(
  unname(bartcoreForestCalibration(internalLogistic, 0L)[
    1L,
    "prior.scale"
  ]),
  pi / sqrt(3.0),
  tolerance = 1e-12
)
expect_error(
  dbarts:::bartcoreBCFSampler(
    dbarts(x, yContinuous, control = seededControlBcfFamily()),
    z,
    family = "aft"
  ),
  "a treatment forest does not support an AFT"
)

# and it takes the same FAMILY-AWARE default the public route does, which is
# what the creation oracle's draw-for-draw comparison needs: a route that
# resolved sd.control to a literal would build a model the public surface
# cannot express. Both halves - the map's own product and the half-Cauchy
# median beside it - read off each route's own reader.
internalRow <- bartcoreForestCalibration(internalProbit, 0L)[1L, ]
expect_equal(unname(internalRow["amplitude.prior.scale"]), 1)
expect_equal(
  unname(internalRow[c("prior.scale", "amplitude.prior.scale")]),
  unname(probit$getCalibration(1L)[
    1L,
    c("prior.scale", "amplitude.prior.scale")
  ])
)
# gaussian keeps 2 on this route as on the public one, which is why the
# gate harness (bcf-equivalence, gaussian throughout) re-records nothing
internalGaussian <- dbarts:::bartcoreBCFSampler(
  dbarts(x, yContinuous, control = seededControlBcfFamily()),
  z
)
expect_equal(
  unname(
    bartcoreForestCalibration(internalGaussian, 0L)[
      1L,
      "amplitude.prior.scale"
    ]
  ),
  2
)

rm(
  anchors,
  balanced,
  basisSampler,
  bases,
  binaryFits,
  counts,
  declaredAll,
  declaredShape,
  doorSpec,
  eta,
  factors,
  family,
  fit,
  forestParams,
  freeBases2,
  freeBases3,
  freeBases4,
  g,
  gBasis,
  hBasis,
  host,
  internalGaussian,
  internalLogistic,
  internalProbit,
  internalRow,
  kForest,
  logistic,
  medianRowNorm,
  n,
  name,
  ones,
  p,
  paramSlot,
  params,
  pinBases,
  pinScales,
  priorScales,
  probit,
  rare,
  response,
  result,
  s,
  scaled,
  scaledBases,
  scales,
  seededControlBcfFamily,
  spec,
  transportParams,
  twoForests,
  unitBases2,
  unitBases3,
  unitBases4,
  x,
  yBalanced,
  yContinuous,
  yRare,
  z,
  zBasis
)
