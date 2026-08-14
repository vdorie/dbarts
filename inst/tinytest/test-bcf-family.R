# The K-forest coupling under a non-gaussian response (docs/design/bcf.md,
# docs/design/multiplier-combiner.md): probit and logistic are built, the wider
# families are doors that refuse by name. What is pinned here is (i) the
# ALL-BASIS shape, in which every forest is fixed-variance and nothing in the
# sampler can absorb a mis-stated prior, (ii) the CALIBRATION ANCHOR itself,
# which is the only number in the slice that a running, mixing sampler cannot
# reveal to be wrong, (iii) the guards that fire on the plainest possible
# binary two-forest call unless they are keyed to the family, and (iv) the
# family-keyed refusals that flip the moment the sampler reports its own
# family. The anchor assertions are exact rather than statistical: under a
# latent family the response transform is the identity and the map's sqrt(m)
# cancels, so $getCalibration()'s prior.scale IS the map's node scale.

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

seededControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    n.samples = 4L,
    updateState = FALSE,
    rngSeed = 29L,
    ...
  )
}

# every entry non-NULL, so both/all forests take the fixed-variance channel:
# the shape in which the anchor is load-bearing, since no half-Cauchy
# amplitude and no drawn sigma stands between the prior and the index
ones <- matrix(1, n, 1L)
zBasis <- cbind(1 - z, z)
gBasis <- unname(model.matrix(~ g - 1L))
unitBases2 <- list(ones, zBasis)
unitBases3 <- list(ones, zBasis, gBasis)
# and one cell whose rows are NOT unit norm, which is what pins the map's
# per-unit-of-row-norm divisor rather than assuming it
scaledBases <- list(ones, 4 * zBasis)

basisSampler <- function(y, bases, family = "auto", ...) {
  dbarts(
    dbartsData(x, y, bases = bases),
    control = seededControl(),
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
      seededControl()
    )$control,
    "bartcore.bcf"
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
  # the default amplitude prior variance and unit row norms, which the map
  # disperses as 1.04912 sqrt(K)
  for (bases in list(unitBases2, unitBases3)) {
    scales <- priorScales(basisSampler(yBalanced, bases, family))
    expect_equal(
      sqrt(sum(scales^2 * 0.5)),
      1.04912 * sqrt(length(bases)) * s,
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
    control = seededControl(),
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
      control = seededControl()
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
      control = seededControl()
    ),
    "a named 'prior.scale'"
  )
  # node.scale is written by the family switch rather than by the caller, so
  # the bridge's own backstop is where a non-default one can be stated at all
  spec <- dbartsSpec(
    dbartsData(x, yBalanced),
    seededControl(),
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
for (fit in list(probit, logistic)) {
  expect_error(fit$setSigma(2), "fixes the residual standard deviation")
  expect_error(
    fit$setWeights(rep(2, n)),
    "cannot be set after creation|do not support case weights"
  )
  # the response-side conduit OPENS: the latents are refreshed against the
  # combined location, which is what used to make this unsafe
  fit$setResponse(1 - yBalanced)
  fit$setOffset(rep(0.1, n), updateScale = FALSE)
  expect_true(all(is.finite(fit$run(2L, 1L)$train)))
  expect_error(fit$setResponse(yContinuous), "must be coded 0 or 1")
  expect_error(
    fit$setOffset(rep(0.1, n), updateScale = TRUE),
    "does not support a BCF sampler"
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
    control = seededControl()
  ),
  "does not support family \"nbinom\""
)
expect_error(
  dbarts(
    x,
    cbind(exp(yContinuous), 1),
    forests = twoForests,
    family = "aft",
    control = seededControl()
  ),
  "does not support family \"aft\""
)
# createHolder's backstop: a hand-built model states the family the R surface
# would have refused
doorSpec <- dbartsSpec(
  dbartsData(x, yContinuous),
  seededControl(),
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
host <- dbarts(x, yBalanced, control = seededControl())
internalProbit <- dbarts:::bartcoreBCFSampler(host, z)
expect_equal(
  unname(dbarts:::bartcoreForestCalibration(internalProbit, 0L)[
    1L,
    "prior.scale"
  ]),
  1.0,
  tolerance = 1e-12
)
internalLogistic <- dbarts:::bartcoreBCFSampler(host, z, family = "logistic")
expect_equal(
  unname(dbarts:::bartcoreForestCalibration(internalLogistic, 0L)[
    1L,
    "prior.scale"
  ]),
  pi / sqrt(3.0),
  tolerance = 1e-12
)
expect_error(
  dbarts:::bartcoreBCFSampler(
    dbarts(x, yContinuous, control = seededControl()),
    z,
    family = "aft"
  ),
  "BCF does not support an AFT"
)

rm(
  anchors,
  balanced,
  basisSampler,
  bases,
  counts,
  doorSpec,
  eta,
  family,
  fit,
  forestParams,
  g,
  gBasis,
  host,
  internalLogistic,
  internalProbit,
  kForest,
  logistic,
  medianRowNorm,
  n,
  ones,
  p,
  params,
  priorScales,
  probit,
  rare,
  result,
  s,
  scaled,
  scaledBases,
  scales,
  seededControl,
  spec,
  twoForests,
  unitBases2,
  unitBases3,
  x,
  yBalanced,
  yContinuous,
  yRare,
  z,
  zBasis
)
