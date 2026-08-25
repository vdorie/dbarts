# dbartsDrawLatents / dbartsWorkingResponse: the engine's per-sweep
# augmentation draws, run from R against R's own stream.
#
# THREE ORACLES, each gating what the others cannot see.
# (1) DISTRIBUTIONAL: the draw laws against CLOSED FORMS - truncated-normal
# moments at several truncation points and BOTH tails, Polya-Gamma mean and
# variance at several psi including the psi -> 0 limit - each a z on the mean
# at a stated Monte Carlo error. The only instrument that sees which SIDE a
# latent is truncated on or which SHAPE a Polya-Gamma draw carries; a shape
# defect is invisible at shape 1, so the weighted-logistic (b = 3) and nbinom
# (b = y + r = 5) cells are load-bearing.
# (2) AGREEMENT: a pure-R Gibbs loop built ONLY from the two helpers reproduces
# the engine's own posterior for that family. The only instrument that sees the
# WORKING RESPONSE, which a host can get wrong with perfectly drawn latents -
# the composition defect these functions exist to prevent.
# (3) OFFSET DISCRIMINATION: fit is the location WITHOUT the offset, so
# fit = f, offset = o and fit = f + o, offset = NULL draw the same latent and
# part company in the working response.
# Plus the RNG CONTRACT, BOTH halves: set.seed reproduces a call AND a call
# advances R's stream.

truncatedNormalMoments <- function(mu, sd, lower = -Inf, upper = Inf) {
  alpha <- (lower - mu) / sd
  beta <- (upper - mu) / sd
  z <- pnorm(beta) - pnorm(alpha)
  ratio <- (dnorm(alpha) - dnorm(beta)) / z
  # the alpha phi(alpha) terms vanish at an infinite bound
  aPhiA <- if (is.finite(alpha)) alpha * dnorm(alpha) else 0
  bPhiB <- if (is.finite(beta)) beta * dnorm(beta) else 0
  c(mean = mu + sd * ratio, var = sd^2 * (1 + (aPhiA - bPhiB) / z - ratio^2))
}

# PG(b, psi) has mean b tanh(psi / 2) / (2 psi) and variance
# b (sinh psi - psi) / (4 psi^3 cosh^2(psi / 2)), both 0/0 at psi = 0 where the
# limits are b/4 and b/24.
pgMoments <- function(b, psi) {
  if (psi == 0) {
    return(c(mean = b / 4, var = b / 24))
  }
  c(
    mean = b * tanh(psi / 2) / (2 * psi),
    var = b * (sinh(psi) - psi) / (4 * psi^3 * cosh(psi / 2)^2)
  )
}

# z on the mean of N iid draws at the closed-form variance: the stated Monte
# Carlo error is sqrt(var / N) and |z| < 4 is every moment cell's threshold
# (the observed maximum at this seed is 3.2).
momentZ <- function(draws, moments) {
  (mean(draws) - moments[["mean"]]) / sqrt(moments[["var"]] / length(draws))
}

# --- oracle 1: the draw laws against their closed forms

N <- 20000L
set.seed(5L)

# probit, BOTH TAILS at five truncation points: y = 1 truncates N(psi, 1) below
# at 0 and y = 0 above it, so a wrong-sided draw moves the mean by hundreds of
# Monte Carlo errors
for (psi in c(-1.5, -0.5, 0, 0.75, 2)) {
  for (yValue in c(1, 0)) {
    draws <- dbartsDrawLatents("probit", rep(psi, N), rep(yValue, N))
    moments <- if (yValue == 1) {
      truncatedNormalMoments(psi, 1, lower = 0)
    } else {
      truncatedNormalMoments(psi, 1, upper = 0)
    }
    expect_true(abs(momentZ(draws, moments)) < 4)
    expect_true(abs(var(draws) / moments[["var"]] - 1) < 0.05)
    expect_true(all(if (yValue == 1) draws > 0 else draws < 0))
  }
}

# aft imputes a censored log survival time from N(fit + offset, sigma^2)
# truncated BELOW at the observed log time, the one arm carrying a sigma
fitN <- rep(0.25, N)
for (sigma in c(0.5, 2)) {
  for (bound in c(-1, 0.5)) {
    boundN <- rep(bound, N)
    draws <- dbartsDrawLatents("aft", fitN, boundN, sigma = sigma)
    moments <- truncatedNormalMoments(0.25, sigma, lower = bound)
    expect_true(abs(momentZ(draws, moments)) < 4)
    expect_true(abs(var(draws) / moments[["var"]] - 1) < 0.05)
    expect_true(all(draws > bound))
  }
}

# Polya-Gamma at three shapes: 1 (unit-weight logistic), 3 (a logistic count
# weight, the sum of 3 independent PG(1, psi) draws) and 5 (nbinom's y + r,
# which is why that fixture is y = 2, r = 3 and not the shape-1 corner)
oneN <- rep(1, N)
twoN <- rep(2, N)
threeN <- rep(3, N)
for (psi in c(0, 0.5, 2, 5)) {
  psiN <- rep(psi, N)
  unit <- dbartsDrawLatents("logistic", psiN, oneN)
  weighted <- dbartsDrawLatents("logistic", psiN, oneN, weights = threeN)
  counts <- dbartsDrawLatents("nbinom", psiN, twoN, dispersion = 3)
  expect_true(abs(momentZ(unit, pgMoments(1, psi))) < 4)
  expect_true(abs(momentZ(weighted, pgMoments(3, psi))) < 4)
  expect_true(abs(momentZ(counts, pgMoments(5, psi))) < 4)
  expect_true(abs(var(unit) / pgMoments(1, psi)[["var"]] - 1) < 0.05)
  expect_true(abs(var(weighted) / pgMoments(3, psi)[["var"]] - 1) < 0.05)
  expect_true(abs(var(counts) / pgMoments(5, psi)[["var"]] - 1) < 0.05)
  expect_true(all(unit > 0) && all(weighted > 0) && all(counts > 0))
}

# Student-t's scale mixer is Gamma((df + 1) / 2, 2 / (df + r^2 / sigma^2)) at
# the residual r = y - fit - offset
studentFit <- rep(0.2, N)
for (df in c(3, 10)) {
  draws <- dbartsDrawLatents("student", studentFit, oneN, sigma = 1.5, df = df)
  shape <- 0.5 * (df + 1)
  scale <- 2 / (df + 0.8^2 / 1.5^2)
  moments <- c(mean = shape * scale, var = shape * scale^2)
  expect_true(abs(momentZ(draws, moments)) < 4)
}

# ordinal draws z inside its own category interval, the boundary categories
# one-sided; cutpoints are the K - 1 = 3 boundaries of K = 4 categories
cutpoints <- c(-0.5, 0.4, 1.1)
ordinalY <- rep(1:4, each = 50L)
ordinalFit <- rep(0.3, length(ordinalY))
ordinalDraws <- dbartsDrawLatents(
  "ordinal",
  ordinalFit,
  ordinalY,
  cutpoints = cutpoints
)
expect_true(all(
  ordinalDraws > c(-Inf, cutpoints)[ordinalY] &
    ordinalDraws <= c(cutpoints, Inf)[ordinalY]
))

# and the ordinal law itself against the SAME closed forms, since containment
# passes on any draw inside the interval: category 2 is a TWO-SIDED truncation
# of N(psi, 1) and category 4 a one-sided one, and psi enters both
ordinalPsi <- rep(0.3, N)
for (case in list(c(2, -0.5, 0.4), c(4, 1.1, Inf))) {
  draws <- dbartsDrawLatents(
    "ordinal",
    ordinalPsi,
    rep(case[1L], N),
    cutpoints = cutpoints
  )
  moments <- truncatedNormalMoments(0.3, 1, case[2L], case[3L])
  expect_true(abs(momentZ(draws, moments)) < 4)
  expect_true(abs(var(draws) / moments[["var"]] - 1) < 0.05)
}

# the reported quantity, per family, matching what $getLatents() reports
quantityOf <- function(family) {
  attr(
    dbartsDrawLatents(
      family,
      0.2,
      1,
      cutpoints = if (family == "ordinal") 0,
      dispersion = if (family == "nbinom") 2,
      df = if (family == "student") 5
    ),
    "quantity"
  )
}
expect_identical(
  vapply(c("probit", "ordinal", "aft"), quantityOf, ""),
  c(probit = "location", ordinal = "location", aft = "location")
)
expect_identical(
  vapply(c("logistic", "nbinom", "student"), quantityOf, ""),
  c(logistic = "precision", nbinom = "precision", student = "precision")
)

# --- oracle 2: a composed Gibbs loop reproduces the engine's own posterior

set.seed(20260815L)
nObs <- 150L
xAug <- matrix(runif(nObs * 3L), nObs, 3L)
fTrue <- 1.5 * xAug[, 1L] - 1.0 + 0.8 * xAug[, 2L]
augOffset <- rep_len(0.3, nObs)
nBurn <- 200L
nDraws <- 500L

augControl <- function() {
  dbartsControl(
    n.threads = 1L,
    n.trees = 50L,
    n.chains = 1L,
    updateState = FALSE,
    seed = 917L
  )
}

# The host is a GAUSSIAN sampler with sigma pinned at 1 and driven one sweep at
# a time - the model the augmentation reduces its family to - regressing on the
# working response the helpers build. Its leaf prior is restated to the native
# sampler's and k is pinned there (a binary sampler draws k from a hyperprior
# by default), so the two target the same posterior and this is not a prior
# comparison.
augAgreement <- function(family, y, coldStart, link) {
  native <- dbarts(
    xAug,
    y,
    offset = augOffset,
    family = family,
    control = augControl(),
    node.prior = normal(k = 2)
  )
  nativeFit <- rowMeans(native$run(nBurn, nDraws)$train) - augOffset
  host <- dbarts(
    xAug,
    coldStart - augOffset,
    control = augControl(),
    resid.prior = fixed(1)
  )
  host$setCalibration(prior.scale = native$getCalibration()[1L, "prior.scale"])
  fits <- matrix(0, nObs, nDraws)
  for (s in seq_len(nBurn + nDraws)) {
    latent <- dbartsDrawLatents(
      family,
      host$getFitsWithoutOffset(),
      y,
      offset = augOffset
    )
    working <- dbartsWorkingResponse(family, latent, y, offset = augOffset)
    if (family == "logistic") {
      host$setWeights(latent)
    }
    host$setResponse(working)
    host$run(0L, 1L)
    if (s > nBurn) fits[, s - nBurn] <- host$getFitsWithoutOffset()
  }
  composed <- rowMeans(fits)
  list(
    cor = cor(composed, nativeFit),
    rms = sqrt(mean((composed - nativeFit)^2)),
    probRms = sqrt(mean(
      (link(composed + augOffset) -
        link(nativeFit + augOffset))^2
    )),
    signal = cor(composed, fTrue),
    sigma = as.numeric(host$getSigmas()),
    scales = c(
      host$getCalibration()[1L, "prior.scale"],
      native$getCalibration()[1L, "prior.scale"]
    )
  )
}

# the tolerance: two independent chains of the same posterior, so agreement is
# stated as correlation and root-mean-square difference between the posterior
# means. Measured here: cor 0.997 / 0.996, latent rms 0.057 / 0.081, prob rms
# 0.016 / 0.012 (worst over three other seeds: 0.995, 0.095, 0.017).
yProbit <- rbinom(nObs, 1L, pnorm(fTrue + augOffset))
probitArm <- augAgreement("probit", yProbit, 2 * yProbit - 1, pnorm)
expect_equal(probitArm$scales[1L], probitArm$scales[2L])
expect_equal(probitArm$sigma, 1)
expect_true(probitArm$cor > 0.98)
expect_true(probitArm$rms < 0.2)
expect_true(probitArm$probRms < 0.05)
expect_true(probitArm$signal > 0.7)

# the logistic cold start is kappa / omega at omega = 1/4
yLogistic <- rbinom(nObs, 1L, plogis(2 * (fTrue + augOffset)))
logisticArm <- augAgreement(
  "logistic",
  yLogistic,
  4 * (yLogistic - 0.5),
  plogis
)
expect_equal(logisticArm$scales[1L], logisticArm$scales[2L])
expect_true(logisticArm$cor > 0.98)
expect_true(logisticArm$rms < 0.2)
expect_true(logisticArm$probRms < 0.05)
expect_true(logisticArm$signal > 0.7)

# --- oracle 3: the offset convention, which only the PAIR discriminates

set.seed(31L)
f0 <- rnorm(nObs)
o0 <- rnorm(nObs, sd = 0.4)
set.seed(101L)
splitZ <- dbartsDrawLatents("probit", f0, yProbit, offset = o0)
set.seed(101L)
foldedZ <- dbartsDrawLatents("probit", f0 + o0, yProbit)
splitW <- dbartsWorkingResponse("probit", splitZ, yProbit, offset = o0)
expect_identical(as.numeric(splitZ), as.numeric(foldedZ))
expect_identical(splitW, as.numeric(splitZ) - o0)
expect_false(isTRUE(all.equal(
  splitW,
  dbartsWorkingResponse("probit", foldedZ, yProbit)
)))

# and again on a precision family, where the working response is a quotient
# rather than a shift
set.seed(102L)
splitOmega <- dbartsDrawLatents("logistic", f0, yProbit, offset = o0)
set.seed(102L)
foldedOmega <- dbartsDrawLatents("logistic", f0 + o0, yProbit)
splitQ <- dbartsWorkingResponse("logistic", splitOmega, yProbit, offset = o0)
expect_identical(as.numeric(splitOmega), as.numeric(foldedOmega))
expect_equal(splitQ, (yProbit - 0.5) / as.numeric(splitOmega) - o0)
expect_false(isTRUE(all.equal(
  splitQ,
  dbartsWorkingResponse("logistic", foldedOmega, yProbit)
)))

# --- the RNG contract, BOTH halves: reproducibility alone would pass on a
# generator that never wrote R's state back, and stream advance alone on one
# that never read it

fitRng <- rep(0.4, 6L)
yRng <- c(1, 0, 1, 1, 0, 1)
set.seed(77L)
firstCall <- dbartsDrawLatents("probit", fitRng, yRng)
set.seed(77L)
expect_identical(dbartsDrawLatents("probit", fitRng, yRng), firstCall)

set.seed(77L)
withoutDraw <- runif(1L)
set.seed(77L)
invisible(dbartsDrawLatents("probit", fitRng, yRng))
expect_false(isTRUE(all.equal(runif(1L), withoutDraw)))

# and a call composes with any other R draw: same seed, same interleaving,
# same everything
set.seed(78L)
interleavedNormal <- rnorm(3L)
interleavedLatent <- dbartsDrawLatents("probit", fitRng, yRng)
set.seed(78L)
expect_identical(rnorm(3L), interleavedNormal)
expect_identical(dbartsDrawLatents("probit", fitRng, yRng), interleavedLatent)

# --- the by-name refusals: each optional argument belongs to a family, and its
# own family requires it

expect_error(
  dbartsDrawLatents("probit", fitRng, yRng, weights = rep(2, 6L)),
  "'weights' applies only to family \"logistic\""
)
expect_error(
  dbartsDrawLatents("aft", fitRng, yRng, dispersion = 2),
  "'dispersion' applies only to family \"nbinom\""
)
expect_error(
  dbartsDrawLatents("probit", fitRng, yRng, cutpoints = 0),
  "'cutpoints' applies only to family \"ordinal\""
)
expect_error(
  dbartsDrawLatents("logistic", fitRng, yRng, df = 5),
  "'df' applies only to family \"student\""
)
expect_error(
  dbartsDrawLatents("probit", fitRng, yRng, sigma = 2),
  "'sigma' applies only to family \"aft\" and \"student\""
)
# sigma's own formal default (NULL) does not refuse when passed explicitly -
# it reproduces the omitted call's draws bitwise, unlike a real value
set.seed(9L)
withDefault <- dbartsDrawLatents("probit", fitRng, yRng, sigma = NULL)
set.seed(9L)
omitted <- dbartsDrawLatents("probit", fitRng, yRng)
expect_identical(withDefault, omitted)
expect_error(
  dbartsWorkingResponse("probit", fitRng, yRng, weights = rep(2, 6L)),
  "'weights' applies only to family \"logistic\""
)
expect_error(
  dbartsDrawLatents("nbinom", fitRng, rep(2, 6L)),
  "family \"nbinom\" requires 'dispersion'"
)
expect_error(
  dbartsDrawLatents("ordinal", fitRng, rep(1, 6L)),
  "family \"ordinal\" requires 'cutpoints'"
)
expect_error(
  dbartsDrawLatents("student", fitRng, fitRng),
  "family \"student\" requires 'df'"
)
expect_error(
  dbartsDrawLatents("gaussian", fitRng, fitRng),
  "'arg' should be one of"
)

# the response support is the sampler's own rule, stated by the same function
# every conduit that swaps a y calls
expect_error(
  dbartsDrawLatents("probit", fitRng, rnorm(6L)),
  "must be coded 0 or 1"
)
expect_error(
  dbartsDrawLatents("nbinom", fitRng, rep(-1, 6L), dispersion = 2),
  "non-negative integer"
)
expect_error(
  dbartsDrawLatents("ordinal", fitRng, rep(4, 6L), cutpoints = c(0, 1)),
  "integer category index in \\[1, 3\\]"
)

# shapes, counts, and the precision a working response divides by
expect_error(
  dbartsDrawLatents("probit", fitRng, yRng[-1L]),
  "'y' must have length 6, that of 'fit'"
)
expect_error(
  dbartsDrawLatents("logistic", fitRng, yRng, weights = rep(1.5, 6L)),
  "positive whole numbers"
)
expect_error(
  dbartsDrawLatents("nbinom", fitRng, rep(2, 6L), dispersion = 1.5),
  "'dispersion' must be a whole number"
)
expect_error(
  dbartsDrawLatents("ordinal", fitRng, rep(1, 6L), cutpoints = c(1, 0)),
  "'cutpoints' must be strictly increasing"
)
expect_error(
  dbartsWorkingResponse("logistic", rep(0, 6L), yRng),
  "it must be positive"
)
expect_error(
  dbartsDrawLatents("probit", c(1, NA), c(1, 0)),
  "'fit' must be a finite"
)
