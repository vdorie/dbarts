# samplePriorPredictive (issue #31): draws from the BART prior on a private
# freshly built sampler, so the caller's own sampler is left exactly as it
# was.

library(dbarts, quietly = TRUE)

set.seed(0)
n <- 30L
p <- 2L
x <- matrix(runif(n * p), n, p)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3)
z <- rbinom(n, 1L, pnorm(x[, 1L] - x[, 2L]))

control.plain <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L
)

sampler.gauss <- dbarts(y ~ x, control = control.plain, n.trees = 15L)
sampler.binary <- dbarts(z ~ x, control = control.plain, n.trees = 15L)
colnames(x) <- colnames(sampler.gauss$data@x)

x.test <- matrix(runif(6L * p), 6L, p)
colnames(x.test) <- colnames(sampler.gauss$data@x)

# (a) shapes: default x.test (training predictors) and an explicit x.test,
# both response families, both types
ev.gauss <- samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ev")
expect_equal(dim(ev.gauss), c(25L, n))
ppd.gauss <-
  samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ppd")
expect_equal(dim(ppd.gauss), c(25L, n))

ev.gauss.test <- samplePriorPredictive(
  sampler.gauss,
  x.test = x.test,
  n.samples = 10L,
  type = "ev"
)
expect_equal(dim(ev.gauss.test), c(10L, nrow(x.test)))

ev.binary <-
  samplePriorPredictive(sampler.binary, n.samples = 25L, type = "ev")
expect_equal(dim(ev.binary), c(25L, n))
expect_true(all(ev.binary >= 0 & ev.binary <= 1))

ppd.binary <-
  samplePriorPredictive(sampler.binary, n.samples = 25L, type = "ppd")
expect_equal(dim(ppd.binary), c(25L, n))
expect_true(all(ppd.binary %in% c(0, 1)))

ppd.binary.test <- samplePriorPredictive(
  sampler.binary,
  x.test = x.test,
  n.samples = 8L,
  type = "ppd"
)
expect_equal(dim(ppd.binary.test), c(8L, nrow(x.test)))

# (b) the caller's sampler is left untouched: fits off its live trees, and
# its residual sd, are bit-identical before and after
fitsBefore <- sampler.gauss$predict(x)
sigmasBefore <- sampler.gauss$getSigmas()
invisible(samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ppd"))
expect_identical(sampler.gauss$predict(x), fitsBefore)
expect_identical(sampler.gauss$getSigmas(), sigmasBefore)

# (c) seeding semantics. Successive calls are independent by default: each
# call constructs its private sampler fresh, seeding the engine RNG from
# R's stream when the control's seed is NA
indep1 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ev")
indep2 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ev")
expect_false(identical(indep1, indep2))

# set.seed() reproduces the whole draw, engine and observation layer alike
set.seed(11)
repeat1 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ppd")
set.seed(11)
repeat2 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ppd")
expect_identical(repeat1, repeat2)

set.seed(12)
repeatBinary1 <-
  samplePriorPredictive(sampler.binary, n.samples = 20L, type = "ppd")
set.seed(12)
repeatBinary2 <-
  samplePriorPredictive(sampler.binary, n.samples = 20L, type = "ppd")
expect_identical(repeatBinary1, repeatBinary2)

# a control-fixed seed pins the engine stream instead, so two ev calls
# repeat identically without any set.seed
control.pinned <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  seed = 99L
)
sampler.pinned <- dbarts(y ~ x, control = control.pinned, n.trees = 15L)
pinned1 <- samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ev")
pinned2 <- samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ev")
expect_identical(pinned1, pinned2)

# (d) statistical sanity: ppd adds observation noise on top of ev, so its
# aggregate variance is at least as large (gaussian). The pinned sampler
# pairs the forests across the two calls, making the comparison structural
# rather than luck across independent forest sets; binary ev draws staying
# valid probabilities is checked above
ppd.pinned <-
  samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ppd")
expect_true(var(as.vector(ppd.pinned)) >= var(as.vector(pinned1)))

# that noise is DRAWN from the sigma prior rather than held at sigest: with
# the forests pinned, the residual scale of draw s is the row sd of the
# difference, and those vary far more than the row estimates' own sampling
# error (which over 30 observations spans well under a factor of 2)
ev.spread <- samplePriorPredictive(sampler.pinned, n.samples = 60L, type = "ev")
ppd.spread <-
  samplePriorPredictive(sampler.pinned, n.samples = 60L, type = "ppd")
sigmaByDraw <- apply(ppd.spread - ev.spread, 1L, sd)
expect_true(max(sigmaByDraw) / min(sigmaByDraw) > 3)

# (e) validation errors
expect_error(samplePriorPredictive(1:5), pattern = "dbartsSampler")
expect_error(
  samplePriorPredictive(sampler.gauss, type = "bogus"),
  pattern = "should be one of"
)
# a fractional n.samples is refused, naming the argument, rather than
# silently truncated (coerceOrError's integer branch)
expect_error(
  samplePriorPredictive(sampler.gauss, n.samples = 5.5, type = "ev"),
  "'n.samples' must be a whole number; got '5.5'",
  fixed = TRUE
)

# (f) a heteroscedastic sampler has no scalar sigma to add as observation
# noise: "ppd" draws the variance forest from its prior too and adds
# s(x) eps at the predicted rows. "ev" never reaches that draw and is a
# mean-surface prior draw either way. The distribution of the added noise is
# pinned in test-heteroscedastic.R, beside the accessor it reads s(x) with.
sampler.variance <- dbarts(
  y ~ x,
  control = control.plain,
  variance = varianceForest(n.trees = 5L),
  n.trees = 15L
)
# a fractional varianceForest(n.trees) is refused by name too, not merely
# by "value" (R/spec.R hoists the if-expression coerceOrError otherwise
# sees into a named local so coercionSubject can quote it)
expect_error(
  dbarts(
    y ~ x,
    control = control.plain,
    variance = varianceForest(n.trees = 5.5),
    n.trees = 15L
  ),
  "'n.trees' must be a whole number; got '5.5'",
  fixed = TRUE
)
ppd.variance <- samplePriorPredictive(
  sampler.variance,
  n.samples = 5L,
  type = "ppd"
)
expect_equal(dim(ppd.variance), c(5L, n))
expect_true(all(is.finite(ppd.variance)))
ev.variance <- samplePriorPredictive(
  sampler.variance,
  n.samples = 5L,
  type = "ev"
)
expect_equal(dim(ev.variance), c(5L, n))
expect_true(all(is.finite(ev.variance)))

# --- the prior the sampler runs ---------------------------------------------
# Each sample draws k, a Student-t df, a negative-binomial shape, the ordinal
# cutpoints and DART's split weights from their priors, then the trees and
# leaves, and each family's observation layer is its own. The checks below
# are exact distributions: under constant leaves each tree puts one leaf at a
# point, so f(x) | k ~ N(prior.mean, (k.scale / k)^2) at any tree structure,
# and with k = scale chi_df, (f(x) - prior.mean) scale sqrt(df) / k.scale is
# t_df. Under a seeded control the "ev" and "ppd" calls share their forests,
# so ppd given ev is the observation layer alone. Every seed derives from
# seedBase; each KS and chi-squared test passes at p > 1e-3.
normal <- dbartsPriors$normal
invchi <- dbartsPriors$invchi
chi <- dbartsPriors$chi
fixed <- dbartsPriors$fixed
dart <- dbartsPriors$dart
seedBase <- 2611L
control.seeded <- dbartsControl(n.chains = 1L, n.threads = 1L, seed = seedBase)
pAbove <- 1e-3
ksP <- function(...) suppressWarnings(ks.test(...))$p.value
# a KS test of the values inside (lower, upper) against the CDF truncated
# there, for an "ev" whose latent reads back only inside the range a double
# probability resolves
truncatedKsP <- function(z, cdf, lower, upper) {
  z <- z[z > lower & z < upper]
  low <- cdf(lower)
  mass <- cdf(upper) - low
  list(n = length(z), p = ksP(z, function(v) (cdf(v) - low) / mass))
}
chiSquaredP <- function(observed, expected) {
  pchisq(
    sum((observed - expected)^2 / expected),
    length(observed) - 1L,
    lower.tail = FALSE
  )
}

set.seed(seedBase)
nPrior <- 60L
xPrior <- matrix(runif(nPrior * 2L), nPrior, 2L)
colnames(xPrior) <- c("a", "b")

# k drawn from its prior: the default probit chi(1.5, 2), 20 trees, and
# normal(sd = invchi(3, 0.5)), which is k ~ chi(3, k.scale / 0.5), the sd
# 0.5 / chi_3. qnorm(ev) recovers the latent where pnorm resolves it, so the
# test truncates at |latent| < 7
kScaleOf <- list(
  default = function(lp) 2,
  invchi = function(lp) lp$k.scale / 0.5
)
kForms <- list(
  default = list(prior = NULL, df = 1.5),
  invchi = list(prior = normal(sd = invchi(3, 0.5)), df = 3)
)
for (name in names(kForms)) {
  set.seed(seedBase + 1L)
  args <- list(
    xPrior,
    rbinom(nPrior, 1L, 0.5),
    family = "probit",
    n.trees = 20L,
    control = control.plain
  )
  if (!is.null(kForms[[name]]$prior)) {
    args$leaf.prior <- kForms[[name]]$prior
  }
  probit <- do.call(dbarts, args)
  lp <- probit$getLeafPrior()
  df <- kForms[[name]]$df
  multiplier <- kScaleOf[[name]](lp) * sqrt(df) / lp$k.scale
  ev <- samplePriorPredictive(probit, n.samples = 2200L, type = "ev")
  fit <- truncatedKsP(
    (qnorm(ev[, 1L]) - lp$prior.mean) * multiplier,
    function(v) pt(v, df),
    (-7 - lp$prior.mean) * multiplier,
    (7 - lp$prior.mean) * multiplier
  )
  expect_true(fit$n >= 2000L, info = paste("k draw", name))
  expect_true(fit$p > pAbove, info = paste("k draw", name))
}

# a fixed-k gaussian sampler draws what a hand loop of the two existing prior
# entries draws under one seed: the hyperparameter entry draws nothing
set.seed(seedBase + 2L)
handEv <- samplePriorPredictive(sampler.gauss, n.samples = 10L, type = "ev")
set.seed(seedBase + 2L)
handControl <- sampler.gauss$control
handControl@keepTrees <- FALSE
handControl@verbose <- FALSE
hand <- dbarts:::dbartsSampler$new(
  handControl,
  sampler.gauss$model,
  sampler.gauss$data
)
handX <- extract(hand, "predictors")
handLoop <- do.call(
  rbind,
  lapply(seq_len(10L), function(i) {
    hand$sampleTreesFromPrior(updateState = FALSE)
    hand$sampleLeafParametersFromPrior(updateState = FALSE)
    dbarts:::withCodedPredictors(hand$predict(handX, NULL, 1L))
  })
)
expect_identical(handEv, handLoop)

# Student-t: at a fixed df, (ppd - ev) / sigma is t_3; with df drawn, the
# mixture of t's over the engine's nine-point grid under its gamma(2, 0.1)
# weights
student <- dbartsFamilies$student
set.seed(seedBase + 3L)
heavy <- dbarts(
  xPrior,
  rnorm(nPrior),
  family = student(df = 3, sigma = fixed(1)),
  n.trees = 10L,
  control = control.seeded
)
set.seed(seedBase + 3L)
ev <- samplePriorPredictive(heavy, n.samples = 40L, type = "ev")
set.seed(seedBase + 3L)
ppd <- samplePriorPredictive(heavy, n.samples = 40L, type = "ppd")
expect_true(
  ksP(as.vector(ppd - ev) / heavy$getSigmas()[1L], "pt", df = 3) > pAbove
)
heavyDrawn <- dbarts(
  xPrior,
  rnorm(nPrior),
  family = student(sigma = fixed(1)),
  n.trees = 10L,
  control = control.seeded
)
dfGrid <- c(3, 4, 5, 6, 8, 10, 12, 15, 20)
dfWeights <- dfGrid * exp(-0.1 * dfGrid)
dfWeights <- dfWeights / sum(dfWeights)
set.seed(seedBase + 4L)
ev <- samplePriorPredictive(heavyDrawn, n.samples = 2000L, type = "ev")
set.seed(seedBase + 4L)
ppd <- samplePriorPredictive(heavyDrawn, n.samples = 2000L, type = "ppd")
expect_true(
  ksP(
    (ppd - ev)[, 1L] / heavyDrawn$getSigmas()[1L],
    function(v) {
      vapply(v, function(u) sum(dfWeights * pt(u, dfGrid)), numeric(1L))
    }
  ) >
    pAbove
)

# nbinom: "ev" is exp(predict), the center c in it, at the drawn k; "ppd" a
# randomized PIT against pnbinom at the fixed shape and, with the shape drawn,
# against the mixture over the engine's 13-point grid
nbinomFamily <- dbartsFamilies$nbinom
set.seed(seedBase + 5L)
counts <- rnbinom(nPrior, mu = 3, size = 2)
countDefault <- dbarts(
  xPrior,
  counts,
  family = nbinomFamily(),
  n.trees = 10L,
  control = control.seeded
)
lp <- countDefault$getLeafPrior()
multiplier <- 2 * sqrt(1.5) / lp$k.scale
ev <- samplePriorPredictive(countDefault, n.samples = 2100L, type = "ev")
expect_true(all(ev >= 0))
fit <- truncatedKsP(
  (log(ev[, 1L]) - lp$prior.mean) * multiplier,
  function(v) pt(v, 1.5),
  (-700 - lp$prior.mean) * multiplier,
  (700 - lp$prior.mean) * multiplier
)
expect_true(fit$n >= 2000L && fit$p > pAbove)

randomizedPit <- function(y, lowerCdf, density) {
  lowerCdf + runif(length(y)) * density
}
countFixed <- dbarts(
  xPrior,
  counts,
  family = nbinomFamily(shape = 2),
  n.trees = 10L,
  control = control.seeded,
  leaf.prior = normal(k = 2)
)
set.seed(seedBase + 6L)
ev <- samplePriorPredictive(countFixed, n.samples = 50L, type = "ev")
set.seed(seedBase + 6L)
ppd <- samplePriorPredictive(countFixed, n.samples = 50L, type = "ppd")
expect_true(all(ppd == round(ppd) & ppd >= 0))
pit <- randomizedPit(
  ppd,
  pnbinom(ppd - 1, mu = ev, size = 2),
  dnbinom(ppd, mu = ev, size = 2)
)
expect_true(ksP(as.vector(pit), "punif") > pAbove)

shapeGrid <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
shapeWeights <- shapeGrid * exp(-0.1 * shapeGrid)
shapeWeights <- shapeWeights / sum(shapeWeights)
countDrawn <- dbarts(
  xPrior,
  counts,
  family = nbinomFamily(),
  n.trees = 10L,
  control = control.seeded,
  leaf.prior = normal(k = 2)
)
set.seed(seedBase + 7L)
ev <- samplePriorPredictive(countDrawn, n.samples = 2000L, type = "ev")[, 1L]
set.seed(seedBase + 7L)
ppd <- samplePriorPredictive(countDrawn, n.samples = 2000L, type = "ppd")[, 1L]
mixture <- function(f, y) {
  rowSums(vapply(
    seq_along(shapeGrid),
    function(j) shapeWeights[j] * f(y, mu = ev, size = shapeGrid[j]),
    numeric(length(y))
  ))
}
pit <- randomizedPit(ppd, mixture(pnbinom, ppd - 1), mixture(dnbinom, ppd))
expect_true(ksP(pit, "punif") > pAbove)

# ordinal (K = 4, fixed k = 2): one x row at a grid of offsets, so that each
# cutpoint reads back at a row where its cumulative probability resolves. The
# latent, f = -qnorm(ev[, , 1]) less the offset, read at each sample's row
# nearest ev[, , 1] = 1 / 2, is N(prior.mean, (k.scale / 2)^2); each
# sample's cutpoints, gamma_s = qnorm(cumulative ev) + f + offset, agree
# across rows, with gamma_1 = 0 and the two log gaps N(0, 1.5^2); the
# categories drawn follow the probabilities
set.seed(seedBase + 8L)
ordinal <- dbarts(
  xPrior,
  factor(sample(1:4, nPrior, TRUE), ordered = TRUE),
  family = "ordinal",
  n.trees = 10L,
  control = control.seeded,
  leaf.prior = normal(k = 2)
)
lp <- ordinal$getLeafPrior()
offsets <- seq(-16, 2400, by = 8)
zeroRow <- which(offsets == 0)
xOffsets <- xPrior[rep(1L, length(offsets)), , drop = FALSE]
set.seed(seedBase + 9L)
ev <- samplePriorPredictive(
  ordinal,
  x.test = xOffsets,
  offset.test = offsets,
  n.samples = 2000L,
  type = "ev"
)
set.seed(seedBase + 9L)
ppd <- samplePriorPredictive(
  ordinal,
  x.test = xOffsets,
  offset.test = offsets,
  n.samples = 2000L,
  type = "ppd"
)
expect_identical(dim(ev), c(2000L, length(offsets), 4L))
centralRow <- apply(abs(ev[,, 1L] - 0.5), 1L, which.min)
f <- -qnorm(ev[cbind(seq_len(2000L), centralRow, 1L)]) - offsets[centralRow]
expect_true(
  ksP(f, "pnorm", mean = lp$prior.mean, sd = lp$k.scale / 2) > pAbove
)
offsetMatrix <- matrix(offsets, 2000L, length(offsets), byrow = TRUE)
cutpoints <- vapply(
  1:3,
  function(s) {
    cumulative <- apply(ev[,, seq_len(s), drop = FALSE], c(1L, 2L), sum)
    gamma <- qnorm(cumulative) + f + offsetMatrix
    gamma[!(cumulative > 1e-6 & cumulative < 1 - 1e-6)] <- NA
    spread <- apply(gamma, 1L, function(g) diff(range(g, na.rm = TRUE)))
    expect_true(all(is.finite(spread) & spread < 1e-8), info = s)
    rowMeans(gamma, na.rm = TRUE)
  },
  numeric(2000L)
)
expect_true(all(abs(cutpoints[, 1L]) < 1e-8))
expect_true(
  ksP(log(as.vector(apply(cutpoints, 1L, diff))), "pnorm", sd = 1.5) > pAbove
)
expect_true(is.integer(ppd) && all(ppd %in% 1:4))
expect_true(
  chiSquaredP(
    tabulate(ppd[, zeroRow], 4L),
    colSums(ev[, zeroRow, ])
  ) >
    pAbove
)

# a seeded ordinal "ev" repeats across calls and leaves R's stream alone: the
# cutpoints are the engine's draw
seedBefore <- .Random.seed
ordinal1 <- samplePriorPredictive(ordinal, n.samples = 5L, type = "ev")
ordinal2 <- samplePriorPredictive(ordinal, n.samples = 5L, type = "ev")
expect_identical(ordinal1, ordinal2)
expect_identical(.Random.seed, seedBefore)

# multinomial (K = 3, fixed k = 2): rows of "ev" sum to 1, the log ratio
# log(ev_k / ev_1) is N(0, 2 (k.scale / k)^2), an offset.test moves it by
# o_k - o_1 exactly, and the categories drawn follow the probabilities
set.seed(seedBase + 10L)
categories <- factor(sample(c("u", "v", "w"), nPrior, TRUE))
multinomial <- dbarts(
  xPrior,
  categories,
  family = "multinomial",
  n.trees = 10L,
  control = control.seeded,
  leaf.prior = normal(k = 2)
)
lp <- multinomial$getLeafPrior()[[1L]]
xTwo <- xPrior[1:2, , drop = FALSE]
set.seed(seedBase + 11L)
ev <- samplePriorPredictive(
  multinomial,
  x.test = xTwo,
  n.samples = 2000L,
  type = "ev"
)
expect_identical(dim(ev), c(2000L, 2L, 3L))
expect_true(all(abs(apply(ev, c(1L, 2L), sum) - 1) < 1e-12))
expect_true(
  ksP(
    log(ev[, 1L, 2L] / ev[, 1L, 1L]),
    "pnorm",
    sd = sqrt(2) * lp$k.scale / 2
  ) >
    pAbove
)
categoryOffset <- matrix(c(0, 0.5, -1.25, 0.25, 0, 2), 2L, 3L, byrow = TRUE)
set.seed(seedBase + 11L)
evOffset <- samplePriorPredictive(
  multinomial,
  x.test = xTwo,
  offset.test = categoryOffset,
  n.samples = 2000L,
  type = "ev"
)
for (k in 2:3) {
  shift <- log(evOffset[,, k] / evOffset[,, 1L]) - log(ev[,, k] / ev[,, 1L])
  expected <- matrix(
    categoryOffset[, k] - categoryOffset[, 1L],
    2000L,
    2L,
    byrow = TRUE
  )
  expect_true(max(abs(shift - expected)) < 1e-10, info = k)
}
set.seed(seedBase + 11L)
ppd <- samplePriorPredictive(
  multinomial,
  x.test = xTwo,
  n.samples = 2000L,
  type = "ppd"
)
expect_true(is.integer(ppd) && identical(dim(ppd), c(2000L, 2L)))
expect_true(chiSquaredP(tabulate(ppd[, 1L], 3L), colSums(ev[, 1L, ])) > pAbove)

# a count-row sampler draws a one-trial count vector per row
countRows <- dbarts(
  xPrior,
  t(rmultinom(nPrior, 2L, rep(1, 3) / 3)),
  family = "multinomial",
  n.trees = 10L,
  control = control.seeded
)
ppd <- samplePriorPredictive(
  countRows,
  x.test = xTwo,
  n.samples = 5L,
  type = "ppd"
)
expect_identical(dim(ppd), c(5L, 2L, 3L))
expect_true(all(apply(ppd, c(1L, 2L), sum) == 1))

# a sampler carrying a category offset needs offset.test, the all-zero matrix
# for the offset-free surface
offsetSampler <- dbarts(
  xPrior,
  categories,
  family = "multinomial",
  n.trees = 10L,
  control = control.seeded,
  offset = matrix(rnorm(nPrior * 3L), nPrior, 3L)
)
expect_error(
  samplePriorPredictive(offsetSampler, n.samples = 2L),
  "'offset.test' must be given on a sampler that carries a category offset",
  fixed = TRUE
)
expect_identical(
  dim(samplePriorPredictive(
    offsetSampler,
    x.test = xTwo,
    offset.test = matrix(0, 2L, 3L),
    n.samples = 2L
  )),
  c(2L, 2L, 3L)
)

# an improper leaf prior has no prior draws, by either spelling
for (prior in list(normal(sd = invchi(3, 0)), normal(k = chi(1.5, Inf)))) {
  improper <- dbarts(
    xPrior,
    rbinom(nPrior, 1L, 0.5),
    family = "probit",
    n.trees = 5L,
    control = control.plain,
    leaf.prior = prior
  )
  expect_error(
    samplePriorPredictive(improper, n.samples = 2L),
    paste0(
      "samplePriorPredictive does not support an improper leaf prior ",
      "(invchi() with scale 0, or chi() with an infinite scale): it has no ",
      "prior draws"
    ),
    fixed = TRUE
  )
}

# under dart() the split weights and the trees are one draw, rejected whole
# while any tree has an empty leaf, so the weights' prior is the stated
# Dirichlet times P(every tree held | s). In a small case that probability is
# exact: x1 binary, x2 on {0, 1, 2} (two equal cuts), four of the six cells
# occupied, 20 trees, alpha fixed at 2 so s_1 is uniform under the stated
# prior. The tree prior grows a node with probability 0.95 / (1 + depth)^2
# while a cut is left in its interval, the column drawn by s over the columns
# with one, the cut uniform; a tree is held when no leaf is empty.
cells <- rbind(c(0, 0), c(0, 1), c(0, 2), c(1, 2))
xCells <- cells[rep(1:4, each = 5L), ]
colnames(xCells) <- c("x1", "x2")
occupied <- function(lo1, hi1, lo2, hi2) {
  any(
    cells[, 1L] >= lo1 &
      cells[, 1L] <= hi1 &
      cells[, 2L] >= lo2 &
      cells[, 2L] <= hi2
  )
}
heldProbability <- function(s1, lo1, hi1, lo2, hi2, depth) {
  if (!occupied(lo1, hi1, lo2, hi2)) {
    return(0)
  }
  open1 <- hi1 > lo1
  open2 <- hi2 > lo2
  if (!open1 && !open2) {
    return(1)
  }
  growth <- 0.95 / (1 + depth)^2
  weights <- c(if (open1) s1 else 0, if (open2) 1 - s1 else 0)
  weights <- weights / sum(weights)
  grown <- 0
  if (open1) {
    grown <- grown +
      weights[1L] *
        heldProbability(s1, lo1, lo1, lo2, hi2, depth + 1) *
        heldProbability(s1, hi1, hi1, lo2, hi2, depth + 1)
  }
  if (open2) {
    split <- 0
    for (cut in lo2:(hi2 - 1L)) {
      split <- split +
        heldProbability(s1, lo1, hi1, lo2, cut, depth + 1) *
          heldProbability(s1, lo1, hi1, cut + 1L, hi2, depth + 1)
    }
    grown <- grown + weights[2L] * split / (hi2 - lo2)
  }
  (1 - growth) + growth * grown
}
grid <- seq(0.0005, 0.9995, by = 0.001)
density <- vapply(
  grid,
  function(s1) heldProbability(s1, 0, 1, 0, 2, 0),
  numeric(1L)
)^20
cdf <- c(0, cumsum(density / sum(density)))
dartOracle <- function(v) approx(c(0, grid + 0.0005), cdf, xout = v, rule = 2)$y
set.seed(seedBase + 12L)
dartSmall <- dbarts(
  xCells,
  rnorm(nrow(xCells)),
  n.trees = 20L,
  control = control.plain,
  split.probs = dart(alpha = 2, update.alpha = FALSE)
)
dartDraws <- lapply(seq_len(2000L), function(i) {
  .Call(
    dbarts:::C_dbarts_bartcore_sampleHyperparametersFromPrior,
    dartSmall$getPointer(),
    1000L
  )
})
expect_true(all(vapply(dartDraws, `[[`, logical(1L), "forests.drawn")))
dartFirst <- vapply(
  dartDraws,
  function(d) d$split.probabilities[1L],
  numeric(1L)
)
expect_true(ksP(dartFirst, dartOracle) > pAbove)
# the stated prior alone is far off here: the draw is not the per-tree retry
expect_true(ksP(dartFirst, "punif") < 1e-6)

# dart() moves the draws: under one seeded control, "ev" differs from the
# same sampler's at uniform split.probs
set.seed(seedBase + 13L)
yPrior <- rnorm(nPrior)
dartEv <- samplePriorPredictive(
  dbarts(
    xPrior,
    yPrior,
    n.trees = 10L,
    control = control.seeded,
    split.probs = dart()
  ),
  n.samples = 10L,
  type = "ev"
)
uniformEv <- samplePriorPredictive(
  dbarts(xPrior, yPrior, n.trees = 10L, control = control.seeded),
  n.samples = 10L,
  type = "ev"
)
expect_false(identical(dartEv, uniformEv))

# where too few joint draws are held (60 rows, 8 columns, 200 trees: about
# one in a hundred), the draws come from the sampler run with every row
# inactive, and a warning says so
set.seed(seedBase + 14L)
xWide <- matrix(runif(nPrior * 8L), nPrior, 8L)
wide <- dbarts(
  xWide,
  rnorm(nPrior),
  n.trees = 200L,
  control = control.plain,
  split.probs = dart()
)
expect_warning(
  fallback <- samplePriorPredictive(wide, n.samples = 5L, type = "ev"),
  "fewer than one in 50 independent draws"
)
expect_identical(dim(fallback), c(5L, nPrior))
expect_true(all(is.finite(fallback)))

rm(
  normal,
  invchi,
  chi,
  fixed,
  dart,
  seedBase,
  control.seeded,
  pAbove,
  ksP,
  truncatedKsP,
  chiSquaredP,
  nPrior,
  xPrior,
  kScaleOf,
  kForms,
  name,
  args,
  probit,
  df,
  lp,
  multiplier,
  ev,
  fit,
  handEv,
  handControl,
  hand,
  handX,
  handLoop,
  student,
  heavy,
  ppd,
  heavyDrawn,
  dfGrid,
  dfWeights,
  nbinomFamily,
  counts,
  countDefault,
  randomizedPit,
  countFixed,
  pit,
  shapeGrid,
  shapeWeights,
  countDrawn,
  mixture,
  ordinal,
  offsets,
  zeroRow,
  centralRow,
  xOffsets,
  f,
  offsetMatrix,
  cutpoints,
  seedBefore,
  ordinal1,
  ordinal2,
  categories,
  multinomial,
  xTwo,
  categoryOffset,
  evOffset,
  k,
  shift,
  expected,
  countRows,
  offsetSampler,
  prior,
  improper,
  cells,
  xCells,
  occupied,
  heldProbability,
  grid,
  density,
  cdf,
  dartOracle,
  dartSmall,
  dartDraws,
  dartFirst,
  yPrior,
  dartEv,
  uniformEv,
  xWide,
  wide,
  fallback
)

rm(
  sampler.variance,
  ppd.variance,
  ev.variance,
  n,
  p,
  x,
  y,
  z,
  control.plain,
  control.pinned,
  sampler.gauss,
  sampler.binary,
  sampler.pinned,
  x.test,
  ev.gauss,
  ppd.gauss,
  ev.gauss.test,
  ev.binary,
  ppd.binary,
  ppd.binary.test,
  fitsBefore,
  sigmasBefore,
  indep1,
  indep2,
  repeat1,
  repeat2,
  repeatBinary1,
  repeatBinary2,
  pinned1,
  pinned2,
  ppd.pinned
)
