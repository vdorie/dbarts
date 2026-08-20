#!/usr/bin/env Rscript

# Exact-posterior gate for the hurdle / two-part family (docs/design/hurdle.md
# sections 0, 1 and 13). The companion reduction gate (hurdle-reduction.R)
# proves the wrapper's two component fits equal standalone bart2() calls DRAW
# FOR DRAW, but it builds those standalone calls from the package's OWN
# response splitter, so it is blind by construction to a defect in the split
# itself: a wrong split still reduces bitwise to a correct probit plus a
# correct gaussian on the wrong rows. This gate closes that hole - the split is
# re-derived here from the response alone, and the target is the composed
# posterior on the ORIGINAL non-negative scale.
#
# THE DERIVATION. The two parts are conditionally independent given the data
# (section 0), so the joint posterior FACTORIZES and each factor is derived
# separately from y:
#
#     occupancy   z_i = 1{y_i > 0} over all n, probit;
#     positive    w_i = log(y_i) over S = {i : y_i > 0}, gaussian.
#
# Nothing else in the response reaches either part: the zero rows carry no
# information about the positive part beyond their membership in the
# complement of S, and the positive VALUES carry none about occupancy. A
# threshold other than exactly zero, or a subset misaligned with its
# covariates, moves both factors at once - which is what arm (b) below shows.
#
# THE CONFIGURATION. One predictor with two cells and n.cuts = 1, so the root
# holds a single available cut and its children none: the single-tree space is
# exactly the two structures (a shared root leaf, prior 1 - base; one leaf per
# cell, prior base), as in negbin-exact.R. The occupancy factor is then a 1-D
# quadrature over the leaf log-probit; the positive factor, with sigma PINNED
# by resid.prior = fixed(), is normal-normal conjugate in CLOSED FORM on the
# engine's internal [-0.5, 0.5] rescaling of log y (aft-exact.R's convention:
# range * mu + shift returns the reported log-scale fit). Both cells hold zeros
# and positives, so both fits see the same two-cell cut grid.
#
# THE FUNCTIONALS, all three reported channels at both cells:
#
#     type = "prob"  E[pi(x)], the occupancy factor alone;
#     type = "bart"  E[f(x)], the positive part's log-scale fit - read at ALL
#                    n rows through the forced full-x test channel, so the
#                    zero rows the positive part never trained on are gated
#                    too (section 13 hardening c);
#     type = "ev"    E[pi exp(f + sigma^2/2)], the composed natural-scale
#                    mean. The wrapper forms it PER DRAW and averages; the two
#                    chains are independent by construction (independently
#                    derived seeds, section 13 hardening b), so its limit is
#                    E[pi] E[exp(f)] exp(sigma^2/2) - a product of the two
#                    exactly-derived factors times the lognormal retransform.
#
# Two deterministic preconditions ride along and fire before the distributional
# arms could: the pinned sigma must really be pinned, and the wrapper's own
# ingested occupancy response must equal the z derived above from y alone.
#
# On these data both factors' structure posteriors are effectively degenerate
# on the split (log Bayes factor about 20 for the occupancy factor and about
# 110 for the positive one), so the priorRoot / priorSplit weights and the two
# logMarginal terms are exercised but NOT discriminated - the gate would read
# the same with a fairly wide error in either. The root/split mixture itself is
# gated in negbin-exact.R, on the same two-structure space.
#
# Sensitivity, measured by poisoning the reference rather than the engine, at
# full size (against a clean gate reading 0.0002 prob / 0.0000 bart / 0.0002
# ev): (a) dropping the exp(sigma^2/2) retransform from the composed reference
# leaves prob and bart clean and reads max ev gap 0.3855; (b) deriving the
# split at the threshold y > 1 rather than y > 0 - a wrong split rule, which
# reclassifies small positives as zeros and drops them from S - reads max prob
# gap 0.1984, max bart gap 0.3228 and max ev gap 0.0982, and trips the
# deterministic split precondition besides. Arm (b) is the one the reduction
# gate cannot see. Tolerances below bound sampler MC error only; the ev
# channel's is wider because it multiplies two MC-noisy factors on a scale
# of about 3.
#
# Usage: Rscript hurdle-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

ndpost <- if (quick) 20000L else 60000L
nburn <- 5000L
nSeeds <- if (quick) 2L else 3L
tolProb <- if (quick) 0.003 else 0.0015 # vs single-chain MC error
tolBart <- if (quick) 0.003 else 0.0015
tolEv <- if (quick) 0.008 else 0.004

# ---- shipped constants ----

nCells <- 2L
base <- 0.5
power <- 2
k <- 2
probitScale <- 3.0 # probit node.scale (R/model.R's defaultNodeScale)
gaussianScale <- 0.5 # gaussian node.scale, in units of the response range
tauProbit <- probitScale / k # one tree, so no sqrt(n.trees) divisor
tauGaussian <- gaussianScale / k # on the internal [-0.5, 0.5] scale
sigmaFixed <- 0.5 # resid.prior = fixed(sigmaFixed^2) on the log scale

# ---- fixed data: two cells, both mixing zeros with positives ----

set.seed(1618L)
nPerCell <- 150L
n <- nCells * nPerCell
cellOf <- rep(seq_len(nCells), each = nPerCell)
cellValue <- c(0.25, 0.75) # one cut at the midpoint separates them
occupancyTrue <- c(0.45, 0.80)
logMeanTrue <- c(0.20, 1.20)
occupied <- runif(n) < occupancyTrue[cellOf]
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), logMeanTrue[cellOf[occupied]], 0.5))
x <- matrix(cellValue[cellOf], ncol = 1L)

# the split, re-derived here from y alone
z <- as.double(y > 0)
positive <- y > 0
logPositive <- log(y[positive])
stopifnot(all(tapply(z, cellOf, function(v) any(v == 1) && any(v == 0))))
cat(sprintf(
  "occupancy %s of %s | positives %d | log y range [%.2f, %.2f]\n",
  paste(as.vector(tapply(z, cellOf, sum)), collapse = " "),
  paste(as.vector(table(cellOf)), collapse = " "),
  sum(positive),
  min(logPositive),
  max(logPositive)
))

# the two structures a single cut admits: a shared root leaf, or one leaf per
# cell, whose children hold no remaining cut and so are forced leaves
priorRoot <- 1 - base
priorSplit <- base

# ---- factor 1: the occupancy posterior (probit, quadrature) ----

muGrid <- seq(-10, 10, by = 0.002)
dmu <- muGrid[2L] - muGrid[1L]
logMuPrior <- dnorm(muGrid, 0, tauProbit, log = TRUE)
logProbitP <- pnorm(muGrid, log.p = TRUE)
logProbit1mP <- pnorm(muGrid, lower.tail = FALSE, log.p = TRUE)

# marginal likelihood and posterior-mean probability for a probit leaf holding
# s successes out of m trials
occupancyLeaf <- function(s, m) {
  logIntegrand <- logMuPrior + s * logProbitP + (m - s) * logProbit1mP
  peak <- max(logIntegrand)
  w <- exp(logIntegrand - peak)
  list(
    logMarginal = log(sum(w) * dmu) + peak,
    meanProbability = sum(pnorm(muGrid) * w) / sum(w)
  )
}

exactOccupancy <- function(z) {
  successes <- as.vector(tapply(z, cellOf, sum))
  counts <- as.vector(table(cellOf))
  root <- occupancyLeaf(sum(successes), sum(counts))
  cells <- lapply(seq_len(nCells), function(a) {
    occupancyLeaf(successes[a], counts[a])
  })
  logSplit <- sum(vapply(cells, function(q) q$logMarginal, 0))
  scale <- max(root$logMarginal, logSplit)
  wRoot <- priorRoot * exp(root$logMarginal - scale)
  wSplit <- priorSplit * exp(logSplit - scale)
  vapply(
    seq_len(nCells),
    function(a) {
      (wRoot * root$meanProbability + wSplit * cells[[a]]$meanProbability) /
        (wRoot + wSplit)
    },
    0
  )
}

# ---- factor 2: the positive-part posterior (gaussian, conjugate) ----

# Posterior for one gaussian leaf on the engine's internal scale, sigma known:
# with prior mu ~ N(0, tau^2) the precision and mean are closed form, and so is
# the marginal likelihood that weights the two structures.
gaussianLeaf <- function(values, sigmaInternal) {
  m <- length(values)
  precision <- 1 / tauGaussian^2 + m / sigmaInternal^2
  mean <- (sum(values) / sigmaInternal^2) / precision
  list(
    logMarginal = -0.5 *
      m *
      log(2 * pi * sigmaInternal^2) -
      sum(values^2) / (2 * sigmaInternal^2) -
      0.5 * log(tauGaussian^2 * precision) +
      0.5 * mean^2 * precision,
    postMean = mean,
    precision = precision
  )
}

# Posterior-mean log-scale fit E[f] and natural-scale factor E[exp(f)] per
# cell. The engine rescales the working response to [-0.5, 0.5], so the
# reported fit is range * mu + shift and sigma enters divided by the range.
exactPositive <- function(positive, logPositive) {
  cells <- cellOf[positive]
  low <- min(logPositive)
  range <- max(logPositive) - low
  shift <- range * 0.5 + low
  internal <- (logPositive - low) / range - 0.5
  sigmaInternal <- sigmaFixed / range

  root <- gaussianLeaf(internal, sigmaInternal)
  perCell <- lapply(seq_len(nCells), function(a) {
    gaussianLeaf(internal[cells == a], sigmaInternal)
  })
  logSplit <- sum(vapply(perCell, function(q) q$logMarginal, 0))
  scale <- max(root$logMarginal, logSplit)
  wRoot <- priorRoot * exp(root$logMarginal - scale)
  wSplit <- priorSplit * exp(logSplit - scale)

  # E[exp(range mu)] under a normal leaf posterior is its lognormal mean
  fitOf <- function(q) range * q$postMean + shift
  expOf <- function(q) {
    exp(shift + range * q$postMean + range^2 / (2 * q$precision))
  }
  list(
    fit = vapply(
      seq_len(nCells),
      function(a) {
        (wRoot * fitOf(root) + wSplit * fitOf(perCell[[a]])) / (wRoot + wSplit)
      },
      0
    ),
    expFit = vapply(
      seq_len(nCells),
      function(a) {
        (wRoot * expOf(root) + wSplit * expOf(perCell[[a]])) / (wRoot + wSplit)
      },
      0
    )
  )
}

# ---- the composed reference ----

exactProb <- exactOccupancy(z)
exactPos <- exactPositive(positive, logPositive)
exactBart <- exactPos$fit
# the wrapper forms pi exp(f + sigma^2 / 2) per draw; the two chains are
# independent, so the draw average converges to the product of the factors
exactEv <- exactProb * exactPos$expFit * exp(sigmaFixed^2 / 2)

# ---- sampler fit ----

representative <- vapply(
  seq_len(nCells),
  function(a) which(cellOf == a)[1L],
  integer(1L)
)

fitSeed <- function(seed) {
  fit <- bart2(
    x,
    y,
    family = "hurdle.lognormal",
    n.trees = 1L,
    n.cuts = nCells - 1L,
    n.burn = nburn,
    n.samples = ndpost,
    n.chains = 1L,
    n.threads = 1L,
    k = k,
    power = power,
    base = base,
    resid.prior = fixed(sigmaFixed^2),
    verbose = FALSE,
    seed = seed
  )
  # deterministic preconditions: the pinned sigma really is pinned, and the
  # wrapper's own split matches the one derived above from y alone
  stopifnot(
    max(abs(fit$positive$sigma - sigmaFixed)) < 1e-10,
    identical(fit$occupancy$y, z),
    identical(length(fit$positive$y), sum(positive))
  )
  c(
    fitted(fit, type = "prob")[representative],
    fitted(fit, type = "bart")[representative],
    fitted(fit, type = "ev")[representative]
  )
}

fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
fitProb <- fit[seq_len(nCells)]
fitBart <- fit[nCells + seq_len(nCells)]
fitEv <- fit[2L * nCells + seq_len(nCells)]

gapProb <- max(abs(fitProb - exactProb))
gapBart <- max(abs(fitBart - exactBart))
gapEv <- max(abs(fitEv - exactEv))

report <- function(label, exact, sampler, gap, tolerance) {
  cat(sprintf(
    "  %-5s exact %s | sampler %s | max gap %.4f (tol %.4f)%s\n",
    label,
    paste(sprintf("%.4f", exact), collapse = " "),
    paste(sprintf("%.4f", sampler), collapse = " "),
    gap,
    tolerance,
    if (gap > tolerance) "  <- FAIL" else ""
  ))
}

cat("Hurdle exact-posterior gate (single tree, two cells, sigma pinned):\n")
report("prob", exactProb, fitProb, gapProb, tolProb)
report("bart", exactBart, fitBart, gapBart, tolBart)
report("ev", exactEv, fitEv, gapEv, tolEv)

if (gapProb > tolProb || gapBart > tolBart || gapEv > tolEv) {
  quit(status = 1L)
}
cat("\nOK: the hurdle fit matches the exact composed posterior\n")
