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
#     zero        z_i = 1{y_i > 0} over all n, probit;
#     positive    w_i = log(y_i) over S = {i : y_i > 0}, gaussian.
#
# Nothing else in the response reaches either part: the zero rows carry no
# information about the positive part beyond their membership in the
# complement of S, and the positive VALUES carry none about the zero part. A
# threshold other than exactly zero, or a subset misaligned with its
# covariates, moves both factors at once - which is what arm (b) below shows.
#
# THE CONFIGURATION. One predictor with two cells and n.cuts = 1, so the root
# holds a single available cut and its children none: the single-tree space is
# exactly the two structures (a shared root leaf, prior 1 - base; one leaf per
# cell, prior base), as in negbin-exact.R. The zero factor is then a 1-D
# quadrature over the leaf log-probit; the positive factor, with sigma PINNED
# by a fixed sigma prior, is normal-normal conjugate in CLOSED FORM on the
# engine's internal [-0.5, 0.5] rescaling of log y (aft-exact.R's convention:
# range * mu + shift returns the reported log-scale fit). Both cells hold zeros
# and positives, so both fits see the same two-cell cut grid.
#
# THE FUNCTIONALS, all three reported channels at both cells:
#
#     type = "prob"  E[pi(x)], the zero factor alone;
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
# ingested zero response must equal the z derived above from y alone.
#
# On these data both factors' structure posteriors are effectively degenerate
# on the split (log Bayes factor about 20 for the zero factor and about
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
# The single-tree fits below run swap at 0.1 rather than the shipped zero:
# swap is the only proposal that rotates a child's rule up the tree, and with
# one tree the change move cannot re-root once the splits beneath the root
# depend on the root's own variable, so nothing else crosses between
# rootings. The kernel's correctness at the shipped mixture is the balance
# gates' job, not this one's.
#
# THE WEIGHT AND OFFSET ARMS. Two further arms run the same configuration with
# the inputs a hurdle fit takes beyond a response, the expected posterior
# derived from the response, the weights and the offsets alone:
#
#     weights   0 rows in each cell, among the zeros and among the positives;
#               a row of weight 0 is out of both likelihoods and still in the
#               design, so the zero factor sums over the weight-1 rows and the
#               positive factor over the weight-1 positive rows;
#     offsets   a different offset on every row for each part, the zero part's
#               added to the probit latent and the positive part's to log y,
#               so the positive working response is log y less its offset and
#               the channels carry the offset of the row they are read at.
#
# Each factor is the same derivation with its rows chosen by the weights and
# its linear predictor shifted by the offsets; the engine's internal rescaling
# of the working response is read over the weight-1 positive rows, the rows in
# the likelihood, as the engine does.
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
probitScale <- 3.0 # probit leaf.scale (R/model.R's defaultLeafScale)
gaussianScale <- 0.5 # gaussian leaf.scale, in units of the response range
tauProbit <- probitScale / k # one tree, so no sqrt(n.trees) divisor
tauGaussian <- gaussianScale / k # on the internal [-0.5, 0.5] scale
sigmaFixed <- 0.5 # sigma = fixed(sigmaFixed^2) on the log scale

# ---- fixed data: two cells, both mixing zeros with positives ----

set.seed(1618L)
nPerCell <- 150L
n <- nCells * nPerCell
cellOf <- rep(seq_len(nCells), each = nPerCell)
cellValue <- c(0.25, 0.75) # one cut at the midpoint separates them
zeroTrue <- c(0.45, 0.80)
logMeanTrue <- c(0.20, 1.20)
occupied <- runif(n) < zeroTrue[cellOf]
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), logMeanTrue[cellOf[occupied]], 0.5))
x <- matrix(cellValue[cellOf], ncol = 1L)

# the split, re-derived here from y alone
z <- as.double(y > 0)
positive <- y > 0
logPositive <- log(y[positive])
stopifnot(all(tapply(z, cellOf, function(v) any(v == 1) && any(v == 0))))
cat(sprintf(
  "zero-part z=1 %s of %s | positives %d | log y range [%.2f, %.2f]\n",
  paste(as.vector(tapply(z, cellOf, sum)), collapse = " "),
  paste(as.vector(table(cellOf)), collapse = " "),
  sum(positive),
  min(logPositive),
  max(logPositive)
))

# weights of 0 and 1: in each cell, 0 rows among the zeros and among the
# positives, in blocks that are not the first rows of the response
weightsArm <- rep(1, n)
for (a in seq_len(nCells)) {
  inCell <- which(cellOf == a)
  weightsArm[tail(inCell[y[inCell] > 0], 25L)] <- 0
  weightsArm[head(inCell[y[inCell] == 0], 15L)] <- 0
}
# a different offset on every row for each part
offsetZero <- round(runif(n, -0.6, 0.6), 3L)
offsetPositive <- round(runif(n, -0.5, 0.5), 3L)

# the two structures a single cut admits: a shared root leaf, or one leaf per
# cell, whose children hold no remaining cut and so are forced leaves
priorRoot <- 1 - base
priorSplit <- base

# ---- factor 1: the zero posterior (probit, quadrature) ----

muGrid <- seq(-10, 10, by = 0.002)
dmu <- muGrid[2L] - muGrid[1L]
logMuPrior <- dnorm(muGrid, 0, tauProbit, log = TRUE)

# The log-likelihood of the weight-1 rows over the grid of leaf values, row i
# contributing z_i log Phi(mu + o_i) + (1 - z_i) log(1 - Phi(mu + o_i)); a leaf
# sums its rows. The mean probability read at a row is E[Phi(mu + o)] over the
# leaf posterior.
zeroLogLik <- function(z, weights, offset) {
  rows <- which(weights == 1)
  logLik <- matrix(0, n, length(muGrid))
  for (i in rows) {
    shifted <- muGrid + offset[i]
    logLik[i, ] <- if (z[i] == 1) {
      pnorm(shifted, log.p = TRUE)
    } else {
      pnorm(shifted, lower.tail = FALSE, log.p = TRUE)
    }
  }
  logLik
}

zeroLeaf <- function(logLik, rows, readOffset) {
  logIntegrand <- logMuPrior + colSums(logLik[rows, , drop = FALSE])
  peak <- max(logIntegrand)
  w <- exp(logIntegrand - peak)
  list(
    logMarginal = log(sum(w) * dmu) + peak,
    meanProbability = vapply(
      readOffset,
      function(o) sum(pnorm(muGrid + o) * w) / sum(w),
      0
    )
  )
}

# E[pi] at each cell's reading row, the structures weighted by their marginals
exactZero <- function(z, weights, offset, readRows) {
  logLik <- zeroLogLik(z, weights, offset)
  root <- zeroLeaf(logLik, seq_len(n), offset[readRows])
  cells <- lapply(seq_len(nCells), function(a) {
    zeroLeaf(logLik, which(cellOf == a), offset[readRows])
  })
  logSplit <- sum(vapply(cells, function(q) q$logMarginal, 0))
  scale <- max(root$logMarginal, logSplit)
  wRoot <- priorRoot * exp(root$logMarginal - scale)
  wSplit <- priorSplit * exp(logSplit - scale)
  vapply(
    seq_len(nCells),
    function(a) {
      (wRoot *
        root$meanProbability[a] +
        wSplit * cells[[a]]$meanProbability[a]) /
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
# cell. The working response is log y less the offset; the engine rescales it
# to [-0.5, 0.5] over the weight-1 positive rows, so the reported fit is
# range * mu + shift (plus the offset, added by the caller) and sigma enters
# divided by the range. Only the weight-1 rows reach the likelihood.
exactPositive <- function(positive, logPositive, weights, offset) {
  cells <- cellOf[positive]
  working <- logPositive - offset[positive]
  included <- weights[positive] == 1
  low <- min(working[included])
  range <- max(working[included]) - low
  shift <- range * 0.5 + low
  internal <- (working - low) / range - 0.5
  sigmaInternal <- sigmaFixed / range

  root <- gaussianLeaf(internal[included], sigmaInternal)
  perCell <- lapply(seq_len(nCells), function(a) {
    gaussianLeaf(internal[included & cells == a], sigmaInternal)
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

# ---- the composed reference and the sampler fit, per arm ----

fitSeed <- function(seed, arm, readRows) {
  weightsGiven <- if (arm$weighted) arm$weights
  offsetGiven <- if (arm$offsets) {
    list(zero = arm$offsetZero, positive = arm$offsetPositive)
  }
  fit <- bart2(
    x,
    y,
    n.trees = 1L,
    n.cuts = nCells - 1L,
    n.burn = nburn,
    n.samples = ndpost,
    n.chains = 1L,
    n.threads = 1L,
    k = k,
    power = power,
    base = base,
    proposal.probs = c(
      birth_death = 0.5,
      swap = 0.1,
      change = 0.4,
      birth = 0.5
    ),
    family = hurdle.lognormal(sigma = fixed(sigmaFixed^2)),
    weights = weightsGiven,
    offset = offsetGiven,
    verbose = FALSE,
    seed = seed
  )
  # deterministic preconditions: the pinned sigma really is pinned (a fixed
  # sigma stores no draws, and its stored value is rescaled, so to rounding),
  # and the wrapper's own split matches the one derived above from y alone
  stopifnot(
    abs(fit$positive$fixed$sigma - sigmaFixed) < 1e-10,
    identical(fit$zero$y, z),
    identical(length(fit$positive$y), sum(positive))
  )
  c(
    fitted(fit, type = "prob")[readRows],
    fitted(fit, type = "bart")[readRows],
    fitted(fit, type = "ev")[readRows]
  )
}

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

arms <- list(
  list(
    name = "no weights, no offsets",
    weighted = FALSE,
    offsets = FALSE,
    weights = rep(1, n),
    offsetZero = rep(0, n),
    offsetPositive = rep(0, n)
  ),
  list(
    name = "0 rows in each cell",
    weighted = TRUE,
    offsets = FALSE,
    weights = weightsArm,
    offsetZero = rep(0, n),
    offsetPositive = rep(0, n)
  ),
  list(
    name = "an offset in each part",
    weighted = FALSE,
    offsets = TRUE,
    weights = rep(1, n),
    offsetZero = offsetZero,
    offsetPositive = offsetPositive
  )
)

failed <- FALSE
cat("Hurdle exact-posterior gate (single tree, two cells, sigma pinned):\n")
for (arm in arms) {
  # each cell is read at its first weight-1 row
  readRows <- vapply(
    seq_len(nCells),
    function(a) which(cellOf == a & arm$weights == 1)[1L],
    integer(1L)
  )
  exactProb <- exactZero(z, arm$weights, arm$offsetZero, readRows)
  exactPos <- exactPositive(
    positive,
    logPositive,
    arm$weights,
    arm$offsetPositive
  )
  exactBart <- exactPos$fit + arm$offsetPositive[readRows]
  # the wrapper forms pi exp(f + o + sigma^2 / 2) per draw; the two chains are
  # independent, so the draw average converges to the product of the factors
  exactEv <- exactProb *
    exactPos$expFit *
    exp(arm$offsetPositive[readRows] + sigmaFixed^2 / 2)

  fit <- colMeans(do.call(
    rbind,
    lapply(seq_len(nSeeds), fitSeed, arm = arm, readRows = readRows)
  ))
  fitProb <- fit[seq_len(nCells)]
  fitBart <- fit[nCells + seq_len(nCells)]
  fitEv <- fit[2L * nCells + seq_len(nCells)]

  gapProb <- max(abs(fitProb - exactProb))
  gapBart <- max(abs(fitBart - exactBart))
  gapEv <- max(abs(fitEv - exactEv))

  cat(sprintf("arm: %s\n", arm$name))
  report("prob", exactProb, fitProb, gapProb, tolProb)
  report("bart", exactBart, fitBart, gapBart, tolBart)
  report("ev", exactEv, fitEv, gapEv, tolEv)
  if (gapProb > tolProb || gapBart > tolBart || gapEv > tolEv) {
    failed <- TRUE
  }
}

if (failed) {
  quit(status = 1L)
}
cat("\nOK: the hurdle fit matches the exact composed posterior\n")
