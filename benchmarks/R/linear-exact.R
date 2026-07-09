#!/usr/bin/env Rscript

# Exact-posterior gate for bartcore's linear leaves. A single-tree problem on
# one ordinal predictor with K = 4 distinct values (3 uniform cuts, each
# separating adjacent cells) admits brute-force enumeration of the tree space
# (15 contiguous-interval trees) while the linear leaf's ridge marginal is
# closed form, so the joint posterior over trees and leaf coefficients is
# computable exactly. The gate reproduces dbarts's range scaling, the
# fixed-sigma internal scale, the leaf standardization (training mean and
# sample sd), and the ridge prior tau * sigma^2 with tau = (k / scale)^2
# (model.hpp LinearGaussianLeaf::logIntegratedLikelihoodForNode), and checks
#   * the posterior leaf-count distribution P(#leaves = 1..4), and
#   * the posterior mean fit at each distinct predictor value,
# against long-run engine frequencies with batch-means z-scores. The
# leaf-count distribution is the sensitive statistic: a marginal-likelihood
# ridge error (e.g. dropping sigma^2 at model.hpp:304) reweights tree sizes
# by ~0.11 in P(1 leaf) here while moving E[fit] by only ~0.011. Failure
# means the linear leaf's marginal likelihood, posterior draw, or
# standardization is wrong.
#
# Usage: Rscript linear-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nKept <- if (quick) 50000L else 200000L
batchSize <- if (quick) 12500L else 25000L
nThin <- 10L
nBurn <- 2000L
engineSeed <- 20260709L
zBound <- 4

# ---- fixed design ----

set.seed(11L)
K <- 4L # distinct predictor values
nPer <- 10L # small n so the posterior spreads across tree sizes
muCell <- c(-0.9, -0.3, 0.3, 0.9)
noiseSd <- 0.6
cell <- rep(seq_len(K), each = nPer)
x <- as.double(cell)
n <- length(x)
y <- muCell[cell] + rnorm(n, sd = noiseSd)

base <- 0.8
power <- 2
kLeaf <- 2
nodeScale <- 0.5 # gaussian node.scale; scale = nodeScale / sqrt(ntree = 1)

# uniform cuts (n.cuts = K - 1) separate the cells exactly
cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
stopifnot(identical(findInterval(x, cuts) + 1L, cell))

# ---- range scaling, fixed sigma, and standardization, reproduced ----

yRange <- max(y) - min(y)
zScaled <- (y - min(y)) / yRange - 0.5
residVar <- (1 / yRange)^2 # resid.prior = fixed(1)
ridge <- (kLeaf / nodeScale)^2 * residVar
meanX <- mean(x)
sdX <- sd(x) # sample sd, matching standardizationMomentsForColumn
u <- (x - meanX) / sdX
uTest <- (as.double(seq_len(K)) - meanX) / sdX

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# ---- tree enumeration (contiguous-cell leaves over one ordinal column) ----

enumerate <- function(loCell, hiCell, loCut, hiCut, depth) {
  growth <- if (hiCut >= loCut) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(loCell, hiCell)),
    logPrior = log(1 - growth)
  ))
  if (hiCut < loCut) {
    return(result)
  }
  for (j in loCut:hiCut) {
    lefts <- enumerate(loCell, j, loCut, j - 1L, depth + 1L)
    rights <- enumerate(j + 1L, hiCell, j + 1L, hiCut, depth + 1L)
    for (left in lefts) {
      for (right in rights) {
        result[[length(result) + 1L]] <- list(
          leaves = c(left$leaves, right$leaves),
          logPrior = log(growth) -
            log(hiCut - loCut + 1) +
            left$logPrior +
            right$logPrior
        )
      }
    }
  }
  result
}

trees <- enumerate(1L, K, 1L, K - 1L, 0L)
cat(sprintf("enumerated %d trees\n", length(trees)))

# ---- exact posterior ----
#
# Per leaf with design U = [1, u] over its members (unit weights), the ridge
# marginal drops the same n- and w-only terms the engine drops:
#   0.5 p log(ridge) - 0.5 log det M - 0.5 (z'z - b' M^-1 b) / sigma^2,
#   M = U'U + ridge I, b = U'z, p = 2,
# and the conditional posterior mean of the coefficients is M^-1 b.

leafFit <- function(cells) {
  idx <- cell >= cells[1L] & cell <= cells[2L]
  design <- cbind(1, u[idx])
  z <- zScaled[idx]
  m <- crossprod(design) + ridge * diag(2L)
  b <- crossprod(design, z)
  ch <- chol(m)
  mInvB <- backsolve(ch, forwardsolve(t(ch), b))
  list(
    beta = mInvB,
    logML = 0.5 *
      2 *
      log(ridge) -
      sum(log(diag(ch))) -
      0.5 * (sum(z * z) - sum(b * mInvB)) / residVar
  )
}

logW <- numeric(length(trees))
numLeaves <- integer(length(trees))
fitScaled <- matrix(0, length(trees), K)
for (t in seq_along(trees)) {
  lw <- trees[[t]]$logPrior
  leafOfCell <- integer(K)
  betas <- vector("list", length(trees[[t]]$leaves))
  for (li in seq_along(trees[[t]]$leaves)) {
    cellRange <- trees[[t]]$leaves[[li]]
    lf <- leafFit(cellRange)
    lw <- lw + lf$logML
    betas[[li]] <- lf$beta
    leafOfCell[cellRange[1L]:cellRange[2L]] <- li
  }
  logW[t] <- lw
  numLeaves[t] <- length(trees[[t]]$leaves)
  for (p in seq_len(K)) {
    beta <- betas[[leafOfCell[p]]]
    fitScaled[t, p] <- beta[1L] + beta[2L] * uTest[p]
  }
}
w <- exp(logW - logSumExp(logW))
exactSizes <- vapply(seq_len(K), function(s) sum(w[numLeaves == s]), 0)
exactFit <- (colSums(w * fitScaled) + 0.5) * yRange + min(y)

# ---- engine arm ----

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  keepTrees = TRUE,
  n.samples = batchSize,
  n.burn = nBurn,
  n.thin = nThin,
  updateState = TRUE,
  rngSeed = engineSeed,
  n.cuts = K - 1L
)
sampler <- dbarts(
  matrix(x, ncol = 1L),
  y,
  test = matrix(as.double(seq_len(K)), ncol = 1L),
  control = ctl,
  tree.prior = cgm(power, base),
  node.prior = linear(1L, k = kLeaf),
  resid.prior = fixed(1)
)
stopifnot(is.null(sampler$data@offset))

nBatch <- as.integer(ceiling(nKept / batchSize))
engineLeaves <- integer(0L)
engineFits <- NULL
first <- TRUE
for (b in seq_len(nBatch)) {
  r <- if (first) sampler$run(nBurn, batchSize) else sampler$run(0L, batchSize)
  first <- FALSE
  tr <- sampler$getTrees()
  engineLeaves <- c(
    engineLeaves,
    as.integer(tapply(tr$var, tr$sample, function(v) sum(v == -1L)))
  )
  engineFits <- cbind(engineFits, r$test)
}
nDraws <- length(engineLeaves)
stopifnot(nDraws == ncol(engineFits))

# ---- batch-means MC error and z-scores ----

batchMeanSE <- function(v, nBatches = 400L) {
  len <- (length(v) %/% nBatches) * nBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
  sd(bm) / sqrt(nBatches)
}

cat(sprintf("\nengine draws: %d (thin %d)\n", nDraws, nThin))
cat(sprintf(
  "%-12s %10s %10s %10s %8s\n",
  "statistic",
  "engine",
  "exact",
  "MCse",
  "z"
))

anyFailure <- FALSE
report <- function(label, engineValue, exactValue, se) {
  z <- (engineValue - exactValue) / se
  failed <- abs(z) > zBound
  cat(sprintf(
    "%-12s %10.4f %10.4f %10.5f %+8.1f%s\n",
    label,
    engineValue,
    exactValue,
    se,
    z,
    if (failed) " <- FAIL" else ""
  ))
  failed
}

# leaf-count distribution (size 4 folded into size 3+ to keep counts high)
for (s in 1:2) {
  anyFailure <- report(
    sprintf("P(%d leaves)", s),
    mean(engineLeaves == s),
    exactSizes[s],
    batchMeanSE(engineLeaves == s)
  ) ||
    anyFailure
}
anyFailure <- report(
  "P(3+ leaves)",
  mean(engineLeaves >= 3L),
  sum(exactSizes[3:4]),
  batchMeanSE(engineLeaves >= 3L)
) ||
  anyFailure

# posterior mean fit at each distinct predictor value
for (p in seq_len(K)) {
  anyFailure <- report(
    sprintf("E[f(x=%d)]", p),
    mean(engineFits[p, ]),
    exactFit[p],
    batchMeanSE(engineFits[p, ])
  ) ||
    anyFailure
}

if (anyFailure) {
  cat("\nFAIL: linear-leaf sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("\nOK: linear-leaf sampler matches the exact posterior\n")
