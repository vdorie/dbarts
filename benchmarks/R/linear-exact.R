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
# A second arm adds 200 rows at weight 0, every one at the first value, so
# the cells and the tree space are unchanged and the likelihood is the first
# arm's: its exact posterior is the first arm's only when the leaf is
# standardized over the rows in the likelihood. Standardized over every row
# (mean 1.25, sd 0.72 against 2.5 and 1.13) P(1 leaf) moves by ~0.095. The
# arm does not see the gp leaf's lengthscale heuristic over those rows, a
# mask, the refusal of rows in at one value or the keep of a re-derivation;
# tests/cpp's testLeafStandardizationRows and test-leaf-standardization-rows.R
# hold those to literals.
#
# The single-tree fits below run swap at 0.1 rather than the shipped zero:
# swap is the only proposal that rotates a child's rule up the tree, and with
# one tree the change move cannot re-root once the splits beneath the root
# depend on the root's own variable, so nothing else crosses between
# rootings. The kernel's correctness at the shipped mixture is the balance
# gates' job, not this one's.
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
nodeScale <- 0.5 # gaussian leaf.scale; scale = nodeScale / sqrt(ntree = 1)

# uniform cuts (n.cuts = K - 1) separate the cells exactly
cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
stopifnot(identical(findInterval(x, cuts) + 1L, cell))

# ---- range scaling, fixed sigma, and standardization, reproduced ----

yRange <- max(y) - min(y)
zScaled <- (y - min(y)) / yRange - 0.5
residVar <- (1 / yRange)^2 # sigma = fixed(1)
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

# ---- engine arms ----

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  keepTrees = TRUE,
  n.samples = batchSize,
  n.burn = nBurn,
  n.thin = nThin,
  updateState = TRUE,
  seed = engineSeed,
  n.cuts = K - 1L
)
engineArm <- function(xs, ys, weights = NULL) {
  # the test rows are read for their fits alone, so the notice that weights
  # do not reach them says nothing here
  sampler <- withCallingHandlers(
    dbarts(
      matrix(xs, ncol = 1L),
      ys,
      test = matrix(as.double(seq_len(K)), ncol = 1L),
      weights = weights,
      control = ctl,
      power = power,
      base = base,
      proposal.probs = c(
        birth_death = 0.5,
        swap = 0.1,
        change = 0.4,
        birth = 0.5
      ),
      leaf.prior = linear(1L, k = kLeaf),
      family = gaussian(sigma = fixed(1))
    ),
    warning = function(w) {
      if (grepl("ignored for test data", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  stopifnot(is.null(sampler$data@offset))

  nBatch <- as.integer(ceiling(nKept / batchSize))
  engineLeaves <- integer(0L)
  engineFits <- NULL
  first <- TRUE
  for (b in seq_len(nBatch)) {
    r <- if (first) {
      sampler$run(nBurn, batchSize)
    } else {
      sampler$run(0L, batchSize)
    }
    first <- FALSE
    tr <- sampler$getTrees()
    engineLeaves <- c(
      engineLeaves,
      as.integer(tapply(tr$var, tr$sample, function(v) sum(v == -1L)))
    )
    engineFits <- cbind(engineFits, r$test)
  }
  stopifnot(length(engineLeaves) == ncol(engineFits))
  list(leaves = engineLeaves, fits = engineFits)
}
outRows <- 200L
arms <- list(
  "every row in" = engineArm(x, y),
  "200 rows out at the first value" = engineArm(
    c(x, rep(1, outRows)),
    c(y, rep(mean(y), outRows)),
    rep(c(1, 0), c(n, outRows))
  )
)

# ---- batch-means MC error and z-scores ----

batchMeanSE <- function(v, nBatches = 400L) {
  len <- (length(v) %/% nBatches) * nBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
  sd(bm) / sqrt(nBatches)
}

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

for (label in names(arms)) {
  engineLeaves <- arms[[label]]$leaves
  engineFits <- arms[[label]]$fits
  cat(sprintf(
    "\n%s: engine draws %d (thin %d)\n",
    label,
    length(engineLeaves),
    nThin
  ))
  cat(sprintf(
    "%-12s %10s %10s %10s %8s\n",
    "statistic",
    "engine",
    "exact",
    "MCse",
    "z"
  ))
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
}

if (anyFailure) {
  cat("\nFAIL: linear-leaf sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("\nOK: linear-leaf sampler matches the exact posterior\n")
