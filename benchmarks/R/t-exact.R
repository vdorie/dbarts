#!/usr/bin/env Rscript

# Validation of the bartcore Student-t (robust) residual sampler against the
# exact posterior (docs/design/robust-errors.md).
#
# The gate: single-tree BART on one predictor with 3 cut points admits
# brute-force enumeration of the tree space, and with the degrees of freedom
# nu FIXED and sigma FIXED (conditioned on) each leaf's marginal and
# posterior-mean fit are 1-D quadratures over the leaf parameter mu. The
# reference integrates the EXACT product of marginal Student-t densities
# (sqrt(w) (z - mu) / sigma ~ t_nu) against the N(0, (nodeScale/k)^2) leaf
# prior; it NEVER augments the per-observation precisions, so agreement with
# the sampler validates the scale-mixture augmentation the engine runs
# (lambda_i ~ Gamma, entering the node statistics as the logistic omega does).
# Data carry genuine t_nu noise so the tails - where lambda downweighting
# matters - are exercised. The recorded fit is the per-cell posterior-mean
# leaf value, compared to range * E[mu] + shift under the engine's internal
# [-0.5, 0.5] rescaling. Failure means the augmentation or the reference is
# wrong (perturbing the lambda rate in the engine moves the fit and this gate
# catches it).
#
# Usage: Rscript t-exact.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

exactNdpost <- if (quick) 30000L else 90000L
exactBurn <- if (quick) 5000L else 12000L
exactSeeds <- if (quick) 2L else 4L
exactTolerance <- if (quick) 0.02 else 0.012 # vs single-chain MC error

# ---- fixed data: four cells, genuine Student-t errors ----

nu <- 4
sigmaTrue <- 0.6
sigmaFixed <- sigmaTrue

set.seed(3141L)
nPerCell <- 12L
cellCenters <- c(0.125, 0.375, 0.625, 0.875)
x1 <- matrix(
  rep(cellCenters, each = nPerCell) + runif(4L * nPerCell, -0.1, 0.1),
  ncol = 1L
)
cell <- rep(1:4, each = nPerCell)
cellMean <- c(-0.8, -0.2, 0.4, 1.0)
y1 <- cellMean[cell] + sigmaFixed * rt(length(cell), df = nu)

cuts <- min(x1) + (1:3) * (max(x1) - min(x1)) / 4
stopifnot(identical(cell, findInterval(x1[, 1L], cuts) + 1L))

base <- 0.5
power <- 2
k <- 2
nodeScale <- 1.5

# the engine's internal [-0.5, 0.5] rescaling of the observed response
fitMin <- min(y1)
fitRange <- max(y1) - fitMin
fitShift <- fitRange * 0.5 + fitMin
zInternal <- (y1 - fitMin) / fitRange - 0.5
sigmaInt <- sigmaFixed / fitRange
zByCell <- split(zInternal, cell)

# ---- exact posterior over the enumerated single-tree space ----

# CGM prior over trees on cut interval [loCut, hiCut]; leaves are contiguous
# cell ranges (identical to logistic-reference.R's enumeration).
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
trees <- enumerate(1L, 4L, 1L, 3L, 0L)

muGrid <- seq(-4, 4, by = 0.001)
dmu <- 0.001
tau <- nodeScale / k
logMuPrior <- dnorm(muGrid, 0, tau, log = TRUE)

# leaf marginal from the EXACT marginal-t likelihood (lambda integrated out
# analytically), never augmented: sqrt(w) (z - mu) / sigma ~ t_nu, w = 1.
leafQuantities <- function(cells) {
  zc <- unlist(zByCell[as.character(cells)])
  logLik <- numeric(length(muGrid))
  for (z in zc) {
    logLik <- logLik +
      dt((z - muGrid) / sigmaInt, df = nu, log = TRUE) -
      log(sigmaInt)
  }
  logIntegrand <- logMuPrior + logLik
  m <- max(logIntegrand)
  w <- exp(logIntegrand - m)
  list(
    logMarginal = log(sum(w) * dmu) + m,
    postMean = sum(muGrid * w) / sum(w)
  )
}

exactCellFits <- function() {
  logWeights <- numeric(length(trees))
  postMeanMu <- matrix(0, length(trees), 4L)
  for (t in seq_along(trees)) {
    logWeight <- trees[[t]]$logPrior
    for (leaf in trees[[t]]$leaves) {
      cells <- leaf[1L]:leaf[2L]
      q <- leafQuantities(cells)
      logWeight <- logWeight + q$logMarginal
      postMeanMu[t, cells] <- q$postMean
    }
    logWeights[t] <- logWeight
  }
  weights <- exp(logWeights - max(logWeights))
  muHat <- colSums(weights * postMeanMu) / sum(weights)
  fitRange * muHat + fitShift # de-scale to the response scale
}

# ---- sampler fit: fixed nu, conditioned (fixed) sigma ----

fitSingleTree <- function(seed) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    n.cuts = 3L,
    updateState = FALSE
  )
  sampler <- dbarts(
    x1,
    y1,
    control = control,
    node.prior = normal(k),
    tree.prior = cgm(power, base),
    resid.prior = fixed(sigmaFixed^2),
    resid.dist = student(df = nu)
  )
  sampler$model@node.scale <- nodeScale
  bc <- dbarts:::bartcoreSampler(sampler)
  r <- bartcoreRun(bc, exactBurn, exactNdpost)
  # every observation in a cell shares its leaf's fit; take one per cell
  reps <- vapply(1:4, function(cc) which(cell == cc)[1L], integer(1L))
  rowMeans(r$train[reps, , drop = FALSE])
}

exact <- exactCellFits()
fit <- colMeans(do.call(rbind, lapply(seq_len(exactSeeds), fitSingleTree)))
gap <- max(abs(fit - exact))

cat(sprintf(
  "t (nu=%g)  exact %s | sampler %s | max gap %.4f%s\n",
  nu,
  paste(sprintf("%.4f", exact), collapse = " "),
  paste(sprintf("%.4f", fit), collapse = " "),
  gap,
  if (gap > exactTolerance) " <- FAIL" else ""
))

if (gap > exactTolerance) {
  quit(status = 1L)
}
cat("\nOK: Student-t sampler matches the exact posterior\n")
