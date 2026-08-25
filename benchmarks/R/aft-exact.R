#!/usr/bin/env Rscript

# Validation of the bartcore AFT log-normal survival sampler against the
# exact posterior (docs/design/survival.md).
#
# The gate: single-tree BART on one predictor with 3 cut points admits
# brute-force enumeration of the tree space, and with sigma FIXED each leaf's
# marginal and posterior-mean linear predictor are 1-D quadratures over the
# leaf parameter. Uncensored observations enter as gaussian data on log T;
# right-censored ones contribute the upper normal tail (the survival
# probability past their log censoring time) - the exact analogue of the
# engine's truncated-normal latent augmentation. The recorded fit is the
# linear predictor E[log T | x] on the log scale, so the sampler's per-cell
# posterior-mean fit is compared to range * E[mu] + shift under the engine's
# internal rescaling. Failure means the truncated-latent mechanism or the
# leaf conjugacy is wrong (the poison test in the header of the design work
# perturbs the truncation bound and this gate catches it).
#
# Usage: Rscript aft-exact.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

exactNdpost <- if (quick) 20000L else 60000L
exactSeeds <- if (quick) 2L else 4L
exactTolerance <- if (quick) 0.02 else 0.012 # vs single-chain MC error

# ---- fixed data: four cells, a mix of events and right-censored ----

set.seed(2718L)
nPerCell <- 12L
cellCenters <- c(0.125, 0.375, 0.625, 0.875)
x1 <- matrix(
  rep(cellCenters, each = nPerCell) + runif(4L * nPerCell, -0.1, 0.1),
  ncol = 1L
)
cell <- rep(1:4, each = nPerCell)

# increasing true log-time means across cells, gaussian log-time errors
sigmaTrue <- 0.6
cellMean <- c(-0.8, -0.2, 0.4, 1.0)
logTtrue <- cellMean[cell] + sigmaTrue * rnorm(length(cell))
# independent right-censoring; censored observations keep their censoring time
logC <- cellMean[cell] + 0.3 + sigmaTrue * rnorm(length(cell))
status <- as.numeric(logTtrue <= logC) # 1 event, 0 censored
obsLogT <- ifelse(status == 1, logTtrue, logC)
cat(sprintf("censoring rate: %.2f\n", mean(status == 0)))

cuts <- min(x1) + (1:3) * (max(x1) - min(x1)) / 4
stopifnot(identical(cell, findInterval(x1[, 1L], cuts) + 1L))

base <- 0.5
power <- 2
k <- 2
nodeScale <- 1.5
sigmaFixed <- sigmaTrue

# the engine's internal [-0.5, 0.5] rescaling of the observed log-times
fitMin <- min(obsLogT)
fitRange <- max(obsLogT) - fitMin
fitShift <- fitRange * 0.5 + fitMin
zInternal <- (obsLogT - fitMin) / fitRange - 0.5
sigmaInt <- sigmaFixed / fitRange
zByCell <- split(zInternal[status == 1], cell[status == 1])
boundByCell <- split(zInternal[status == 0], cell[status == 0])
getCellList <- function(lst, cells) {
  unlist(lapply(cells, function(cc) {
    v <- lst[[as.character(cc)]]
    if (is.null(v)) numeric() else v
  }))
}

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

leafQuantities <- function(cells) {
  zc <- getCellList(zByCell, cells)
  bc <- getCellList(boundByCell, cells)
  logLik <- numeric(length(muGrid))
  for (z in zc) {
    logLik <- logLik + dnorm(z, muGrid, sigmaInt, log = TRUE)
  }
  for (b in bc) {
    # log P(log T > log C) = upper normal tail on the internal scale
    logLik <- logLik +
      pnorm(b, muGrid, sigmaInt, lower.tail = FALSE, log.p = TRUE)
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
  fitRange * muHat + fitShift # de-scale to the log-time linear predictor
}

# ---- sampler fit ----

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
    obsLogT,
    control = control,
    node.prior = normal(k),
    tree.prior = cgm(power, base),
    resid.prior = fixed(sigmaFixed^2)
  )
  sampler$model@node.scale <- nodeScale
  ctrl <- sampler$control
  attr(ctrl, "bartcore.survival") <- status
  sampler$control <- ctrl
  bc <- dbarts:::bartcoreSampler(sampler, family = "aft")
  r <- bartcoreRun(bc, 5000L, exactNdpost)
  # every observation in a cell shares its leaf's fit; take one per cell
  reps <- vapply(1:4, function(cc) which(cell == cc)[1L], integer(1L))
  rowMeans(r$train[reps, , drop = FALSE])
}

exact <- exactCellFits()
fit <- colMeans(do.call(rbind, lapply(seq_len(exactSeeds), fitSingleTree)))
gap <- max(abs(fit - exact))

cat(sprintf(
  "aft exact  %s | sampler %s | max gap %.4f%s\n",
  paste(sprintf("%.4f", exact), collapse = " "),
  paste(sprintf("%.4f", fit), collapse = " "),
  gap,
  if (gap > exactTolerance) " <- FAIL" else ""
))

if (gap > exactTolerance) {
  quit(status = 1L)
}
cat("\nOK: aft sampler matches the exact posterior\n")
