#!/usr/bin/env Rscript

# Restricted-moderator exact-posterior gate for the bartcore BCF two-forest
# sampler (docs/design/bcf.md). The treatment forest tau(x) is confined to a
# strict subset of the design's columns (its moderators); the prognostic
# forest mu(x) is unrestricted. This gate builds a two-predictor problem -
# x1 over K1 = 3 ordinal cells and x2 over K2 = 2 ordinal cells, 6 crossed
# cells, a fixed 0/1 treatment balanced within each crossed cell, gaussian
# response - with mu depending on both predictors and tau on x1 only, and
# restricts tau to x1. The restricted posterior is well defined and differs
# from the unrestricted one, so a mask that leaked (tau reading x2) or a
# backfit that ignored the restriction would show up as an exact-vs-sampler
# gap.
#
# The exact side enumerates each single-tree forest's tree space and takes
# the product. tau's space is the 1-D contiguous-cell enumeration over x1's
# 3 cells (the restriction makes x1 its only available variable, so the
# variable-selection factor is trivial). mu's space is a bounded 2-D
# axis-aligned enumeration over the 6-cell grid: at each node the split-
# variable prior is uniform over the available variables and the cut prior
# is uniform over that variable's available cuts, exactly the engine's
# CGMTreePrior. Conditional on the glue and sigma the leaf parameters
# integrate in closed form (a Gaussian block marginal) with sigma vectorized
# through the whitened eigendecomposition; the range scaling, the chisq sigma
# calibration, and the leaf-scale map are reproduced as in the full-store
# gate. Only mode 1 (fixed glue a = 1, b0 = 0, b1 = 1) is run - the free-glue
# modes are already validated full-store - matching E[mu] over the 6 crossed
# cells and E[tau] over the 3 x1 cells.
#
# Usage: Rscript bcf-exact-restricted.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nSigma <- if (quick) 50L else 90L
nburn <- if (quick) 1000L else 3000L
nSeeds <- if (quick) 1L else 3L
tolerance <- if (quick) 0.05 else 0.015
ndpost <- if (quick) 1200L else 2500L
thin <- if (quick) 6L else 15L

# ---- fixed design: two ordinal predictors, tau restricted to x1 ----

set.seed(20260711L)
K1 <- 3L # x1 cells
K2 <- 2L # x2 cells
nCross <- K1 * K2 # crossed cells
nPerCell <- 40L # even, so z is balanced within each crossed cell

crossedCell <- rep(seq_len(nCross), each = nPerCell)
x1cell <- ((crossedCell - 1L) %/% K2) + 1L
x2cell <- ((crossedCell - 1L) %% K2) + 1L
x <- cbind(as.double(x1cell), as.double(x2cell))
z <- as.double(rep(rep(c(0, 1), nPerCell / 2L), nCross))

# mu depends on both predictors, tau on x1 only
muX1 <- c(-2.0, 0.0, 2.0)
muX2 <- c(-0.5, 0.5)
muTrue <- muX1[x1cell] + muX2[x2cell]
tauFun <- c(0.8, 1.0, 1.2)
noiseSd <- 0.2
# generated under the mode-1 glue (a = 1, b0 = 0, b1 = 1)
y <- muTrue + z * tauFun[x1cell] + rnorm(nCross * nPerCell, sd = noiseSd)
n <- length(y)

sigEst <- 0.3
sigDf <- 3
sigQuant <- 0.9

muBase <- 0.8
muPower <- 2
tauBase <- 0.25
tauPower <- 3
sdControl <- 2
sdModerate <- 1
bVar <- 0.5

# per-column n.cuts so the uniform grid places exactly one cut between each
# pair of adjacent cells (x2's two values admit a single cut; the shared
# default would place two degenerate ones)
nCutsX1 <- K1 - 1L
nCutsX2 <- K2 - 1L
cuts1 <- min(x[, 1L]) + seq_len(nCutsX1) * (max(x[, 1L]) - min(x[, 1L])) / K1
cuts2 <- min(x[, 2L]) + seq_len(nCutsX2) * (max(x[, 2L]) - min(x[, 2L])) / K2
stopifnot(identical(findInterval(x[, 1L], cuts1) + 1L, x1cell))
stopifnot(identical(findInterval(x[, 2L], cuts2) + 1L, x2cell))

# ---- range scaling and the calibration map, reproduced ----

yRange <- max(y) - min(y)
yScaled <- (y - min(y)) / yRange - 0.5
s <- sd(yScaled)
scaleMu <- s # mu leaf prior sd (single tree, k = 1)
scaleTau <- sdModerate * s / 0.674 # tau leaf prior sd (half-normal median)

# chisq sigma prior on the internal scale: sigma^2 ~ InvGamma(df/2, df*lambda/2)
initialSigma <- sigEst / yRange
sigmaRawScale <- qchisq(1 - sigQuant, sigDf) / sigDf
lambda <- initialSigma^2 * sigmaRawScale

# ---- tau tree enumeration (1-D contiguous-cell leaves over x1) ----

enumerate <- function(loCell, hiCell, loCut, hiCut, depth, base, power) {
  growth <- if (hiCut >= loCut) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(loCell, hiCell)),
    logPrior = log(1 - growth)
  ))
  if (hiCut < loCut) {
    return(result)
  }
  for (j in loCut:hiCut) {
    lefts <- enumerate(loCell, j, loCut, j - 1L, depth + 1L, base, power)
    rights <- enumerate(j + 1L, hiCell, j + 1L, hiCut, depth + 1L, base, power)
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

tauLeafAssignment <- function(tree) {
  leafOf <- integer(K1)
  for (li in seq_along(tree$leaves)) {
    rng <- tree$leaves[[li]]
    leafOf[rng[1L]:rng[2L]] <- li
  }
  leafOf
}

# ---- mu tree enumeration (2-D axis-aligned leaves over the crossed grid) ----
#
# A node covers a rectangle of cells [lo1, hi1] x [lo2, hi2]; x1 offers
# hi1 - lo1 cuts and x2 offers hi2 - lo2. The split-variable factor is
# 1 / numAvailable and the cut factor is 1 / (cuts of the chosen variable),
# matching CGMTreePrior with uniform split probabilities.

crossedIndex <- function(i1, i2) (i1 - 1L) * K2 + i2
x1OfCross <- ((seq_len(nCross) - 1L) %/% K2) + 1L

enumerate2D <- function(lo1, hi1, lo2, hi2, depth, base, power) {
  n1 <- hi1 - lo1 # x1 cuts available at this node
  n2 <- hi2 - lo2 # x2 cuts available at this node
  numAvail <- (n1 > 0L) + (n2 > 0L)
  growth <- if (numAvail > 0L) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(lo1, hi1, lo2, hi2)),
    logPrior = log(1 - growth)
  ))
  if (numAvail == 0L) {
    return(result)
  }
  logSplit <- log(growth) - log(numAvail)
  combine <- function(lefts, rights, logCut) {
    for (left in lefts) {
      for (right in rights) {
        result[[length(result) + 1L]] <<- list(
          leaves = c(left$leaves, right$leaves),
          logPrior = logSplit - logCut + left$logPrior + right$logPrior
        )
      }
    }
  }
  if (n1 > 0L) {
    for (j in lo1:(hi1 - 1L)) {
      combine(
        enumerate2D(lo1, j, lo2, hi2, depth + 1L, base, power),
        enumerate2D(j + 1L, hi1, lo2, hi2, depth + 1L, base, power),
        log(n1)
      )
    }
  }
  if (n2 > 0L) {
    for (j in lo2:(hi2 - 1L)) {
      combine(
        enumerate2D(lo1, hi1, lo2, j, depth + 1L, base, power),
        enumerate2D(lo1, hi1, j + 1L, hi2, depth + 1L, base, power),
        log(n2)
      )
    }
  }
  result
}

muLeafAssignment <- function(tree) {
  leafOf <- integer(nCross)
  for (li in seq_along(tree$leaves)) {
    r <- tree$leaves[[li]]
    for (i1 in r[1L]:r[2L]) {
      for (i2 in r[3L]:r[4L]) {
        leafOf[crossedIndex(i1, i2)] <- li
      }
    }
  }
  leafOf
}

tauTrees <- enumerate(1L, K1, 1L, K1 - 1L, 0L, tauBase, tauPower)
muTrees <- enumerate2D(1L, K1, 1L, K2, 0L, muBase, muPower)
tauAssign <- lapply(tauTrees, tauLeafAssignment)
muAssign <- lapply(muTrees, muLeafAssignment)
tauP <- vapply(tauTrees, function(t) length(t$leaves), integer(1L))
muP <- vapply(muTrees, function(t) length(t$leaves), integer(1L))
cat(sprintf(
  "enumerated %d mu-trees, %d tau-trees\n",
  length(muTrees),
  length(tauTrees)
))
# each tree prior is proper: the enumerated weights sum to 1
stopifnot(
  abs(sum(exp(vapply(muTrees, function(t) t$logPrior, 0))) - 1) < 1e-9,
  abs(sum(exp(vapply(tauTrees, function(t) t$logPrior, 0))) - 1) < 1e-9
)

# ---- cell sufficient statistics keyed on (crossed cell, z) ----

nCz <- matrix(0, nCross, 2L)
syCz <- matrix(0, nCross, 2L)
for (i in seq_len(n)) {
  zg <- if (z[i] != 0) 2L else 1L
  nCz[crossedCell[i], zg] <- nCz[crossedCell[i], zg] + 1
  syCz[crossedCell[i], zg] <- syCz[crossedCell[i], zg] + yScaled[i]
}
syy <- sum(yScaled^2)

# ---- sigma quadrature grid and its log prior density ----

sigmaGrid <- seq(0.03 * s, 1.5 * s, length.out = nSigma)
invSig2 <- 1 / sigmaGrid^2
logSigmaPrior <- (sigDf / 2) *
  log(sigDf * lambda / 2) -
  lgamma(sigDf / 2) -
  (sigDf / 2 + 1) * log(sigmaGrid^2) -
  sigDf * lambda / (2 * sigmaGrid^2) +
  log(2 * sigmaGrid)

# ---- exact posterior over the enumerated joint space ----
#
# Conditional on a glue row (a, b0, b1) and sigma the leaf block integrates in
# closed form; the sigma dimension is vectorized through the whitened
# eigendecomposition of the design crossproduct. The mu block spans the 6
# crossed cells and the tau block the 3 x1 cells, coupled per crossed cell
# through its x1 component.

logDet2pi <- n * log(2 * pi)
priorPrecMu <- 1 / scaleMu^2
priorPrecTau <- 1 / scaleTau^2

exactBCF <- function(glueGrid) {
  maxLW <- -Inf
  wSum <- 0
  accMu <- numeric(nCross)
  accTau <- numeric(K1)
  for (gi in seq_len(nrow(glueGrid))) {
    a <- glueGrid[gi, 1L]
    bz <- c(glueGrid[gi, 2L], glueGrid[gi, 3L])
    gluePrior <- glueGrid[gi, 4L]
    for (mt in seq_along(muTrees)) {
      muLeaf <- muAssign[[mt]]
      pmu <- muP[mt]
      muPrior <- muTrees[[mt]]$logPrior
      for (tt in seq_along(tauTrees)) {
        tauLeaf <- tauAssign[[tt]]
        qtau <- tauP[tt]
        k <- pmu + qtau
        dtd <- matrix(0, k, k)
        dty <- numeric(k)
        for (cc in seq_len(nCross)) {
          l <- muLeaf[cc]
          m <- pmu + tauLeaf[x1OfCross[cc]]
          for (zg in 1:2) {
            ncz <- nCz[cc, zg]
            if (ncz == 0) {
              next
            }
            sycz <- syCz[cc, zg]
            dtd[l, l] <- dtd[l, l] + a * a * ncz
            dtd[m, m] <- dtd[m, m] + bz[zg]^2 * ncz
            dtd[l, m] <- dtd[l, m] + a * bz[zg] * ncz
            dtd[m, l] <- dtd[l, m]
            dty[l] <- dty[l] + a * sycz
            dty[m] <- dty[m] + bz[zg] * sycz
          }
        }
        lVec <- sqrt(c(rep(priorPrecMu, pmu), rep(priorPrecTau, qtau)))
        logDetLambda0 <- 2 * sum(log(lVec))
        mMat <- dtd / outer(lVec, lVec)
        eg <- eigen(mMat, symmetric = TRUE)
        eVal <- eg$values
        eVal[eVal < 0] <- 0
        c0 <- as.vector(crossprod(eg$vectors, dty / lVec))
        den <- 1 + outer(eVal, invSig2) # k x nSigma
        logDetA <- logDetLambda0 + colSums(log(den))
        quad <- (invSig2^2) * colSums((c0^2) / den)
        logMarg <- 0.5 *
          logDetLambda0 -
          0.5 * logDetA -
          0.5 * syy * invSig2 +
          0.5 * quad -
          0.5 * (logDet2pi - n * log(invSig2))
        gMat <- (c0 %o% invSig2) / den
        theta <- (eg$vectors %*% gMat) / lVec # k x nSigma
        logW <- muPrior +
          tauTrees[[tt]]$logPrior +
          gluePrior +
          logSigmaPrior +
          logMarg
        muCellByS <- theta[muLeaf, , drop = FALSE]
        tauCellByS <- theta[pmu + tauLeaf, , drop = FALSE]
        batchMax <- max(logW)
        if (batchMax > maxLW) {
          f <- exp(maxLW - batchMax)
          wSum <- wSum * f
          accMu <- accMu * f
          accTau <- accTau * f
          maxLW <- batchMax
        }
        wS <- exp(logW - maxLW)
        wSum <- wSum + sum(wS)
        accMu <- accMu + as.vector(muCellByS %*% wS)
        accTau <- accTau + as.vector(tauCellByS %*% wS)
      }
    }
  }
  list(mu = accMu / wSum, tau = accTau / wSum)
}

# ---- sampler collection ----

repMu <- vapply(seq_len(nCross), function(cc) which(crossedCell == cc)[1L], 1L)
repTau <- vapply(seq_len(K1), function(c1) which(x1cell == c1)[1L], 1L)

samplerFit <- function(seed) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    n.cuts = c(nCutsX1, nCutsX2),
    updateState = FALSE
  )
  host <- dbarts(
    x,
    y,
    control = control,
    sigma = sigEst,
    tree.prior = cgm(muPower, muBase),
    node.prior = normal(2)
  )
  bc <- dbarts:::bartcoreBCFSampler(
    host,
    z,
    n.trees.treatment = 1L,
    treatment.base = tauBase,
    treatment.power = tauPower,
    sd.control = sdControl,
    sd.moderate = sdModerate,
    b.prior.variance = bVar,
    update.a = FALSE,
    update.b = FALSE,
    moderators = 1L
  )
  bartcoreRun(bc, nburn, 1L)
  # containment: a restricted tau forest splits on x1 only, so its fit is
  # constant across x2 within each x1 cell
  tauFull <- bartcoreForestFits(bc, 1L)[, 1L]
  stopifnot(all(tapply(tauFull, x1cell, function(v) diff(range(v))) < 1e-9))
  muM <- matrix(0, ndpost, nCross)
  tauM <- matrix(0, ndpost, K1)
  for (d in seq_len(ndpost)) {
    bartcoreRun(bc, 0L, thin)
    muM[d, ] <- bartcoreForestFits(bc, 0L)[repMu, 1L]
    tauM[d, ] <- bartcoreForestFits(bc, 1L)[repTau, 1L]
  }
  list(mu = colMeans(muM), tau = colMeans(tauM))
}

fits <- lapply(seq_len(nSeeds), samplerFit)
fit <- list(
  mu = rowMeans(vapply(fits, function(f) f$mu, numeric(nCross))),
  tau = rowMeans(vapply(fits, function(f) f$tau, numeric(K1)))
)

# ---- mode 1: glue fixed at (a, b0, b1) = (1, 0, 1) ----

exact <- exactBCF(matrix(c(1, 0, 1, 0), nrow = 1L))

report <- function(label, exactValue, fitValue) {
  gap <- max(abs(exactValue - fitValue))
  failed <- !is.finite(gap) || gap > tolerance
  cat(sprintf(
    "  %-6s exact %s | sampler %s (gap %.4f%s)\n",
    label,
    paste(sprintf("%+.4f", exactValue), collapse = " "),
    paste(sprintf("%+.4f", fitValue), collapse = " "),
    gap,
    if (failed) " <- FAIL" else ""
  ))
  failed
}

cat("restricted moderators (tau on x1 only), fixed glue:\n")
anyFailure <- report("E[mu]", exact$mu, fit$mu)
anyFailure <- report("E[tau]", exact$tau, fit$tau) || anyFailure

if (anyFailure) {
  cat("FAIL: restricted BCF sampler deviates from exact\n")
  quit(status = 1L)
}
cat("OK: restricted BCF sampler matches the exact posterior\n")
