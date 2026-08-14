#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore BCF two-forest sampler
# (docs/design/bcf.md). One ordinal predictor over a handful of distinct
# values, a fixed 0/1 treatment balanced within each cell, and a gaussian
# response admit brute-force enumeration of each single-tree forest's tree
# space; the joint space is the product. Conditional on the glue (a, b0, b1)
# and sigma the leaf parameters integrate in closed form (a Gaussian block
# marginal), leaving a low-dimensional quadrature over the glue and a 1-D
# quadrature over sigma against the sampler's actual prior on the range-scaled
# response. The gate reproduces dbarts's range scaling, the chisq sigma
# calibration, and the deliverable-1 leaf-scale map, so it validates the map
# and the two-forest backfit end to end. Failure means the backfit, the glue
# draw, or the calibration map is wrong.
#
# Modes: (1) glue fixed a = 1, b0 = 0, b1 = 1 - match E[mu], E[tau];
# (2a) a free (Cauchy(0, aPriorScale)), b fixed - match E[a mu], E[tau];
# (2b) b0/b1 free (N(0, bPriorVariance)), a fixed - match E[mu], E[(b1-b0)tau].
#
# Usage: Rscript bcf-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nSigma <- if (quick) 50L else 90L
aN <- if (quick) 121L else 201L
bN <- if (quick) 45L else 71L
bRange <- if (quick) 3.5 else 4.0
nburn <- if (quick) 1000L else 3000L
nSeeds <- if (quick) 1L else 3L
tolerance <- if (quick) 0.05 else 0.015
# modes 1 and 2a mix well. Mode 2b frees b0/b1, whose scale trades off with
# tau's (the bscale ridge); mu couples to the poorly-mixing ratio
# (b0 + b1) / (b1 - b0), so its E[mu] is a high-variance MCMC target needing
# many chains and a looser bound. The identified E[(b1-b0)tau] stays tight.
toleranceMu2b <- if (quick) 0.07 else 0.035
ndLight <- if (quick) 1200L else 2500L
thinLight <- if (quick) 6L else 15L
nSeeds2b <- if (quick) 4L else 8L
nd2b <- if (quick) 2000L else 3000L
thin2b <- if (quick) 150L else 200L

# ---- fixed design ----

set.seed(20260707L)
K <- 3L # distinct predictor cells
nPerCell <- 50L # even, so z is balanced within each cell
cell <- rep(seq_len(K), each = nPerCell)
x <- matrix(as.double(cell), ncol = 1L)
z <- as.double(rep(rep(c(0, 1), nPerCell / 2L), K))
# controls carry a treatment-scaled signal too (b0True != 0), so both b
# coefficients are identifiable and mode 2b's mu is well conditioned
muTrue <- c(-2.0, 0.0, 2.0)
tauFun <- c(0.8, 1.0, 1.2)
b0True <- 0.5
b1True <- 1.5
noiseSd <- 0.2
y <- muTrue[cell] +
  (b0True + (b1True - b0True) * z) * tauFun[cell] +
  rnorm(K * nPerCell, sd = noiseSd)
n <- length(y)

sigEst <- 0.3 # the sigma estimate anchoring the chisq prior
sigDf <- 3
sigQuant <- 0.9

muBase <- 0.8
muPower <- 2
tauBase <- 0.25
tauPower <- 3
sdControl <- 2
sdModerate <- 1
bVar <- 0.5

# uniform cuts separate the cells exactly (dbarts n.cuts = K - 1)
cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
stopifnot(identical(findInterval(x[, 1L], cuts) + 1L, cell))

# ---- range scaling and the calibration map, reproduced ----

yRange <- max(y) - min(y)
yScaled <- (y - min(y)) / yRange - 0.5
s <- sd(yScaled)
scaleMu <- s # mu leaf prior sd (single tree, k = 1)
scaleTau <- sdModerate * s / 0.674 # tau leaf prior sd (half-normal median)
aPriorScale <- sdControl

# chisq sigma prior on the internal scale (ChiSquaredScalePrior): sigma^2 ~
# df * lambda / chisq(df), i.e. InvGamma(df/2, df*lambda/2)
initialSigma <- sigEst / yRange
sigmaRawScale <- qchisq(1 - sigQuant, sigDf) / sigDf
lambda <- initialSigma^2 * sigmaRawScale

# ---- tree enumeration (ordinal contiguous-cell leaves) ----

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

leafAssignment <- function(tree) {
  leafOf <- integer(K)
  for (li in seq_along(tree$leaves)) {
    rng <- tree$leaves[[li]]
    leafOf[rng[1L]:rng[2L]] <- li
  }
  leafOf
}

muTrees <- enumerate(1L, K, 1L, K - 1L, 0L, muBase, muPower)
tauTrees <- enumerate(1L, K, 1L, K - 1L, 0L, tauBase, tauPower)
muAssign <- lapply(muTrees, leafAssignment)
tauAssign <- lapply(tauTrees, leafAssignment)
muP <- vapply(muTrees, function(t) length(t$leaves), integer(1L))
tauP <- vapply(tauTrees, function(t) length(t$leaves), integer(1L))
cat(sprintf(
  "enumerated %d mu-trees, %d tau-trees\n",
  length(muTrees),
  length(tauTrees)
))

# ---- cell sufficient statistics ----

nCz <- matrix(0, K, 2L)
syCz <- matrix(0, K, 2L)
for (i in seq_len(n)) {
  zg <- if (z[i] != 0) 2L else 1L
  nCz[cell[i], zg] <- nCz[cell[i], zg] + 1
  syCz[cell[i], zg] <- syCz[cell[i], zg] + yScaled[i]
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

# glueGrid columns: a, b0, b1, logPrior. Conditional on a glue row and sigma
# the leaf block integrates in closed form; the sigma dimension is vectorized
# through the whitened eigendecomposition of the design crossproduct.
logDet2pi <- n * log(2 * pi)
priorPrecMu <- 1 / scaleMu^2
priorPrecTau <- 1 / scaleTau^2

exactBCF <- function(glueGrid) {
  maxLW <- -Inf
  wSum <- 0
  accMu <- numeric(K)
  accTau <- numeric(K)
  accAMu <- numeric(K)
  accBTau <- numeric(K)
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
        for (c in seq_len(K)) {
          l <- muLeaf[c]
          m <- pmu + tauLeaf[c]
          for (zg in 1:2) {
            ncz <- nCz[c, zg]
            if (ncz == 0) {
              next
            }
            sycz <- syCz[c, zg]
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
          accAMu <- accAMu * f
          accBTau <- accBTau * f
          maxLW <- batchMax
        }
        wS <- exp(logW - maxLW)
        wSum <- wSum + sum(wS)
        muW <- as.vector(muCellByS %*% wS)
        tauW <- as.vector(tauCellByS %*% wS)
        accMu <- accMu + muW
        accTau <- accTau + tauW
        accAMu <- accAMu + a * muW
        accBTau <- accBTau + (bz[2L] - bz[1L]) * tauW
      }
    }
  }
  list(
    mu = accMu / wSum,
    tau = accTau / wSum,
    aMu = accAMu / wSum,
    bTau = accBTau / wSum
  )
}

# ---- sampler collection ----

repObs <- vapply(seq_len(K), function(c) which(cell == c)[1L], integer(1L))

samplerFit <- function(seed, updateA, updateB, ndpost, thin) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    n.cuts = K - 1L,
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
    update.a = updateA,
    update.b = updateB
  )
  dbarts:::bartcoreRun(bc, nburn, 1L)
  muM <- matrix(0, ndpost, K)
  tauM <- matrix(0, ndpost, K)
  aVec <- numeric(ndpost)
  bDiff <- numeric(ndpost)
  for (d in seq_len(ndpost)) {
    dbarts:::bartcoreRun(bc, 0L, thin)
    muM[d, ] <- dbarts:::bartcoreForestFits(bc, 0L)[repObs, 1L]
    tauM[d, ] <- dbarts:::bartcoreForestFits(bc, 1L)[repObs, 1L]
    g <- dbarts:::bartcoreForestAmplitudes(bc)[, 1L]
    aVec[d] <- g[1L]
    bDiff[d] <- g[3L] - g[2L]
  }
  list(
    mu = colMeans(muM),
    tau = colMeans(tauM),
    aMu = colMeans(muM * aVec),
    bTau = colMeans(tauM * bDiff)
  )
}

samplerAverage <- function(updateA, updateB, ndpost, thin, seeds) {
  fits <- lapply(
    seq_len(seeds),
    samplerFit,
    updateA = updateA,
    updateB = updateB,
    ndpost = ndpost,
    thin = thin
  )
  list(
    mu = rowMeans(vapply(fits, function(f) f$mu, numeric(K))),
    tau = rowMeans(vapply(fits, function(f) f$tau, numeric(K))),
    aMu = rowMeans(vapply(fits, function(f) f$aMu, numeric(K))),
    bTau = rowMeans(vapply(fits, function(f) f$bTau, numeric(K)))
  )
}

# ---- run the three modes ----

reportMode <- function(
  name,
  exLeft,
  exRight,
  fitLeft,
  fitRight,
  labLeft,
  labRight,
  tolLeft,
  tolRight
) {
  gapLeft <- max(abs(exLeft - fitLeft))
  gapRight <- max(abs(exRight - fitRight))
  cat(sprintf("mode %s\n", name))
  cat(sprintf(
    "  %-10s exact %s | sampler %s (gap %.4f)\n",
    labLeft,
    paste(sprintf("%+.4f", exLeft), collapse = " "),
    paste(sprintf("%+.4f", fitLeft), collapse = " "),
    gapLeft
  ))
  cat(sprintf(
    "  %-10s exact %s | sampler %s (gap %.4f)\n",
    labRight,
    paste(sprintf("%+.4f", exRight), collapse = " "),
    paste(sprintf("%+.4f", fitRight), collapse = " "),
    gapRight
  ))
  failed <- gapLeft > tolLeft || gapRight > tolRight
  if (failed) {
    cat("  <- FAIL\n")
  }
  failed
}

anyFailure <- FALSE

# mode 1: glue fixed at (1, 0, 1)
ex1 <- exactBCF(matrix(c(1, 0, 1, 0), nrow = 1L))
fit1 <- samplerAverage(FALSE, FALSE, ndLight, thinLight, nSeeds)
anyFailure <- reportMode(
  "1 (fixed glue)",
  ex1$mu,
  ex1$tau,
  fit1$mu,
  fit1$tau,
  "E[mu]",
  "E[tau]",
  tolerance,
  tolerance
) ||
  anyFailure

# mode 2a: a free ~ Cauchy(0, aPriorScale); b fixed at (0, 1)
aGrid <- seq(-6, 6, length.out = aN)
glue2a <- cbind(aGrid, 0, 1, dcauchy(aGrid, 0, aPriorScale, log = TRUE))
ex2a <- exactBCF(glue2a)
fit2a <- samplerAverage(TRUE, FALSE, ndLight, thinLight, nSeeds)
anyFailure <- reportMode(
  "2a (a free)",
  ex2a$aMu,
  ex2a$tau,
  fit2a$aMu,
  fit2a$tau,
  "E[a mu]",
  "E[tau]",
  tolerance,
  tolerance
) ||
  anyFailure

# mode 2b: b0/b1 free ~ N(0, bPriorVariance); a fixed at 1
bGrid <- seq(-bRange, bRange, length.out = bN)
bPair <- expand.grid(b0 = bGrid, b1 = bGrid)
glue2b <- cbind(
  1,
  bPair$b0,
  bPair$b1,
  dnorm(bPair$b0, 0, sqrt(bVar), log = TRUE) +
    dnorm(bPair$b1, 0, sqrt(bVar), log = TRUE)
)
ex2b <- exactBCF(glue2b)
fit2b <- samplerAverage(FALSE, TRUE, nd2b, thin2b, nSeeds2b)
anyFailure <- reportMode(
  "2b (b free)",
  ex2b$mu,
  ex2b$bTau,
  fit2b$mu,
  fit2b$bTau,
  "E[mu]",
  "E[(b1-b0)tau]",
  toleranceMu2b,
  tolerance
) ||
  anyFailure

if (anyFailure) {
  quit(status = 1L)
}
cat("OK: BCF sampler matches the exact posterior\n")
