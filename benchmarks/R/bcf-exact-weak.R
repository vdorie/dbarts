#!/usr/bin/env Rscript

# Prior-dominated exact-posterior gate for the BCF a-glue draw
# (gate-hardening-1.0 sub-item 1; poison-sweep finding: dropping the prior
# precision 1/aVariance from the a full conditional at bcf-exact's data size
# moves E[a mu] by 1e-4 and NO gate sees it). Here n = 40 with a WEAK mu
# signal and a SMALL Cauchy scale (sd.control = 0.5), so the a posterior is
# prior-dominated: the 1/aVariance term carries a large share of the
# conditional precision and dropping it visibly spreads the a posterior.
# Because the ridge-regime a posterior keeps Cauchy tails (|a| has no finite
# mean), the spread channel is gated through P(|a| <= q) at fixed
# thresholds - the variance-channel check a mean-only gate is blind to -
# alongside the identified E[a mu] and E[tau].
#
# The exact side reuses bcf-exact.R's machinery (range scaling, chisq sigma
# calibration, leaf-scale map, tree enumeration, closed-form leaf block with
# sigma vectorized through the whitened eigendecomposition), with the a
# integral substituted a = sd.control * tan(t): the Cauchy prior becomes
# uniform in t, so the quadrature covers the whole real line and the heavy
# tails cost nothing. An SBC experiment (bcf-weak) validated that this
# prior-dominant configuration is well-behaved; per its caution about the
# (a, mu) scale ridge, burn-in is generous and draws are thinned (the |a|
# autocorrelation at this n is reported).
#
# Usage: Rscript bcf-exact-weak.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nSigma <- if (quick) 50L else 90L
tN <- if (quick) 301L else 601L # t-grid points for a = scale * tan(t)
nburn <- if (quick) 2000L else 4000L
nSeeds <- if (quick) 2L else 4L
ndpost <- if (quick) 1500L else 3000L
thin <- if (quick) 8L else 15L
toleranceMean <- if (quick) 0.05 else 0.02 # E[a mu], E[tau] (internal x range)
toleranceProb <- if (quick) 0.05 else 0.025 # P(|a| <= q)

# ---- fixed design: small n, weak mu, real tau ----

set.seed(20260711L)
K <- 2L
nPerCell <- 20L # n = 40, the SBC-validated prior-dominant size
cell <- rep(seq_len(K), each = nPerCell)
x <- matrix(as.double(cell), ncol = 1L)
z <- as.double(rep(rep(c(0, 1), nPerCell / 2L), K))
muTrue <- c(0.2, -0.2) # weak: the a data precision stays small
tauFun <- c(0.8, 1.2)
noiseSd <- 1.0
y <- muTrue[cell] + z * tauFun[cell] + rnorm(K * nPerCell, sd = noiseSd)
n <- length(y)

sigEst <- 1.0
sigDf <- 3
sigQuant <- 0.9

muBase <- 0.8
muPower <- 2
tauBase <- 0.25
tauPower <- 3
sdControl <- 0.5 # SMALL Cauchy scale: prior precision dominates
sdModerate <- 1
bVar <- 0.5 # inert (b fixed), passed through to the sampler

aThresholds <- c(0.25, 0.5, 1.0)

# uniform cuts separate the cells exactly (dbarts n.cuts = K - 1)
cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
stopifnot(identical(findInterval(x[, 1L], cuts) + 1L, cell))

# ---- range scaling and the calibration map, reproduced (bcf-exact.R) ----

yRange <- max(y) - min(y)
yScaled <- (y - min(y)) / yRange - 0.5
s <- sd(yScaled)
scaleMu <- s # mu leaf prior sd (single tree, k = 1)
scaleTau <- sdModerate * s / 0.674 # tau leaf prior sd (half-normal median)
aPriorScale <- sdControl

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
#
# glueGrid columns: a, b0, b1, logPrior. With the tan substitution the
# Cauchy prior on a is uniform in t, so every row carries logPrior = 0 and
# the normalizing constant absorbs the (dt / pi) factor. The rest is
# bcf-exact.R's exactBCF, with P(|a| <= q) accumulators added.

logDet2pi <- n * log(2 * pi)
priorPrecMu <- 1 / scaleMu^2
priorPrecTau <- 1 / scaleTau^2

exactWeak <- function(glueGrid, thresholds) {
  maxLW <- -Inf
  wSum <- 0
  accTau <- numeric(K)
  accAMu <- numeric(K)
  accProbA <- numeric(length(thresholds))
  for (gi in seq_len(nrow(glueGrid))) {
    a <- glueGrid[gi, 1L]
    bz <- c(glueGrid[gi, 2L], glueGrid[gi, 3L])
    gluePrior <- glueGrid[gi, 4L]
    inBall <- abs(a) <= thresholds
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
          accTau <- accTau * f
          accAMu <- accAMu * f
          accProbA <- accProbA * f
          maxLW <- batchMax
        }
        wS <- exp(logW - maxLW)
        wTotal <- sum(wS)
        wSum <- wSum + wTotal
        muW <- as.vector(muCellByS %*% wS)
        tauW <- as.vector(tauCellByS %*% wS)
        accTau <- accTau + tauW
        accAMu <- accAMu + a * muW
        accProbA <- accProbA + inBall * wTotal
      }
    }
  }
  list(
    tau = accTau / wSum,
    aMu = accAMu / wSum,
    probA = accProbA / wSum
  )
}

# ---- sampler collection ----

repObs <- vapply(seq_len(K), function(c) which(cell == c)[1L], integer(1L))

samplerFit <- function(seed) {
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
    update.a = TRUE,
    update.b = FALSE
  )
  dbarts:::bartcoreRun(bc, nburn, 1L)
  muM <- matrix(0, ndpost, K)
  tauM <- matrix(0, ndpost, K)
  aVec <- numeric(ndpost)
  for (d in seq_len(ndpost)) {
    dbarts:::bartcoreRun(bc, 0L, thin)
    muM[d, ] <- dbarts:::bartcoreForestFits(bc, 0L)[repObs, 1L]
    tauM[d, ] <- dbarts:::bartcoreForestFits(bc, 1L)[repObs, 1L]
    aVec[d] <- dbarts:::bartcoreForestAmplitudes(bc)[1L, 1L]
  }
  list(
    tau = colMeans(tauM),
    aMu = colMeans(muM * aVec),
    probA = vapply(aThresholds, function(q) mean(abs(aVec) <= q), 0),
    aAcf1 = acf(abs(aVec), plot = FALSE, lag.max = 1L)$acf[2L]
  )
}

fits <- lapply(seq_len(nSeeds), samplerFit)
fit <- list(
  tau = rowMeans(vapply(fits, function(f) f$tau, numeric(K))),
  aMu = rowMeans(vapply(fits, function(f) f$aMu, numeric(K))),
  probA = rowMeans(vapply(
    fits,
    function(f) f$probA,
    numeric(length(aThresholds))
  ))
)
cat(sprintf(
  "|a| lag-1 autocorrelation of thinned draws, per seed: %s\n",
  paste(sprintf("%.2f", vapply(fits, function(f) f$aAcf1, 0)), collapse = " ")
))

# ---- exact side: uniform t-grid, a = aPriorScale * tan(t) ----

tGrid <- seq(-pi / 2, pi / 2, length.out = tN + 2L)
tGrid <- tGrid[-c(1L, tN + 2L)] # drop the poles
glueGrid <- cbind(aPriorScale * tan(tGrid), 0, 1, 0)
exact <- exactWeak(glueGrid, aThresholds)

# ---- verdict ----

anyFailure <- FALSE
report <- function(label, exactValue, fitValue, tolerance) {
  gap <- max(abs(exactValue - fitValue))
  failed <- !is.finite(gap) || gap > tolerance
  cat(sprintf(
    "  %-14s exact %s | sampler %s (gap %.4f%s)\n",
    label,
    paste(sprintf("%+.4f", exactValue), collapse = " "),
    paste(sprintf("%+.4f", fitValue), collapse = " "),
    gap,
    if (failed) " <- FAIL" else ""
  ))
  failed
}

cat("prior-dominated a (Cauchy scale 0.5, n = 40):\n")
anyFailure <- report("E[a mu]", exact$aMu, fit$aMu, toleranceMean) || anyFailure
anyFailure <- report("E[tau]", exact$tau, fit$tau, toleranceMean) || anyFailure
anyFailure <- report(
  sprintf("P(|a|<=%s)", paste(aThresholds, collapse = ",")),
  exact$probA,
  fit$probA,
  toleranceProb
) ||
  anyFailure

if (anyFailure) {
  cat("FAIL: BCF a-glue posterior deviates from exact\n")
  quit(status = 1L)
}
cat("OK: prior-dominated BCF a posterior matches exact\n")
