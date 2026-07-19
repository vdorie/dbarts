#!/usr/bin/env Rscript

# Exact-posterior gate for the monotone (mBART) constrained constant leaf
# (docs/design/monotone.md section 9). Two parts, each matched to a direct
# constrained quadrature to MC error; failure means the constrained draw, the
# B' structure-move marginal, the neighbor geometry, or the c-inflation is
# wrong. Both parts fit a single-tree, one-constrained-predictor gaussian model
# with fixed sigma, so the constrained integrals are 1-D/2-D (part a) or a 3-D
# monotone-cone integral over three ordered leaves (part b).
#
# Usage: Rscript monotone-reference.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
zBound <- 4.0

cScale <- sqrt(pi / (pi - 1.0)) # the c-inflation constant (design section 6)

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}
batchMeanSE <- function(v, nBatches = 200L) {
  len <- (length(v) %/% nBatches) * nBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
  sd(bm) / sqrt(nBatches)
}

# One leaf's unconstrained gaussian posterior N(mean, sd^2) from its scaled
# sufficient statistic (sum z, count) under a leaf prior sd priorSd.
leafPosterior <- function(sumZ, n, residVar, priorSd) {
  priorPrec <- 1 / priorSd^2
  postPrec <- n / residVar
  variance <- 1 / (priorPrec + postPrec)
  list(mean = (sumZ / residVar) * variance, sd = sqrt(variance))
}

# ---- shared design helpers -------------------------------------------------

runSampler <- function(
  y,
  x,
  testX,
  nodeScale,
  priorSigma,
  nBurn,
  nKept,
  batch,
  seed
) {
  ctl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    keepTrees = TRUE,
    n.samples = batch,
    updateState = TRUE,
    rngSeed = seed,
    n.cuts = length(unique(x)) - 1L
  )
  monotone <- c(1L)
  names(monotone) <- colnames(x)
  sampler <- dbarts(
    x,
    y,
    test = testX,
    control = ctl,
    tree.prior = cgm(power, base),
    node.prior = normal(k),
    # fixed() takes the residual VARIANCE on the original scale, so the engine's
    # internal (range-scaled) residual variance is priorSigma^2 / yRange^2 =
    # (priorSigma / yRange)^2, matching residVar in the quadrature
    resid.prior = fixed(priorSigma^2),
    monotone = monotone
  )
  stopifnot(sampler$model@p.birth_death == 1)
  stopifnot(is.null(sampler$data@offset))

  nBatch <- as.integer(ceiling(nKept / batch))
  leaves <- integer(0L)
  fits <- NULL
  first <- TRUE
  for (b in seq_len(nBatch)) {
    r <- if (first) sampler$run(nBurn, batch) else sampler$run(0L, batch)
    first <- FALSE
    tr <- sampler$getTrees()
    leaves <- c(
      leaves,
      as.integer(tapply(tr$var, tr$sample, function(v) sum(v == -1L)))
    )
    fits <- cbind(fits, r$test)
  }
  list(leaves = leaves, fitsScaled = (fits - min(y)) / (max(y) - min(y)) - 0.5)
}

base <- 0.5
power <- 2
k <- 2

anyFailure <- FALSE
report <- function(label, engineValue, exactValue, se) {
  z <- (engineValue - exactValue) / se
  failed <- abs(z) > zBound
  anyFailure <<- anyFailure || failed
  cat(sprintf(
    "%-28s engine %9.4f  exact %9.4f  MCse %7.4f  z %6.2f%s\n",
    label,
    engineValue,
    exactValue,
    se,
    z,
    if (failed) "  <- FAIL" else ""
  ))
}

# ---- part (a): one-cut enumerable gate (root vs single ordered split) -------

cat("part (a): one-cut structure posterior and constrained fit\n")
{
  ndpost <- if (quick) 20000L else 120000L
  batch <- 2000L
  nBurn <- 2000L

  set.seed(207L)
  nPer <- 20L
  cell <- rep(1:2, each = nPer)
  x <- matrix(as.double(cell), ncol = 1L, dimnames = list(NULL, "x1"))
  # a mild increasing signal: the +1 constraint agrees with the data, so the
  # cone probability stays well-conditioned (away from the deep tail) and the
  # structure posterior spreads across root and split
  muCell <- c(0.12, 0.32)
  noiseSd <- 0.4
  y <- muCell[cell] + rnorm(length(cell), sd = noiseSd)

  yRange <- max(y) - min(y)
  z <- (y - min(y)) / yRange - 0.5
  residVar <- (noiseSd / yRange)^2
  nodeScale <- 0.5
  priorSd <- nodeScale / k
  priorSdC <- cScale * priorSd

  sumZ <- tapply(z, cell, sum)
  nCell <- as.integer(table(cell))

  gLo <- -3
  gHi <- 3
  d <- if (quick) 0.001 else 0.0005
  grid <- seq(gLo, gHi, by = d)
  reduced <- function(mu, sumZ, n) {
    exp(sumZ * mu / residVar - 0.5 * n * mu^2 / residVar)
  }

  # root: one unconstrained leaf over both cells
  f0 <- reduced(grid, sum(sumZ), sum(nCell)) * dnorm(grid, 0, priorSd)
  logMLroot <- log(sum(f0) * d)
  meanRoot <- sum(grid * f0) / sum(f0)

  # split: {mu_L <= mu_R}, both c-inflated; the constrained prior normalizer is
  # P(X <= Y) = 1/2 for X, Y iid N(0, priorSdC^2)
  fL <- reduced(grid, sumZ[[1L]], nCell[1L]) * dnorm(grid, 0, priorSdC)
  fR <- reduced(grid, sumZ[[2L]], nCell[2L]) * dnorm(grid, 0, priorSdC)
  cumL <- cumsum(fL)
  cumGL <- cumsum(grid * fL)
  coneZ <- sum(fR * cumL)
  logMLsplit <- log(2 * coneZ * d^2)
  meanL <- sum(fR * cumGL) / coneZ
  meanR <- sum(grid * fR * cumL) / coneZ

  # the CGM structure prior over {root, split}: with a single available cut the
  # two children admit no further split, so their depth-1 growth is forced to 0
  # (each a leaf w.p. 1) and the split's prior is just base
  priorRoot <- log(1 - base)
  priorSplit <- log(base)
  logW <- c(priorRoot + logMLroot, priorSplit + logMLsplit)
  w <- exp(logW - logSumExp(logW))
  pSplitExact <- w[2L]
  fitExact <- c(
    w[1L] * meanRoot + w[2L] * meanL,
    w[1L] * meanRoot + w[2L] * meanR
  )

  runs <- lapply(seq_len(if (quick) 1L else 2L), function(s) {
    runSampler(
      y,
      x,
      x[c(1L, nPer + 1L), , drop = FALSE],
      nodeScale,
      noiseSd,
      nBurn,
      ndpost,
      batch,
      300L + s
    )
  })
  leaves <- do.call(c, lapply(runs, `[[`, "leaves"))
  fitsScaled <- do.call(cbind, lapply(runs, `[[`, "fitsScaled"))

  report(
    "P(split)",
    mean(leaves == 2L),
    pSplitExact,
    batchMeanSE(as.double(leaves == 2L))
  )
  report(
    "fit[bin 1]",
    mean(fitsScaled[1L, ]),
    fitExact[1L],
    batchMeanSE(fitsScaled[1L, ])
  )
  report(
    "fit[bin 2]",
    mean(fitsScaled[2L, ]),
    fitExact[2L],
    batchMeanSE(fitsScaled[2L, ])
  )
}

# ---- part (b): fixed 3-leaf double-bounded interior leaf posterior ----------

cat("\npart (b): 3-leaf monotone-cone leaf posterior (means and covariances)\n")
{
  ndpost <- if (quick) 40000L else 200000L
  batch <- 2000L
  nBurn <- 4000L

  set.seed(913L)
  nPer <- 25L
  cell <- rep(1:3, each = nPer)
  x <- matrix(as.double(cell), ncol = 1L, dimnames = list(NULL, "x1"))
  # separated enough that 3 leaves carry substantial posterior mass, but with a
  # non-monotone middle so the interior leaf is genuinely double-bounded
  muCell <- c(0.25, 0.1, 0.6)
  noiseSd <- 0.35
  y <- muCell[cell] + rnorm(length(cell), sd = noiseSd)

  yRange <- max(y) - min(y)
  z <- (y - min(y)) / yRange - 0.5
  residVar <- (noiseSd / yRange)^2
  nodeScale <- 0.5
  priorSdC <- cScale * nodeScale / k # all three leaves are constrained

  sumZ <- tapply(z, cell, sum)
  nCell <- as.integer(table(cell))
  post <- lapply(1:3, function(b) {
    leafPosterior(sumZ[[b]], nCell[b], residVar, priorSdC)
  })

  # 3-D constrained posterior moments over {mu1 <= mu2 <= mu3} for independent
  # per-leaf normals, by cumulative sums (O(grid)): each leaf's unconstrained
  # posterior times the ordering indicator. The interior mean's O(d) cumsum
  # discretization must sit well below the sampler MCse (~2e-4): at d = 0.0015 it
  # biased E[mu2] by ~3e-4 and spuriously inflated its z to ~3.2; 1e-4 converges
  # it (E[mu2] stable past the fifth decimal), so the z reflects sampler-vs-truth
  # agreement, not grid error.
  d <- if (quick) 0.001 else 1e-04
  grid <- seq(-2, 2.5, by = d)
  p1 <- dnorm(grid, post[[1L]]$mean, post[[1L]]$sd)
  p2 <- dnorm(grid, post[[2L]]$mean, post[[2L]]$sd)
  p3 <- dnorm(grid, post[[3L]]$mean, post[[3L]]$sd)
  # E[mu1^a mu2^b mu3^c] over the cone {mu1<=mu2<=mu3} by nested cumulative sums:
  #   sum_k g_k^c p3_k * cumsum_j( g_j^b p2_j * cumsum_i(g_i^a p1_i) )_k
  coneMoment <- function(a, b, cc) {
    cA <- cumsum(grid^a * p1)
    Bk <- cumsum(grid^b * p2 * cA)
    sum(grid^cc * p3 * Bk) * d^3
  }
  Z <- coneMoment(0, 0, 0)
  m1 <- coneMoment(1, 0, 0) / Z
  m2 <- coneMoment(0, 1, 0) / Z
  m3 <- coneMoment(0, 0, 1) / Z
  m11 <- coneMoment(2, 0, 0) / Z
  m22 <- coneMoment(0, 2, 0) / Z
  m33 <- coneMoment(0, 0, 2) / Z
  m12 <- coneMoment(1, 1, 0) / Z
  m13 <- coneMoment(1, 0, 1) / Z
  m23 <- coneMoment(0, 1, 1) / Z
  meanExact <- c(m1, m2, m3)
  covExact <- c(
    m11 - m1^2,
    m22 - m2^2,
    m33 - m3^2,
    m12 - m1 * m2,
    m13 - m1 * m3,
    m23 - m2 * m3
  )

  runs <- lapply(seq_len(if (quick) 1L else 3L), function(s) {
    runSampler(
      y,
      x,
      x[c(1L, nPer + 1L, 2L * nPer + 1L), , drop = FALSE],
      nodeScale,
      noiseSd,
      nBurn,
      ndpost,
      batch,
      900L + s
    )
  })
  leaves <- do.call(c, lapply(runs, `[[`, "leaves"))
  fitsScaled <- do.call(cbind, lapply(runs, `[[`, "fitsScaled"))
  keep <- which(leaves == 3L)
  cat(sprintf(
    "3-leaf sweeps: %d of %d (%.1f%%)\n",
    length(keep),
    length(leaves),
    100 * length(keep) / length(leaves)
  ))
  mu <- fitsScaled[, keep, drop = FALSE]

  report("E[mu1]", mean(mu[1L, ]), meanExact[1L], batchMeanSE(mu[1L, ]))
  report(
    "E[mu2] (interior)",
    mean(mu[2L, ]),
    meanExact[2L],
    batchMeanSE(mu[2L, ])
  )
  report("E[mu3]", mean(mu[3L, ]), meanExact[3L], batchMeanSE(mu[3L, ]))
  report(
    "Var[mu1]",
    var(mu[1L, ]),
    covExact[1L],
    batchMeanSE((mu[1L, ] - mean(mu[1L, ]))^2)
  )
  report(
    "Var[mu2]",
    var(mu[2L, ]),
    covExact[2L],
    batchMeanSE((mu[2L, ] - mean(mu[2L, ]))^2)
  )
  report(
    "Var[mu3]",
    var(mu[3L, ]),
    covExact[3L],
    batchMeanSE((mu[3L, ] - mean(mu[3L, ]))^2)
  )
  cp <- (mu[1L, ] - mean(mu[1L, ])) * (mu[2L, ] - mean(mu[2L, ]))
  report("Cov[mu1,mu2]", mean(cp), covExact[4L], batchMeanSE(cp))
  cp <- (mu[2L, ] - mean(mu[2L, ])) * (mu[3L, ] - mean(mu[3L, ]))
  report("Cov[mu2,mu3]", mean(cp), covExact[6L], batchMeanSE(cp))
}

if (anyFailure) {
  quit(status = 1L)
}
cat("\nOK: monotone sampler matches the constrained posterior (both parts)\n")
