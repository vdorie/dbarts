#!/usr/bin/env Rscript

# Exact-posterior gate for heteroscedastic BART (docs/design/heteroscedastic.md).
# Two parts, both matching the sampler's posterior to a computable reference to
# Monte Carlo error - never widen a tolerance to pass, a failure is a real model
# error.
#
#   (a) m' = 1 homoscedastic reduction. A constant predictor stops both forests
#       from splitting, so a single-tree variance forest is one chi^-2 leaf whose
#       posterior IS the homoscedastic sigma^2 posterior (the Section-3.4
#       calibration collapses to (nu, scale) exactly). The heteroscedastic s^2
#       posterior mean must match a plain homoscedastic sigma^2 posterior mean.
#
#   (b) m' = 2 CLOSING exact gate - the definitive check of the multiplicative
#       roll. Two variance trees on ONE binary predictor (two cells) plus a mean
#       forest, sigma fixed at 1 (the variance forest IS the variance). The joint
#       tree space (each of the two variance trees and the mean tree a root or a
#       two-cell split) is enumerated; per combination the variance leaves are
#       integrated by quadrature over their PRODUCT s^2(cell) = h_1 h_2 and the
#       mean leaves analytically (Gaussian-conjugate given the covering variance).
#       The sampler's posterior mean of the IDENTIFIED s^2(cell) and f(cell) match
#       the quadrature. A self-consistently wrong divisor (s^2 instead of s^2_{-j})
#       fails here.
#
# Usage: Rscript heteroscedastic-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

anyFailure <- FALSE
report <- function(label, sampler, exact, tol) {
  gap <- max(abs(sampler - exact))
  fail <- gap > tol
  cat(sprintf(
    "%-26s gap %.4f  (tol %.3f)%s\n",
    label,
    gap,
    tol,
    if (fail) "  <- FAIL" else ""
  ))
  if (fail) anyFailure <<- TRUE
}

# ---------------------------------------------------------------------------
# Part (a): m' = 1 homoscedastic reduction
# ---------------------------------------------------------------------------

partA <- function() {
  ndpost <- if (quick) 4000L else 20000L
  nburn <- 2000L
  tol <- if (quick) 0.08 else 0.05

  set.seed(51L)
  n <- 300L
  x <- matrix(0.5, n, 1L) # constant: neither forest can split
  y <- 2 + 1.5 * rnorm(n)

  ctl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    rngSeed = 51L
  )

  het <- dbarts(x, y, control = ctl, variance = TRUE, n.trees.variance = 1L)
  rHet <- het$run(nburn, ndpost)
  s2Post <- mean(rHet$variance)

  homo <- dbarts(x, y, control = ctl)
  rHomo <- homo$run(nburn, ndpost)
  sigma2Post <- mean(rHomo$sigma^2)

  cat("Part (a) m'=1 reduction:\n")
  cat(sprintf(
    "  heteroscedastic E[s^2] %.4f   homoscedastic E[sigma^2] %.4f\n",
    s2Post,
    sigma2Post
  ))
  report("  (a) m'=1 reduction", s2Post, sigma2Post, tol)
}

# ---------------------------------------------------------------------------
# Part (b): m' = 2 closing exact gate
# ---------------------------------------------------------------------------

partB <- function() {
  ndpost <- if (quick) 30000L else 120000L
  nburn <- 5000L
  nSeeds <- if (quick) 2L else 4L
  # tolerances bound MC + quadrature error a several-fold above the observed
  # gaps (full-mode s^2 gap ~0.003, f gap ~0.001; quick-mode s^2 gap ~0.025 at
  # the shorter run). NEVER widen to pass - a mismatch signals a wrong roll.
  s2Tol <- if (quick) 0.11 else 0.045
  fTol <- if (quick) 0.035 else 0.015

  set.seed(707L)
  base <- 0.6
  power <- 2.0
  sigdf <- 3.0
  sigquant <- 0.9
  nPerCell <- 30L
  cell <- rep(0:1, each = nPerCell)
  n <- length(cell)
  x <- matrix(as.double(cell), ncol = 1L)
  # a distinct mean and variance per cell
  fTrue <- c(-0.4, 0.6)
  sTrue <- c(0.5, 1.3)
  y <- fTrue[cell + 1L] + sTrue[cell + 1L] * rnorm(n)

  # ---- reproduce the sampler's transform and priors on the working scale ----
  yRange <- max(y) - min(y)
  z <- (y - min(y)) / yRange - 0.5
  zA <- z[cell == 0L]
  zB <- z[cell == 1L]
  nA <- length(zA)
  nB <- length(zB)
  sumZA <- sum(zA)
  sumZB <- sum(zB)
  sumZ2A <- sum(zA * zA)
  sumZ2B <- sum(zB * zB)

  # mean leaf prior N(0, tauf2), one mean tree: scale = nodeScale / sqrt(1)
  nodeScale <- 0.5
  kLeaf <- 2.0
  tauf2 <- (nodeScale / kLeaf)^2

  # variance leaf chi^-2(nu', lambda'^2) on the working scale (m' = 2). The
  # sampler anchors the level at the resolved sigest; extract it from the fit.
  sigmaRawScale <- qchisq(1 - sigquant, sigdf) / sigdf

  # ---- exact posterior by joint tree enumeration + variance quadrature ----
  # log chi^-2(h; nu, tau2) prior density
  logChiInv <- function(h, nu, tau2) {
    0.5 *
      nu *
      log(0.5 * nu * tau2) -
      lgamma(0.5 * nu) -
      (0.5 * nu + 1) * log(h) -
      0.5 * nu * tau2 / h
  }
  # full mean-leaf log marginal of one cell (n iid N(mu, s2), mu ~ N(0, tauf2)),
  # keeping the s2-dependent normalizer since s2 varies over the integral
  meanMargCell <- function(s2, nn, sumZ, sumZ2) {
    A <- nn / s2 + 1 / tauf2
    b <- sumZ / s2
    -0.5 *
      nn *
      log(2 * pi * s2) -
      0.5 * log(tauf2 * A) -
      0.5 * sumZ2 / s2 +
      0.5 * b * b / A
  }
  # mean root: a shared mu couples the two cells
  meanMargRoot <- function(s2A, s2B) {
    A <- nA / s2A + nB / s2B + 1 / tauf2
    b <- sumZA / s2A + sumZB / s2B
    -0.5 *
      nA *
      log(2 * pi * s2A) -
      0.5 * nB * log(2 * pi * s2B) -
      0.5 * log(tauf2 * A) -
      0.5 * (sumZ2A / s2A + sumZ2B / s2B) +
      0.5 * b * b / A
  }

  computeExact <- function(sigest) {
    initialSigma <- sigest / yRange
    priorScale <- initialSigma^2 * sigmaRawScale
    leafScale <- priorScale^(1 / 2) # lambda'^2, m' = 2
    leafDf <- 2 / (1 - sqrt(1 - 2 / sigdf))

    # a log-spaced quadrature grid for one chi^-2 variance factor; 40 nodes puts
    # the discretization error well below the sampler MCse (a both-split combo is
    # 4-D, so the node count trades against nNode^4 grid size)
    nNode <- if (quick) 36L else 44L
    center <- log(leafDf * leafScale / (leafDf))
    lo <- center - 9
    hi <- center + 12
    t <- seq(lo, hi, length.out = nNode)
    dt <- t[2L] - t[1L]
    h <- exp(t)
    logW <- logChiInv(h, leafDf, leafScale) + t + log(dt) # p(h) dh, h = e^t

    # each variance tree is a root (both cells share h, prob 1 - base) or a
    # two-cell split (independent h per cell, prob base). Enumerate the joint
    # (var1, var2, mean) structures and integrate. Numerators accumulate in
    # linear space (relative to each combo's max) since f(cell) can be negative.
    comboLogScale <- c() # treePrior + m0 per combo
    comboZ <- c() # marginal integral (linear, post-shift)
    comboS2A <- c()
    comboS2B <- c()
    comboFA <- c()
    comboFB <- c()

    for (r1 in c(FALSE, TRUE)) {
      for (r2 in c(FALSE, TRUE)) {
        for (mroot in c(FALSE, TRUE)) {
          # number of root variance trees determines the shared factor
          roots <- c(r1, r2)
          treePrior <- sum(ifelse(
            roots,
            log(1 - base),
            log(base)
          )) +
            ifelse(mroot, log(1 - base), log(base))

          # build (s2A, s2B, weight) as vectors over the variance grid. Root
          # trees share one h across cells; split trees draw independently.
          # We accumulate on the log scale via a small grid product.
          # index sets: shared roots vary once; split trees vary per cell.
          gridList <- list()
          # represent the joint over needed h-dims as a data.frame of columns
          # then form s2A, s2B and logw
          dims <- list()
          nmShared <- integer(0)
          nmA <- integer(0)
          nmB <- integer(0)
          treeIsRoot <- roots
          # assign a grid column per required h
          colIndex <- 0L
          treeCols <- vector("list", 2L)
          for (j in 1:2) {
            if (treeIsRoot[j]) {
              colIndex <- colIndex + 1L
              treeCols[[j]] <- list(shared = colIndex)
            } else {
              colIndex <- colIndex + 1L
              aCol <- colIndex
              colIndex <- colIndex + 1L
              bCol <- colIndex
              treeCols[[j]] <- list(a = aCol, b = bCol)
            }
          }
          nCols <- colIndex
          grids <- rep(list(seq_len(nNode)), nCols)
          idx <- as.matrix(expand.grid(grids))
          # log weight = sum of node log-weights over every column
          lw <- rowSums(matrix(logW[idx], nrow = nrow(idx)))
          s2A <- rep(1, nrow(idx))
          s2B <- rep(1, nrow(idx))
          for (j in 1:2) {
            tc <- treeCols[[j]]
            if (!is.null(tc$shared)) {
              hv <- h[idx[, tc$shared]]
              s2A <- s2A * hv
              s2B <- s2B * hv
            } else {
              s2A <- s2A * h[idx[, tc$a]]
              s2B <- s2B * h[idx[, tc$b]]
            }
          }
          if (mroot) {
            mm <- meanMargRoot(s2A, s2B)
          } else {
            mm <- meanMargCell(s2A, nA, sumZA, sumZ2A) +
              meanMargCell(s2B, nB, sumZB, sumZ2B)
          }
          integrand <- lw + mm
          m0 <- max(integrand)
          w <- exp(integrand - m0) # linear integrand, shifted by m0
          # posterior mean of f(cell): E[mu_cell | s2] = b/A. Mean root shares mu.
          if (mroot) {
            A <- nA / s2A + nB / s2B + 1 / tauf2
            b <- sumZA / s2A + sumZB / s2B
            fA <- b / A
            fB <- fA
          } else {
            fA <- (sumZA / s2A) / (nA / s2A + 1 / tauf2)
            fB <- (sumZB / s2B) / (nB / s2B + 1 / tauf2)
          }
          comboLogScale <- c(comboLogScale, treePrior + m0)
          comboZ <- c(comboZ, sum(w))
          comboS2A <- c(comboS2A, sum(s2A * w))
          comboS2B <- c(comboS2B, sum(s2B * w))
          comboFA <- c(comboFA, sum(fA * w))
          comboFB <- c(comboFB, sum(fB * w))
        }
      }
    }

    gmax <- max(comboLogScale + log(comboZ))
    wt <- exp(comboLogScale - gmax) # bounded per-combo weight
    den <- sum(wt * comboZ)
    s2ExactA <- sum(wt * comboS2A) / den
    s2ExactB <- sum(wt * comboS2B) / den
    fExactAWorking <- sum(wt * comboFA) / den
    fExactBWorking <- sum(wt * comboFB) / den

    # map back to the original scale: s^2 by yRange^2, f by yRange (+ shift)
    list(
      s2A = s2ExactA * yRange^2,
      s2B = s2ExactB * yRange^2,
      fA = (fExactAWorking + 0.5) * yRange + min(y),
      fB = (fExactBWorking + 0.5) * yRange + min(y)
    )
  }

  # ---- engine arm ----
  ctl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    updateState = FALSE
  )

  fitSeed <- function(seed) {
    set.seed(seed)
    sampler <- dbarts(
      x,
      y,
      control = ctl,
      tree.prior = cgm(power, base),
      resid.prior = chisq(sigdf, sigquant),
      variance = TRUE,
      n.trees.variance = 2L,
      power.variance = power,
      base.variance = base
    )
    sigest <- sampler$data@sigma
    r <- sampler$run(nburn, ndpost)
    s2A <- mean(r$variance[cell == 0L, ])
    s2B <- mean(r$variance[cell == 1L, ])
    fA <- mean(r$train[cell == 0L, ])
    fB <- mean(r$train[cell == 1L, ])
    c(s2A = s2A, s2B = s2B, fA = fA, fB = fB, sigest = sigest)
  }

  fits <- do.call(rbind, lapply(seq_len(nSeeds), fitSeed))
  fit <- colMeans(fits)
  exact <- computeExact(fit["sigest"])

  cat("\nPart (b) m'=2 closing exact gate:\n")
  cat(sprintf(
    "  s^2(A) exact %.4f sampler %.4f | s^2(B) exact %.4f sampler %.4f\n",
    exact$s2A,
    fit["s2A"],
    exact$s2B,
    fit["s2B"]
  ))
  cat(sprintf(
    "  f(A)   exact %.4f sampler %.4f | f(B)   exact %.4f sampler %.4f\n",
    exact$fA,
    fit["fA"],
    exact$fB,
    fit["fB"]
  ))
  report(
    "  (b) m'=2 s^2(cell)",
    c(fit["s2A"], fit["s2B"]),
    c(exact$s2A, exact$s2B),
    s2Tol
  )
  report(
    "  (b) m'=2 f(cell)",
    c(fit["fA"], fit["fB"]),
    c(exact$fA, exact$fB),
    fTol
  )
}

partA()
cat("\n")
partB()
cat("\n")

if (anyFailure) {
  cat("FAIL: heteroscedastic sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("OK: heteroscedastic sampler matches the exact posterior on both parts\n")
