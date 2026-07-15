#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore multinomial (softmax) sampler
# (docs/design/multinomial.md). Three arms, each matching the sampler's
# posterior mean of the IDENTIFIED softmax probabilities to a closed-form
# quadrature, to Monte Carlo error. A failure means the softmax coupling, the
# interleaved one-vs-rest Polya-Gamma draw, the level-centering move, or the
# log-sum-exp margin formation is wrong - fix the model, never loosen the gate.
#
# The sampler's per-forest total-fit prior is symmetric N(0, tau^2)^K with
#   tau = nodeScale / k,  nodeScale = pi*sqrt(3)/sqrt(2),  k = 2
# (the K = 2 pairwise-log-odds anchor: f_ik - f_ij has sd sqrt(2)*tau =
# pi*sqrt(3)/2, matching the shipped logistic single-forest calibration). The
# raw f_ik are non-identified (softmax is invariant to a per-observation common
# shift), so every arm compares the identified softmax probabilities, never the
# raw f_ik, and every quadrature marginalizes the level analytically rather than
# fixing a reference forest to zero.
#
#   (1) intercept-only K = 3: a constant predictor gives root-only trees; the
#       induced prior on the differences d_k = f_k - f_K is the CORRELATED
#       Gaussian N(0, tau^2 (I + 11')), integrated by 2-D quadrature.
#   (2) K = 2 == logistic: the two-category multinomial matches the shipped
#       logistic sampler DISTRIBUTIONALLY (different rng consumption, NOT
#       draw-for-draw) - the strongest internal-consistency check.
#   (3) covariate-dependent K = 3: one binary categorical predictor (two cells),
#       one tree per forest, the joint tree space enumerated and integrated by
#       per-cell quadrature over the difference-space Gaussian - the only exact
#       gate on tree growth under softmax.
#
# Usage: Rscript multinomial-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

tau <- pi * sqrt(3) / sqrt(2) / 2 # per-forest total-fit prior sd
tau2 <- tau * tau

# Each arm's tolerance is a fixed constant bounding sampler MC error plus
# quadrature discretization, sized several-fold above the cross-seed spread
# observed when the gate was recorded (full-mode max gaps 0.0000 / 0.0008 /
# 0.0012 against tolerances 0.008 / 0.008 / 0.015). Never widen a tolerance to
# make a failing arm pass - a failure means the model is wrong.
anyFailure <- FALSE
report <- function(name, gap, tol) {
  fail <- gap > tol
  cat(sprintf(
    "%-22s max gap %.4f  (tol %.3f)%s\n",
    name,
    gap,
    tol,
    if (fail) "  <- FAIL" else ""
  ))
  if (fail) anyFailure <<- TRUE
}

# ---------------------------------------------------------------------------
# Arm 1: intercept-only K = 3 softmax
# ---------------------------------------------------------------------------

arm1 <- function() {
  ndpost <- if (quick) 10000L else 40000L
  nburn <- 5000L
  nSeeds <- if (quick) 1L else 3L
  tolerance <- if (quick) 0.012 else 0.008

  set.seed(2026L)
  K <- 3L
  n <- 240L
  # a modest, deliberately imbalanced multinomial sample
  trueP <- c(0.5, 0.3, 0.2)
  labels <- sample.int(K, n, replace = TRUE, prob = trueP) - 1L
  counts <- tabulate(labels + 1L, K)
  # a constant predictor: no valid cut points, so every tree stays a root
  x <- matrix(0.5, n, 1L)

  # exact posterior mean of the identified softmax probabilities. Marginalizing
  # the level, the prior on d = (f0 - f2, f1 - f2) is N(0, tau^2 (I + 11')); the
  # softmax depends only on d (f2 = 0 reference), so 2-D quadrature suffices.
  Sigma <- tau2 * matrix(c(2, 1, 1, 2), 2L, 2L)
  Prec <- solve(Sigma)
  M <- 12
  g <- seq(-M, M, length.out = 301L)
  grid <- as.matrix(expand.grid(d0 = g, d1 = g))
  d0 <- grid[, 1L]
  d1 <- grid[, 2L]
  lse <- log1p(exp(d0) + exp(d1)) # log(exp(d0) + exp(d1) + 1); f2 = 0
  logP <- cbind(d0 - lse, d1 - lse, -lse) # log softmax, categories 0,1,2
  loglik <- logP %*% counts
  q <- Prec[1, 1] * d0^2 + 2 * Prec[1, 2] * d0 * d1 + Prec[2, 2] * d1^2
  w <- as.vector(loglik) - 0.5 * q
  w <- exp(w - max(w))
  exact <- colSums(exp(logP) * w) / sum(w)

  fitSeed <- function(seed) {
    set.seed(seed)
    control <- dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 50L,
      updateState = FALSE
    )
    host <- dbarts(x, as.double(labels), control = control)
    bc <- dbarts:::bartcoreMultinomialSampler(host, labels, K = K)
    r <- dbarts:::bartcoreRun(bc, nburn, ndpost)
    # every observation shares the intercept-only probabilities; average them
    apply(r$train, 2L, mean)
  }
  fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))

  cat("Arm 1 (intercept-only K = 3):\n")
  cat(sprintf("  exact   %s\n", paste(sprintf("%.4f", exact), collapse = " ")))
  cat(sprintf("  sampler %s\n", paste(sprintf("%.4f", fit), collapse = " ")))
  report("  arm1 intercept K=3", max(abs(fit - exact)), tolerance)
}

# ---------------------------------------------------------------------------
# Arm 2: K = 2 multinomial == shipped logistic (distributional)
# ---------------------------------------------------------------------------

arm2 <- function() {
  # Distributional match, NOT draw-for-draw: the two samplers consume the rng
  # differently (K forests + interleaved PG + centering vs one forest). INTERCEPT
  # ONLY (a constant predictor, root-only trees): there the K = 2 multinomial's
  # log-odds is the single scalar f1 - f0 with prior sd sqrt(2)*tau = pi*sqrt(3)/2,
  # exactly the shipped logistic single-forest sd, so the two are the SAME model
  # on P(y = 1) and match to MC error. (With a covariate the multinomial log-odds
  # is a DIFFERENCE of two ensembles, not a single ensemble, so the functional
  # priors differ - covariate K = 2 is covered by arm 3's K = 3 tree-growth gate,
  # not by a logistic equivalence.)
  ndpost <- if (quick) 10000L else 30000L
  nburn <- 3000L
  nSeeds <- if (quick) 3L else 5L
  tolerance <- if (quick) 0.015 else 0.008

  set.seed(3131L)
  n <- 200L
  x <- matrix(0.5, n, 1L) # constant predictor: intercept-only
  y <- rbinom(n, 1L, 0.35)

  control <- function() {
    dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 50L,
      updateState = FALSE
    )
  }

  fitLogistic <- function(seed) {
    set.seed(seed)
    # family = "logistic" installs the pi*sqrt(3) node scale; offset 0 centers
    # the log-odds prior at zero, as the multinomial's is
    host <- dbarts(
      x,
      y,
      offset = 0,
      family = "logistic",
      control = control(),
      verbose = FALSE
    )
    bc <- dbarts:::bartcoreSampler(host, family = "logistic")
    r <- dbarts:::bartcoreRun(bc, nburn, ndpost)
    rowMeans(plogis(r$train)) # r$train is the latent log-odds
  }
  fitMultinomial <- function(seed) {
    set.seed(seed)
    host <- dbarts(x, as.double(y), control = control())
    bc <- dbarts:::bartcoreMultinomialSampler(host, as.integer(y), K = 2L)
    r <- dbarts:::bartcoreRun(bc, nburn, ndpost)
    apply(r$train[, 2L, , drop = FALSE], 1L, mean) # P(y = 1)
  }

  pLog <- rowMeans(vapply(seq_len(nSeeds), fitLogistic, numeric(n)))
  pMult <- rowMeans(vapply(1000L + seq_len(nSeeds), fitMultinomial, numeric(n)))

  cat("Arm 2 (K = 2 multinomial vs logistic, per-observation P(y = 1)):\n")
  cat(sprintf(
    "  mean |logistic - multinomial| = %.4f, max = %.4f\n",
    mean(abs(pLog - pMult)),
    max(abs(pLog - pMult))
  ))
  report("  arm2 K=2 == logistic", max(abs(pLog - pMult)), tolerance)
}

# ---------------------------------------------------------------------------
# Arm 3: covariate-dependent K = 3 (two cells, one tree per forest)
# ---------------------------------------------------------------------------

arm3 <- function() {
  ndpost <- if (quick) 20000L else 60000L
  nburn <- 4000L
  nSeeds <- if (quick) 2L else 3L
  tolerance <- if (quick) 0.02 else 0.015

  set.seed(7007L)
  K <- 3L
  base <- 0.5
  power <- 2.0
  nPerCell <- 40L
  # two cells (a binary categorical predictor); a distinct softmax per cell
  cellProb <- list(c(0.6, 0.25, 0.15), c(0.2, 0.35, 0.45))
  cell <- rep(0:1, each = nPerCell)
  labels <- integer(length(cell))
  for (c in 0:1) {
    idx <- which(cell == c)
    labels[idx] <- sample.int(
      K,
      length(idx),
      replace = TRUE,
      prob = cellProb[[c + 1L]]
    ) -
      1L
  }
  x <- matrix(as.double(cell), ncol = 1L)
  countsA <- tabulate(labels[cell == 0L] + 1L, K)
  countsB <- tabulate(labels[cell == 1L] + 1L, K)

  # ---- exact posterior by joint tree enumeration + nested quadrature ----
  # One tree per forest over a 2-level categorical: each forest is a root
  # (prob 1 - base, one leaf shared by both cells) or splits into the two
  # single-cell leaves (prob base), so 2^K joint tree combinations. Every leaf
  # is N(0, tau^2). Given a combination, the two cells are conditionally
  # INDEPENDENT once the shared (root) forests' leaves are fixed - a root forest
  # ties f_kA = f_kB, so the two cells couple only through it. So marginalize by
  # nesting: an OUTER quadrature over the root leaves, and, for each cell, an
  # INNER quadrature over that cell's split leaves. Both cells share the root
  # grid, so no difference-space degeneracy and no matrix inversion arises.
  nT <- if (quick) 31L else 41L
  M <- 10
  tnode <- seq(-M, M, length.out = nT)
  dt <- tnode[2L] - tnode[1L]
  wnode <- dnorm(tnode, 0, tau) * dt # normalized N(0, tau^2) node weight

  logSoftmax3vec <- function(vals) {
    m <- pmax(vals[, 1L], vals[, 2L], vals[, 3L])
    s <- log(exp(vals[, 1L] - m) + exp(vals[, 2L] - m) + exp(vals[, 3L] - m))
    (vals - m) - s
  }

  # Integrate one cell over its split leaves for every outer root-grid point:
  # returns g[root] = E_split[softmax^counts] and h[root, k] = E_split[softmax_k
  # softmax^counts], the marginal likelihood and its k-th predictive moment.
  cellIntegral <- function(rootIdx, splitIdx, counts) {
    nS <- length(splitIdx)
    nSplit <- nT^nS
    nRoot <- nT^length(rootIdx)
    orderF <- c(splitIdx, rootIdx) # split forests vary fastest
    gr <- as.matrix(expand.grid(rep(list(tnode), K)))
    vals <- matrix(0, nrow(gr), K)
    for (j in seq_along(orderF)) {
      vals[, orderF[j]] <- gr[, j]
    }
    wSplit <- rep(1, nrow(gr))
    for (j in seq_len(nS)) {
      wSplit <- wSplit * (dnorm(gr[, j], 0, tau) * dt)
    }
    ls <- logSoftmax3vec(vals)
    sm <- exp(ls)
    lik <- exp(as.vector(ls %*% counts))
    wl <- wSplit * lik
    g <- colSums(matrix(wl, nrow = nSplit)) # length nRoot
    h <- vapply(
      seq_len(K),
      function(k) colSums(matrix(wl * sm[, k], nrow = nSplit)),
      numeric(nRoot)
    )
    if (nRoot == 1L) {
      h <- matrix(h, 1L, K)
    }
    list(g = g, h = h)
  }

  combos <- as.matrix(expand.grid(
    r0 = c(FALSE, TRUE),
    r1 = c(FALSE, TRUE),
    r2 = c(FALSE, TRUE)
  ))
  logDen <- numeric(nrow(combos))
  logNumA <- matrix(0, nrow(combos), K)
  logNumB <- matrix(0, nrow(combos), K)
  for (ci in seq_len(nrow(combos))) {
    roots <- combos[ci, ]
    rootIdx <- which(roots)
    splitIdx <- which(!roots)
    nR <- length(rootIdx)
    if (nR == 0L) {
      wRoot <- 1
    } else {
      rg <- as.matrix(expand.grid(rep(list(wnode), nR)))
      wRoot <- rep(1, nrow(rg))
      for (j in seq_len(nR)) {
        wRoot <- wRoot * rg[, j]
      }
    }
    cA <- cellIntegral(rootIdx, splitIdx, countsA)
    cB <- cellIntegral(rootIdx, splitIdx, countsB)
    Zc <- sum(wRoot * cA$g * cB$g)
    predA <- colSums(cA$h * (wRoot * cB$g))
    predB <- colSums(cB$h * (wRoot * cA$g))
    comboLogPrior <- sum(ifelse(roots, log(1 - base), log(base)))
    logDen[ci] <- comboLogPrior + log(Zc)
    logNumA[ci, ] <- comboLogPrior + log(predA)
    logNumB[ci, ] <- comboLogPrior + log(predB)
  }
  gmax <- max(logDen)
  den <- sum(exp(logDen - gmax))
  exactA <- colSums(exp(logNumA - gmax)) / den
  exactB <- colSums(exp(logNumB - gmax)) / den

  fitSeed <- function(seed) {
    set.seed(seed)
    control <- dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 1L,
      updateState = FALSE
    )
    host <- dbarts(
      x,
      as.double(labels),
      control = control,
      tree.prior = cgm(power, base)
    )
    host$data@varTypes[1L] <- 1L # mark the predictor categorical
    bc <- dbarts:::bartcoreMultinomialSampler(host, labels, K = K)
    r <- dbarts:::bartcoreRun(bc, nburn, ndpost)
    pA <- apply(r$train[cell == 0L, , , drop = FALSE], 2L, mean)
    pB <- apply(r$train[cell == 1L, , , drop = FALSE], 2L, mean)
    c(pA, pB)
  }
  fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
  fitA <- fit[1:K]
  fitB <- fit[(K + 1L):(2L * K)]

  cat("Arm 3 (covariate K = 3, cell A then cell B):\n")
  cat(sprintf(
    "  exact A   %s\n",
    paste(sprintf("%.4f", exactA), collapse = " ")
  ))
  cat(sprintf("  sampler A %s\n", paste(sprintf("%.4f", fitA), collapse = " ")))
  cat(sprintf(
    "  exact B   %s\n",
    paste(sprintf("%.4f", exactB), collapse = " ")
  ))
  cat(sprintf("  sampler B %s\n", paste(sprintf("%.4f", fitB), collapse = " ")))
  report(
    "  arm3 covariate K=3",
    max(abs(c(fitA - exactA, fitB - exactB))),
    tolerance
  )
}

arm1()
cat("\n")
arm2()
cat("\n")
arm3()
cat("\n")

if (anyFailure) {
  quit(status = 1L)
}
cat("OK: multinomial sampler matches the exact posterior on all three arms\n")
