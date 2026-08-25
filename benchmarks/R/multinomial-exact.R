#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore multinomial (softmax) sampler
# (docs/design/multinomial.md). Arms 1-6 each match the sampler's posterior
# mean of the IDENTIFIED softmax probabilities to a closed-form quadrature, to
# Monte Carlo error; arm 7 gates the one quantity they cannot see, the
# NON-identified level the centering move draws. A failure means the softmax
# coupling, the interleaved one-vs-rest Polya-Gamma draw, the level-centering
# move, or the log-sum-exp margin formation is wrong - fix the model, never
# loosen the gate.
# Arm 4 reuses arm 3's fit with two test rows duplicating the two cells and
# gates that the TEST channel (combinedTestFits) equals the same quadrature
# target as the train channel - the softmax-invariance correctness fact. Arm 5
# is the grouped-count analog of arm 1 (n_i > 1), the only exact gate on the
# PG(n_i) summing draw and the (y - n_i/2) working response. Arm 6 is arm 1
# under a fixed category offset, the only oracle that catches a consistently
# wrong sign or placement for it: a create-vs-swap parity cannot, since both of
# its arms would make the same error. Arm 7 reads the grand level itself
# against the marginal N(0, tau^2/K) the move is an exact independence sampler
# from at the intercept-only configuration - the only arm that sees the scale
# of the shift rather than the direction it leaves alone.
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

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

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
    r <- bartcoreRun(bc, nburn, ndpost)
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
    r <- bartcoreRun(bc, nburn, ndpost)
    rowMeans(plogis(r$train)) # r$train is the latent log-odds
  }
  fitMultinomial <- function(seed) {
    set.seed(seed)
    host <- dbarts(x, as.double(y), control = control())
    bc <- dbarts:::bartcoreMultinomialSampler(host, as.integer(y), K = 2L)
    r <- bartcoreRun(bc, nburn, ndpost)
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

  # test rows duplicating the two cells: the test channel is the same softmax
  # blend on held-out rows, so a cell's test-set probabilities must match its
  # train/quadrature target - the direct gate on combinedTestFits equalling the
  # (already-gated) combinedFits train blend
  xTest <- matrix(c(0, 1), 2L, 1L)

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
      test = xTest,
      control = control,
      tree.prior = cgm(power, base)
    )
    host$data@varTypes[1L] <- 1L # mark the predictor categorical
    bc <- dbarts:::bartcoreMultinomialSampler(host, labels, K = K)
    r <- bartcoreRun(bc, nburn, ndpost)
    pA <- apply(r$train[cell == 0L, , , drop = FALSE], 2L, mean)
    pB <- apply(r$train[cell == 1L, , , drop = FALSE], 2L, mean)
    tA <- apply(r$test[1L, , , drop = FALSE], 2L, mean)
    tB <- apply(r$test[2L, , , drop = FALSE], 2L, mean)
    c(pA, pB, tA, tB)
  }
  fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
  fitA <- fit[1:K]
  fitB <- fit[(K + 1L):(2L * K)]
  testA <- fit[(2L * K + 1L):(3L * K)]
  testB <- fit[(3L * K + 1L):(4L * K)]

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

  cat("Arm 4 (test-at-creation K = 3, held-out cell A then cell B):\n")
  cat(sprintf(
    "  test A    %s\n",
    paste(sprintf("%.4f", testA), collapse = " ")
  ))
  cat(sprintf(
    "  test B    %s\n",
    paste(sprintf("%.4f", testB), collapse = " ")
  ))
  report(
    "  arm4 test K=3",
    max(abs(c(testA - exactA, testB - exactB))),
    tolerance
  )
}

# ---------------------------------------------------------------------------
# Arm 5: intercept-only K = 3 GROUPED COUNTS (n_i > 1)
# ---------------------------------------------------------------------------

armCount <- function() {
  # The count path's only exact gate: the PG(n_i) summing draw and the
  # (y - n_i/2) working response. Intercept-only, so every row shares the softmax
  # p and the joint likelihood collapses to multinomial(sum_i n_i, p) on the
  # COLUMN SUMS - the same correlated-difference-Gaussian quadrature arm 1 uses,
  # with the aggregate counts as its sufficient statistic. (K = 2 counts would
  # reduce to binomial(n_i, p); this K = 3 arm gates the softmax coupling.)
  ndpost <- if (quick) 10000L else 40000L
  nburn <- 5000L
  nSeeds <- if (quick) 1L else 3L
  tolerance <- if (quick) 0.012 else 0.008

  set.seed(4242L)
  K <- 3L
  n <- 120L
  trueP <- c(0.5, 0.3, 0.2)
  # grouped rows, n_i in 2..6 trials each: a genuine multi-trial count matrix
  trials <- sample(2:6, n, replace = TRUE)
  counts <- t(vapply(
    seq_len(n),
    function(i) rmultinom(1L, trials[i], trueP)[, 1L],
    integer(K)
  ))
  totalCounts <- colSums(counts)
  x <- matrix(0.5, n, 1L) # constant predictor: no cut points, root-only trees

  # exact posterior mean of the identified softmax probabilities, the aggregate-
  # count analog of arm 1: the prior on d = (f0 - f2, f1 - f2) is N(0,
  # tau^2 (I + 11')), the softmax depends only on d, so 2-D quadrature over the
  # total counts suffices.
  Sigma <- tau2 * matrix(c(2, 1, 1, 2), 2L, 2L)
  Prec <- solve(Sigma)
  M <- 12
  g <- seq(-M, M, length.out = 301L)
  grid <- as.matrix(expand.grid(d0 = g, d1 = g))
  d0 <- grid[, 1L]
  d1 <- grid[, 2L]
  lse <- log1p(exp(d0) + exp(d1)) # log(exp(d0) + exp(d1) + 1); f2 = 0
  logP <- cbind(d0 - lse, d1 - lse, -lse) # log softmax, categories 0,1,2
  loglik <- logP %*% totalCounts
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
    # the host response is ignored by the multinomial sampler; pass a varying
    # dummy so the gaussian host builds without a degenerate scale
    host <- dbarts(x, as.double(counts[, 1L]), control = control)
    bc <- dbarts:::bartcoreMultinomialCountSampler(host, counts, K = K)
    r <- bartcoreRun(bc, nburn, ndpost)
    # every observation shares the intercept-only probabilities; average them
    apply(r$train, 2L, mean)
  }
  fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))

  cat("Arm 5 (intercept-only K = 3 grouped counts):\n")
  cat(sprintf(
    "  total counts %s\n",
    paste(totalCounts, collapse = " ")
  ))
  cat(sprintf("  exact   %s\n", paste(sprintf("%.4f", exact), collapse = " ")))
  cat(sprintf("  sampler %s\n", paste(sprintf("%.4f", fit), collapse = " ")))
  report("  arm5 counts K=3", max(abs(fit - exact)), tolerance)
}

# ---------------------------------------------------------------------------
# Arm 6: intercept-only K = 3 with a fixed CATEGORY OFFSET
# ---------------------------------------------------------------------------

armOffset <- function() {
  # The category offset's only exact gate, and the only oracle in the arc that
  # can see a consistently wrong SIGN or PLACEMENT: a create-vs-swap parity
  # agrees happily with an offset that is added to the margin but not
  # subtracted from the working response, or vice versa, since both of its arms
  # make the same mistake. Here the offset is IN the likelihood the quadrature
  # integrates, so any of those errors moves the sampler off the target.
  #
  # Intercept-only (a constant predictor, root-only trees), so the K forests
  # carry one level each and the identified quantity is the difference vector
  # d = (f0 - f2, f1 - f2) with the same correlated N(0, tau^2 (I + 11')) prior
  # arm 1 uses. The offset varies by category and by GROUP - two blocks of rows
  # with different per-category shifts - so the softmax differs between the
  # blocks while the forests do not, which no common-shift null direction can
  # produce. Within a group the likelihood collapses to multinomial on that
  # group's column sums, so the quadrature is arm 1's with one likelihood
  # factor per group, each evaluated at d + delta_g.
  ndpost <- if (quick) 10000L else 40000L
  nburn <- 5000L
  nSeeds <- if (quick) 1L else 3L
  tolerance <- if (quick) 0.012 else 0.008

  set.seed(5150L)
  K <- 3L
  nPerGroup <- 120L
  n <- 2L * nPerGroup
  group <- rep(1:2, each = nPerGroup)
  # per-group, per-category offsets; neither row is constant, so neither lies
  # along the softmax's null direction
  offsetByGroup <- rbind(c(0.0, 0.8, -0.5), c(0.4, -0.6, 0.9))
  offset <- offsetByGroup[group, ]
  trueP <- t(apply(offsetByGroup, 1L, function(o) exp(o) / sum(exp(o))))
  labels <- vapply(
    seq_len(n),
    function(i) sample.int(K, 1L, prob = trueP[group[i], ]) - 1L,
    integer(1L)
  )
  counts <- vapply(
    1:2,
    function(g) tabulate(labels[group == g] + 1L, K),
    numeric(K)
  )
  x <- matrix(0.5, n, 1L) # constant predictor: no cut points, root-only trees

  # exact posterior mean of each group's identified softmax probabilities:
  # arm 1's 2-D quadrature over d, with the group's offset differences added
  # inside the softmax and one likelihood factor per group
  Sigma <- tau2 * matrix(c(2, 1, 1, 2), 2L, 2L)
  Prec <- solve(Sigma)
  M <- 12
  g <- seq(-M, M, length.out = 301L)
  grid <- as.matrix(expand.grid(d0 = g, d1 = g))
  d0 <- grid[, 1L]
  d1 <- grid[, 2L]
  logPOf <- function(delta) {
    e0 <- d0 + delta[1L]
    e1 <- d1 + delta[2L]
    lse <- log1p(exp(e0) + exp(e1)) # log(exp(e0) + exp(e1) + 1); f2 = 0
    cbind(e0 - lse, e1 - lse, -lse)
  }
  # the reference category's offset drops out of the differences, exactly as
  # the reference forest's level does
  logP <- lapply(1:2, function(gr) {
    logPOf(offsetByGroup[gr, 1:2] - offsetByGroup[gr, 3L])
  })
  q <- Prec[1, 1] * d0^2 + 2 * Prec[1, 2] * d0 * d1 + Prec[2, 2] * d1^2
  w <- as.vector(logP[[1L]] %*% counts[, 1L] + logP[[2L]] %*% counts[, 2L]) -
    0.5 * q
  w <- exp(w - max(w))
  exact <- vapply(
    1:2,
    function(gr) {
      colSums(exp(logP[[gr]]) * w) / sum(w)
    },
    numeric(K)
  )

  fitSeed <- function(seed) {
    set.seed(seed)
    control <- dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 50L,
      updateState = FALSE
    )
    host <- dbarts(x, as.double(labels), control = control)
    bc <- dbarts:::bartcoreMultinomialSampler(
      host,
      labels,
      K = K,
      offset = offset
    )
    r <- bartcoreRun(bc, nburn, ndpost)
    # every row of a group shares that group's probabilities; average them
    c(
      apply(r$train[group == 1L, , , drop = FALSE], 2L, mean),
      apply(r$train[group == 2L, , , drop = FALSE], 2L, mean)
    )
  }
  fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
  fit <- matrix(fit, K, 2L)

  cat("Arm 6 (intercept-only K = 3 with a category offset):\n")
  for (gr in 1:2) {
    cat(sprintf(
      "  group %d exact   %s\n",
      gr,
      paste(sprintf("%.4f", exact[, gr]), collapse = " ")
    ))
    cat(sprintf(
      "  group %d sampler %s\n",
      gr,
      paste(sprintf("%.4f", fit[, gr]), collapse = " ")
    ))
  }
  report("  arm6 category offset", max(abs(fit - exact)), tolerance)
}

# ---------------------------------------------------------------------------
# Arm 7: the level-centering move's own conditional (the NON-identified level)
# ---------------------------------------------------------------------------

# Arms 1-6 all compare IDENTIFIED softmax probabilities, so none of them sees
# the scale of the level-centering shift itself - only its likelihood-invariant
# direction. This arm reads the level directly. At the intercept-only
# configuration every tree is a root, so the move reduces to an exact
# INDEPENDENCE sampler drawing the grand level from its marginal N(0, tau^2/K)
# (docs/design/multinomial.md): each sweep's level is an independent draw with
# sd tau/sqrt(K), whatever n.trees is - the identity that fixes the per-leaf
# prior sd the conditional is written on. Run one sweep at a time and read the
# raw per-category fits, which the query returns AFTER the move.
#
# The gate is a Monte Carlo band: the sample sd of N draws has se ~ sd/sqrt(2N),
# and the tolerance is 3 se. The second line is the independence half of the
# claim, at the 3/sqrt(N) band for a lag-1 autocorrelation of white noise.
armLevel <- function() {
  ndpost <- if (quick) 1000L else 3000L
  nburn <- 200L
  K <- 3L
  n <- 60L
  target <- tau / sqrt(K)

  set.seed(2027L)
  labels <- sample.int(K, n, replace = TRUE) - 1L
  x <- matrix(0.5, n, 1L) # constant: no valid cut points, every tree a root

  set.seed(5027L)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 75L,
    updateState = FALSE
  )
  host <- dbarts(x, as.double(labels), control = control)
  bc <- dbarts:::bartcoreMultinomialSampler(host, labels, K = K)
  bartcoreRun(bc, nburn, 1L)
  level <- vapply(
    seq_len(ndpost),
    function(i) {
      bartcoreRun(bc, 0L, 1L)
      mean(vapply(
        seq_len(K) - 1L,
        function(k) mean(bartcoreForestFits(bc, k)),
        numeric(1L)
      ))
    },
    numeric(1L)
  )
  observed <- sd(level)
  acf1 <- cor(level[-ndpost], level[-1L])

  cat("Arm 7 (level-centering conditional, intercept-only K = 3):\n")
  cat(sprintf(
    "  marginal sd  %.6f  target tau/sqrt(K) %.6f  (N = %d)\n",
    observed,
    target,
    ndpost
  ))
  cat(sprintf("  lag-1 acf    %+.4f  (band %.4f)\n", acf1, 3 / sqrt(ndpost)))
  report(
    "  arm7 level sd",
    abs(observed - target),
    3 * target / sqrt(2 * ndpost)
  )
  report("  arm7 level acf1", abs(acf1), 3 / sqrt(ndpost))
}

arm1()
cat("\n")
arm2()
cat("\n")
arm3()
cat("\n")
armCount()
cat("\n")
armOffset()
cat("\n")
armLevel()
cat("\n")

if (anyFailure) {
  quit(status = 1L)
}
cat("OK: multinomial sampler matches the exact posterior on all arms\n")
