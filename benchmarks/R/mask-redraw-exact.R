#!/usr/bin/env Rscript

# Exact-posterior gate for an active-row mask that a larger sampler redraws
# every sweep. The model is a two-part mixture: row i is in the BART component
# with probability p (a_i = 1) and otherwise in a fixed background component.
# The combined sampler alternates one BART sweep given the mask with a redraw
# of every a_i from its conditional given the fit, the fit's own latents
# integrated out, and installs the new mask with $setActiveRows. That second
# step is the right conditional only if the prior over trees does not read the
# mask - a leaf is empty only when no row at all reaches it, so a leaf of only
# switched-off rows is legal and keeps its prior - and, on a latent family,
# only if a row switched back in has its latent redrawn against the current
# fit before the trees next move (docs/design/empty-leaf-veto.md,
# docs/design/active-rows-mask.md).
#
# Ten rows, one predictor, one tree, so the joint posterior over (tree, mask)
# is enumerable: every tree x all 1024 masks, leaf values integrated out in
# closed form (gaussian) or by 1-D quadrature (probit). Three arms:
#   gaussian  one ordinal predictor with 5 values, residual sd fixed;
#   factor    one unordered factor with 4 levels, split by subsets, one level
#             holding two rows that look like background, so the mask often
#             leaves it no active row;
#   probit    the ordinal design with a binary response and a Bernoulli
#             background.
# Compared, per arm: the long-run membership probability P(a_i = 1 | y) of
# each row (Rao-Blackwellized: the mean of the conditional the mask is drawn
# from) and the posterior mean of the fit at each cell, by batch-means z.
#
# Sizes and what they detect. Measured on this design, a sampler that refuses
# a leaf of only switched-off rows misses the membership probabilities by 0.066
# (gaussian), 0.21 (factor) and 0.058 (probit) and the fits by 0.19, 0.67 and
# 0.45; a probit sampler that keeps a reactivated row's stale latent misses
# them by 0.004 and 0.055. Batch-means standard errors run near 0.8 / sqrt(N)
# and 2.4 / sqrt(N) (gaussian), 0.5 / sqrt(N) and 2.9 / sqrt(N) (probit) over
# N sweeps. Quick mode keeps 1e6 sweeps on the gaussian and factor arms and
# 2e6 on the probit arm, which puts the first defect at 35 standard errors or
# more and the second at about 27 (fit) and 12 (membership); full mode keeps
# 8e6 each. The bound is |z| <= 4.5 on each of 44 statistics, their standard
# errors from 100 or 200 batch means (800 in full mode): by a union bound with
# Student-t tails a correct sampler fails a fresh seed once in 1400 quick runs
# and once in 2900 full ones. The seed is fixed, so a given build either
# passes or does not.
#
# Usage: Rscript mask-redraw-exact.R [quick] [gaussian] [factor] [probit]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
armNames <- intersect(c("gaussian", "factor", "probit"), args)
if (length(armNames) == 0L) {
  armNames <- c("gaussian", "factor", "probit")
}

numChains <- 4L
batchSize <- 10000L
numBatches <- if (quick) {
  c(gaussian = 25L, factor = 25L, probit = 50L)
} else {
  c(gaussian = 200L, factor = 200L, probit = 200L)
}
engineSeed <- 20261005L
zBound <- 4.5

# ---- the designs ----

designs <- list(
  gaussian = list(
    cell = rep(1:5, each = 2L),
    y = c(-0.2, 0.3, 0.4, 1.1, 1.9, 2.3, 0.2, 0.9, 2.6, 3.1),
    family = "gaussian",
    factor = FALSE,
    nodeScale = 0.5
  ),
  factor = list(
    cell = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L),
    y = c(-0.2, 0.3, 0.6, 2.1, 2.7, 0.2, 0.9, 1.9, 2.6, 3.1),
    family = "gaussian",
    factor = TRUE,
    nodeScale = 0.5
  ),
  probit = list(
    cell = rep(1:5, each = 2L),
    y = c(0, 0, 0, 1, 1, 1, 1, 0, 1, 1),
    family = "probit",
    factor = FALSE,
    nodeScale = 3 # the binary families' node scale, over one tree
  )
)
# shared by the three: P(a_i = 1), the gaussian background N(mu0, sigma0^2)
# and BART residual sd, the probit background Bernoulli(q0), and the priors
pActive <- 0.5
mu0 <- 0.5
sigma0 <- 1
sigma <- 0.5
q0 <- 0.25
kLeaf <- 2
base <- 0.95
power <- 2

# ---- tree enumeration ----
#
# Every tree as its leaves, each a bit set over the cells, and its log prior.
# Ordinal: bd-balance.R's recursion over runs of adjacent cells. Factor:
# categorical-exact.R's, uniform over the 2^R - 2 direction assignments of the
# R levels reachable at a node.

cellBits <- function(cells) sum(bitwShiftL(1L, cells - 1L))
bitCells <- function(bits, K) {
  which(bitwAnd(bits, bitwShiftL(1L, 0:(K - 1L))) != 0L)
}

enumerateOrdinalTrees <- function(K) {
  enumerate <- function(loCell, hiCell, depth) {
    growth <- if (hiCell > loCell) base / (1 + depth)^power else 0
    result <- list(list(
      leaves = cellBits(loCell:hiCell),
      logPrior = log(1 - growth)
    ))
    for (j in seq_len(hiCell - loCell) + loCell - 1L) {
      for (left in enumerate(loCell, j, depth + 1L)) {
        for (right in enumerate(j + 1L, hiCell, depth + 1L)) {
          result[[length(result) + 1L]] <- list(
            leaves = c(left$leaves, right$leaves),
            logPrior = log(growth) -
              log(hiCell - loCell) +
              left$logPrior +
              right$logPrior
          )
        }
      }
    }
    result
  }
  enumerate(1L, K, 0L)
}

enumerateFactorTrees <- function(K) {
  enumerate <- function(mask, depth) {
    numReachable <- length(bitCells(mask, K))
    growth <- if (numReachable >= 2L) base / (1 + depth)^power else 0
    result <- list(list(leaves = mask, logPrior = log(1 - growth)))
    if (numReachable < 2L) {
      return(result)
    }
    subsets <- Filter(
      function(s) s != 0L && s != mask && bitwAnd(s, bitwNot(mask)) == 0L,
      seq_len(bitwShiftL(1L, K) - 1L)
    )
    for (directions in subsets) {
      for (left in enumerate(bitwAnd(mask, bitwNot(directions)), depth + 1L)) {
        for (right in enumerate(directions, depth + 1L)) {
          result[[length(result) + 1L]] <- list(
            leaves = c(left$leaves, right$leaves),
            logPrior = log(growth) -
              log(2^numReachable - 2) +
              left$logPrior +
              right$logPrior
          )
        }
      }
    }
    result
  }
  enumerate(bitwShiftL(1L, K) - 1L, 0L)
}

# ---- the exact posterior of the joint over (tree, mask) ----
#
# The tree prior is the engine's CGM prior and does not read the mask. Given a
# tree and a mask the leaves are independent: each integrates its active rows
# against its N(0, tau^2) prior, and a leaf with no active row contributes 1
# and keeps that prior.
exactPosterior <- function(d) {
  K <- max(d$cell)
  n <- length(d$y)
  trees <- if (d$factor) enumerateFactorTrees(K) else enumerateOrdinalTrees(K)
  # all 2^n masks; row r has code r - 1, row i's bit being 2^(i - 1)
  codes <- 0:(2^n - 1)
  masks <- vapply(
    seq_len(n),
    function(i) bitwAnd(bitwShiftR(codes, i - 1L), 1L),
    codes
  )

  if (d$family == "gaussian") {
    yRange <- max(d$y) - min(d$y)
    shift <- min(d$y) + 0.5 * yRange
    tau <- yRange * d$nodeScale / kLeaf
    y0 <- d$y - shift
    backgroundLogDensity <- dnorm(d$y, mu0, sigma0, log = TRUE)
  } else {
    shift <- 0
    tau <- d$nodeScale / kLeaf
    step <- 0.004
    muGrid <- seq(-12, 12, by = step)
    muDensity <- dnorm(muGrid, 0, tau)
    logPhi <- pnorm(muGrid, log.p = TRUE)
    logPhiC <- pnorm(-muGrid, log.p = TRUE)
    backgroundLogDensity <- d$y * log(q0) + (1 - d$y) * log(1 - q0)
  }

  # per distinct leaf and mask: the log marginal of the leaf's active rows and
  # the posterior mean of the leaf value
  leafBits <- sort(unique(unlist(lapply(trees, function(tree) tree$leaves))))
  leafLogMarginal <- matrix(0, length(leafBits), nrow(masks))
  leafMean <- matrix(0, length(leafBits), nrow(masks))
  for (l in seq_along(leafBits)) {
    inLeaf <- as.double(d$cell %in% bitCells(leafBits[l], K))
    count <- as.vector(masks %*% inLeaf)
    if (d$family == "gaussian") {
      s <- as.vector(masks %*% (inLeaf * y0))
      ss <- as.vector(masks %*% (inLeaf * y0^2))
      precision <- 1 / tau^2 + count / sigma^2
      leafLogMarginal[l, ] <- -0.5 *
        count *
        log(2 * pi * sigma^2) -
        0.5 * ss / sigma^2 +
        0.5 * log((1 / tau^2) / precision) +
        0.5 * (s / sigma^2)^2 / precision
      leafMean[l, ] <- (s / sigma^2) / precision
    } else {
      successes <- as.vector(masks %*% (inLeaf * d$y))
      key <- successes * (n + 1L) + count
      for (k in unique(key)) {
        hit <- key == k
        w <- muDensity *
          exp(
            successes[hit][1L] * logPhi + (count - successes)[hit][1L] * logPhiC
          )
        leafLogMarginal[l, hit] <- log(sum(w) * step)
        leafMean[l, hit] <- sum(muGrid * w) / sum(w)
      }
    }
  }

  numActive <- rowSums(masks)
  maskLogPrior <- numActive *
    log(pActive) +
    (n - numActive) * log(1 - pActive) +
    as.vector((1 - masks) %*% backgroundLogDensity)

  logWeight <- matrix(0, length(trees), nrow(masks))
  cellMean <- array(0, c(length(trees), nrow(masks), K))
  for (t in seq_along(trees)) {
    kinds <- match(trees[[t]]$leaves, leafBits)
    logWeight[t, ] <- trees[[t]]$logPrior +
      colSums(leafLogMarginal[kinds, , drop = FALSE]) +
      maskLogPrior
    for (kind in kinds) {
      for (cell in bitCells(leafBits[kind], K)) {
        cellMean[t, , cell] <- leafMean[kind, ]
      }
    }
  }
  w <- exp(logWeight - max(logWeight))
  w <- w / sum(w)
  list(
    numTrees = length(trees),
    membership = as.vector(colSums(w) %*% masks),
    fit = shift + vapply(seq_len(K), function(c) sum(w * cellMean[,, c]), 0)
  )
}

# ---- the combined sampler ----

makeSampler <- function(d, seed) {
  K <- max(d$cell)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    keepTrees = FALSE,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = 1L,
    updateState = FALSE,
    seed = seed,
    n.cuts = K - 1L
  )
  if (d$factor) {
    # as categorical-exact.R: with one tree and one variable only swap rotates
    # a child's rule up the tree
    control@proposal.probs <- c(
      birth_death = 0.5,
      swap = 0.1,
      change = 0.4,
      perturb = 0,
      rule_gibbs = 0,
      birth = 0.5
    )
    predictors <- data.frame(x1 = factor(d$cell, levels = seq_len(K)))
  } else {
    predictors <- matrix(as.double(d$cell), ncol = 1L)
  }
  dbarts(
    predictors,
    d$y,
    control = control,
    tree.prior = cgm(power, base),
    leaf.prior = normal(kLeaf),
    family = if (d$family == "probit") {
      "probit"
    } else {
      gaussian(sigma = fixed(sigma^2))
    }
  )
}

# log odds of a_i = 1 given the fit, the mask's full conditional with any
# latent of the fit integrated out
logOddsFor <- function(d) {
  if (d$family == "probit") {
    backgroundLogDensity <- d$y * log(q0) + (1 - d$y) * log(1 - q0)
    function(f) {
      log(pActive) -
        log(1 - pActive) +
        pnorm(ifelse(d$y == 1, f, -f), log.p = TRUE) -
        backgroundLogDensity
    }
  } else {
    constant <- log(pActive) -
      log(1 - pActive) -
      dnorm(d$y, mu0, sigma0, log = TRUE)
    function(f) constant + dnorm(d$y, f, sigma, log = TRUE)
  }
}

# batch means of the membership conditional and of the fit at each cell, one
# row per batch; the first batch of each chain is warm-up and dropped
runChain <- function(d, seed, numBatches) {
  sampler <- makeSampler(d, seed)
  set.seed(seed + 1000000L)
  n <- length(d$y)
  K <- max(d$cell)
  firstRow <- match(seq_len(K), d$cell)
  logOdds <- logOddsFor(d)
  # the bridge entries $run and $setActiveRows call, without the per-call
  # argument handling, which at one sweep a call is most of the run time
  pointer <- sampler$getPointer()
  keepFits <- sampler$control@keepFits
  runEntry <- dbarts:::C_dbarts_bartcore_run
  maskEntry <- dbarts:::C_dbarts_bartcore_setActiveRows

  membership <- matrix(0, numBatches, n)
  fit <- matrix(0, numBatches, K)
  for (b in 0:numBatches) {
    membershipSum <- numeric(n)
    fitSum <- numeric(K)
    for (i in seq_len(batchSize)) {
      f <- .Call(runEntry, pointer, 0L, 1L, NULL, NULL, keepFits)$train
      p <- 1 / (1 + exp(-logOdds(f)))
      .Call(maskEntry, pointer, as.double(runif(n) < p))
      membershipSum <- membershipSum + p
      fitSum <- fitSum + f[firstRow]
    }
    if (b > 0L) {
      membership[b, ] <- membershipSum / batchSize
      fit[b, ] <- fitSum / batchSize
    }
  }
  list(membership = membership, fit = fit)
}

# ---- verdict ----

compare <- function(label, batches, exact) {
  estimate <- colMeans(batches)
  se <- apply(batches, 2L, sd) / sqrt(nrow(batches))
  z <- (estimate - exact) / se
  cat(sprintf("  %s\n", label))
  cat(sprintf("    %9s %9s %9s %7s\n", "exact", "sampler", "MCse", "z"))
  for (j in seq_along(exact)) {
    cat(sprintf(
      "    %9.5f %9.5f %9.5f %+7.1f%s\n",
      exact[j],
      estimate[j],
      se[j],
      z[j],
      if (is.na(z[j]) || abs(z[j]) > zBound) " <- FAIL" else ""
    ))
  }
  worst <- which.max(abs(estimate - exact))
  cat(sprintf(
    "    largest difference %.5f, largest |z| %.1f\n",
    abs(estimate - exact)[worst],
    max(abs(z))
  ))
  anyNA(z) || any(abs(z) > zBound)
}

anyFailure <- FALSE
for (armIndex in seq_along(armNames)) {
  arm <- armNames[armIndex]
  d <- designs[[arm]]
  exact <- exactPosterior(d)
  started <- proc.time()[3L]
  chains <- lapply(
    engineSeed + 100L * armIndex + seq_len(numChains),
    function(seed) runChain(d, seed, numBatches[[arm]])
  )
  cat(sprintf(
    "\n%s arm: %d trees x 1024 masks enumerated; %d chains x %d batches x %d sweeps in %.0f s\n",
    arm,
    exact$numTrees,
    numChains,
    numBatches[[arm]],
    batchSize,
    proc.time()[3L] - started
  ))
  failedMembership <- compare(
    "membership probability P(a_i = 1 | y), by row",
    do.call(rbind, lapply(chains, function(chain) chain$membership)),
    exact$membership
  )
  failedFit <- compare(
    "posterior mean of the fit, by cell",
    do.call(rbind, lapply(chains, function(chain) chain$fit)),
    exact$fit
  )
  anyFailure <- anyFailure || failedMembership || failedFit
}

if (anyFailure) {
  cat("\nFAIL: the combined sampler deviates from the exact joint posterior\n")
  quit(status = 1L)
}
cat("\nOK: the combined sampler matches the exact joint posterior\n")
