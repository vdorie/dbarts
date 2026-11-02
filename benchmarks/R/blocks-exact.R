#!/usr/bin/env Rscript

# Exact-posterior gate for the blocks() prior: a forest partitioned into G
# groups gives each group 1 / G of its prior variance whatever the groups'
# tree counts, so a tree of a group of m_g trees has leaf variance
# T^2 / (G m_g), T the forest's prior sd. With equal counts that is every
# tree's T^2 / m; with unequal counts the groups' leaf scales differ, and that
# is what this gate holds the sampler to.
#
# One forest of three trees in two groups: x1, a 0/1 column, takes trees 1 and
# 2 and x2, another 0/1 column, takes tree 3. Each column has one cut, so each
# tree is a root or its one split and the forest has 8 structures; the tree
# prior is the engine's CGM prior, a root splitting with probability base and
# a child of a split having nothing left to split on. Twelve rows, three in
# each of the four cells; gaussian with sigma fixed. Given a structure the
# leaves are jointly normal, trees 1 and 2's at T^2 / 4 and tree 3's at
# T^2 / 2, so the marginal likelihood and the posterior mean of the fit are in
# closed form, and the posterior is enumerated exactly. Two arms:
#   fixed  k = 2;
#   chi    k drawn under chi(1.5, 2), integrated by an outer quadrature on a
#          grid of k.
# Compared, per arm, by batch-means z: the posterior probability that 0, 1 or
# 2 of the x1 trees split (trees 1 and 2 are exchangeable, so this is their
# joint), that tree 3 splits, the posterior mean of the fit at each cell, and
# in the chi arm the posterior mean of k.
#
# Sizes and what they detect. The per-tree prior, every leaf at T^2 / 3 (the
# blocks() prior before each group took an equal part), moves the fits at the
# x1-only and x2-only cells by 0.033 in the fixed arm and 0.032 in the chi
# arm, the split probabilities by 0.002 to 0.020 and k by 0.018. Batch-means
# standard errors in quick mode run near 1e-4 on the fits, 2.5e-4 on the
# split probabilities and 1.3e-3 on k, so the per-tree prior, run against a
# build carrying it, misses the fits by more than 300 of them and fails both
# arms. A k statistic that sums the leaves' squares without each group's
# factor misses the chi arm's k by 0.016 and two of its fits by 0.002, failing
# them at |z| near 13. The bound is |z| <= 4.5 on each of 8 statistics in the
# fixed arm and 9 in the chi arm, their standard errors from 4 chains of 100
# batches of 1e4 sweeps (400 batches in full mode): by a union bound with
# Student-t tails a correct sampler fails a fresh seed once in about 6600
# quick runs and 8000 full ones. The seed is fixed, so a given build either
# passes or does not. Quick mode is about 10 seconds.
#
# Usage: Rscript blocks-exact.R [quick] [fixed] [chi]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
armNames <- intersect(c("fixed", "chi"), args)
if (length(armNames) == 0L) {
  armNames <- c("fixed", "chi")
}

numChains <- 4L
batchSize <- 10000L
numBatches <- if (quick) 100L else 400L
engineSeed <- 20261029L
zBound <- 4.5

# ---- the design ----

x1 <- rep(c(0, 1, 0, 1), each = 3L)
x2 <- rep(c(0, 0, 1, 1), each = 3L)
y <- c(-0.24, 0.01, -0.61, 0.15, 1.17, 0.33, 0.88, 0.6, 0.33, 0.65, 0.72, 0.91)
cellOfRow <- 1L + x1 + 2L * x2
sigma <- 0.4
base <- 0.5
power <- 2
kFixed <- 2
kDegreesOfFreedom <- 1.5
kScale <- 2
# the trees' groups and each tree's share of the forest's prior variance:
# group x1's two trees each 1 / (2 * 2), group x2's one tree 1 / (2 * 1)
treeColumn <- c(1L, 1L, 2L)
treeShare <- c(1 / 4, 1 / 4, 1 / 2)

# ---- the exact posterior ----
#
# Every structure, its log prior, and per value of k its log marginal
# likelihood and posterior mean of the fit at each cell; the k arm weighs a
# grid of k by its chi prior and the k = 2 arm puts all its mass on 2.
exactPosterior <- function(kGrid, kLogPrior) {
  n <- length(y)
  yRange <- max(y) - min(y)
  shift <- min(y) + 0.5 * yRange
  y0 <- y - shift
  # the forest's prior sd at k is T / k, T the leaf prior's k.scale
  forestScale <- 0.5 * yRange
  columns <- cbind(x1, x2)
  cellColumns <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
  structures <- as.matrix(expand.grid(0:1, 0:1, 0:1))
  numStructures <- nrow(structures)
  logWeight <- matrix(0, numStructures, length(kGrid))
  cellFit <- array(0, c(numStructures, length(kGrid), 4L))
  for (s in seq_len(numStructures)) {
    design <- NULL
    cellDesign <- NULL
    share <- NULL
    for (t in 1:3) {
      if (structures[s, t] == 0) {
        design <- cbind(design, rep(1, n))
        cellDesign <- cbind(cellDesign, rep(1, 4L))
      } else {
        column <- columns[, treeColumn[t]]
        cellColumn <- cellColumns[, treeColumn[t]]
        design <- cbind(design, column == 0, column == 1)
        cellDesign <- cbind(cellDesign, cellColumn == 0, cellColumn == 1)
      }
      share <- c(share, rep(treeShare[t], 1L + structures[s, t]))
    }
    logPrior <- sum(ifelse(structures[s, ] == 1, log(base), log(1 - base)))
    for (g in seq_along(kGrid)) {
      leafVariance <- diag((forestScale / kGrid[g])^2 * share, length(share))
      covariance <- sigma^2 *
        diag(n) +
        design %*% leafVariance %*% t(design)
      root <- chol(covariance)
      solved <- backsolve(root, forwardsolve(t(root), y0))
      logWeight[s, g] <- logPrior +
        kLogPrior[g] -
        sum(log(diag(root))) -
        0.5 * sum(y0 * solved)
      cellFit[s, g, ] <- shift +
        cellDesign %*% leafVariance %*% t(design) %*% solved
    }
  }
  w <- exp(logWeight - max(logWeight))
  w <- w / sum(w)
  numSplitX1 <- rowSums(structures[, 1:2])
  list(
    splits = c(
      vapply(0:2, function(j) sum(w[numSplitX1 == j, ]), 0),
      sum(w[structures[, 3] == 1, ])
    ),
    fit = vapply(seq_len(4L), function(c) sum(w * cellFit[,, c]), 0),
    k = sum(w * rep(kGrid, each = numStructures))
  )
}

# ---- the sampler ----

makeSampler <- function(leafPrior, seed) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = FALSE,
    n.samples = batchSize,
    n.burn = 0L,
    n.thin = 1L,
    updateState = FALSE,
    seed = seed
  )
  dbarts(
    data.frame(x1 = x1, x2 = x2),
    y,
    control = control,
    n.trees = 3L,
    n.cuts = 1L,
    tree.prior = cgm(power, base),
    leaf.prior = leafPrior,
    family = gaussian(sigma = fixed(sigma^2)),
    blocks = dbartsForests$blocks(
      groups = list("x1", "x2"),
      trees.per.group = c(2L, 1L)
    )
  )
}

# batch means of the split indicators, the fit at each cell and k, one row
# per batch; the first batch of each chain is warm-up and dropped
runChain <- function(leafPrior, seed) {
  sampler <- makeSampler(leafPrior, seed)
  # the calibration the reference assumes: T = k.scale, half the range
  calibration <- sampler$getLeafPrior()
  stopifnot(isTRUE(all.equal(calibration$k.scale, 0.5 * diff(range(y)))))
  firstRow <- match(seq_len(4L), cellOfRow)
  splits <- matrix(0, numBatches, 4L)
  fit <- matrix(0, numBatches, 4L)
  k <- numeric(numBatches)
  for (b in 0:numBatches) {
    result <- sampler$run(0L, batchSize)
    if (b == 0L) {
      next
    }
    splitX1 <- result$varcount[1L, ]
    splits[b, ] <- c(
      mean(splitX1 == 0L),
      mean(splitX1 == 1L),
      mean(splitX1 == 2L),
      mean(result$varcount[2L, ] == 1L)
    )
    fit[b, ] <- rowMeans(result$train[firstRow, , drop = FALSE])
    # a fixed k is not reported
    k[b] <- if (is.null(result$k)) NA_real_ else mean(result$k)
  }
  list(splits = splits, fit = fit, k = k)
}

# ---- verdict ----

compare <- function(label, batches, exact) {
  batches <- as.matrix(batches)
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
  anyNA(z) || any(abs(z) > zBound)
}

arms <- list(
  fixed = list(
    leafPrior = dbartsPriors$normal(kFixed),
    kGrid = kFixed,
    kLogPrior = 0
  ),
  chi = local({
    grid <- seq(0.005, 15, by = 0.005)
    list(
      leafPrior = dbartsPriors$normal(dbartsPriors$chi(
        kDegreesOfFreedom,
        kScale
      )),
      kGrid = grid,
      kLogPrior = (kDegreesOfFreedom - 1) * log(grid) - grid^2 / (2 * kScale^2)
    )
  })
)

anyFailure <- FALSE
for (armIndex in seq_along(armNames)) {
  arm <- armNames[armIndex]
  spec <- arms[[arm]]
  exact <- exactPosterior(spec$kGrid, spec$kLogPrior)
  started <- proc.time()[3L]
  chains <- lapply(
    engineSeed + 100L * armIndex + seq_len(numChains),
    function(seed) runChain(spec$leafPrior, seed)
  )
  cat(sprintf(
    "\n%s arm: 8 structures enumerated; %d chains x %d batches x %d sweeps in %.0f s\n",
    arm,
    numChains,
    numBatches,
    batchSize,
    proc.time()[3L] - started
  ))
  failed <- compare(
    "P(0, 1, 2 of the x1 trees split), P(the x2 tree splits)",
    do.call(rbind, lapply(chains, function(chain) chain$splits)),
    exact$splits
  )
  failed <- compare(
    "posterior mean of the fit, by cell",
    do.call(rbind, lapply(chains, function(chain) chain$fit)),
    exact$fit
  ) ||
    failed
  if (arm == "chi") {
    failed <- compare(
      "posterior mean of k",
      unlist(lapply(chains, function(chain) chain$k)),
      exact$k
    ) ||
      failed
  }
  anyFailure <- anyFailure || failed
}

if (anyFailure) {
  cat("\nFAIL: the blocks() sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("\nOK: the blocks() sampler matches the exact posterior\n")
