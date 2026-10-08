#!/usr/bin/env Rscript

# Exact-posterior gate for the birth/death move's detailed balance,
# specifically the reverse-move node-selection counts (the death
# acceptance's P(select the collapsed node for birth) must be counted on
# the POST-death tree, and the birth acceptance's P(select for death) on
# the POST-birth tree; counting on the wrong tree biases the stationary
# distribution over tree sizes, a defect no other gate isolates - the
# change-balance verdicts condition on the root being split, which cancels
# size effects). A single-tree, one-ordinal-predictor (K = 4 distinct
# values, 3 uniform cuts separating the cells), constant-leaf, fixed-sigma
# problem restricted to a birth/death-dominated kernel (proposal.probs
# birth_death = 0.99; the bridge requires < 1, and the 1% change moves
# target the same posterior - change-balance validates that - so the
# exact arm is unchanged) has a fully enumerable tree space:
# contiguous-cell partitions, no depth truncation (cuts exhaust). Because leaf marginals
# depend only on the partition and single-cell leaves are not birthable,
# death transitions here cross birthable-count boundaries in both
# directions, so a wrong-tree reverse count shifts the partition
# distribution detectably in both tails. The exact posterior over the 8
# cut-set partitions (prior x integrated likelihood summed over split
# orders) is compared against long-run engine frequencies with
# batch-means z-scores.
#
# Calibration matched to the engine exactly (as change-balance.R):
# uniform cuts, range scaling, fixed(1) residual variance -> internal
# sigma = 1 / range, constant-leaf conjugate marginal with
# priorPrecision = (k / scale)^2, scale = leaf.scale / sqrt(ntree) = 0.5.
#
# The zero-weight arm (Rscript bd-balance.R zeroweight) replaces the exact arm
# with the same gate under a weight vector, installed between samples on a
# grown tree, that zeroes two adjacent cells' rows. A zero-weight row is in the
# design and not in the likelihood: a partition isolating those cells holds a
# leaf no likelihood term reaches, which is legal and scores nothing
# (docs/design/empty-leaf-veto.md). The exact target is then the WHOLE
# enumeration, the tree prior unchanged, with the leaf marginals taken over the
# positive-weight rows alone - not the enumeration restricted to the partitions
# whose every leaf holds a weighted row, which the arm reports its distance
# from.
#
# Two arms put the gate on a column whose cut grid would repeat a point if
# equal neighbours were kept, where a grid holds each point once
# (docs/design/cut-grid.md); each first asserts the grid.
# - narrow (Rscript bd-balance.R narrow): the four cells sit on four adjacent
#   doubles and n.cuts is 5. Five evenly spaced points over that range round
#   to three distinct ones, which separate the cells, so the target is the
#   enumeration above unchanged; with the two repeats kept the prior over
#   thresholds is another one.
# - constmissing (Rscript bd-balance.R constmissing): cells 1 and 2 hold one
#   value and cells 3 and 4 are missing, n.cuts 5. The grid is one point, so
#   two trees have mass: the stump, and the split of present from missing,
#   whose rule is the one cut, at half the mass for the direction missing
#   values take, and whose children hold no further cut. With five copies of
#   the point the children would count as splittable and the odds move by
#   the two factors of not splitting them.
#
# Usage: Rscript bd-balance.R [quick] [zeroweight | narrow | constmissing]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
zeroWeightArm <- "zeroweight" %in% args
narrowArm <- "narrow" %in% args
constMissingArm <- "constmissing" %in% args
stopifnot(zeroWeightArm + narrowArm + constMissingArm <= 1L)

nKept <- if (quick) 100000L else 300000L
batchSize <- if (quick) 25000L else 50000L
nThin <- 10L
nBurn <- 2000L
engineSeed <- 20260710L
zBound <- 4

# ---- fixed design ----

set.seed(17L)
K <- 4L
nPer <- 10L
muCell <- c(-0.5, -0.2, 0.2, 0.5)
noiseSd <- 1.0
cell <- rep(seq_len(K), each = nPer)
x <- as.double(cell)
n <- length(x)
y <- muCell[cell] + rnorm(n, sd = noiseSd)

base <- 0.8
power <- 2
kLeaf <- 2
nodeScale <- 0.5

cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
stopifnot(identical(findInterval(x, cuts) + 1L, cell))
numCutsAsked <- K - 1L
if (narrowArm) {
  # a row goes left of a cut it does not exceed, so the cuts sit on the
  # first three cells' values
  x <- 1 + (cell - 1L) * .Machine$double.eps
  cuts <- 1 + (seq_len(K - 1L) - 1L) * .Machine$double.eps
  stopifnot(identical(findInterval(x, cuts, left.open = TRUE) + 1L, cell))
  numCutsAsked <- 5L
} else if (constMissingArm) {
  x <- ifelse(cell <= 2L, 1, NA_real_)
  cuts <- 1
  numCutsAsked <- 5L
}

yRange <- max(y) - min(y)
zScaled <- (y - min(y)) / yRange - 0.5
residVar <- (1 / yRange)^2

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# ---- integrated likelihood, verbatim from ConstantGaussianLeaf ----

priorPrecision <- (kLeaf / nodeScale)^2
# a zero weight leaves the row out of the leaf's sufficient statistics and in
# its member count. Two ADJACENT cells are zeroed rather than one so that the
# trees include a split with no weighted row on either side, whose two leaves
# and their parent all score nothing.
kept <- if (zeroWeightArm) cell == 1L | cell == K else rep(TRUE, n)
logIL <- function(cells) {
  idx <- cell >= cells[1L] & cell <= cells[2L] & kept
  z <- zScaled[idx]
  nLeaf <- length(z)
  if (nLeaf == 0L) {
    return(0)
  }
  posteriorPrecision <- nLeaf / residVar
  mean <- sum(z) / nLeaf
  centeredSumOfSquares <- sum(z * z) - sum(z) * mean
  0.5 *
    log(priorPrecision / (priorPrecision + posteriorPrecision)) -
    0.5 * centeredSumOfSquares / residVar -
    0.5 *
      ((priorPrecision * mean) * (posteriorPrecision * mean)) /
      (priorPrecision + posteriorPrecision)
}

# ---- exact posterior over partitions ----
#
# Enumerate every tree (split orders distinct, prior depth-dependent) and
# aggregate prior x likelihood onto the partition = the set of cut indices
# used, which is what the engine arm can observe per sample.

enumerate <- function(loCell, hiCell, loCut, hiCut, depth) {
  growth <- if (hiCut >= loCut) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(loCell, hiCell)),
    cutsUsed = integer(0L),
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
          cutsUsed = c(j, left$cutsUsed, right$cutsUsed),
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

trees <- enumerate(1L, K, 1L, K - 1L, 0L)
if (constMissingArm) {
  # the stump, and present against missing: one cut to choose, half the mass
  # for sending missing values right, and neither child with a cut left
  trees <- list(
    list(
      leaves = list(c(1L, K)),
      cutsUsed = integer(0L),
      logPrior = log(1 - base)
    ),
    list(
      leaves = list(c(1L, 2L), c(3L, K)),
      cutsUsed = 1L,
      logPrior = log(base) - log(2)
    )
  )
}
cat(sprintf("enumerated %d trees\n", length(trees)))

signatureOf <- function(cutIndices) {
  if (length(cutIndices) == 0L) {
    return("(none)")
  }
  paste(sort(cutIndices), collapse = "+")
}

# whether every leaf of a tree holds a positive-weight row: the trees a rule
# counting weight rather than members would keep
holdsWeightThroughout <- function(leaves) {
  all(vapply(
    leaves,
    function(range) any(kept & cell >= range[1L] & cell <= range[2L]),
    TRUE
  ))
}

logW <- numeric(length(trees))
signatures <- character(length(trees))
weighted <- logical(length(trees))
for (t in seq_along(trees)) {
  signatures[t] <- signatureOf(trees[[t]]$cutsUsed)
  weighted[t] <- holdsWeightThroughout(trees[[t]]$leaves)
  lw <- trees[[t]]$logPrior
  for (cellRange in trees[[t]]$leaves) {
    lw <- lw + logIL(cellRange)
  }
  logW[t] <- lw
}
w <- exp(logW - logSumExp(logW))
exactPartition <- vapply(split(w, signatures), sum, 0)
partitionNames <- names(sort(exactPartition, decreasing = TRUE))
# the same enumeration restricted to the trees holding weight in every leaf
restrictedPartition <- vapply(
  split(w * weighted / sum(w * weighted), signatures),
  sum,
  0
)

# ---- engine arm: pure birth/death kernel ----

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  keepTrees = TRUE,
  n.samples = batchSize,
  n.burn = nBurn,
  n.thin = nThin,
  updateState = TRUE,
  seed = engineSeed,
  n.cuts = numCutsAsked
)
sampler <- dbarts(
  matrix(x, ncol = 1L),
  y,
  control = ctl,
  tree.prior = cgm(power, base),
  leaf.prior = normal(kLeaf),
  family = gaussian(sigma = fixed(1)),
  proposal.probs = c(birth_death = 0.99, swap = 0, change = 0.01, birth = 0.5)
)
stopifnot(is.null(sampler$data@offset))
if (narrowArm || constMissingArm) {
  grid <- attr(sampler$state, "cutPoints")[[1L]]
  if (!identical(grid, cuts)) {
    cat(sprintf(
      "FAIL: the cut grid holds %d points, %d distinct, where %d were expected\n",
      length(grid),
      length(unique(grid)),
      length(cuts)
    ))
    quit(status = 1L)
  }
  cat(sprintf(
    "cut grid: %d distinct points of %d asked\n",
    length(grid),
    numCutsAsked
  ))
}

# ---- the weights go in on a grown tree ----
#
# Weights do not ride the tree, so the install lands on whatever partition the
# chain holds; the first batch's burn-in follows it.
if (zeroWeightArm) {
  invisible(sampler$run(nBurn, 1L))
  sampler$setWeights(as.double(kept))
  cat(sprintf(
    "cells %s zeroed; %.4f of the exact posterior sits on partitions with a leaf of only zero-weight rows\n",
    paste(which(!kept[seq_len(K) * nPer]), collapse = "+"),
    1 - sum(w * weighted)
  ))
}

nBatch <- as.integer(ceiling(nKept / batchSize))
engineSignatures <- character(0L)
first <- TRUE
for (b in seq_len(nBatch)) {
  if (first) {
    invisible(sampler$run(nBurn, batchSize))
    first <- FALSE
  } else {
    invisible(sampler$run(0L, batchSize))
  }
  tr <- sampler$getTrees()
  splits <- tr[tr$var == 1L, c("sample", "value")]
  splits$cut <- match(splits$value, cuts)
  stopifnot(!anyNA(splits$cut))
  bySample <- split(splits$cut, splits$sample)
  sig <- rep("(none)", batchSize)
  present <- as.integer(names(bySample))
  sig[present] <- vapply(bySample, signatureOf, "")
  engineSignatures <- c(engineSignatures, sig)
}
nDraws <- length(engineSignatures)

# ---- batch-means MC error and z-scores ----

batchMeanSE <- function(v, nBatches = 400L) {
  len <- (length(v) %/% nBatches) * nBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
  sd(bm) / sqrt(nBatches)
}

cat(sprintf("\nengine draws: %d (thin %d)\n", nDraws, nThin))
cat(sprintf(
  "%-10s %10s %10s %10s %8s\n",
  "partition",
  "engine",
  "exact",
  "MCse",
  "z"
))

anyFailure <- FALSE
for (name in partitionNames) {
  engineProb <- mean(engineSignatures == name)
  se <- batchMeanSE(engineSignatures == name)
  z <- (engineProb - exactPartition[name]) / se
  failed <- is.na(z) || abs(z) > zBound
  anyFailure <- anyFailure || failed
  cat(sprintf(
    "%-10s %10.4f %10.4f %10.5f %+8.1f%s\n",
    name,
    engineProb,
    exactPartition[name],
    se,
    z,
    if (failed) " <- FAIL" else ""
  ))
}

if (zeroWeightArm) {
  # non-vacuity: the target that counts weight instead of members is a
  # different distribution, and the chain is not on it
  engineProbs <- vapply(
    partitionNames,
    function(name) mean(engineSignatures == name),
    0
  )
  cat(sprintf(
    "total variation from the exact posterior %.4f, from the one restricted to weighted leaves %.4f\n",
    0.5 * sum(abs(engineProbs - exactPartition[partitionNames])),
    0.5 * sum(abs(engineProbs - restrictedPartition[partitionNames]))
  ))
}

if (anyFailure) {
  cat("\nFAIL: birth/death stationary distribution deviates from exact\n")
  quit(status = 1L)
}
cat("\nOK: birth/death chain matches the exact posterior over partitions\n")
