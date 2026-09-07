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
# priorPrecision = (k / scale)^2, scale = node.scale / sqrt(ntree) = 0.5.
#
# The veto arm (Rscript bd-balance.R veto) replaces the exact arm with the
# same gate run from OUTSIDE the admissible set: a weight vector installed
# between samples zeroes one cell's rows, so every partition isolating that
# cell holds a leaf no likelihood term reaches and is vetoed. The exact target
# is then this enumeration restricted to the admissible partitions and
# renormalized, with the leaf marginals taken over the positive-weight rows
# alone (docs/design/empty-leaf-veto.md). The chain is started STRANDED - grown
# under positive weights until it holds the offending cut, then handed the
# zeroing vector - so the arm also measures the absorption time and checks that
# the admissible set is absorbing thereafter.
#
# Usage: Rscript bd-balance.R [quick] [veto]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
vetoArm <- "veto" %in% args

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

yRange <- max(y) - min(y)
zScaled <- (y - min(y)) / yRange - 0.5
residVar <- (1 / yRange)^2

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# ---- integrated likelihood, verbatim from ConstantGaussianLeaf ----

priorPrecision <- (kLeaf / nodeScale)^2
# a zero weight is absence: the row leaves the leaf's sufficient statistics
# outright, exactly as it leaves the veto's count. Two ADJACENT cells are
# zeroed rather than one so that the stranded start below is vetoed on both
# sides of a split - the comparison the veto rank exists for. One zeroed cell
# would only ever strand a leaf beside a live sibling, whose collapse is a
# rank DECREASE and was already accepted before the rank existed.
kept <- if (vetoArm) cell == 1L | cell == K else rep(TRUE, n)
logIL <- function(cells) {
  idx <- cell >= cells[1L] & cell <= cells[2L] & kept
  z <- zScaled[idx]
  nLeaf <- length(z)
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
cat(sprintf("enumerated %d trees\n", length(trees)))

signatureOf <- function(cutIndices) {
  if (length(cutIndices) == 0L) {
    return("(none)")
  }
  paste(sort(cutIndices), collapse = "+")
}

# the veto's own definition of the target: a tree is admissible when every one
# of its leaves holds a positive-weight row, and the target is the enumeration
# restricted to those and renormalized
isAdmissible <- function(leaves) {
  all(vapply(
    leaves,
    function(range) any(kept & cell >= range[1L] & cell <= range[2L]),
    TRUE
  ))
}

logW <- numeric(length(trees))
signatures <- character(length(trees))
for (t in seq_along(trees)) {
  signatures[t] <- signatureOf(trees[[t]]$cutsUsed)
  if (!isAdmissible(trees[[t]]$leaves)) {
    logW[t] <- -Inf
    next
  }
  lw <- trees[[t]]$logPrior
  for (cellRange in trees[[t]]$leaves) {
    lw <- lw + logIL(cellRange)
  }
  logW[t] <- lw
}
w <- exp(logW - logSumExp(logW))
exactPartition <- vapply(split(w, signatures), sum, 0)
vetoedNames <- names(exactPartition)[exactPartition == 0]
partitionNames <- names(sort(
  exactPartition[exactPartition > 0],
  decreasing = TRUE
))

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
  n.cuts = K - 1L
)
sampler <- dbarts(
  matrix(x, ncol = 1L),
  y,
  control = ctl,
  tree.prior = cgm(power, base),
  node.prior = normal(kLeaf),
  resid.prior = fixed(1),
  proposal.probs = c(birth_death = 0.99, swap = 0, change = 0.01, birth = 0.5)
)
stopifnot(is.null(sampler$data@offset))

# ---- the stranded start ----
#
# Grow under positive weights until the live tree splits the cells that are
# about to be zeroed apart from each other, then install the zeroing vector.
# The current state is now outside the admissible set - the state no
# install-time gate can prevent, since weights do not ride the tree - and it is
# stranded on BOTH sides of a split, so the collapse that repairs it is an
# equal-rank comparison. That comparison is NaN without the rank, and a chain
# that reaches it never leaves.
absorbed <- NA_integer_
if (vetoArm) {
  cutsHeld <- function() {
    live <- sampler$getTrees(current = TRUE)
    match(live$value[live$var == 1L], cuts)
  }
  partitionHeld <- function() signatureOf(cutsHeld())
  # The one shape that makes the two zeroed cells SIBLINGS: the cuts in
  # depth-first order are 1, then 3, then 2, so the leaves holding cells 2 and
  # 3 are the children of one parent whose own range is entirely zero-weight.
  # Their collapse is then an EQUAL-rank comparison - the one the veto rank
  # decides and the one that is NaN without it - and it is the tree's only
  # death, so a chain that cannot take it cannot leave this partition at all.
  strandedShape <- c(1L, 3L, 2L)
  stranded <- FALSE
  for (sweep in seq_len(20000L)) {
    invisible(sampler$run(1L, 1L))
    if (identical(cutsHeld(), strandedShape)) {
      stranded <- TRUE
      break
    }
  }
  if (!stranded) {
    cat("\nFAIL: could not grow the tree the arm strands\n")
    quit(status = 1L)
  }
  sampler$setWeights(as.double(kept))
  for (sweep in seq_len(20000L)) {
    invisible(sampler$run(0L, 1L))
    if (partitionHeld() %in% partitionNames) {
      absorbed <- sweep
      break
    }
  }
  if (is.na(absorbed)) {
    cat("\nFAIL: the stranded chain never re-entered the admissible set\n")
    quit(status = 1L)
  }
  cat(sprintf(
    "stranded with cells %s zeroed as siblings, absorbed after %d sweeps\n",
    paste(which(!kept[seq_len(K) * nPer]), collapse = "+"),
    absorbed
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
if (vetoArm) {
  # the admissible set is absorbing: from inside it, a vetoed partition has
  # acceptance exactly zero, so the run must contain none at all
  visits <- sum(engineSignatures %in% vetoedNames)
  cat(sprintf(
    "vetoed partitions (%s): %d visits after absorption\n",
    paste(vetoedNames, collapse = ", "),
    visits
  ))
  anyFailure <- visits > 0L
}
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

if (anyFailure) {
  cat("\nFAIL: birth/death stationary distribution deviates from exact\n")
  quit(status = 1L)
}
cat("\nOK: birth/death chain matches the exact posterior over partitions\n")
