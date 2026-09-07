#!/usr/bin/env Rscript

# Permanent detailed-balance gate for the PERTURB move (docs/design/
# perturb-move.md section 4). The move keeps an interior node's split variable
# and displaces only its cut, drawing uniformly from the window
# W(c) = {j in [lo, hi] : 0 < |j - c| <= w} at w = 1 and correcting by
# log|W(c)| - log|W(c')|. The node's own prior factors cancel exactly rather
# than against a proposal density, so the whole Hastings term IS that window
# ratio: it fires only at the ends of the descendant-valid interval, and a
# kernel that dropped it would still look correct anywhere in the middle. This
# gate is the one that reads the WITHIN-VARIABLE cut distribution; the shipped
# balance gates read a root-split-VARIABLE marginal (change-balance.R), a
# partition distribution (bd-balance.R) and a two-state prior ratio
# (swap-balance.R), none of which a cut displacement moves.
#
# ---- the prior-only arm ----
#
# An all-zero weight vector installed with $setWeights turns the likelihood
# off: every occupied leaf takes veto rank 1, every member-empty leaf rank 2,
# no leaf enters a likelihood term, and both branches of every acceptance score
# exactly 0. The kernel is then reversible with respect to the CGM tree prior
# truncated to the member-occupied trees, and the design makes that truncation
# vacuous: two ordinal columns of 6 and 4 distinct values as a FULL FACTORIAL
# with useQuantiles = TRUE, which puts 5 and 3 cuts at the midpoints between
# consecutive distinct values, so every leaf of every reachable tree holds a
# cell and the target is the plain CGM prior. Every chain starts at a bare
# root, so the run needs a burn-in it would not otherwise.
#
# The tree space holds 33,610,060,775 trees, so no chi-square over trees.
# Three statistics instead, all closed form or a dynamic program:
#   (1) the root's (variable, cut) marginal plus the stump, NINE states,
#       P(grow) x P(v) x 1/|SI_v|;
#   (2) the leaf count, from a dynamic program whose state is (remaining x1
#       cuts, remaining x2 cuts, DEPTH) - growth is base/(1 + depth)^power, so
#       remaining cuts alone does not define the recursion - pooled at >= 6;
#   (3) the (root cut, left-child cut) joint on the SAME variable, thirteen
#       states, closed form: the descendant-valid interval and the clipped
#       window, which nothing else gates.
# Both closed forms and the program draw the split VARIABLE uniformly among
# the variables that still have a cut at the node, which is what the engine
# does with no split probabilities set: a variable exhausted below its
# ancestors drops out of the draw rather than wasting it, and the leaf-count
# vector reproduces only under that convention.
#
# States of prior mass below 0.004 are dropped by this PRE-STATED rule, their
# counts being degenerate at any feasible run length: statistic 1 keeps all
# nine, statistic 2 one to five leaves plus the pooled >= 6 bin, statistic 3
# the six states at root cut 2 and 3, dropping the seven at 4 and 5. Family
# size m = 21, Holm at family alpha 0.05, so the thresholds run from |z| = 1.96
# on the least extreme test to |z| = 3.04 on the most extreme.
#
# The arm's mixture is birth_death 0.10, change 0.10, perturb 0.80: change is
# retained because it moves the root's VARIABLE directly (a perturb-only chain
# would have to pass through a stump, 5 percent of the prior mass) and
# birth/death because statistic 2 moves through nothing else, both at the
# smallest share keeping their own statistic non-degenerate while leaving
# perturb dominant on statistics 1 and 3.
#
# The z is a pooled batch-means z: a state's indicator is averaged within a
# chain over 500 batches of consecutive kept draws, s_c is the sd of chain c's
# batch means, and the four independent chains pool as the mean of the four
# estimates over sqrt(sum_c s_c^2 / 500) / 4. Batching rather than a binomial
# formula is what absorbs the within-chain autocorrelation. The burn-in is not
# assumed: the script reads the first lag at which each scored series'
# kept-draw autocorrelation falls under 0.1 and REFUSES to score a run whose
# burn-in is under fifty of those lags, doubling and re-running instead.
#
# ---- the confirmation arm ----
#
# The prior-only arm never exercises the likelihood term. The confirmation arm
# runs the same grid and the same mixture with positive weights and NO mask -
# the mask is what makes the likelihood constant, so the arms cannot share a
# configuration - and scores the WITHIN-VARIABLE root cut law against an exact
# region dynamic program in change-balance.R's shape: a node's reachable set is
# a rectangle in (x1 value index, x2 value index), splitInterval is
# ancestor-only, and no depth truncation is needed because 5 + 3 cuts cap every
# path at depth 8. Calibration is matched to the engine exactly, else the
# comparison is void: quantile cuts, range scaling, fixed(1) residual variance
# -> internal sigma = 1 / range, and the constant-leaf conjugate marginal with
# priorPrecision = (k / scale)^2, scale = node.scale / sqrt(ntree) = 0.5. The
# same Holm machinery scores its eight states as their own family, and the
# distance to the PRIOR conditional (uniform over the variable's cuts) is
# reported as the arm's power: a kernel that dropped the likelihood would sit
# there.
#
# ---- the poison arms ----
#
# Both poisons of the design are engine defects, so each is run here the way
# change-balance.R and swap-balance.R run their wrong targets: the flag mutates
# the SCRIPT's expectation to the law the defective kernel would target, and
# the gate must FAIL against it. Both are leading-order, undiluted predictions
# - the arm's change and birth/death shares pull a really-poisoned chain back
# toward the true law - and both mutate statistic 1 alone, which is the
# statistic the design requires each to fail.
#   poison1 : logProposalCorrection dropped. The uncorrected chain is
#             reversible for p(c) proportional to pi(c)|W(c)|, so moves out of
#             an end are over-accepted and moves into one under-accepted and
#             the boundary cuts starve: x1's five root cuts go from 0.2 each to
#             (1,2,2,2,1)/8, a 37 percent relative shift at the ends.
#   poison2 : a one-sided +w window. |W| is then 1 everywhere but at hi, where
#             it is 0, so every proposal is c -> c+1 and the cut is absorbed at
#             the top of the interval.
# The matching engine mutations are carried by benchmarks/R/mutation-battery.R,
# which builds them for real; the interval invariance the reverse count rests
# on is asserted in tests/cpp instead of poisoned here.
#
# Usage: Rscript perturb-balance.R [quick] [poison1|poison2]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
poison <- if ("poison1" %in% args) {
  1L
} else if ("poison2" %in% args) {
  2L
} else {
  0L
}

# ---- fixed design ----

K1 <- 6L # x1 distinct values -> 5 quantile cuts
K2 <- 4L # x2 distinct values -> 3 quantile cuts
base <- 0.95
power <- 2
kLeaf <- 2
nodeScale <- 0.5 # gaussian node.scale default; scale = nodeScale / sqrt(1)
windowWidth <- 1L # moves.hpp perturbWidth

proposalProbs <- c(
  birth_death = 0.10,
  swap = 0,
  change = 0.10,
  perturb = 0.80,
  birth = 0.5
)

nChains <- 4L
nThin <- 20L
nBatches <- 500L
massFloor <- 0.004
familyAlpha <- 0.05
acfLagMax <- 200L
acfBurnRatio <- 50L
acfDraws <- 25000L # kept draws per chain the acf ladder reads
maxBurnDoublings <- 3L

priorBatchLen <- if (quick) 100L else 500L
priorBlock <- if (quick) 25000L else 50000L
priorBurn <- 20000L
priorSeed <- 20260907L

confirmBatchLen <- if (quick) 50L else 200L
confirmBlock <- if (quick) 25000L else 50000L
confirmBurn <- 100000L
confirmSeed <- 20260908L
confirmPer <- 2L # replicates per grid cell

cuts1 <- (seq_len(K1)[-K1] + seq_len(K1)[-1L]) / 2
cuts2 <- (seq_len(K2)[-K2] + seq_len(K2)[-1L]) / 2

logSumExp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(v - m)))
}

# ---- statistic 1: the root (variable, cut) marginal, closed form ----
#
# P(grow at the root) = base, the variable is uniform over the two available
# columns, and the rule is uniform over that variable's cuts. Nine states,
# summing to one exactly.

rootLabels <- c(
  "stump",
  paste0("x1c", seq_len(K1 - 1L)),
  paste0("x2c", seq_len(K2 - 1L))
)
rootTarget <- c(
  1 - base,
  rep(base * 0.5 / (K1 - 1L), K1 - 1L),
  rep(base * 0.5 / (K2 - 1L), K2 - 1L)
)
names(rootTarget) <- rootLabels
stopifnot(abs(sum(rootTarget) - 1) < 1e-12)

# The window count |W(c)| = min(hi, c + w) - max(lo, c - w) at a root, whose
# own interval is the whole grid. Poison 1's target reweights the cut law by
# it within each variable; poison 2's absorbs the cut at the top.
windowCount <- function(cut, nCut) {
  min(nCut, cut + windowWidth) - max(1L, cut - windowWidth)
}

poisonRootTarget <- function(which) {
  out <- rootTarget
  for (v in 1:2) {
    nCut <- if (v == 1L) K1 - 1L else K2 - 1L
    idx <- match(paste0("x", v, "c", seq_len(nCut)), rootLabels)
    shape <- if (which == 1L) {
      vapply(seq_len(nCut), windowCount, 0L, nCut = nCut)
    } else {
      c(rep(0L, nCut - 1L), 1L)
    }
    out[idx] <- sum(rootTarget[idx]) * shape / sum(shape)
  }
  out
}

# ---- statistic 2: the leaf count, by dynamic program ----
#
# State (remaining x1 cuts, remaining x2 cuts, depth). A node with no available
# variable cannot grow (growthProbability returns 0 there), the split variable
# is uniform over the variables that still have a cut, and splitting a variable
# with a available cuts at its j-th sends j - 1 of them left and a - j right.
# The leaf count of a split node is the convolution of its children's.

leafCountTarget <- function() {
  memo <- new.env(parent = emptyenv())
  dist <- function(avail1, avail2, depth) {
    key <- paste(avail1, avail2, depth, sep = ",")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    numAvail <- (avail1 > 0L) + (avail2 > 0L)
    growth <- if (numAvail == 0L) 0 else base / (1 + depth)^power
    out <- numeric((avail1 + 1L) * (avail2 + 1L))
    out[1L] <- 1 - growth
    if (growth > 0) {
      for (v in 1:2) {
        avail <- if (v == 1L) avail1 else avail2
        if (avail == 0L) {
          next
        }
        weight <- growth / (numAvail * avail)
        for (j in seq_len(avail)) {
          if (v == 1L) {
            left <- dist(j - 1L, avail2, depth + 1L)
            right <- dist(avail - j, avail2, depth + 1L)
          } else {
            left <- dist(avail1, j - 1L, depth + 1L)
            right <- dist(avail1, avail - j, depth + 1L)
          }
          joint <- convolve(left, rev(right), type = "open")
          out[2L:(length(joint) + 1L)] <- out[2L:(length(joint) + 1L)] +
            weight * joint
        }
      }
    }
    memo[[key]] <- out
    out
  }
  dist(K1 - 1L, K2 - 1L, 0L)
}

leafFull <- leafCountTarget()
stopifnot(abs(sum(leafFull) - 1) < 1e-9)
leafLabels <- c(paste0("leaves", 1:5), "leaves6+")
leafTarget <- c(leafFull[1:5], sum(leafFull[6:length(leafFull)]))
names(leafTarget) <- leafLabels

# ---- statistic 3: the (root cut, left-child cut) joint on one variable ----
#
# The root splits variable v at cut c, and its LEFT child splits the same
# variable: the child's ancestor-constrained interval is 1..c-1, so the pair
# exists only for c >= 2, and the other variable is always available there, so
# the child's variable draw is uniform over two. Thirteen states across the two
# columns; the depth-1 growth factor is base / 2^power.

pairTable <- do.call(
  rbind,
  lapply(1:2, function(v) {
    nCut <- if (v == 1L) K1 - 1L else K2 - 1L
    do.call(
      rbind,
      lapply(2L:nCut, function(c) {
        data.frame(
          label = paste0("x", v, "c", c, ".", seq_len(c - 1L)),
          variable = v,
          rootCut = c,
          childCut = seq_len(c - 1L),
          prob = base *
            0.5 /
            nCut *
            (base / 2^power) *
            0.5 /
            (c - 1L),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
stopifnot(nrow(pairTable) == 13L)
pairCodeOf <- function(variable, rootCut, childCut) {
  variable * 100L + rootCut * 10L + childCut
}
pairTable$code <- pairCodeOf(
  pairTable$variable,
  pairTable$rootCut,
  pairTable$childCut
)

# ---- the retained family ----

retainedRoot <- rootTarget[rootTarget >= massFloor]
retainedLeaf <- leafTarget[leafTarget >= massFloor]
retainedPair <- pairTable[pairTable$prob >= massFloor, ]
stopifnot(
  length(retainedRoot) == 9L,
  length(retainedLeaf) == 6L,
  nrow(retainedPair) == 6L
)

# ---- Holm over a family of two-sided batch-means z ----

holmThresholds <- function(m) {
  qnorm(1 - (familyAlpha / (m - seq_len(m) + 1L)) / 2)
}

holmVerdict <- function(z) {
  m <- length(z)
  # na.last = FALSE puts a degenerate state (an empty subsequence, a state the
  # run never visits) at the most extreme rank, where it is rejected rather
  # than silently stopping the step-down short of itself
  ordered <- order(abs(z), decreasing = TRUE, na.last = FALSE)
  bound <- holmThresholds(m)
  rejected <- logical(m)
  threshold <- numeric(m)
  threshold[ordered] <- bound
  exceeds <- abs(z[ordered]) > bound | !is.finite(z[ordered])
  firstKept <- which(!exceeds)[1L]
  numRejected <- if (is.na(firstKept)) m else firstKept - 1L
  if (numRejected > 0L) {
    rejected[ordered[seq_len(numRejected)]] <- TRUE
  }
  list(rejected = rejected, threshold = threshold)
}

# ---- batch-means MC error, pooled over independent chains ----
#
# Same estimator as change-balance.R's batchMeanSE at one chain; the four
# chains pool as the mean of their estimates over sqrt(sum_c s_c^2 / n) / 4.

batchMeanSE <- function(v, numBatches) {
  len <- (length(v) %/% numBatches) * numBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = numBatches))
  sd(bm) / sqrt(numBatches)
}

batchMean <- function(v, numBatches) {
  len <- (length(v) %/% numBatches) * numBatches
  mean(v[seq_len(len)])
}

poolChains <- function(series, numBatches) {
  stopifnot(all(lengths(series) >= numBatches))
  estimates <- vapply(series, batchMean, 0, numBatches = numBatches)
  errors <- vapply(series, batchMeanSE, 0, numBatches = numBatches)
  list(
    estimate = mean(estimates),
    se = sqrt(sum(errors^2)) / length(errors)
  )
}

# The ladder sbc.R uses: the first lag at which the autocorrelation falls
# under 0.1, NA when it does not within acfLagMax.
firstUnder <- function(v) {
  if (sd(v) <= 0) {
    return(0L)
  }
  a <- acf(v, lag.max = acfLagMax, plot = FALSE)$acf[,, 1]
  hit <- which(a < 0.1)[1L]
  if (is.na(hit)) NA_integer_ else hit - 1L
}

acfLadder <- function(seriesByChain) {
  lags <- vapply(
    seriesByChain,
    function(series) {
      vapply(
        series,
        function(v) firstUnder(v[seq_len(min(length(v), acfDraws))]),
        0L
      )
    },
    numeric(length(seriesByChain[[1L]]))
  )
  apply(lags, 1L, function(row) if (anyNA(row)) NA_integer_ else max(row))
}

# ---- engine arms ----
#
# One sampler of four chains at one thread, so the chains run sequentially off
# independent streams. getTrees rows are chain-major, sample-ordered and
# preorder within a sample, so the root is a sample's first row and the left
# child its second.

runChains <- function(x, y, seed, numBurn, numKept, blockSize, mask) {
  ctl <- dbartsControl(
    n.chains = nChains,
    n.threads = 1L,
    n.trees = 1L,
    useQuantiles = TRUE,
    keepTrees = TRUE,
    n.samples = blockSize,
    n.burn = numBurn,
    n.thin = nThin,
    updateState = TRUE,
    seed = seed,
    n.cuts = 100L
  )
  sampler <- dbarts(
    x,
    y,
    control = ctl,
    tree.prior = cgm(power, base),
    node.prior = normal(kLeaf),
    resid.prior = fixed(1),
    proposal.probs = proposalProbs
  )
  stopifnot(is.null(sampler$data@offset))
  if (mask) {
    sampler$setWeights(rep(0, nrow(x)))
  }
  numBlocks <- as.integer(numKept %/% blockSize)
  stopifnot(numBlocks * blockSize == numKept)
  rootCode <- matrix(0L, numKept, nChains)
  leafCode <- matrix(0L, numKept, nChains)
  pairCode <- matrix(0L, numKept, nChains)
  for (b in seq_len(numBlocks)) {
    if (b == 1L) {
      invisible(sampler$run(numBurn, blockSize))
    } else {
      invisible(sampler$run(0L, blockSize))
    }
    block <- extractBlock(sampler$getTrees(), blockSize)
    rows <- ((b - 1L) * blockSize + 1L):(b * blockSize)
    rootCode[rows, ] <- block$root
    leafCode[rows, ] <- block$leaf
    pairCode[rows, ] <- block$pair
  }
  list(root = rootCode, leaf = leafCode, pair = pairCode)
}

extractBlock <- function(tr, blockSize) {
  total <- nChains * blockSize
  index <- (tr$chain - 1L) * blockSize + tr$sample
  stopifnot(!is.unsorted(index))
  counts <- tabulate(index, nbins = total)
  stopifnot(all(counts > 0L))
  start <- cumsum(c(1L, counts[-total]))

  rootVar <- tr$var[start]
  rootValue <- tr$value[start]
  isSplit <- rootVar != -1L
  rootCut <- integer(total)
  rootCut[rootVar == 1L] <- match(rootValue[rootVar == 1L], cuts1)
  rootCut[rootVar == 2L] <- match(rootValue[rootVar == 2L], cuts2)
  stopifnot(all(rootCut[isSplit] > 0L), !anyNA(rootCut))

  root <- integer(total)
  root[!isSplit] <- 1L
  root[rootVar == 1L] <- 1L + rootCut[rootVar == 1L]
  root[rootVar == 2L] <- 1L + (K1 - 1L) + rootCut[rootVar == 2L]

  leaf <- pmin(tabulate(index[tr$var == -1L], nbins = total), 6L)

  childVar <- integer(total)
  childValue <- rep(NA_real_, total)
  childVar[isSplit] <- tr$var[start[isSplit] + 1L]
  childValue[isSplit] <- tr$value[start[isSplit] + 1L]
  same <- isSplit & childVar == rootVar
  childCut <- integer(total)
  childCut[same & rootVar == 1L] <- match(
    childValue[same & rootVar == 1L],
    cuts1
  )
  childCut[same & rootVar == 2L] <- match(
    childValue[same & rootVar == 2L],
    cuts2
  )
  stopifnot(all(childCut[same] > 0L), !anyNA(childCut))
  pair <- integer(total)
  pair[same] <- pairCodeOf(rootVar[same], rootCut[same], childCut[same])

  list(
    root = matrix(root, blockSize, nChains),
    leaf = matrix(leaf, blockSize, nChains),
    pair = matrix(pair, blockSize, nChains)
  )
}

indicatorSeries <- function(codes, value) {
  lapply(seq_len(ncol(codes)), function(c) as.numeric(codes[, c] == value))
}

# ---- the prior-only arm ----

set.seed(priorSeed)
priorGrid <- expand.grid(x1 = seq_len(K1), x2 = seq_len(K2))
priorX <- cbind(x1 = as.double(priorGrid$x1), x2 = as.double(priorGrid$x2))
priorY <- rnorm(nrow(priorX))

priorKept <- nBatches * priorBatchLen
rootScored <- if (poison == 0L) rootTarget else poisonRootTarget(poison)

cat(sprintf(
  "perturb-balance gate  (%s%s)\n",
  if (quick) "quick" else "full",
  if (poison == 0L) "" else sprintf(", POISON %d target", poison)
))
cat(sprintf(
  paste0(
    "  prior-only arm: %d rows, x1 %d cuts, x2 %d cuts,",
    " %d chains x %d kept at thin %d\n"
  ),
  nrow(priorX),
  K1 - 1L,
  K2 - 1L,
  nChains,
  priorKept,
  nThin
))
cat(sprintf(
  "  mixture: birth_death %.2f, change %.2f, perturb %.2f; window w = %d\n",
  proposalProbs[["birth_death"]],
  proposalProbs[["change"]],
  proposalProbs[["perturb"]],
  windowWidth
))

priorDraws <- NULL
priorLags <- NULL
burnSweeps <- priorBurn
scored <- FALSE
for (attempt in seq_len(maxBurnDoublings + 1L)) {
  priorDraws <- runChains(
    priorX,
    priorY,
    priorSeed,
    burnSweeps,
    priorKept,
    priorBlock,
    mask = TRUE
  )
  ladderSeries <- lapply(seq_len(nChains), function(c) {
    series <- list(leafCount = as.numeric(priorDraws$leaf[, c]))
    for (name in names(retainedRoot)) {
      series[[paste0("root:", name)]] <-
        as.numeric(priorDraws$root[, c] == match(name, rootLabels))
    }
    for (name in names(retainedLeaf)) {
      series[[paste0("leaf:", name)]] <-
        as.numeric(priorDraws$leaf[, c] == match(name, leafLabels))
    }
    for (i in seq_len(nrow(retainedPair))) {
      series[[paste0("pair:", retainedPair$label[i])]] <-
        as.numeric(priorDraws$pair[, c] == retainedPair$code[i])
    }
    series
  })
  priorLags <- acfLadder(ladderSeries)
  burnKept <- burnSweeps %/% nThin
  worst <- if (anyNA(priorLags)) NA_integer_ else max(priorLags)
  cat(sprintf(
    "  burn-in %d sweeps (%d kept); slowest acf<0.1 lag %s (%s), need <= %d\n",
    burnSweeps,
    burnKept,
    if (is.na(worst)) "NONE within lag.max" else format(worst),
    if (is.na(worst)) "-" else names(priorLags)[which.max(priorLags)],
    burnKept %/% acfBurnRatio
  ))
  if (!is.na(worst) && burnKept >= acfBurnRatio * worst) {
    scored <- TRUE
    break
  }
  burnSweeps <- 2L * burnSweeps
}
if (!scored) {
  cat("\nFAIL: burn-in never reached fifty autocorrelation times; not scored\n")
  quit(status = 1L, save = "no")
}

# ---- score the family ----

scoreFamily <- function(labels, targets, seriesByState, numBatches) {
  m <- length(labels)
  estimate <- numeric(m)
  se <- numeric(m)
  for (i in seq_len(m)) {
    pooled <- poolChains(seriesByState[[i]], numBatches)
    estimate[i] <- pooled$estimate
    se[i] <- pooled$se
  }
  z <- (estimate - targets) / se
  verdict <- holmVerdict(z)
  data.frame(
    label = labels,
    engine = estimate,
    target = targets,
    se = se,
    z = z,
    threshold = verdict$threshold,
    rejected = verdict$rejected,
    stringsAsFactors = FALSE
  )
}

priorLabels <- unname(c(
  paste0("root ", names(retainedRoot)),
  paste0("leaf ", names(retainedLeaf)),
  paste0("pair ", retainedPair$label)
))
priorTargets <- unname(c(
  rootScored[names(retainedRoot)],
  retainedLeaf,
  retainedPair$prob
))
priorSeries <- c(
  lapply(names(retainedRoot), function(name) {
    indicatorSeries(priorDraws$root, match(name, rootLabels))
  }),
  lapply(names(retainedLeaf), function(name) {
    indicatorSeries(priorDraws$leaf, match(name, leafLabels))
  }),
  lapply(retainedPair$code, function(code) {
    indicatorSeries(priorDraws$pair, code)
  })
)
priorResult <- scoreFamily(priorLabels, priorTargets, priorSeries, nBatches)

reportFamily <- function(result, title) {
  cat(sprintf("\n%s\n", title))
  cat(sprintf(
    "%-14s %10s %10s %10s %8s %8s\n",
    "state",
    "engine",
    "target",
    "MCse",
    "z",
    "holm"
  ))
  for (i in seq_len(nrow(result))) {
    cat(sprintf(
      "%-14s %10.5f %10.5f %10.5f %+8.2f %8.2f%s\n",
      result$label[i],
      result$engine[i],
      result$target[i],
      result$se[i],
      result$z[i],
      result$threshold[i],
      if (result$rejected[i]) " <- REJECT" else ""
    ))
  }
  cat(sprintf(
    "  worst |z| = %.2f at %s; %d of %d rejected\n",
    max(abs(result$z)),
    result$label[which.max(abs(result$z))],
    sum(result$rejected),
    nrow(result)
  ))
}

reportFamily(
  priorResult,
  sprintf(
    "prior-only arm (m = %d, Holm at alpha %.2f, |z| bound %.2f to %.2f)",
    nrow(priorResult),
    familyAlpha,
    min(priorResult$threshold),
    max(priorResult$threshold)
  )
)
priorOk <- !any(priorResult$rejected)

# ---- the confirmation arm: exact posterior on the same grid ----
#
# Regions are rectangles in (x1 value index, x2 value index); the recursion is
# exact rather than depth-truncated because 5 + 3 cuts cap every path at depth
# 8. A full factorial leaves every rectangle occupied, so no split is refused
# by the empty-leaf veto and none has to be skipped here.

makeLogIL <- function(residVar) {
  priorPrecision <- (kLeaf / nodeScale)^2
  function(n, S, SS) {
    if (n == 0) {
      return(0)
    }
    posteriorPrecision <- n / residVar
    mean <- S / n
    centeredSumOfSquares <- SS - S * mean
    0.5 *
      log(priorPrecision / (priorPrecision + posteriorPrecision)) -
      0.5 * centeredSumOfSquares / residVar -
      0.5 *
        ((priorPrecision * mean) * (posteriorPrecision * mean)) /
        (priorPrecision + posteriorPrecision)
  }
}

exactRootMarginal <- function(x1, x2, y) {
  yRange <- max(y) - min(y)
  yScaled <- (y - min(y)) / yRange - 0.5
  logIL <- makeLogIL((1 / yRange)^2)
  i1 <- match(x1, seq_len(K1))
  i2 <- match(x2, seq_len(K2))
  cellN <- matrix(0, K1, K2)
  cellS <- matrix(0, K1, K2)
  cellSS <- matrix(0, K1, K2)
  for (o in seq_along(y)) {
    cellN[i1[o], i2[o]] <- cellN[i1[o], i2[o]] + 1
    cellS[i1[o], i2[o]] <- cellS[i1[o], i2[o]] + yScaled[o]
    cellSS[i1[o], i2[o]] <- cellSS[i1[o], i2[o]] + yScaled[o]^2
  }
  stopifnot(all(cellN > 0))
  prefix <- function(m) {
    p <- matrix(0, K1 + 1L, K2 + 1L)
    p[-1L, -1L] <- m
    p <- t(apply(p, 1L, cumsum))
    apply(p, 2L, cumsum)
  }
  pn <- prefix(cellN)
  ps <- prefix(cellS)
  pss <- prefix(cellSS)
  rect <- function(P, a1, b1, a2, b2) {
    P[b1 + 1L, b2 + 1L] - P[a1, b2 + 1L] - P[b1 + 1L, a2] + P[a1, a2]
  }
  regionIL <- function(a1, b1, a2, b2) {
    logIL(
      rect(pn, a1, b1, a2, b2),
      rect(ps, a1, b1, a2, b2),
      rect(pss, a1, b1, a2, b2)
    )
  }
  memo <- new.env(parent = emptyenv())
  logM <- function(a1, b1, a2, b2, depth) {
    key <- paste(a1, b1, a2, b2, depth, sep = ",")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    avail1 <- b1 - a1
    avail2 <- b2 - a2
    numAvail <- (avail1 > 0L) + (avail2 > 0L)
    growth <- if (numAvail == 0L) 0 else base / (1 + depth)^power
    terms <- log(1 - growth) + regionIL(a1, b1, a2, b2)
    if (avail1 > 0L) {
      for (c in a1:(b1 - 1L)) {
        terms <- c(
          terms,
          log(growth) -
            log(numAvail) -
            log(avail1) +
            logM(a1, c, a2, b2, depth + 1L) +
            logM(c + 1L, b1, a2, b2, depth + 1L)
        )
      }
    }
    if (avail2 > 0L) {
      for (c in a2:(b2 - 1L)) {
        terms <- c(
          terms,
          log(growth) -
            log(numAvail) -
            log(avail2) +
            logM(a1, b1, a2, c, depth + 1L) +
            logM(a1, b1, c + 1L, b2, depth + 1L)
        )
      }
    }
    result <- logSumExp(terms)
    memo[[key]] <- result
    result
  }
  logWeights <- c(stump = log(1 - base) + regionIL(1L, K1, 1L, K2))
  for (c in seq_len(K1 - 1L)) {
    logWeights[paste0("x1c", c)] <- log(base) -
      log(2) -
      log(K1 - 1L) +
      logM(1L, c, 1L, K2, 1L) +
      logM(c + 1L, K1, 1L, K2, 1L)
  }
  for (c in seq_len(K2 - 1L)) {
    logWeights[paste0("x2c", c)] <- log(base) -
      log(2) -
      log(K2 - 1L) +
      logM(1L, K1, 1L, c, 1L) +
      logM(1L, K1, c + 1L, K2, 1L)
  }
  exp(logWeights - logSumExp(logWeights))
}

confirmOk <- NA
if (poison == 0L) {
  set.seed(confirmSeed)
  confirmGrid <- expand.grid(x1 = seq_len(K1), x2 = seq_len(K2))
  confirmGrid <- confirmGrid[
    rep(seq_len(nrow(confirmGrid)), each = confirmPer),
  ]
  confirmX1 <- as.double(confirmGrid$x1)
  confirmX2 <- as.double(confirmGrid$x2)
  confirmSignal <- 0.5 * (confirmX1 > 3) + 0.3 * (confirmX2 > 2)
  confirmY <- confirmSignal + rnorm(length(confirmX1))
  confirmX <- cbind(x1 = confirmX1, x2 = confirmX2)
  confirmKept <- nBatches * confirmBatchLen

  exactRoot <- exactRootMarginal(confirmX1, confirmX2, confirmY)
  cat(sprintf(
    "\n  confirmation arm: %d rows, positive weights, %d chains x %d kept\n",
    nrow(confirmX),
    nChains,
    confirmKept
  ))
  confirmDraws <- runChains(
    confirmX,
    confirmY,
    confirmSeed,
    confirmBurn,
    confirmKept,
    confirmBlock,
    mask = FALSE
  )

  # the WITHIN-VARIABLE cut law, the quantity the cut move governs: the
  # subsequence of draws whose root splits that variable, exactly as
  # swap-balance.R scores its A/B subsequence
  confirmLabels <- character(0)
  confirmTargets <- numeric(0)
  confirmSeries <- list()
  confirmPrior <- numeric(0)
  for (v in 1:2) {
    nCut <- if (v == 1L) K1 - 1L else K2 - 1L
    states <- match(paste0("x", v, "c", seq_len(nCut)), rootLabels)
    conditional <- exactRoot[states] / sum(exactRoot[states])
    for (k in seq_len(nCut)) {
      if (conditional[k] < massFloor) {
        next
      }
      confirmLabels <- c(confirmLabels, sprintf("x%dc%d|x%d", v, k, v))
      confirmTargets <- c(confirmTargets, conditional[k])
      confirmPrior <- c(confirmPrior, 1 / nCut)
      confirmSeries[[length(confirmSeries) + 1L]] <- lapply(
        seq_len(nChains),
        function(c) {
          sub <- confirmDraws$root[, c]
          sub <- sub[sub %in% states]
          as.numeric(sub == states[k])
        }
      )
    }
  }
  confirmResult <- scoreFamily(
    confirmLabels,
    confirmTargets,
    confirmSeries,
    nBatches
  )
  reportFamily(
    confirmResult,
    sprintf(
      "confirmation arm, within-variable cut law (m = %d, Holm at alpha %.2f)",
      nrow(confirmResult),
      familyAlpha
    )
  )
  powerZ <- (confirmResult$engine - confirmPrior) / confirmResult$se
  cat(sprintf(
    "  power: |z| vs the PRIOR conditional (a dropped likelihood) up to %.1f\n",
    max(abs(powerZ))
  ))
  confirmOk <- !any(confirmResult$rejected)
} else {
  cat("\n  confirmation arm skipped: the poison arms score statistic 1 alone\n")
}

# ---- verdict ----

cat("\n================ VERDICT ================\n")
cat(sprintf(
  "  prior-only arm      : %s (worst |z| = %.2f, %d of %d Holm rejections)\n",
  if (priorOk) "PASS" else "FAIL",
  max(abs(priorResult$z)),
  sum(priorResult$rejected),
  nrow(priorResult)
))
if (poison == 0L) {
  cat(sprintf(
    "  confirmation arm    : %s (worst |z| = %.2f, %d of %d Holm rejections)\n",
    if (confirmOk) "PASS" else "FAIL",
    max(abs(confirmResult$z)),
    sum(confirmResult$rejected),
    nrow(confirmResult)
  ))
}
pass <- priorOk && (poison != 0L || confirmOk)
if (poison != 0L) {
  cat(sprintf(
    "  poison %d target     : the gate must FAIL here, and does: %s\n",
    poison,
    if (pass) "NO" else "yes"
  ))
}
cat(sprintf("\nPERTURB BALANCE GATE: %s\n", if (pass) "PASS" else "FAIL"))

# Exit nonzero on failure so a runner/CI can catch it (matches bd-balance.R).
if (!pass) {
  quit(status = 1L, save = "no")
}
