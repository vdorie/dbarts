#!/usr/bin/env Rscript

# Permanent detailed-balance gate for the RULE_GIBBS move (docs/design/
# nog-gibbs.md section 5). The move replaces the split rule at a nog node - an
# interior node whose two children are both leaves - with a draw from that
# rule's own full conditional. The neighbourhood is closed (it reads ANCESTORS
# only, so it is identical from every state in it), the acceptance is
# identically one, and no proposal count survives; what the kernel has to get
# right is the WEIGHT it draws from, and every term of that weight is a term
# this gate reads:
#
#     log w = S                             the scan's rank-admitted marginal
#           + log P(split variable)         constant off DART
#           - log |SI|                      the node's own rule prior
#           + log(1 - growth(left))
#           + log(1 - growth(right))        the prior strictly below the node
#
# The shipped balance gates read a root-split-VARIABLE marginal
# (change-balance.R), a partition distribution (bd-balance.R), a two-state
# prior ratio (swap-balance.R) and a within-variable cut displacement
# (perturb-balance.R). None of them reads an enumerated rule conditional, and
# none of them would move if the two `1 - growth` factors below the node were
# dropped.
#
# ---- the prior-only arm ----
#
# An all-zero weight vector installed with $setWeights turns the likelihood
# off. Under it every leaf holding rows takes veto rank 1 and every leaf
# holding none rank 2, so the smallest rank the candidates carry is 1, the
# stratum the kernel draws over is exactly the member-occupied candidate set,
# and the scan's rank-admitted marginal S is 0 across all of it: EVERY
# candidate's occupancy rank is the same, the stratum is the whole
# neighbourhood, and the draw is the prior conditional on the rule at that
# node. So this arm exercises the node selection, the enumeration, the rank
# stratum, the rule prior and the two `1 - growth` factors, and it does NOT
# exercise one scan entry - which is why the confirmation arm below is
# mandatory here rather than optional as it is for a Metropolis move.
#
# The kernel is then reversible with respect to the CGM tree prior truncated
# to the member-occupied trees, and the design makes that truncation vacuous:
# two ordinal columns of 6 and 2 distinct values as a FULL FACTORIAL with
# useQuantiles = TRUE, which puts 5 and 1 cuts at the midpoints between
# consecutive distinct values, so every leaf of every reachable tree holds a
# cell and the target is the plain CGM prior. tree.prior = cgm(power = 0.5,
# base = 0.95): the LOPSIDED cut counts, 5 against 1, are what give poison 2
# its size, and the low power is what gives poison 1 its size, deepening trees
# so that a candidate can exhaust a child's last variable. Every chain starts
# at a bare root, so the run needs a burn-in it would not otherwise.
#
# Three statistics, all closed form or a dynamic program:
#   (1) the root's (variable, cut) marginal plus the stump, SEVEN states,
#       P(grow) x P(v) x 1/|SI_v|: 0.095 for each of x1's five cuts, 0.475 for
#       x2's one, 0.05 for the stump;
#   (2) the leaf count, from a dynamic program whose state is (remaining x1
#       cuts, remaining x2 cuts, DEPTH) - growth is base/(1 + depth)^power, so
#       remaining cuts alone does not define the recursion - support 1 to 12
#       with the tail pooled into the >= j bin the program fixes below;
#   (3) the left child's cut CONDITIONAL on the root splitting x2 and that
#       child being itself a nog node, FIVE states. It is the one statistic
#       where the two `1 - growth` factors VARY: the root's only x2 cut
#       exhausts x2, so the child holds x1 alone, and cutting x1 at either end
#       leaves a grandchild with no available variable and a factor of 1 while
#       a middle cut leaves both splittable at 1 - 0.95/3^0.5 = 0.4515. That
#       is why the grid is 6 x 2 rather than perturb-balance.R's 6 x 4, where
#       a depth-1 child still holds cuts of both columns and the two factors
#       are constant across its whole neighbourhood.
# All three draw the split VARIABLE uniformly among the variables that STILL
# have a cut at the node, which is what the engine does with no split
# probabilities set: a variable exhausted below its ancestors drops out of the
# draw rather than wasting it, and the leaf-count vector reproduces only under
# that convention.
#
# States of prior mass below 0.004 are dropped by this PRE-STATED rule, their
# counts being degenerate at any feasible run length: statistic 1 keeps all
# seven, statistic 3 all five, and statistic 2 keeps one to ten leaves plus
# the pooled >= 11 bin (the twelve-leaf state alone is under the floor). Family
# size m = 23, Holm at family alpha 0.05, so the thresholds run from
# |z| = 1.96 on the least extreme test to |z| = 3.07 on the most extreme.
#
# The arm's mixture is birth_death 0.10, change 0.10, rule_gibbs 0.80: change
# is retained because statistic 1 needs a mechanism at a non-nog root and
# birth/death because statistic 2 moves through nothing else, both at the
# smallest share keeping their own statistic non-degenerate while leaving
# rule_gibbs dominant.
#
# The z is a pooled batch-means z: a state's indicator is averaged within a
# chain over B batches of consecutive kept draws, s_c is the sd of chain c's
# batch means, and the four independent chains pool as the mean of the four
# estimates over sqrt(sum_c s_c^2 / B) / 4. Batching rather than a binomial
# formula is what absorbs the within-chain autocorrelation, and it is the
# batch LENGTH that has to clear that autocorrelation, so both modes batch at
# 400 kept draws or more. The burn-in is not assumed either: the script reads
# the first lag at which each scored series' kept-draw autocorrelation falls
# under 0.1 and REFUSES to score a run whose burn-in is under fifty of those
# lags, doubling and re-running instead.
# Statistic 3's draws are CONDITIONAL, so the same refuse-and-double rule
# covers them from the other side: the conditioning event's mass is closed
# form here, and the script refuses to score statistic 3 unless the REALIZED
# conditional effective count clears ten times its own detection floor.
#
# ---- the confirmation arm ----
#
# The prior-only arm never exercises the likelihood term S, which is the whole
# of what the rank-aware scan contributes. The confirmation arm runs the same
# grid and the same mixture with positive weights and NO mask - the mask is
# what makes the likelihood constant, so the arms cannot share a configuration
# - and scores the root (variable, cut) marginal against an exact region
# dynamic program in change-balance.R's shape: a node's reachable set is a
# rectangle in (x1 value index, x2 value index), splitInterval is
# ancestor-only, and no depth truncation is needed because 5 + 1 cuts cap
# every path at depth 6. Calibration is matched to the engine exactly, else
# the comparison is void: quantile cuts, range scaling, fixed(1) residual
# variance -> internal sigma = 1 / range, and the constant-leaf conjugate
# marginal with priorPrecision = (k / scale)^2, scale = node.scale /
# sqrt(ntree) = 0.5. The same Holm machinery scores its states as their own
# family, and the distance to the PRIOR marginal is reported as the arm's
# power: a kernel that dropped S would sit there.
#
# ---- the poison arms ----
#
# Both poisons of the design are engine defects, so each is run here the way
# change-balance.R and swap-balance.R run their wrong targets: the flag
# mutates the SCRIPT's expectation to the law the defective kernel would
# target, and the gate must FAIL against it. Both are leading-order, undiluted
# predictions - the arm's change and birth/death shares pull a really-poisoned
# chain back toward the true law - and each mutates the ONE statistic the
# design requires it to fail.
#   poison1 : the two log(1 - growth(child)) terms dropped. Statistic 3's five
#             states go from (0.2981, 0.1346, 0.1346, 0.1346, 0.2981) to
#             uniform 0.2 - the ends fall 33 percent, the middles rise 49.
#             Statistic 1 is BLIND to this poison by construction, the root's
#             children always retaining the other variable, which is why
#             statistic 3 exists.
#   poison2 : the 1/|SI_v| rule factor dropped. The root's six candidates then
#             carry equal weight: x1's cuts go 0.095 to 0.158333 (+67 percent)
#             and x2's single cut 0.475 to 0.158333 (-67 percent), the same
#             low-cardinality bias in mirror image that change-balance.R's
#             gate repaired. Statistic 3 is blind to it, the depth-1 child
#             holding one variable so that 1/|SI| is constant over its whole
#             neighbourhood.
# The matching engine mutations are carried by benchmarks/R/mutation-battery.R
# (m26, m27), which builds them for real; the neighbourhood's ancestor-only
# closure, which is what removes the reverse count, is asserted in tests/cpp
# instead of poisoned here.
#
# Usage: Rscript rule-gibbs-balance.R [quick] [poison1|poison2]

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
K2 <- 2L # x2 distinct values -> 1 quantile cut
base <- 0.95
power <- 0.5
kLeaf <- 2
nodeScale <- 0.5 # gaussian node.scale default; scale = nodeScale / sqrt(1)

proposalProbs <- c(
  birth_death = 0.10,
  swap = 0,
  change = 0.10,
  perturb = 0,
  rule_gibbs = 0.80,
  birth = 0.5
)

nChains <- 4L
nThin <- 20L
massFloor <- 0.004
familyAlpha <- 0.05
acfLagMax <- 200L
acfBurnRatio <- 50L
acfDraws <- 25000L # kept draws per chain the acf ladder reads
maxBurnDoublings <- 3L
conditionalFloorRatio <- 10 # statistic 3's own refuse-and-double rule

# It is the batch LENGTH that has to clear the series' autocorrelation time,
# not the batch count, so quick trades batches for length rather than
# shortening them: both modes batch at 400 kept draws or more against a
# slowest integrated autocorrelation of about 30.
nBatches <- if (quick) 125L else 500L
priorBatchLen <- if (quick) 400L else 500L
priorBlock <- if (quick) 25000L else 50000L
# Sized off the ladder rather than assumed: the slowest series here is the
# root's x2 cut indicator, not the leaf count, and it reads about 55 kept
# draws to an autocorrelation of 0.1, so fifty of those is 2,750 kept draws
# and the burn-in starts at 5,000. The ladder below still governs.
priorBurn <- 100000L
priorSeed <- 20260907L

confirmBatchLen <- 200L
confirmBlock <- if (quick) 25000L else 50000L
confirmBurn <- 100000L
confirmSeed <- 20260908L
confirmPer <- 2L # replicates per grid cell
confirmSignal1 <- 0.8
confirmSignal2 <- 0.5

cuts1 <- (seq_len(K1)[-K1] + seq_len(K1)[-1L]) / 2
cuts2 <- (seq_len(K2)[-K2] + seq_len(K2)[-1L]) / 2

logSumExp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(v - m)))
}

growthProbability <- function(depth, numAvailable) {
  if (numAvailable == 0L) 0 else base / (1 + depth)^power
}

# ---- statistic 1: the root (variable, cut) marginal, closed form ----
#
# P(grow at the root) = base, the variable is uniform over the two available
# columns, and the rule is uniform over that variable's cuts. Seven states,
# summing to one exactly. The two `1 - growth` factors are CONSTANT over this
# neighbourhood - every root split leaves both children holding the other
# column, so both children grow with the depth-1 probability whatever the rule
# is - which is what makes statistic 1 blind to poison 1 and statistic 3
# necessary.

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

# Poison 2 drops the rule prior, so the six split candidates carry equal
# weight; the stump comes from birth/death and is untouched.
poisonRootTarget <- function() {
  out <- rootTarget
  splits <- setdiff(rootLabels, "stump")
  out[splits] <- base / length(splits)
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
    growth <- growthProbability(depth, numAvail)
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

# The pooling point the program fixes ahead of the run: the FINEST resolution
# at which every retained singleton bin and the pooled tail both clear the
# mass floor.
poolAt <- max(which(vapply(
  seq_along(leafFull),
  function(j) {
    j > 1L &&
      all(leafFull[seq_len(j - 1L)] >= massFloor) &&
      sum(leafFull[j:length(leafFull)]) >= massFloor
  },
  TRUE
)))
leafLabels <- c(
  paste0("leaves", seq_len(poolAt - 1L)),
  paste0("leaves", poolAt, "+")
)
leafTarget <- c(
  leafFull[seq_len(poolAt - 1L)],
  sum(leafFull[poolAt:length(leafFull)])
)
names(leafTarget) <- leafLabels

# ---- statistic 3: the depth-1 child's cut, conditional ----
#
# The root splits x2 at its ONLY cut, so the left child's ancestor-constrained
# rectangle holds all six x1 values and one x2 value: x1 is the child's only
# available variable, the variable draw is forced and 1/|SI| = 1/5 is constant
# over the child's whole neighbourhood. Conditional on that child being a nog
# node - it splits, and both grandchildren are leaves - the only term left
# varying is the pair of `1 - growth` factors, so the law over the five cuts is
# proportional to (1 - growth(left grandchild)) x (1 - growth(right)), the
# factor being 1 where a cut leaves that grandchild no available variable at
# all.

childCuts <- seq_len(K1 - 1L)
childWeights <- vapply(
  childCuts,
  function(c) {
    (1 - growthProbability(2L, min(c - 1L, 1L))) *
      (1 - growthProbability(2L, min(K1 - 1L - c, 1L)))
  },
  0
)
childTarget <- childWeights / sum(childWeights)
childLabels <- paste0("x1c", childCuts, "|nog")
names(childTarget) <- childLabels

# The conditioning event's own mass, which the effective-count refusal reads:
# the root takes its x2 cut, the depth-1 child grows, and its cut leaves both
# grandchildren leaves.
childCondMass <- rootTarget[["x2c1"]] *
  growthProbability(1L, 2L) *
  sum(childWeights) /
  (K1 - 1L)

# Poison 1 drops both factors, leaving the child's five cuts uniform.
poisonChildTarget <- function() {
  out <- rep(1 / length(childTarget), length(childTarget))
  names(out) <- childLabels
  out
}

# ---- the retained family ----

retainedRoot <- rootTarget[rootTarget >= massFloor]
retainedLeaf <- leafTarget[leafTarget >= massFloor]
retainedChild <- childTarget[childTarget >= massFloor]
stopifnot(
  length(retainedRoot) == 7L,
  length(retainedLeaf) == length(leafTarget),
  length(retainedChild) == 5L
)
familySize <- length(retainedRoot) +
  length(retainedLeaf) +
  length(retainedChild)

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

# The undiluted detection floor: the effective draw count at which a shift
# from p to q clears the family's STRICTEST Holm threshold.
detectionFloor <- function(p, q, m) {
  ceiling(max(holmThresholds(m))^2 * p * (1 - p) / (p - q)^2)
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
# preorder within a sample, so a sample's first row is the root, its second the
# root's left child, and - where that child is itself a split whose own left
# child is a leaf - its third and fourth are the two grandchildren.

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
    tree.prior = cgm(power = power, base = base),
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
  childCode <- matrix(0L, numKept, nChains)
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
    childCode[rows, ] <- block$child
  }
  list(root = rootCode, leaf = leafCode, child = childCode)
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

  leaf <- pmin(tabulate(index[tr$var == -1L], nbins = total), poolAt)

  # statistic 3's conditioning event, read straight off the preorder: the root
  # takes x2, its left child splits, and both of that child's children are
  # leaves. A root split puts at least three nodes in the sample and a split
  # left child at least five, so the offsets below are always in range.
  onX2 <- rootVar == 2L
  childVar <- integer(total)
  childVar[onX2] <- tr$var[start[onX2] + 1L]
  childSplits <- onX2 & childVar != -1L
  stopifnot(all(childVar[childSplits] == 1L))
  isNog <- childSplits
  isNog[childSplits] <- tr$var[start[childSplits] + 2L] == -1L &
    tr$var[start[childSplits] + 3L] == -1L
  child <- integer(total)
  child[isNog] <- match(tr$value[start[isNog] + 1L], cuts1)
  stopifnot(all(child[isNog] > 0L), !anyNA(child))

  list(
    root = matrix(root, blockSize, nChains),
    leaf = matrix(leaf, blockSize, nChains),
    child = matrix(child, blockSize, nChains)
  )
}

indicatorSeries <- function(codes, value) {
  lapply(seq_len(ncol(codes)), function(c) as.numeric(codes[, c] == value))
}

# the conditional subsequence, scored as swap-balance.R scores its A/B one
conditionalSeries <- function(codes, value, keep) {
  lapply(seq_len(ncol(codes)), function(c) {
    sub <- codes[, c]
    as.numeric(sub[sub %in% keep] == value)
  })
}

# ---- the prior-only arm ----

set.seed(priorSeed)
priorGrid <- expand.grid(x1 = seq_len(K1), x2 = seq_len(K2))
priorX <- cbind(x1 = as.double(priorGrid$x1), x2 = as.double(priorGrid$x2))
priorY <- rnorm(nrow(priorX))

priorKept <- nBatches * priorBatchLen
rootScored <- if (poison == 2L) poisonRootTarget() else rootTarget
childScored <- if (poison == 1L) poisonChildTarget() else childTarget

cat(sprintf(
  "rule-gibbs-balance gate  (%s%s)\n",
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
  "  mixture: birth_death %.2f, change %.2f, rule_gibbs %.2f; cgm(%.2f, %.2f)\n",
  proposalProbs[["birth_death"]],
  proposalProbs[["change"]],
  proposalProbs[["rule_gibbs"]],
  power,
  base
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
    for (name in names(retainedChild)) {
      series[[paste0("child:", name)]] <-
        as.numeric(priorDraws$child[, c] == match(name, childLabels))
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
  n <- length(labels)
  estimate <- numeric(n)
  se <- numeric(n)
  for (i in seq_len(n)) {
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
  paste0("child ", names(retainedChild))
))
priorTargets <- unname(c(
  rootScored[names(retainedRoot)],
  retainedLeaf,
  childScored[names(retainedChild)]
))
childStates <- match(names(retainedChild), childLabels)
priorSeries <- c(
  lapply(names(retainedRoot), function(name) {
    indicatorSeries(priorDraws$root, match(name, rootLabels))
  }),
  lapply(names(retainedLeaf), function(name) {
    indicatorSeries(priorDraws$leaf, match(name, leafLabels))
  }),
  lapply(childStates, function(state) {
    conditionalSeries(priorDraws$child, state, childStates)
  })
)
priorResult <- scoreFamily(priorLabels, priorTargets, priorSeries, nBatches)

reportFamily <- function(result, title) {
  cat(sprintf("\n%s\n", title))
  cat(sprintf(
    "%-16s %10s %10s %10s %8s %8s\n",
    "state",
    "engine",
    "target",
    "MCse",
    "z",
    "holm"
  ))
  for (i in seq_len(nrow(result))) {
    cat(sprintf(
      "%-16s %10.5f %10.5f %10.5f %+8.2f %8.2f%s\n",
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

# ---- statistic 3's own adequacy rule ----
#
# The conditional draws are the only ones the run does not control directly, so
# the script refuses to score them unless the realized effective count clears
# ten times poison 1's undiluted detection floor on that state.

childRows <- match(paste0("child ", names(retainedChild)), priorResult$label)
childFloor <- detectionFloor(
  childTarget[names(retainedChild)],
  poisonChildTarget()[names(retainedChild)],
  familySize
)
childEffective <- priorResult$engine[childRows] *
  (1 - priorResult$engine[childRows]) /
  priorResult$se[childRows]^2
cat(sprintf(
  "\n  statistic 3: conditioning mass %.4f, realized %.4f of kept draws\n",
  childCondMass,
  mean(priorDraws$child > 0L)
))
cat(sprintf(
  "  effective count %.0f to %.0f, floor %d to %d, need >= %.0f\n",
  min(childEffective),
  max(childEffective),
  min(childFloor),
  max(childFloor),
  conditionalFloorRatio * max(childFloor)
))
childAdequate <- all(childEffective >= conditionalFloorRatio * childFloor)
if (!childAdequate) {
  cat("\nFAIL: statistic 3's conditional draws are too few to score\n")
  quit(status = 1L, save = "no")
}
priorOk <- !any(priorResult$rejected)

# ---- the confirmation arm: exact posterior on the same grid ----
#
# Regions are rectangles in (x1 value index, x2 value index); the recursion is
# exact rather than depth-truncated because 5 + 1 cuts cap every path at depth
# 6. A full factorial leaves every rectangle occupied, so no split is refused
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
    growth <- growthProbability(depth, numAvail)
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
  confirmSignal <- confirmSignal1 *
    (confirmX1 > 3) +
    confirmSignal2 * (confirmX2 > 1)
  confirmY <- confirmSignal + rnorm(length(confirmX1))
  confirmX <- cbind(x1 = confirmX1, x2 = confirmX2)
  confirmKept <- nBatches * confirmBatchLen

  exactRoot <- exactRootMarginal(confirmX1, confirmX2, confirmY)
  stopifnot(
    identical(names(exactRoot), rootLabels),
    abs(sum(exactRoot) - 1) < 1e-12
  )
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

  keep <- exactRoot >= massFloor
  confirmLabels <- rootLabels[keep]
  confirmTargets <- unname(exactRoot[keep])
  confirmPrior <- unname(rootTarget[keep])
  confirmSeries <- lapply(which(keep), function(state) {
    indicatorSeries(confirmDraws$root, state)
  })
  confirmResult <- scoreFamily(
    confirmLabels,
    confirmTargets,
    confirmSeries,
    nBatches
  )
  reportFamily(
    confirmResult,
    sprintf(
      "confirmation arm, exact root marginal (m = %d, Holm at alpha %.2f)",
      nrow(confirmResult),
      familyAlpha
    )
  )
  powerZ <- (confirmResult$engine - confirmPrior) / confirmResult$se
  cat(sprintf(
    "  power: |z| vs the PRIOR marginal (a dropped likelihood) up to %.1f\n",
    max(abs(powerZ))
  ))
  confirmOk <- !any(confirmResult$rejected)
} else {
  cat("\n  confirmation arm skipped: the poison arms score the prior law\n")
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
cat(sprintf("\nRULE_GIBBS BALANCE GATE: %s\n", if (pass) "PASS" else "FAIL"))

# Exit nonzero on failure so a runner/CI can catch it (matches bd-balance.R).
if (!pass) {
  quit(status = 1L, save = "no")
}
