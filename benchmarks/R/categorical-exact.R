#!/usr/bin/env Rscript

# Exact-posterior gate for bartcore's factor splits, in two arms over one
# single-tree, one-factor-predictor (4 levels) probit problem whose tree
# space admits brute-force enumeration, with 1-D quadrature for the leaf
# marginals.
#
# Typed CATEGORICAL, the rules are direction assignments over the categories
# reachable at each node, uniform over the 2^R - 2 assignments that leave
# neither side empty, both orientations counted. Failure means the
# categorical rule prior, draw, or partition logic is wrong.
#
# Typed ORDERED FACTOR, the rules are thresholds on the level order, uniform
# over the boundaries the ancestors leave available, so a node holding levels
# [low, high] splits the interval and its children inherit the two halves.
# The enumeration is over TREES, not over the partitions they induce: the
# reachable leaf partitions are the contiguous ones (the subsets of the K - 1
# interior boundaries) and several trees realize the same one, so each tree
# carries its own CGM mass and the marginalization is the sum over trees. The
# arm runs twice, at a cut cap above and below the level count, and both must
# match the same posterior - the grid follows the level table, so the cap does
# not apply. Failure means the ordered-factor grid, its cap treatment, or the
# ordinal rule prior is wrong.
#
# Both arms always run and report; the script exits non-zero if either failed.
#
# Usage: Rscript categorical-exact.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
ndpost <- if (quick) 10000L else 50000L
n.seeds <- if (quick) 1L else 3L
tolerance <- if (quick) 0.01 else 0.005

set.seed(1001L)
numCategories <- 4L
nPerCategory <- 30L
category <- rep(0:(numCategories - 1L), each = nPerCategory)
x <- matrix(as.double(category), ncol = 1L)
categoryProbability <- c(0.25, 0.7, 0.45, 0.6)
y <- rbinom(length(category), 1L, categoryProbability[category + 1L])
n <- length(y)
offset <- qnorm(mean(y))
x.test <- matrix(as.double(0:(numCategories - 1L)), ncol = 1L)

base <- 0.5
power <- 2
k <- 2
tau <- 3.0 / k # dbarts binary node scale over a single tree

successesByCategory <- as.vector(tapply(y, category, sum))
countsByCategory <- as.vector(table(category))

popcount <- function(mask) sum(as.integer(intToBits(mask)))

# trees over a reachable-category mask: leaf w.p. 1 - base / (1 + d)^power
# (forced when fewer than two categories reach the node), else an assignment
# sending a proper nonempty subset right, uniform over the 2^R - 2 of them
enumerate <- function(mask, depth) {
  numReachable <- popcount(mask)
  growth <- if (numReachable >= 2L) base / (1 + depth)^power else 0
  result <- list(list(leaves = list(mask), logPrior = log(1 - growth)))
  if (numReachable < 2L) {
    return(result)
  }

  subsets <- Filter(
    function(d) d != 0L && d != mask && bitwAnd(d, bitwNot(mask)) == 0L,
    seq_len(15L)
  )
  for (directions in subsets) {
    lefts <- enumerate(bitwAnd(mask, bitwNot(directions)), depth + 1L)
    rights <- enumerate(directions, depth + 1L)
    for (left in lefts) {
      for (right in rights) {
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

trees <- enumerate(bitwShiftL(1L, numCategories) - 1L, 0L)
cat(sprintf("enumerated %d trees\n", length(trees)))

muGrid <- seq(-10, 10, by = 0.005)
muDensity <- dnorm(muGrid, 0, tau)
leafQuantities <- function(s, m) {
  w <- muDensity *
    exp(
      s *
        pnorm(muGrid + offset, log.p = TRUE) +
        (m - s) * pnorm(-(muGrid + offset), log.p = TRUE)
    )
  list(
    logMarginal = log(sum(w) * 0.005),
    meanProbability = sum(pnorm(muGrid + offset) * w) / sum(w)
  )
}

logWeights <- numeric(length(trees))
predictions <- matrix(0, length(trees), numCategories)
for (t in seq_along(trees)) {
  logWeight <- trees[[t]]$logPrior
  for (leafMask in trees[[t]]$leaves) {
    inLeaf <- which(
      bitwAnd(bitwShiftL(1L, 0:(numCategories - 1L)), leafMask) != 0L
    )
    q <- leafQuantities(
      sum(successesByCategory[inLeaf]),
      sum(countsByCategory[inLeaf])
    )
    logWeight <- logWeight + q$logMarginal
    predictions[t, inLeaf] <- q$meanProbability
  }
  logWeights[t] <- logWeight
}
weights <- exp(logWeights - max(logWeights))
exact <- colSums(weights * predictions) / sum(weights)

fitBartcore <- function(seed) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    updateState = FALSE
  )
  host <- dbarts(
    x,
    y,
    test = x.test,
    offset = offset,
    control = control,
    node.prior = normal(k),
    tree.prior = cgm(power, base)
  )
  # a bare code matrix carries no level table; type the column by hand
  host$data@varTypes[1L] <- 1L
  host$data@offset.test <- NULL
  bc <- dbarts:::bartcoreSampler(host)
  r <- bartcoreRun(bc, 5000L, ndpost)
  colMeans(pnorm(t(r$test) + offset))
}

fit <- colMeans(do.call(rbind, lapply(seq_len(n.seeds), fitBartcore)))

cat(sprintf("%9s %s\n", "exact", paste(sprintf("%.4f", exact), collapse = " ")))
cat(sprintf("%9s %s\n", "sampler", paste(sprintf("%.4f", fit), collapse = " ")))
gap <- max(abs(fit - exact))
cat(sprintf("max gap %.4f%s\n", gap, if (gap > tolerance) " <- FAIL" else ""))

failures <- character()
if (gap > tolerance) {
  failures <- c(failures, "categorical")
} else {
  cat("OK: categorical sampler matches the exact posterior\n")
}

# ------------------------------------------------------------------------
# Ordered-factor arm. The same design typed as an ORDERED factor, whose grid
# is the K - 1 midpoints between consecutive declared level codes, so its
# rules are thresholds on the level order rather than direction assignments
# over the levels. A node holding levels [low, high] draws uniformly over the
# high - low interior boundaries the ancestors leave available, which
# enumerates exactly as the categorical space above: over TREES, each with its
# own CGM mass. Several trees induce the same contiguous leaf partition and
# the leaf marginals depend only on that partition, so summing per tree is
# what marginalizes correctly - the partitions are NOT equally likely, and a
# singleton child's empty cut interval (a forced leaf) is why. Both n.cuts
# arms below must match the SAME exact posterior: the grid follows the level
# table, so n.cuts does not apply to the column, including when it sits
# under K - 1.

enumerateOrdered <- function(low, high, depth) {
  numAvailable <- high - low
  growth <- if (numAvailable >= 1L) base / (1 + depth)^power else 0
  result <- list(list(leaves = list(low:high), logPrior = log(1 - growth)))
  if (numAvailable < 1L) {
    return(result)
  }

  for (cut in low:(high - 1L)) {
    lefts <- enumerateOrdered(low, cut, depth + 1L)
    rights <- enumerateOrdered(cut + 1L, high, depth + 1L)
    for (left in lefts) {
      for (right in rights) {
        result[[length(result) + 1L]] <- list(
          leaves = c(left$leaves, right$leaves),
          logPrior = log(growth) -
            log(numAvailable) +
            left$logPrior +
            right$logPrior
        )
      }
    }
  }
  result
}

orderedTrees <- enumerateOrdered(0L, numCategories - 1L, 0L)
cat(sprintf("enumerated %d ordered-factor trees\n", length(orderedTrees)))

orderedLogWeights <- numeric(length(orderedTrees))
orderedPredictions <- matrix(0, length(orderedTrees), numCategories)
for (t in seq_along(orderedTrees)) {
  logWeight <- orderedTrees[[t]]$logPrior
  for (levelSet in orderedTrees[[t]]$leaves) {
    inLeaf <- levelSet + 1L
    q <- leafQuantities(
      sum(successesByCategory[inLeaf]),
      sum(countsByCategory[inLeaf])
    )
    logWeight <- logWeight + q$logMarginal
    orderedPredictions[t, inLeaf] <- q$meanProbability
  }
  orderedLogWeights[t] <- logWeight
}
orderedWeights <- exp(orderedLogWeights - max(orderedLogWeights))
orderedExact <-
  colSums(orderedWeights * orderedPredictions) / sum(orderedWeights)

fitOrdered <- function(seed, n.cuts) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    updateState = FALSE
  )
  host <- dbarts(
    x,
    y,
    test = x.test,
    offset = offset,
    control = control,
    node.prior = normal(k),
    tree.prior = cgm(power, base)
  )
  # a bare code matrix carries no level table; type the column by hand
  host$data@varTypes[1L] <- 2L
  host$data@n.cuts[1L] <- n.cuts
  host$data@offset.test <- NULL
  bc <- dbarts:::bartcoreSampler(host)
  r <- bartcoreRun(bc, 5000L, ndpost)
  colMeans(pnorm(t(r$test) + offset))
}

cat(sprintf(
  "%9s %s\n",
  "exact",
  paste(sprintf("%.4f", orderedExact), collapse = " ")
))
orderedGap <- 0
for (numCuts in c(100L, 2L)) {
  orderedFit <- colMeans(do.call(
    rbind,
    lapply(seq_len(n.seeds), fitOrdered, n.cuts = numCuts)
  ))
  gap <- max(abs(orderedFit - orderedExact))
  orderedGap <- max(orderedGap, gap)
  cat(sprintf(
    "%9s %s  (n.cuts %3d, max gap %.4f%s)\n",
    "sampler",
    paste(sprintf("%.4f", orderedFit), collapse = " "),
    numCuts,
    gap,
    if (gap > tolerance) " <- FAIL" else ""
  ))
}

if (orderedGap > tolerance) {
  failures <- c(failures, "ordered factor")
} else {
  cat("OK: ordered-factor sampler matches the exact posterior\n")
}

if (length(failures) > 0L) {
  cat(sprintf(
    "FAIL: %s arm(s) miss the exact posterior\n",
    paste(failures, collapse = ", ")
  ))
  quit(status = 1L)
}
