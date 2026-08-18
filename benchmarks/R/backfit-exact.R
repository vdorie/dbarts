#!/usr/bin/env Rscript

# Conditional-exactness gate for the gaussian backfit at ensemble scale. The
# Gibbs step that draws tree j's leaves conditions on every other tree as an
# offset: with the partition settled, each leaf value is normal with precision
#   (k sqrt(m) / node.scale)^2 + sum_leaf w_i / sigma^2
# and mean (sum_leaf w_i r_i / sigma^2) divided by that precision, r being the
# internally scaled response minus every OTHER tree's current fit
# (model.hpp, ConstantGaussianLeaf::drawFromPosterior).
#
# Nothing is frozen by hand - the sweep's own ordering supplies the
# conditioning. Entering tree j the engine's running residual holds the
# response minus trees 1..j-1 at THIS sweep and trees j+1..m at the LAST, and
# consecutive keepTrees records carry both, so every tree in the forest is
# checkable rather than one. The check is the standardized draw
#   t = (recorded leaf value - posterior mean) / posterior sd,
# recomputed in R from the recorded trees: t is the standard normal the engine
# consumed for that leaf. Everything the reference reads precedes the draw, so
# each t is N(0, 1) given all earlier ones - the t's are iid N(0, 1) across
# leaves, trees and sweeps, and exact moment and KS tests apply. Conditional
# validity does not need the chain to have converged, and stratifying on leaf
# size keeps it (the count is known before the draw), which is what gives the
# leaf prior its own arm: the prior carries ~2% of the posterior precision in
# a 700-observation leaf and about half in a 15-observation one, so only small
# leaves see it. A rare strong indicator in the design guarantees they exist.
#
# Two deterministic preconditions ride along, and fire before the
# distributional arms could: the R-side routing must reproduce every node's
# recorded observation count, and the per-tree fits must sum to the recorded
# training fit.
#
# Sensitivity, measured by poisoning the reference rather than the engine, on
# arm A at full size (193k draws): a residual that misses the current sweep's
# already-drawn trees reads z(var) = +434 and KS p = 0; dropping sqrt(m) from
# the leaf prior precision reads z(mean) = -9.3 pooled and -19.3 on
# prior-informed leaves, which is why that stratum is reported apart.
#
# Not seen here: the tree-structure Metropolis step (the realized partition is
# conditioned on, never checked - that is the exact-posterior gates' job), the
# sigma and k draws themselves, and a one-sweep lag error in which sigma the
# leaf draw reads (adjacent sigmas differ by well under a percent, below this
# oracle's resolution). Gaussian constant leaves only.
#
# Usage: Rscript backfit-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

kLeaf <- 2
nodeScale <- 0.5 # gaussian node.scale; leaf scale = node.scale / sqrt(m)
zBound <- 4
ksBound <- 1e-4
fitTolerance <- 1e-9
# prior share of the posterior precision splitting the two leaf strata; below
# it a leaf is data-dominated and says almost nothing about the leaf prior
priorInformed <- 0.25
thinStratum <- 1000L # a stratum below this reports rather than gates

anyFailure <- FALSE

reportValue <- function(label, value, bound, failed, format = "%12.4f") {
  cat(sprintf(
    paste0("%-30s ", format, " %12s%s\n"),
    label,
    value,
    bound,
    if (failed) "  <- FAIL" else ""
  ))
  failed
}

reportZ <- function(label, z) {
  reportValue(label, z, sprintf("|z| < %g", zBound), abs(z) > zBound, "%+12.2f")
}

reportGap <- function(label, gap, tolerance) {
  reportValue(
    label,
    gap,
    sprintf("< %g", tolerance),
    !(gap <= tolerance),
    "%12.3g"
  )
}

## Standardized leaf draws for one seeded configuration, plus the two
## deterministic preconditions. Returns the pivots with the prior's share of
## each draw's posterior precision and its tree's position in the sweep.
backfitArm <- function(
  n,
  p,
  m,
  nSweep,
  nBurn,
  dataSeed,
  engineSeed,
  useWeights,
  useOffset
) {
  set.seed(dataSeed)
  x <- matrix(runif(n * p), n, p)
  # rare strong indicator: the splits it earns make leaves small enough for the
  # leaf prior to carry real weight, which the ensemble's typical leaf does not
  x[, 5L] <- as.double(runif(n) < 0.04)
  signal <- 3 * x[, 1L] * x[, 2L] - 2 * cos(pi * x[, 3L]) + x[, 4L]
  y <- signal + 4 * x[, 5L] + rnorm(n, sd = 0.8)
  weights <- if (useWeights) rgamma(n, 4, 4) else rep(1, n)
  offset <- if (useOffset) runif(n, -1, 1) else NULL

  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = m,
    keepTrees = TRUE,
    n.samples = nSweep,
    n.burn = nBurn,
    n.thin = 1L,
    updateState = FALSE,
    seed = engineSeed
  )
  sampler <- dbarts(
    x,
    y,
    weights = weights,
    offset = offset,
    control = control,
    node.prior = normal(k = kLeaf)
  )
  result <- sampler$run(nBurn, nSweep)
  trees <- sampler$getTrees()

  # the internal response scale, re-derived rather than read
  modelled <- if (useOffset) y - offset else y
  yRange <- max(modelled) - min(modelled)
  shift <- min(modelled) + 0.5 * yRange
  z <- (modelled - min(modelled)) / yRange - 0.5
  priorPrecision <- (kLeaf * sqrt(m) / nodeScale)^2

  # the calibration in force must be the one the reference assumes: fixed k,
  # this response transform, and a forest prior sd of node.scale at k = 1
  calibration <- sampler$getCalibration()
  reported <- function(field) unname(calibration[1L, field])
  stopifnot(
    reported("k.has.hyperprior") == 0,
    reported("k") == kLeaf,
    isTRUE(all.equal(reported("response.scale"), yRange)),
    isTRUE(all.equal(reported("response.shift"), shift)),
    isTRUE(all.equal(reported("prior.scale"), nodeScale * yRange))
  )

  # pre-order node rows, one block per (sweep, tree)
  blockStart <- which(c(
    TRUE,
    diff(trees$sample) != 0L | diff(trees$tree) != 0L
  ))
  blockEnd <- c(blockStart[-1L] - 1L, nrow(trees))
  stopifnot(length(blockStart) == m * nSweep)
  vars <- trees$var
  values <- trees$value

  # descend one tree's block, returning each row's leaf position within it
  routeTree <- function(first, last) {
    splitVar <- vars[first:last]
    splitValue <- values[first:last]
    leafAt <- integer(n)
    descend <- function(pos, rows) {
      if (splitVar[pos] == -1L) {
        leafAt[rows] <<- pos
        return(pos + 1L)
      }
      left <- x[rows, splitVar[pos]] <= splitValue[pos]
      nextPos <- descend(pos + 1L, rows[left])
      descend(nextPos, rows[!left])
    }
    stopifnot(descend(1L, seq_len(n)) == last - first + 2L)
    leafAt
  }

  pivot <- vector("list", nSweep)
  share <- vector("list", nSweep)
  position <- vector("list", nSweep)
  previousFits <- NULL
  countGap <- 0
  fitGap <- 0
  for (sweep in seq_len(nSweep)) {
    currentFits <- matrix(0, n, m)
    leafRow <- matrix(0L, n, m)
    for (j in seq_len(m)) {
      block <- (sweep - 1L) * m + j
      first <- blockStart[block]
      last <- blockEnd[block]
      at <- routeTree(first, last)
      leafRow[, j] <- at + first - 1L
      currentFits[, j] <- values[first:last][at]
      isLeaf <- vars[first:last] == -1L
      routed <- tabulate(at, last - first + 1L)
      countGap <- max(
        countGap,
        max(abs(routed[isLeaf] - trees$n[first:last][isLeaf]))
      )
    }
    fitted <- yRange * rowSums(currentFits) + shift
    if (useOffset) {
      fitted <- fitted + offset
    }
    fitGap <- max(fitGap, max(abs(fitted - result$train[, sweep])))

    if (!is.null(previousFits)) {
      # sigma is drawn after the trees, so the sweep's draws read the sigma
      # recorded at the end of the previous one
      residualVariance <- (result$sigma[sweep - 1L] / yRange)^2
      totalPrevious <- rowSums(previousFits)
      drawnCurrent <- numeric(n) # trees already redrawn this sweep
      drawnPrevious <- numeric(n) # the same trees as they stood before
      armPivot <- vector("list", m)
      armShare <- vector("list", m)
      for (j in seq_len(m)) {
        # the response minus every OTHER tree: those already redrawn this
        # sweep at their new fits, the rest at the fits they held last sweep
        residual <- z -
          drawnCurrent -
          totalPrevious +
          drawnPrevious +
          previousFits[, j]
        at <- leafRow[, j]
        weighted <- rowsum(weights * residual, at)
        node <- as.integer(rownames(weighted))
        sumWeights <- as.vector(rowsum(weights, at))
        posteriorPrecision <- priorPrecision + sumWeights / residualVariance
        posteriorMean <- (as.vector(weighted) / residualVariance) /
          posteriorPrecision
        # (value - mean) / sd, the standard normal the draw consumed
        armPivot[[j]] <- (values[node] - posteriorMean) *
          sqrt(posteriorPrecision)
        armShare[[j]] <- priorPrecision / posteriorPrecision
        drawnCurrent <- drawnCurrent + currentFits[, j]
        drawnPrevious <- drawnPrevious + previousFits[, j]
      }
      pivot[[sweep]] <- unlist(armPivot)
      share[[sweep]] <- unlist(armShare)
      position[[sweep]] <- rep(seq_len(m), lengths(armPivot))
    }
    previousFits <- currentFits
  }

  list(
    pivot = unlist(pivot),
    share = unlist(share),
    position = unlist(position),
    countGap = countGap,
    fitGap = fitGap,
    sigma = mean(result$sigma)
  )
}

## Exact moment and KS tests on a set of iid standard normals.
checkMoments <- function(label, t) {
  n <- length(t)
  if (n < thinStratum) {
    cat(sprintf("%-30s %12d %12s\n", paste(label, "draws"), n, "(thin)"))
    return(FALSE)
  }
  failed <- reportZ(paste(label, "mean"), mean(t) * sqrt(n))
  reportZ(paste(label, "variance"), (var(t) - 1) * sqrt((n - 1) / 2)) || failed
}

checkArm <- function(name, arm, m) {
  cat(sprintf(
    "\n%s: %d leaf draws, %d trees, mean sigma %.3f\n",
    name,
    length(arm$pivot),
    m,
    arm$sigma
  ))
  failed <- reportGap("routed vs recorded counts", arm$countGap, 0)
  failed <- reportGap("tree fits vs training fit", arm$fitGap, fitTolerance) ||
    failed

  t <- arm$pivot
  failed <- checkMoments("pooled", t) || failed
  ks <- stats::ks.test(t, "pnorm")$p.value
  failed <- reportValue(
    "pooled KS vs N(0, 1)",
    ks,
    sprintf("p > %g", ksBound),
    ks <= ksBound,
    "%12.3g"
  ) ||
    failed

  informed <- arm$share >= priorInformed
  failed <- checkMoments("prior-informed leaf", t[informed]) || failed
  failed <- checkMoments("data-dominated leaf", t[!informed]) || failed

  # an order error in the running residual would land on one end of the sweep
  early <- arm$position <= m %/% 3L
  late <- arm$position > m - m %/% 3L
  failed <- reportZ(
    "first-third tree mean",
    mean(t[early]) * sqrt(sum(early))
  ) ||
    failed
  reportZ("last-third tree mean", mean(t[late]) * sqrt(sum(late))) || failed
}

cat(sprintf(
  "%-30s %12s %12s\n",
  "check",
  "statistic",
  "bound"
))

armA <- backfitArm(
  n = 1000L,
  p = 10L,
  m = 200L,
  nSweep = if (quick) 100L else 400L,
  nBurn = 200L,
  dataSeed = 23L,
  engineSeed = 314L,
  useWeights = FALSE,
  useOffset = FALSE
)
anyFailure <- checkArm("ensemble scale", armA, 200L) || anyFailure

armB <- backfitArm(
  n = 500L,
  p = 6L,
  m = 50L,
  nSweep = if (quick) 100L else 300L,
  nBurn = 200L,
  dataSeed = 91L,
  engineSeed = 271L,
  useWeights = TRUE,
  useOffset = TRUE
)
anyFailure <- checkArm("case weights and offset", armB, 50L) || anyFailure

if (anyFailure) {
  cat("\nFAIL: leaf draws deviate from the conditional closed form\n")
  quit(status = 1L)
}
cat("\nOK: leaf draws match the conditional closed form\n")
