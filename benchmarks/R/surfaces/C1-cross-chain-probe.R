#!/usr/bin/env Rscript

# C1, the same-temperature one-tree exchange between chains: a GENERATOR-ONLY
# pre-check of the acceptance rate that move would have. Nothing here is a
# kernel. The sampler runs its shipped moves untouched, the exchange is scored
# offline from the states the sampler leaves behind, and no draw is altered or
# proposed; the script only reads tree structure, fitted values and sigma back
# out of a fitted sampler.
#
# The move priced here swaps ONE tree's structure between two chains of one
# sampler that target the same posterior. Leaf values are integrated out and
# redrawn after, and the two tree priors cancel: the prior reads rule indices
# and grid shape and never a predictor value, and one sampler shares its cut
# grid and column codes across every chain, so a donor structure installs with
# the identical row partition. The Metropolis acceptance on the product target
# is then symmetric-proposal ordinary:
#
#   alpha = L_A(T_{B,k}) L_B(T_{A,j}) / (L_A(T_{A,j}) L_B(T_{B,k}))
#
# the four collapsed marginals of each tree's structure against each chain's
# OWN partial residual - that chain's residual with its own tree removed - and
# its own sigma.
#
# The collapsed marginal. For a tree inducing leaves l with member counts n_l
# and partial-residual sums s_l, under leaf prior mu ~ N(0, tau^2) and
# residual variance sigma^2 the standard conjugate integral gives
#
#   log L = sum_l [ -n_l/2 log(2 pi sigma^2) - SS_l / (2 sigma^2)
#                   + 1/2 log(sigma^2 / (sigma^2 + n_l tau^2))
#                   + 1/2 s_l^2 tau^2 / (sigma^2 (sigma^2 + n_l tau^2)) ]
#
# The first two terms sum over leaves to a quantity that depends on the
# residual alone and not on the partition, so they cancel inside each chain's
# own ratio; the last two are what this script computes. `probeLogMarginalFull`
# keeps the constants and is checked against a brute-force numerical integral,
# `probeLogMarginal` drops them.
#
# Scale. Leaf values come off getTrees on the engine's internal response scale,
# where a gaussian response is mapped to [-0.5, 0.5] by (y - min)/range - 0.5,
# so residuals and sigma are put on that scale too and tau is the shipped node
# prior there: node.scale / (k sqrt(m)) = 0.5 / (2 sqrt(75)). The reduced log
# marginal is invariant to a common rescaling of (residual, sigma, tau), so the
# choice of scale cannot move alpha; working internally just avoids carrying
# the response range through every term.
#
# Configuration: the He and Hahn cell at 10.4 of docs/design/benchmark-surfaces
# - independent design, n = 10000, p = 30, Trig+poly, kappa = 1, 75 trees -
# with eight chains in one sampler, 500 burn-in sweeps, then 100 states spaced
# five sweeps apart. At each state 50 exchanges are drawn (chains A != B and
# tree indices j, k uniform) and scored. A second seed and the low-noise
# Friedman census cell (n = 5000, p = 10, 75 trees, sigma^2 = 0.1) run for
# contrast.
#
# Readouts: the distribution of log alpha and of min(1, alpha), the implied
# acceptance rate overall and with identical structures excluded, the same by
# leaf-count pair, and the share of a chain's tree structures that also appear
# in some other chain.
#
# Usage: Rscript C1-cross-chain-probe.R [outputDir] [quick] [cell ...]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nTrees <- 75L
nChains <- 8L
spacing <- 5L
nBurn <- if (quick) 50L else 500L
nStates <- if (quick) 5L else 100L
nExchanges <- if (quick) 10L else 50L
nMatched <- if (quick) 5L else 25L

# The shipped constant-gaussian leaf prior: mu ~ N(0, (node.scale / (k sqrt m))^2)
# on the internal response scale, at the shipped k = 2 and node.scale = 0.5.
leafK <- 2
nodeScale <- 0.5
tau <- nodeScale / (leafK * sqrt(nTrees))

# Cells. `c1` is 10.4's independent-design He-Hahn arm at its own data seeds;
# `lownoise` is the move census's low-noise cell (n = 5000, p = 10,
# sigma^2 = 0.1) drawn through this battery's own Friedman generator, which
# differs from the census's only in the order it consumes the RNG stream.
probeCells <- list(
  c1seed1 = list(cell = "c1", replicate = 1L),
  c1seed2 = list(cell = "c1", replicate = 2L),
  lownoise = list(cell = "lownoise", replicate = 1L)
)

probeData <- function(cell, replicate) {
  if (cell == "c1") {
    set.seed(surfacesDataSeed("C1", "independent", replicate))
    return(surfacesHeHahn(
      n = if (quick) 2000L else 10000L,
      nTest = 1000L,
      p = 30L,
      which = "trigpoly",
      kappa = 1,
      design = "independent"
    ))
  }
  if (cell == "lownoise") {
    set.seed(surfacesDataSeed("P1", "lownoise", replicate))
    return(surfacesFriedman(
      n = if (quick) 2000L else 5000L,
      nTest = 1000L,
      p = 10L,
      sigma = sqrt(0.1)
    ))
  }
  stop("unknown cell '", cell, "'")
}

outputDir <- surfacesOutputDir(args, flags = c("quick", names(probeCells)))
requested <- intersect(args, names(probeCells))
if (length(requested) == 0L) {
  requested <- names(probeCells)
}

# ------------------------------------------------------- tree readout helpers

# Route every row of `x` through one tree and return its leaf assignment, the
# leaf member counts and the leaf values, all in the tree's own preorder. Trees
# arrive from getTrees depth-first, root first, left subtree before right, with
# var = -1 marking a leaf; splits send x <= cut left and x > cut right.
probeAssignLeaves <- function(varCol, valueCol, x) {
  nObs <- nrow(x)
  leaf <- integer(nObs)
  value <- numeric(0)
  count <- integer(0)
  cursor <- 1L
  nLeaves <- 0L
  descend <- function(rows) {
    node <- cursor
    cursor <<- cursor + 1L
    variable <- varCol[node]
    if (variable < 0L) {
      nLeaves <<- nLeaves + 1L
      leaf[rows] <<- nLeaves
      value[nLeaves] <<- valueCol[node]
      count[nLeaves] <<- length(rows)
      return(invisible(NULL))
    }
    goesLeft <- x[rows, variable] <= valueCol[node]
    descend(rows[goesLeft])
    descend(rows[!goesLeft])
  }
  descend(seq_len(nObs))
  list(leaf = leaf, value = value, count = count)
}

# Sum `values` within each leaf, in leaf order and with absent leaves at zero.
probeLeafSums <- function(values, leaf, nLeaves) {
  out <- numeric(nLeaves)
  grouped <- rowsum(values, leaf, reorder = TRUE)
  out[as.integer(rownames(grouped))] <- grouped[, 1L]
  out
}

# The partition-dependent part of the collapsed marginal: the two terms of the
# header's expression that a repartition can move.
probeLogMarginal <- function(fit, residual, sigmaSq, tauSq) {
  sums <- probeLeafSums(residual, fit$leaf, length(fit$count))
  denominator <- sigmaSq + fit$count * tauSq
  0.5 *
    sum(
      log(sigmaSq / denominator) +
        sums^2 * tauSq / (sigmaSq * denominator)
    )
}

# One leaf's collapsed marginal with every constant kept, for the check below.
probeLogMarginalFull <- function(residual, sigmaSq, tauSq) {
  nObs <- length(residual)
  denominator <- sigmaSq + nObs * tauSq
  -0.5 *
    nObs *
    log(2 * pi * sigmaSq) -
    0.5 * sum(residual^2) / sigmaSq +
    0.5 * log(sigmaSq / denominator) +
    0.5 * sum(residual)^2 * tauSq / (sigmaSq * denominator)
}

# The same integral done numerically, on the exact posterior support: the
# integrand is written relative to its own peak so the quadrature never sees
# an underflowed density.
probeLogMarginalNumeric <- function(residual, sigmaSq, tauSq) {
  nObs <- length(residual)
  posteriorVar <- 1 / (nObs / sigmaSq + 1 / tauSq)
  posteriorMean <- (sum(residual) / sigmaSq) * posteriorVar
  posteriorSd <- sqrt(posteriorVar)
  logKernel <- function(mu) {
    vapply(
      mu,
      function(m) {
        sum(dnorm(residual, m, sqrt(sigmaSq), log = TRUE)) +
          dnorm(m, 0, sqrt(tauSq), log = TRUE)
      },
      numeric(1L)
    )
  }
  peak <- logKernel(posteriorMean)
  area <- integrate(
    function(mu) exp(logKernel(mu) - peak),
    lower = posteriorMean - 12 * posteriorSd,
    upper = posteriorMean + 12 * posteriorSd,
    rel.tol = 1e-10
  )
  peak + log(area$value)
}

# A structure key that ignores leaf VALUES: interior nodes contribute their
# split variable and cut, leaves contribute a placeholder. Two trees sharing a
# key induce the identical partition, so their exchange has alpha exactly 1.
probeStructureKeys <- function(trees, group) {
  token <- ifelse(
    trees$var < 0L,
    "L",
    paste0(trees$var, "@", sprintf("%.17g", trees$value))
  )
  vapply(split(token, group), paste, character(1L), collapse = "|")
}

# ------------------------------------------------------------------ the probe

probeRunCell <- function(name, spec) {
  cat(sprintf(
    "\n=== %s (%s, replicate %d) ===\n",
    name,
    spec$cell,
    spec$replicate
  ))
  data <- probeData(spec$cell, spec$replicate)
  x <- data$x
  y <- data$y
  nObs <- length(y)

  # the internal response scale the engine reports leaf values on
  yMin <- min(y)
  yRange <- diff(range(y))
  yInternal <- (y - yMin) / yRange - 0.5
  fitShift <- yRange * 0.5 + yMin

  control <- dbartsControl(
    verbose = FALSE,
    n.trees = nTrees,
    n.chains = nChains,
    n.threads = 1L,
    n.samples = spacing,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = surfacesSamplerSeed(spec$replicate)
  )
  sampler <- dbarts(x, y, control = control)
  invisible(sampler$run(nBurn, spacing))

  set.seed(1000L + spec$replicate)
  records <- vector("list", nStates)
  shared <- vector("list", nStates)
  checkRows <- NULL
  reconstructionError <- NA_real_

  for (state in seq_len(nStates)) {
    samples <- sampler$run(0L, spacing)
    trees <- sampler$getTrees(sampleNums = spacing)
    fitInternal <- (samples$train[, spacing, ] - fitShift) / yRange
    sigmaInternal <- samples$sigma[spacing, ] / yRange

    group <- (trees$chain - 1L) * nTrees + trees$tree
    nodeRows <- split(seq_len(nrow(trees)), group)
    keys <- probeStructureKeys(trees, group)
    leafCounts <- as.integer(tapply(trees$var < 0L, group, sum))

    cache <- new.env(parent = emptyenv())
    treeFit <- function(chain, tree) {
      slot <- (chain - 1L) * nTrees + tree
      id <- as.character(slot)
      if (exists(id, envir = cache, inherits = FALSE)) {
        return(get(id, envir = cache, inherits = FALSE))
      }
      rows <- nodeRows[[id]]
      fit <- probeAssignLeaves(trees$var[rows], trees$value[rows], x)
      assign(id, fit, envir = cache)
      fit
    }

    # once, at the first state: the leaf-count and fit reconstruction checks
    if (state == 1L) {
      total <- numeric(nObs)
      for (tree in seq_len(nTrees)) {
        fit <- treeFit(1L, tree)
        stopifnot(identical(
          as.integer(trees$n[nodeRows[[as.character(tree)]]][
            trees$var[nodeRows[[as.character(tree)]]] < 0L
          ]),
          fit$count
        ))
        total <- total + fit$value[fit$leaf]
      }
      reconstructionError <- max(abs(total - fitInternal[, 1L]))
      checkRows <- probeCheckMarginal(
        treeFit,
        trees,
        nodeRows,
        yInternal,
        fitInternal,
        sigmaInternal
      )
    }

    # Two draw designs. `uniform` is the move as proposed - chains A != B and
    # tree indices j, k all uniform. `matchedIndex` holds k = j, the cheaper
    # variant in which a chain trades a tree for its own counterpart in
    # another chain; the two are recorded side by side because they are
    # different moves and the second is the one an implementation would reach
    # for first.
    drawExchange <- function(count, matched) {
      chainA <- sample.int(nChains, count, replace = TRUE)
      offset <- sample.int(nChains - 1L, count, replace = TRUE)
      chainB <- 1L + ((chainA - 1L + offset) %% nChains)
      treeJ <- sample.int(nTrees, count, replace = TRUE)
      treeK <- if (matched) treeJ else sample.int(nTrees, count, replace = TRUE)
      logAlpha <- numeric(count)
      for (e in seq_len(count)) {
        fitJ <- treeFit(chainA[e], treeJ[e])
        fitK <- treeFit(chainB[e], treeK[e])
        residualA <- yInternal -
          fitInternal[, chainA[e]] +
          fitJ$value[fitJ$leaf]
        residualB <- yInternal -
          fitInternal[, chainB[e]] +
          fitK$value[fitK$leaf]
        sigmaSqA <- sigmaInternal[chainA[e]]^2
        sigmaSqB <- sigmaInternal[chainB[e]]^2
        logAlpha[e] <- probeLogMarginal(fitK, residualA, sigmaSqA, tau^2) -
          probeLogMarginal(fitJ, residualA, sigmaSqA, tau^2) +
          probeLogMarginal(fitJ, residualB, sigmaSqB, tau^2) -
          probeLogMarginal(fitK, residualB, sigmaSqB, tau^2)
      }
      slotA <- (chainA - 1L) * nTrees + treeJ
      slotB <- (chainB - 1L) * nTrees + treeK
      data.frame(
        state = state,
        design = if (matched) "matchedIndex" else "uniform",
        chainA = chainA,
        chainB = chainB,
        treeJ = treeJ,
        treeK = treeK,
        leavesJ = leafCounts[slotA],
        leavesK = leafCounts[slotB],
        logAlpha = logAlpha,
        sameStructure = keys[as.character(slotA)] == keys[as.character(slotB)],
        stringsAsFactors = FALSE
      )
    }

    records[[state]] <- rbind(
      drawExchange(nExchanges, FALSE),
      drawExchange(nMatched, TRUE)
    )

    keyChain <- rep(seq_len(nChains), each = nTrees)
    elsewhere <- vapply(
      seq_len(nChains),
      function(chain) {
        mine <- keys[keyChain == chain]
        others <- keys[keyChain != chain]
        mean(mine %in% others)
      },
      numeric(1L)
    )
    stumps <- leafCounts == 1L
    elsewhereNonStump <- vapply(
      seq_len(nChains),
      function(chain) {
        mine <- keys[keyChain == chain & !stumps]
        others <- keys[keyChain != chain]
        if (length(mine) == 0L) NA_real_ else mean(mine %in% others)
      },
      numeric(1L)
    )
    shared[[state]] <- data.frame(
      state = state,
      chain = seq_len(nChains),
      sharedAll = elsewhere,
      sharedNonStump = elsewhereNonStump,
      stumpShare = vapply(
        seq_len(nChains),
        function(chain) mean(stumps[keyChain == chain]),
        numeric(1L)
      ),
      meanLeaves = vapply(
        seq_len(nChains),
        function(chain) mean(leafCounts[keyChain == chain]),
        numeric(1L)
      ),
      stringsAsFactors = FALSE
    )
  }

  list(
    cell = spec$cell,
    replicate = spec$replicate,
    exchanges = do.call(rbind, records),
    shared = do.call(rbind, shared),
    marginalCheck = checkRows,
    reconstructionError = reconstructionError,
    tau = tau,
    sigmaInternalLast = sigmaInternal
  )
}

# The closed form against the numerical integral, on the smallest leaves the
# first state offers so the quadrature is on a leaf a person can check.
probeCheckMarginal <- function(
  treeFit,
  trees,
  nodeRows,
  yInternal,
  fitInternal,
  sigmaInternal
) {
  rows <- NULL
  sigmaSq <- sigmaInternal[1L]^2
  for (tree in seq_len(min(20L, nTrees))) {
    fit <- treeFit(1L, tree)
    if (length(fit$count) < 2L) {
      next
    }
    residual <- yInternal - fitInternal[, 1L] + fit$value[fit$leaf]
    smallest <- which.min(fit$count)
    leafResidual <- residual[fit$leaf == smallest]
    rows <- rbind(
      rows,
      data.frame(
        tree = tree,
        leafSize = length(leafResidual),
        closedForm = probeLogMarginalFull(leafResidual, sigmaSq, tau^2),
        numeric = probeLogMarginalNumeric(leafResidual, sigmaSq, tau^2),
        stringsAsFactors = FALSE
      )
    )
    if (!is.null(rows) && nrow(rows) >= 5L) {
      break
    }
  }
  rows$difference <- rows$closedForm - rows$numeric
  rows
}

# --------------------------------------------------------------- the readouts

probeAcceptance <- function(logAlpha) {
  mean(exp(pmin(0, logAlpha)))
}

probeSummary <- function(exchanges) {
  quantiles <- quantile(
    exchanges$logAlpha,
    c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
    names = FALSE
  )
  data.frame(
    n = nrow(exchanges),
    acceptance = probeAcceptance(exchanges$logAlpha),
    pGeOne = mean(exchanges$logAlpha >= 0),
    pGePoint1 = mean(exchanges$logAlpha >= log(0.1)),
    min = quantiles[1L],
    q01 = quantiles[2L],
    q05 = quantiles[3L],
    q25 = quantiles[4L],
    median = quantiles[5L],
    q75 = quantiles[6L],
    q95 = quantiles[7L],
    q99 = quantiles[8L],
    max = quantiles[9L],
    stringsAsFactors = FALSE
  )
}

# Under the product target the exchange is a symmetric-proposal involution, so
# E[alpha] = 1 exactly and Markov bounds the positive tail by
# P(log alpha > t) <= exp(-t). Measuring that tail therefore says whether the
# eight chains are draws from a common posterior at all: an observed tail far
# above the bound is the chains sitting in different places, read off the
# likelihood rather than off a between-chain variance.
probeTailBound <- function(logAlpha) {
  thresholds <- c(0, 1, 2, 5, 10, 20, 50)
  data.frame(
    t = thresholds,
    observed = vapply(
      thresholds,
      function(u) mean(logAlpha > u),
      numeric(1L)
    ),
    stationaryBound = exp(-thresholds),
    stringsAsFactors = FALSE
  )
}

# Acceptance against sweep number: if the tail above is a burn-in artifact it
# falls away as the chains run on, and if it is the standing between-chain
# spread it does not.
probeByStateBlock <- function(exchanges, nBlocks = 4L) {
  edges <- seq(0, max(exchanges$state), length.out = nBlocks + 1L)
  block <- cut(exchanges$state, breaks = edges, include.lowest = TRUE)
  parts <- split(exchanges, block)
  rows <- lapply(names(parts), function(label) {
    part <- parts[[label]]
    data.frame(
      states = label,
      n = nrow(part),
      acceptance = probeAcceptance(part$logAlpha),
      pGeOne = mean(part$logAlpha >= 0),
      medianLogAlpha = median(part$logAlpha),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

probeLeafBin <- function(leaves) {
  ifelse(leaves >= 5L, "5+", as.character(leaves))
}

probeByLeafPair <- function(exchanges) {
  low <- probeLeafBin(pmin(exchanges$leavesJ, exchanges$leavesK))
  high <- probeLeafBin(pmax(exchanges$leavesJ, exchanges$leavesK))
  parts <- split(exchanges, paste0(low, "-", high))
  rows <- lapply(names(parts), function(label) {
    part <- parts[[label]]
    data.frame(
      pair = label,
      n = nrow(part),
      sameStructureShare = mean(part$sameStructure),
      acceptance = probeAcceptance(part$logAlpha),
      acceptanceDistinct = if (any(!part$sameStructure)) {
        probeAcceptance(part$logAlpha[!part$sameStructure])
      } else {
        NA_real_
      },
      medianLogAlpha = median(part$logAlpha),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

probeReport <- function(result) {
  cat(sprintf(
    "\n-- %s replicate %d --\n",
    result$cell,
    result$replicate
  ))
  cat(sprintf(
    "fit reconstruction max abs error %.3g; tau %.6f\n",
    result$reconstructionError,
    result$tau
  ))
  cat("\nclosed form against the numerical integral, one leaf each:\n")
  print(result$marginalCheck, row.names = FALSE, digits = 10)
  for (design in c("uniform", "matchedIndex")) {
    exchanges <- result$exchanges[
      result$exchanges$design == design,
      ,
      drop = FALSE
    ]
    distinct <- exchanges[!exchanges$sameStructure, , drop = FALSE]
    cat(sprintf("\n[%s] log alpha, all exchanges:\n", design))
    print(probeSummary(exchanges), row.names = FALSE, digits = 4)
    cat(sprintf(
      "[%s] identical structures in both slots: %.4f of exchanges\n",
      design,
      mean(exchanges$sameStructure)
    ))
    cat(sprintf("[%s] log alpha, identical structures excluded:\n", design))
    print(probeSummary(distinct), row.names = FALSE, digits = 4)
    cat(sprintf(
      "[%s] by leaf-count pair (unordered, capped at 5+):\n",
      design
    ))
    print(probeByLeafPair(exchanges), row.names = FALSE, digits = 4)
    cat(sprintf("[%s] positive tail against the stationary bound:\n", design))
    print(probeTailBound(exchanges$logAlpha), row.names = FALSE, digits = 4)
    cat(sprintf("[%s] acceptance by sweep block:\n", design))
    print(probeByStateBlock(exchanges), row.names = FALSE, digits = 4)
  }
  cat("\nstructure sharing across chains, per state and chain:\n")
  print(
    data.frame(
      statistic = c(
        "shared with another chain",
        "shared, non-stump trees only",
        "stump share",
        "mean leaves per tree"
      ),
      mean = c(
        mean(result$shared$sharedAll),
        mean(result$shared$sharedNonStump, na.rm = TRUE),
        mean(result$shared$stumpShare),
        mean(result$shared$meanLeaves)
      ),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE,
    digits = 4
  )
  invisible(NULL)
}

results <- list()
for (name in requested) {
  results[[name]] <- probeRunCell(name, probeCells[[name]])
  probeReport(results[[name]])
}

surfacesSave(
  results,
  outputDir,
  if (quick) "C1-cross-chain-probe-quick" else "C1-cross-chain-probe"
)
