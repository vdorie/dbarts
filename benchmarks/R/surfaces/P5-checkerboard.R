#!/usr/bin/env Rscript

# P5, the checkerboard on an autocorrelated design: Zhu, Zeng and Kosorok
# (2015) scenario 3. A pure two-way interaction surface,
# f = 2 x5 x10 + 2 x15 x20, on forty predictors whose covariance is
# 0.9^|j-k|. Inclusion has an exactly right answer, {x5, x10, x15, x20},
# and every true column sits between two decoys correlated 0.9 with it.
#
# There is no published dbarts number for this cell. Its oracle is the
# inclusion truth, and its primary statistic is structural: the
# between-chain standard deviation of the time-averaged inclusion
# proportion on the four true columns and on their immediate neighbours.
# A structural readout is not label invariant, so it is read BETWEEN chains
# and never within one.
#
# The between-chain spread is reported against a mixing null computed from
# the same draws: if the chains were exploring one mode, the spread of their
# time averages would be the Monte Carlo standard error of a chain mean,
# sd(draws) / sqrt(ESS). A ratio near one says the chains agree to within
# their own resolution; a ratio far above one says they are reporting
# different structures for the same fit.
#
# Secondary: pointwise 95% coverage and RMSE of the true f on a held-out
# thousand rows, the summed inclusion share on the true columns, and whether
# any decoy outranks the weakest true column.
#
# Three move-set arms, on the same matched seeds and differing only in
# `proposal.probs`, so that every contrast is paired:
#
#   default     proposal.probs unset, which is the shipped mixture
#               (birth_death 0.6, swap 0, change 0.4, perturb 0)
#   birthdeath  birth_death 1, swap 0, change 0
#   swap        birth_death 0.5, swap 0.1, change 0.4, the former default
#
# This is the cell Tan et al.'s Theorem 5.2 is about: the theorem bounds the
# hitting time for a pure interaction when the change move is disallowed,
# which is what `birthdeath` is and what neither of the other two is.
# `default` is re-run beside the other two rather than read off an earlier
# session, so the contrast is paired within one run; on matched seeds it
# reproduces the single-arm run exactly.
#
# Usage: Rscript P5-checkerboard.R [outputDir] [quick] [arm ...]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nReplicates <- if (quick) 2L else 20L
n <- 1600L
p <- 40L
nTest <- 1000L
nChains <- 8L
nBurn <- 1000L
nSamples <- if (quick) 400L else 2000L

# NULL leaves proposal.probs unset, which is the shipped mixture.
arms <- list(
  default = NULL,
  birthdeath = c(
    birth_death = 1,
    swap = 0,
    change = 0,
    perturb = 0,
    birth = 0.5
  ),
  swap = c(
    birth_death = 0.5,
    swap = 0.1,
    change = 0.4,
    perturb = 0,
    birth = 0.5
  )
)
armNames <- names(arms)
selectedArms <- intersect(armNames, args)
if (length(selectedArms) > 0L) {
  arms <- arms[selectedArms]
}
# The control every other arm is read against.
armControl <- "default"

outputDir <- surfacesOutputDir(args, flags = c("quick", armNames))

signal <- surfacesCheckerboardSignal
neighbours <- surfacesCheckerboardNeighbours(p)
decoys <- setdiff(seq_len(p), signal)

surfacesUptime("uptime before")

rows <- list()
perColumn <- list()
for (armName in names(arms)) {
  probs <- arms[[armName]]
  for (replicate in seq_len(nReplicates)) {
    set.seed(surfacesDataSeed("P5", "checkerboard", replicate))
    data <- surfacesCheckerboard(n = n, p = p, nTest = nTest)
    call <- list(
      data$x,
      data$y,
      test = data$xTest,
      n.chains = nChains,
      n.burn = nBurn,
      n.samples = nSamples,
      n.thin = 1L,
      n.threads = 1L,
      verbose = FALSE,
      seed = surfacesSamplerSeed(replicate)
    )
    if (!is.null(probs)) {
      call$proposal.probs <- probs
    }
    startedAt <- proc.time()
    fit <- do.call(bart2, call)
    elapsed <- (proc.time() - startedAt)[["elapsed"]]

    inclusion <- surfacesInclusion(fit$varcount)
    chainMeans <- matrix(0, nChains, p)
    chainMcse <- matrix(0, nChains, p)
    for (chain in seq_len(nChains)) {
      draws <- inclusion[surfacesChainRows(chain, nSamples), , drop = FALSE]
      chainMeans[chain, ] <- colMeans(draws)
      chainMcse[chain, ] <- vapply(
        seq_len(p),
        function(j) sd(draws[, j]) / sqrt(surfacesEss(draws[, j])),
        numeric(1L)
      )
    }
    betweenSd <- apply(chainMeans, 2L, sd)
    mixingNull <- sqrt(colMeans(chainMcse^2))

    testDraws <- extract(fit, type = "ev", sample = "test")
    pooledInclusion <- colMeans(chainMeans)
    weakestSignal <- min(pooledInclusion[signal])

    perColumn[[length(perColumn) + 1L]] <- data.frame(
      arm = armName,
      replicate = replicate,
      column = seq_len(p),
      inclusion = pooledInclusion,
      betweenChainSd = betweenSd,
      mixingNull = mixingNull,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      arm = armName,
      replicate = replicate,
      signalBetweenSd = mean(betweenSd[signal]),
      signalMixingNull = mean(mixingNull[signal]),
      signalRatio = mean(betweenSd[signal]) / mean(mixingNull[signal]),
      neighbourBetweenSd = mean(betweenSd[neighbours]),
      neighbourRatio = mean(betweenSd[neighbours]) /
        mean(mixingNull[neighbours]),
      inclusionSignal = sum(pooledInclusion[signal]),
      inclusionNeighbour = sum(pooledInclusion[neighbours]),
      maxDecoy = max(pooledInclusion[decoys]),
      decoysAboveWeakest = sum(pooledInclusion[decoys] > weakestSignal),
      coverage = surfacesCoverage(testDraws, data$fTest),
      rmse = surfacesRmse(testDraws, data$fTest),
      wall = elapsed,
      stringsAsFactors = FALSE
    )
    cat(sprintf(
      "%-11s rep %2d  incl(true) %.3f  sd/null %.1f  decoys above weakest %d  cover %.3f  %.0fs\n",
      armName,
      replicate,
      rows[[length(rows)]]$inclusionSignal,
      rows[[length(rows)]]$signalRatio,
      rows[[length(rows)]]$decoysAboveWeakest,
      rows[[length(rows)]]$coverage,
      elapsed
    ))
    rm(fit, testDraws, inclusion)
    invisible(gc(verbose = FALSE))
  }
}
results <- do.call(rbind, rows)
columns <- do.call(rbind, perColumn)

# One statistic across the arms, mean over seeds with the min-max range.
armSummaryRow <- function(label, column, digits = 3L) {
  cat(sprintf("%-28s", label))
  for (armName in names(arms)) {
    keep <- results$arm == armName
    cat(sprintf(
      " %-22s",
      surfacesRange(results[[column]][keep], digits = digits)
    ))
  }
  cat("\n")
}

surfacesHeader("P5 checkerboard: mean over seeds (min-max)")
cat(sprintf("%-28s", "statistic"))
for (armName in names(arms)) {
  cat(sprintf(" %-22s", armName))
}
cat("\n")
armSummaryRow("between-chain sd, true cols", "signalBetweenSd", digits = 4L)
armSummaryRow("mixing null, true cols", "signalMixingNull", digits = 4L)
armSummaryRow("ratio to null, true cols", "signalRatio", digits = 2L)
armSummaryRow("between-chain sd, decoys", "neighbourBetweenSd", digits = 4L)
armSummaryRow("ratio to null, decoys", "neighbourRatio", digits = 2L)
armSummaryRow("inclusion share, true cols", "inclusionSignal")
armSummaryRow("inclusion share, decoys", "inclusionNeighbour")
armSummaryRow("largest non-true inclusion", "maxDecoy")
armSummaryRow("non-true cols above weakest", "decoysAboveWeakest", digits = 1L)
armSummaryRow("95% coverage of true f", "coverage")
armSummaryRow("held-out RMSE", "rmse")
armSummaryRow("wall s per fit", "wall", digits = 1L)

surfacesHeader("pooled inclusion by column, mean over seeds")
cat(sprintf("%-6s", "column"))
for (armName in names(arms)) {
  cat(sprintf(" %-8s", armName))
}
cat("\n")
for (j in seq_len(p)) {
  tag <- if (j %in% signal) {
    "TRUE"
  } else if (j %in% neighbours) {
    "decoy"
  } else {
    ""
  }
  cat(sprintf("x%-5d", j))
  for (armName in names(arms)) {
    keep <- columns$arm == armName & columns$column == j
    cat(sprintf(" %-8.4f", mean(columns$inclusion[keep])))
  }
  cat(sprintf(" %s\n", tag))
}

# The paired differences against the control arm, seed by seed. The primary
# is the first of them; the three that carry a frozen margin are in the
# second table.
contrastArms <- setdiff(names(arms), armControl)
if (armControl %in% names(arms) && length(contrastArms) > 0L) {
  pairedColumn <- function(contrast, base, column) {
    contrast[[column]] - base[[column]]
  }
  pairedArm <- function(armName) {
    base <- results[results$arm == armControl, ]
    contrast <- results[results$arm == armName, ]
    paired <- match(base$replicate, contrast$replicate)
    keepBase <- !is.na(paired)
    list(base = base[keepBase, ], contrast = contrast[paired[keepBase], ])
  }

  surfacesHeader(sprintf(
    "P5 move sets: paired difference against %s, seed by seed (structure)",
    armControl
  ))
  cat(sprintf(
    "%-11s %-35s %-35s %s\n",
    "arm",
    "d between-chain sd, true cols",
    "d ratio to null, true cols",
    "d between-chain sd, decoys"
  ))
  for (armName in contrastArms) {
    pair <- pairedArm(armName)
    if (nrow(pair$contrast) == 0L) {
      next
    }
    cat(sprintf(
      "%-11s %-35s %-35s %s\n",
      armName,
      surfacesPairedDifference(
        pairedColumn(pair$contrast, pair$base, "signalBetweenSd"),
        digits = 4L
      ),
      surfacesPairedDifference(
        pairedColumn(pair$contrast, pair$base, "signalRatio"),
        digits = 2L
      ),
      surfacesPairedDifference(
        pairedColumn(pair$contrast, pair$base, "neighbourBetweenSd"),
        digits = 4L
      )
    ))
  }

  surfacesHeader(sprintf(
    "P5 move sets: paired difference against %s, seed by seed (gated)",
    armControl
  ))
  cat(sprintf(
    "%-11s %-35s %-35s %s\n",
    "arm",
    "d inclusion share, true cols",
    "d 95% coverage",
    "d held-out RMSE"
  ))
  for (armName in contrastArms) {
    pair <- pairedArm(armName)
    if (nrow(pair$contrast) == 0L) {
      next
    }
    cat(sprintf(
      "%-11s %-35s %-35s %s\n",
      armName,
      surfacesPairedDifference(
        pairedColumn(pair$contrast, pair$base, "inclusionSignal")
      ),
      surfacesPairedDifference(
        pairedColumn(pair$contrast, pair$base, "coverage")
      ),
      surfacesPairedDifference(pairedColumn(pair$contrast, pair$base, "rmse"))
    ))
  }

  surfacesHeader(
    "P5 move sets: margins (coverage -0.010, inclusion -0.010, RMSE ratio 1.02)"
  )
  cat(sprintf(
    "%-11s %-28s %-28s %-12s %s\n",
    "arm",
    "coverage",
    "inclusion share, true cols",
    "rmse ratio",
    "rmse verdict"
  ))
  verdicts <- list()
  for (armName in contrastArms) {
    pair <- pairedArm(armName)
    if (nrow(pair$contrast) == 0L) {
      next
    }
    dCoverage <- pairedColumn(pair$contrast, pair$base, "coverage")
    dInclusion <- pairedColumn(pair$contrast, pair$base, "inclusionSignal")
    dRmse <- pairedColumn(pair$contrast, pair$base, "rmse")
    dBetween <- pairedColumn(pair$contrast, pair$base, "signalBetweenSd")
    rmseRatio <- mean(pair$contrast$rmse) / mean(pair$base$rmse)
    # RMSE's margin is a ratio, so the point estimate carries it and the
    # paired difference underneath carries the separation condition.
    rmseVerdict <- if (rmseRatio <= 1.02) {
      "within margin"
    } else {
      surfacesMarginVerdict(dRmse, 0, 1)
    }
    verdicts[[length(verdicts) + 1L]] <- data.frame(
      arm = armName,
      coverageVerdict = surfacesMarginVerdict(dCoverage, -0.010, -1),
      inclusionVerdict = surfacesMarginVerdict(dInclusion, -0.010, -1),
      rmseRatio = rmseRatio,
      rmseVerdict = rmseVerdict,
      # The primary is a spread between chains, so smaller is better.
      primaryImprovement = surfacesImprovementVerdict(dBetween, -1),
      stringsAsFactors = FALSE
    )
    cat(sprintf(
      "%-11s %-28s %-28s %-12.3f %s\n",
      armName,
      verdicts[[length(verdicts)]]$coverageVerdict,
      verdicts[[length(verdicts)]]$inclusionVerdict,
      rmseRatio,
      rmseVerdict
    ))
  }
  verdicts <- do.call(rbind, verdicts)

  surfacesHeader("P5 move sets: the primary against the 4x-SE improvement bar")
  cat(sprintf("%-11s %s\n", "arm", "between-chain sd, true cols (lower wins)"))
  for (i in seq_len(nrow(verdicts))) {
    cat(sprintf(
      "%-11s %s\n",
      verdicts$arm[i],
      verdicts$primaryImprovement[i]
    ))
  }
} else {
  verdicts <- NULL
}

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    columns = columns,
    verdicts = verdicts,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      p = p,
      nTest = nTest,
      nChains = nChains,
      nBurn = nBurn,
      nSamples = nSamples,
      signal = signal,
      neighbours = neighbours,
      arms = arms,
      armControl = armControl
    )
  ),
  outputDir,
  "P5-checkerboard"
)
