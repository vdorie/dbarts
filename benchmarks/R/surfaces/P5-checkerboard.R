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
# Usage: Rscript P5-checkerboard.R [outputDir] [quick]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
outputDir <- surfacesOutputDir(args, flags = "quick")

nReplicates <- if (quick) 2L else 20L
n <- 1600L
p <- 40L
nTest <- 1000L
nChains <- 8L
nBurn <- 1000L
nSamples <- if (quick) 400L else 2000L

signal <- surfacesCheckerboardSignal
neighbours <- surfacesCheckerboardNeighbours(p)
decoys <- setdiff(seq_len(p), signal)

surfacesUptime("uptime before")

rows <- list()
perColumn <- list()
for (replicate in seq_len(nReplicates)) {
  set.seed(surfacesDataSeed("P5", "checkerboard", replicate))
  data <- surfacesCheckerboard(n = n, p = p, nTest = nTest)
  startedAt <- proc.time()
  fit <- bart2(
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

  perColumn[[replicate]] <- data.frame(
    replicate = replicate,
    column = seq_len(p),
    inclusion = pooledInclusion,
    betweenChainSd = betweenSd,
    mixingNull = mixingNull,
    stringsAsFactors = FALSE
  )
  rows[[replicate]] <- data.frame(
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
    "rep %2d  incl(true) %.3f  sd/null %.1f  decoys above weakest %d  cover %.3f  %.0fs\n",
    replicate,
    rows[[replicate]]$inclusionSignal,
    rows[[replicate]]$signalRatio,
    rows[[replicate]]$decoysAboveWeakest,
    rows[[replicate]]$coverage,
    elapsed
  ))
}
results <- do.call(rbind, rows)
columns <- do.call(rbind, perColumn)

surfacesHeader("P5 checkerboard: mean over seeds (min-max)")
cat(sprintf(
  "%-28s %s\n",
  "between-chain sd, true cols",
  surfacesRange(results$signalBetweenSd, digits = 4L)
))
cat(sprintf(
  "%-28s %s\n",
  "mixing null, true cols",
  surfacesRange(results$signalMixingNull, digits = 4L)
))
cat(sprintf(
  "%-28s %s\n",
  "ratio to null, true cols",
  surfacesRange(results$signalRatio, digits = 2L)
))
cat(sprintf(
  "%-28s %s\n",
  "between-chain sd, decoys",
  surfacesRange(results$neighbourBetweenSd, digits = 4L)
))
cat(sprintf(
  "%-28s %s\n",
  "ratio to null, decoys",
  surfacesRange(results$neighbourRatio, digits = 2L)
))
cat(sprintf(
  "%-28s %s\n",
  "inclusion share, true cols",
  surfacesRange(results$inclusionSignal)
))
cat(sprintf(
  "%-28s %s\n",
  "inclusion share, decoys",
  surfacesRange(results$inclusionNeighbour)
))
cat(sprintf(
  "%-28s %s\n",
  "largest non-true inclusion",
  surfacesRange(results$maxDecoy)
))
cat(sprintf(
  "%-28s %s\n",
  "non-true cols above weakest",
  surfacesRange(results$decoysAboveWeakest, digits = 1L)
))
cat(sprintf(
  "%-28s %s\n",
  "95% coverage of true f",
  surfacesRange(results$coverage)
))
cat(sprintf("%-28s %s\n", "held-out RMSE", surfacesRange(results$rmse)))

surfacesHeader("pooled inclusion by column, mean over seeds")
pooled <- tapply(columns$inclusion, columns$column, mean)
for (j in seq_len(p)) {
  tag <- if (j %in% signal) {
    "TRUE"
  } else if (j %in% neighbours) {
    "decoy"
  } else {
    ""
  }
  cat(sprintf("x%-3d %.4f %s\n", j, pooled[[j]], tag))
}

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    columns = columns,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      p = p,
      nTest = nTest,
      nChains = nChains,
      nBurn = nBurn,
      nSamples = nSamples,
      signal = signal,
      neighbours = neighbours
    )
  ),
  outputDir,
  "P5-checkerboard"
)
