#!/usr/bin/env Rscript

# How much of C1's autocorrelation lives in the leaf values, given fixed
# structures?
#
# The recorded C1 arm (independent design, 75 trees, one chain of 1000 burn-in
# and 2500 kept) reads a minimum ESS over 25 fixed test points of about 2 out
# of 2500 kept draws, in a regime whose pooled structural acceptance is a loose
# 25 percent. Two readings fit that: the structural half of the kernel is what
# fails to move, or the leaf Gibbs is slow on its own and no structural kernel
# would rescue it.
#
# The frozen mixture separates them. Fit the recorded arm, then freeze the tree
# structures - all four structural probabilities zero, so no move is proposed -
# and run the same number of draws again. The frozen chain is the leaf Gibbs
# alone, conditional on one fixed forest:
#
#   frozen ESS near 2500  the leaf half mixes freely and the deficit is the
#                         structural half's
#   frozen ESS near 2     the leaf Gibbs itself is the slow part
#
# Two frozen chains are run per seed, one from the last kept draw and one from
# the 1250th, so a frozen number cannot be an artifact of the particular forest
# it was frozen at.
#
# Usage: Rscript C1-frozen-ess.R [outputDir] [quick]

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

# The recorded independent75 arm's cell, at its first five seeds.
nReplicates <- if (quick) 2L else 5L
n <- if (quick) 2000L else 10000L
nTest <- 1000L
p <- 30L
kappa <- 1
nTrees <- 75L
nBurn <- if (quick) 200L else 1000L
nSamples <- if (quick) 300L else 2500L
meanFunction <- "trigpoly"
design <- "independent"

# 25 evenly spaced held-out rows, as the recorded arm reads them.
essPoints <- as.integer(round(seq(1, nTest, length.out = 25L)))

# The frozen mixture: no structural proposal is made, so the trees stand where
# the structural chain left them and only the leaf values and sigma move.
freeze <- function(model) {
  model@p.birth_death <- 0
  model@p.swap <- 0
  model@p.change <- 0
  model@p.perturb <- 0
  model
}

# run() reports test fits as points x draws; the readouts want draws x points.
essDraws <- function(samples) {
  t(samples$test[essPoints, , drop = FALSE])
}

lag1 <- function(draws) {
  vapply(
    seq_len(ncol(draws)),
    function(j) cor(draws[-1L, j], draws[-nrow(draws), j]),
    numeric(1L)
  )
}

# One chain's numbers: the worst-point ESS the battery ranks on, the median
# point's ESS, and the median point's lag-1 autocorrelation. The worst point's
# posterior spread rides along, against the median point's: an ESS read on a
# coordinate that barely moves is a numerical statistic rather than a mixing
# one, and the ratio says which this is.
chainReadout <- function(draws) {
  ess <- surfacesPointEss(draws, seq_len(ncol(draws)))
  acf1 <- lag1(draws)
  spread <- apply(draws, 2L, sd)
  list(
    minEss = min(ess),
    medianEss = median(ess),
    lag1 = median(acf1),
    maxLag1 = max(acf1),
    sdRatio = spread[[which.min(ess)]] / median(spread)
  )
}

surfacesUptime("uptime before")

rows <- list()
for (replicate in seq_len(nReplicates)) {
  set.seed(surfacesDataSeed("C1", paste0(meanFunction, design), replicate))
  data <- surfacesHeHahn(n, nTest, p, meanFunction, kappa, design = design)

  startedAt <- proc.time()
  sampler <- bart2(
    data$x,
    data$y,
    test = data$xTest,
    n.trees = nTrees,
    n.chains = 1L,
    n.burn = nBurn,
    n.samples = nSamples,
    n.thin = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = surfacesSamplerSeed(replicate),
    samplerOnly = TRUE
  )
  # samplerOnly hands back the sampler before bart2 draws its initial forest
  sampler$sampleTreesFromPrior(updateState = FALSE)

  # the structural chain, split so the halfway state can be branched from
  half <- nSamples %/% 2L
  first <- essDraws(sampler$run(nBurn, half))
  sampler$storeState()
  midpoint <- sampler$copy()
  second <- essDraws(sampler$run(0L, nSamples - half))
  structural <- rbind(first, second)
  rm(first, second)

  frozenModel <- freeze(sampler$model)
  sampler$setModel(frozenModel)
  frozenLast <- essDraws(sampler$run(0L, nSamples))

  midpoint$setModel(frozenModel)
  frozenMid <- essDraws(midpoint$run(0L, nSamples))
  elapsed <- (proc.time() - startedAt)[["elapsed"]]

  readouts <- list(
    structural = chainReadout(structural),
    frozenLast = chainReadout(frozenLast),
    frozenMid = chainReadout(frozenMid)
  )
  for (chain in names(readouts)) {
    rows[[length(rows) + 1L]] <- data.frame(
      replicate = replicate,
      chain = chain,
      minEss = readouts[[chain]]$minEss,
      medianEss = readouts[[chain]]$medianEss,
      lag1 = readouts[[chain]]$lag1,
      maxLag1 = readouts[[chain]]$maxLag1,
      sdRatio = readouts[[chain]]$sdRatio,
      nSamples = nSamples,
      wall = elapsed,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf(
    "rep %d  structural min ESS %6.1f  frozen(last) %7.1f  frozen(mid) %7.1f  %.0fs\n",
    replicate,
    readouts$structural$minEss,
    readouts$frozenLast$minEss,
    readouts$frozenMid$minEss,
    elapsed
  ))
  rm(sampler, midpoint, structural, frozenLast, frozenMid)
  invisible(gc(verbose = FALSE))
}
results <- do.call(rbind, rows)

chains <- c("structural", "frozenLast", "frozenMid")
labels <- c(
  structural = "structural",
  frozenLast = "frozen at 2500",
  frozenMid = "frozen at 1250"
)

surfacesHeader(sprintf(
  "C1 frozen-forest ESS: %d kept draws, 25 held-out points",
  nSamples
))
cat(sprintf(
  "%-6s %-16s %-12s %-12s %-12s %-12s %s\n",
  "seed",
  "chain",
  "min ESS",
  "median ESS",
  "lag-1 med",
  "lag-1 max",
  "sd at min/med"
))
for (replicate in seq_len(nReplicates)) {
  for (chain in chains) {
    keep <- results$replicate == replicate & results$chain == chain
    cat(sprintf(
      "%-6d %-16s %-12.1f %-12.1f %-12.3f %-12.3f %.2f\n",
      replicate,
      labels[[chain]],
      results$minEss[keep],
      results$medianEss[keep],
      results$lag1[keep],
      results$maxLag1[keep],
      results$sdRatio[keep]
    ))
  }
}

surfacesHeader("median over seeds")
cat(sprintf(
  "%-16s %-12s %-12s %-12s %-12s %s\n",
  "chain",
  "min ESS",
  "median ESS",
  "lag-1 med",
  "lag-1 max",
  "sd at min/med"
))
for (chain in chains) {
  keep <- results$chain == chain
  cat(sprintf(
    "%-16s %-12.1f %-12.1f %-12.3f %-12.3f %.2f\n",
    labels[[chain]],
    median(results$minEss[keep]),
    median(results$medianEss[keep]),
    median(results$lag1[keep]),
    median(results$maxLag1[keep]),
    median(results$sdRatio[keep])
  ))
}

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      nTest = nTest,
      p = p,
      kappa = kappa,
      nTrees = nTrees,
      nBurn = nBurn,
      nSamples = nSamples,
      meanFunction = meanFunction,
      design = design,
      essPoints = essPoints
    )
  ),
  outputDir,
  "C1-frozen-ess"
)
