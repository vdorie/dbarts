#!/usr/bin/env Rscript

# C1, the He and Hahn factorial: the average-case core cell, and the only
# one with a published BART number. Thirty correlated continuous
# predictors, ten thousand rows, moderate noise, and two of the paper's four
# mean functions - Trig+poly for the interaction and Single index for the
# rotated ridge.
#
# Published numbers to reproduce (arXiv 2002.03375v4 Table 4, BART column,
# kappa = 1, averaged over 100 replications):
#
#   Trig+Poly     95% pointwise coverage 0.74, interval length 2.89, RMSE 1.27
#   Single Index  95% pointwise coverage 0.73, interval length 4.62, RMSE 2.08
#
# Primary statistic: 95% pointwise coverage of the true mean function. The
# paper's own reading of the deficit is that it "may indicate inadequate
# chain length of BART (that is, poor mixing)". Interval length and RMSE
# come with it because they distinguish a coverage miss caused by a
# mis-specified generating process from one caused by the sampler.
#
# Secondary: minimum ESS over 25 fixed held-out points, which is the
# worst-coordinate statistic the correlated-design literature ranks on.
#
# The paper's section 5, where Table 4 lives, fixes the sample size, the
# noise level and the chain length but not which of section 4.1's two
# predictor arms it used, and not the tree count. Both are therefore run as
# diagnostic arms alongside the pre-registered one:
#
#   correlated75      correlated factor design, shipped default tree count
#   correlated75grow  correlated75 plus n.grow.sweeps = 5 (XBART grow-from-root
#                     warm start; k = 5 is stochtree's own default, num_gfr,
#                     the only published ecosystem default for this count -
#                     docs/design/grow-from-root-default.md - and is used here
#                     since neither man/bart2.Rd nor docs/design/grow-from-
#                     root.md names a study value of its own)
#   correlated200     correlated factor design, 200 trees
#   independent75     independent standard normal design, shipped default
#   independent200    independent standard normal design, 200 trees
#
# Chain length follows the paper: one chain, 1000 burn-in, 2500 kept.
#
# Usage: Rscript C1-he-hahn.R [outputDir] [quick]

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
n <- if (quick) 2000L else 10000L
nTest <- 1000L
p <- 30L
kappa <- 1
nBurn <- 1000L
nSamples <- if (quick) 500L else 2500L

meanFunctions <- c("trigpoly", "singleindex")
published <- data.frame(
  meanFunction = meanFunctions,
  coverage = c(0.74, 0.73),
  length = c(2.89, 4.62),
  rmse = c(1.27, 2.08),
  stringsAsFactors = FALSE
)

arms <- list(
  correlated75 = list(design = "correlated", nTrees = 75L, growSweeps = 0L),
  correlated75grow = list(design = "correlated", nTrees = 75L, growSweeps = 5L),
  correlated200 = list(design = "correlated", nTrees = 200L, growSweeps = 0L),
  independent75 = list(design = "independent", nTrees = 75L, growSweeps = 0L),
  independent200 = list(design = "independent", nTrees = 200L, growSweeps = 0L)
)

# 25 evenly spaced held-out rows carry the ESS, as the move-set grid does.
essPoints <- as.integer(round(seq(1, nTest, length.out = 25L)))

surfacesUptime("uptime before")

rows <- list()
for (which in meanFunctions) {
  for (armName in names(arms)) {
    arm <- arms[[armName]]
    for (replicate in seq_len(nReplicates)) {
      set.seed(surfacesDataSeed("C1", paste0(which, arm$design), replicate))
      data <- surfacesHeHahn(n, nTest, p, which, kappa, design = arm$design)
      startedAt <- proc.time()
      fit <- bart2(
        data$x,
        data$y,
        test = data$xTest,
        n.trees = arm$nTrees,
        n.chains = 1L,
        n.burn = nBurn,
        n.samples = nSamples,
        n.thin = 1L,
        n.threads = 1L,
        n.grow.sweeps = arm$growSweeps,
        verbose = FALSE,
        seed = surfacesSamplerSeed(replicate)
      )
      elapsed <- (proc.time() - startedAt)[["elapsed"]]
      trainDraws <- extract(fit, type = "ev", sample = "train")
      testDraws <- extract(fit, type = "ev", sample = "test")
      trainInterval <- apply(
        trainDraws,
        2L,
        quantile,
        probs = c(0.025, 0.975),
        names = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        meanFunction = which,
        arm = armName,
        replicate = replicate,
        coverageTrain = mean(
          data$f >= trainInterval[1L, ] & data$f <= trainInterval[2L, ]
        ),
        lengthTrain = mean(trainInterval[2L, ] - trainInterval[1L, ]),
        rmseTrain = surfacesRmse(trainDraws, data$f),
        coverageTest = surfacesCoverage(testDraws, data$fTest),
        rmseTest = surfacesRmse(testDraws, data$fTest),
        minEss = min(surfacesPointEss(testDraws, essPoints)),
        sigmaTruth = data$sigma,
        sigmaPosterior = mean(fit$sigma),
        wall = elapsed,
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "%-12s %-14s rep %2d  cover %.3f  len %.2f  rmse %.2f  %.0fs\n",
        which,
        armName,
        replicate,
        rows[[length(rows)]]$coverageTrain,
        rows[[length(rows)]]$lengthTrain,
        rows[[length(rows)]]$rmseTrain,
        elapsed
      ))
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader("C1 He and Hahn: mean over seeds (min-max), in-sample f")
cat(sprintf(
  "%-12s %-14s %-20s %-20s %-20s %-16s %s\n",
  "mean fn",
  "arm",
  "95% coverage",
  "interval length",
  "RMSE",
  "held-out cover",
  "min ESS"
))
for (which in meanFunctions) {
  for (armName in names(arms)) {
    keep <- results$meanFunction == which & results$arm == armName
    cat(sprintf(
      "%-12s %-14s %-20s %-20s %-20s %-16s %s\n",
      which,
      armName,
      surfacesRange(results$coverageTrain[keep]),
      surfacesRange(results$lengthTrain[keep], digits = 2L),
      surfacesRange(results$rmseTrain[keep], digits = 2L),
      surfacesRange(results$coverageTest[keep]),
      surfacesRange(results$minEss[keep], digits = 0L)
    ))
  }
}

surfacesHeader("published reference (He and Hahn 2023 Table 4, BART, kappa 1)")
for (i in seq_len(nrow(published))) {
  cat(sprintf(
    "%-12s coverage %.2f  length %.2f  rmse %.2f\n",
    published$meanFunction[i],
    published$coverage[i],
    published$length[i],
    published$rmse[i]
  ))
}

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    published = published,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      nTest = nTest,
      p = p,
      kappa = kappa,
      nBurn = nBurn,
      nSamples = nSamples,
      arms = arms,
      essPoints = essPoints
    )
  ),
  outputDir,
  "C1-he-hahn"
)
