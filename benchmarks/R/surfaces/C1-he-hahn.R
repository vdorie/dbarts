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
#   independent75grow independent75 plus n.grow.sweeps = 5, same rationale
#                     as correlated75grow
#   independent200    independent standard normal design, 200 trees
#
# Chain length follows the paper: one chain, 1000 burn-in, 2500 kept. Three
# further arms vary only that, on the independent 75-tree cell, to tell a
# coverage deficit caused by chain exploration apart from one caused by the
# tree count:
#
#   independent75pool4      four chains pooled at 500 burn-in and 500 kept
#                           each, which is bart2's shipped chain default
#   independent75pool4long  four chains pooled at 1000 burn-in, 2500 kept
#   independent75long       one chain, 1000 burn-in, 25000 kept, with coverage
#                           also read off the first 2500 of those draws so the
#                           length effect is visible inside a single fit
#
# Every arm's chain settings live in the arms list, so the six above keep the
# paper's single chain. Pooled arms are fit with combineChains = FALSE and
# report, beside the pooled coverage: minimum ESS summed over chains, the
# median over chains of each chain's own minimum, and a between-chain ratio -
# the across-chain standard deviation of a chain's posterior mean of f over
# the pooled posterior standard deviation, median over the 25 ESS points.
# Near 0 the chains agree; near 1 each sits in its own place and pooling is
# what widens the interval.
#
# Usage: Rscript C1-he-hahn.R [outputDir] [quick] [arm ...]

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

# Arms carry their own chain settings; omitting them takes the paper's single
# chain at the top-level burn-in and sample counts. `prefix` asks for a second
# set of readouts from the first that many kept draws of each chain.
armSpec <- function(
  design,
  nTrees,
  growSweeps = 0L,
  nChains = 1L,
  armBurn = nBurn,
  armSamples = nSamples,
  prefix = NA_integer_
) {
  list(
    design = design,
    nTrees = nTrees,
    growSweeps = growSweeps,
    nChains = nChains,
    nBurn = armBurn,
    nSamples = armSamples,
    prefix = prefix
  )
}

arms <- list(
  correlated75 = armSpec("correlated", 75L),
  correlated75grow = armSpec("correlated", 75L, growSweeps = 5L),
  correlated200 = armSpec("correlated", 200L),
  independent75 = armSpec("independent", 75L),
  independent75grow = armSpec("independent", 75L, growSweeps = 5L),
  independent200 = armSpec("independent", 200L),
  independent75pool4 = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = if (quick) 100L else 500L
  ),
  independent75pool4long = armSpec("independent", 75L, nChains = 4L),
  independent75long = armSpec(
    "independent",
    75L,
    armSamples = if (quick) 2000L else 25000L,
    prefix = if (quick) 500L else 2500L
  )
)
selectedArms <- intersect(names(arms), args)
if (length(selectedArms) > 0L) {
  arms <- arms[selectedArms]
}

outputDir <- surfacesOutputDir(args, flags = c("quick", names(arms)))

# 25 evenly spaced held-out rows carry the ESS, as the move-set grid does.
essPoints <- as.integer(round(seq(1, nTest, length.out = 25L)))

# --- chain-aware readouts --------------------------------------------------

# Pooled draw matrices are chain-major, so a chain is a contiguous row block.
armChainDraws <- function(draws, chain, nSamples) {
  draws[surfacesChainRows(chain, nSamples), , drop = FALSE]
}

# The first `prefix` kept draws of every chain, as row indices into a pooled
# matrix: the same fit read at a shorter chain length.
armPrefixRows <- function(nChains, nSamples, prefix) {
  unlist(lapply(
    seq_len(nChains),
    function(chain) surfacesChainRows(chain, nSamples)[seq_len(prefix)]
  ))
}

# ESS at each ESS point of each chain on its own, as points x chains. One
# chain reduces to surfacesPointEss on the pooled matrix.
armChainEss <- function(draws, nChains, nSamples, points) {
  vapply(
    seq_len(nChains),
    function(chain) {
      surfacesPointEss(armChainDraws(draws, chain, nSamples), points)
    },
    numeric(length(points))
  )
}

# Summing ESS across chains before taking the worst point is the pooled
# statistic; the median of the per-chain worst points says whether that sum is
# four mixing chains or one chain and three bystanders.
armEssSummedMin <- function(ess) {
  min(rowSums(ess))
}

armEssChainMedianMin <- function(ess) {
  median(apply(ess, 2L, min))
}

# Across-chain spread of the posterior mean of f over the pooled posterior
# spread, median over the ESS points. Near 0 the chains agree and pooling adds
# nothing; near 1 each chain sits in its own place and pooling is what widens
# the interval.
armBetweenChainRatio <- function(draws, nChains, nSamples, points) {
  if (nChains < 2L) {
    return(NA_real_)
  }
  median(vapply(
    points,
    function(j) {
      chainMeans <- vapply(
        seq_len(nChains),
        function(chain) mean(draws[surfacesChainRows(chain, nSamples), j]),
        numeric(1L)
      )
      sd(chainMeans) / sd(draws[, j])
    },
    numeric(1L)
  ))
}

armInterval <- function(draws) {
  apply(draws, 2L, quantile, probs = c(0.025, 0.975), names = FALSE)
}

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
        n.chains = arm$nChains,
        n.burn = arm$nBurn,
        n.samples = arm$nSamples,
        n.thin = 1L,
        n.threads = 1L,
        n.grow.sweeps = arm$growSweeps,
        combineChains = arm$nChains == 1L,
        verbose = FALSE,
        seed = surfacesSamplerSeed(replicate)
      )
      elapsed <- (proc.time() - startedAt)[["elapsed"]]
      # Both extractions pool; the chains are recovered by row block below.
      trainDraws <- extract(fit, type = "ev", sample = "train")
      testDraws <- extract(fit, type = "ev", sample = "test")
      trainInterval <- armInterval(trainDraws)
      ess <- armChainEss(testDraws, arm$nChains, arm$nSamples, essPoints)
      prefixRows <- if (is.na(arm$prefix)) {
        integer(0)
      } else {
        armPrefixRows(arm$nChains, arm$nSamples, arm$prefix)
      }
      prefixTrain <- trainDraws[prefixRows, , drop = FALSE]
      prefixInterval <- if (length(prefixRows) > 0L) {
        armInterval(prefixTrain)
      } else {
        NULL
      }
      rows[[length(rows) + 1L]] <- data.frame(
        meanFunction = which,
        arm = armName,
        replicate = replicate,
        nChains = arm$nChains,
        nBurn = arm$nBurn,
        nSamples = arm$nSamples,
        coverageTrain = mean(
          data$f >= trainInterval[1L, ] & data$f <= trainInterval[2L, ]
        ),
        lengthTrain = mean(trainInterval[2L, ] - trainInterval[1L, ]),
        rmseTrain = surfacesRmse(trainDraws, data$f),
        coverageTest = surfacesCoverage(testDraws, data$fTest),
        rmseTest = surfacesRmse(testDraws, data$fTest),
        minEss = armEssSummedMin(ess),
        minEssChainMedian = armEssChainMedianMin(ess),
        betweenChain = armBetweenChainRatio(
          testDraws,
          arm$nChains,
          arm$nSamples,
          essPoints
        ),
        prefixSamples = arm$prefix,
        coveragePrefixTrain = if (is.null(prefixInterval)) {
          NA_real_
        } else {
          mean(
            data$f >= prefixInterval[1L, ] & data$f <= prefixInterval[2L, ]
          )
        },
        lengthPrefixTrain = if (is.null(prefixInterval)) {
          NA_real_
        } else {
          mean(prefixInterval[2L, ] - prefixInterval[1L, ])
        },
        rmsePrefixTrain = if (is.null(prefixInterval)) {
          NA_real_
        } else {
          surfacesRmse(prefixTrain, data$f)
        },
        coveragePrefixTest = if (is.null(prefixInterval)) {
          NA_real_
        } else {
          surfacesCoverage(
            testDraws[prefixRows, , drop = FALSE],
            data$fTest
          )
        },
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
      rm(fit, trainDraws, testDraws, prefixTrain)
      invisible(gc(verbose = FALSE))
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader("C1 He and Hahn: mean over seeds (min-max), in-sample f")
cat(sprintf(
  "%-12s %-22s %-20s %-20s %-20s %-16s %s\n",
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
      "%-12s %-22s %-20s %-20s %-20s %-16s %s\n",
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

surfacesHeader("C1 chain diagnostics: mean over seeds (min-max)")
cat(sprintf(
  "%-12s %-22s %-9s %-16s %-16s %-16s %s\n",
  "mean fn",
  "arm",
  "chains",
  "min ESS summed",
  "min ESS/chain",
  "between-chain",
  "wall s"
))
for (which in meanFunctions) {
  for (armName in names(arms)) {
    keep <- results$meanFunction == which & results$arm == armName
    cat(sprintf(
      "%-12s %-22s %-9s %-16s %-16s %-16s %.1f\n",
      which,
      armName,
      sprintf("%d x %d", arms[[armName]]$nChains, arms[[armName]]$nSamples),
      surfacesRange(results$minEss[keep], digits = 0L),
      surfacesRange(results$minEssChainMedian[keep], digits = 0L),
      surfacesRange(results$betweenChain[keep], digits = 2L),
      mean(results$wall[keep])
    ))
  }
}

prefixArms <- names(arms)[vapply(arms, function(a) !is.na(a$prefix), TRUE)]
if (length(prefixArms) > 0L) {
  surfacesHeader("C1 chain length within one fit: first draws against all kept")
  cat(sprintf(
    "%-12s %-22s %-9s %-20s %-20s %-20s %s\n",
    "mean fn",
    "arm",
    "draws",
    "95% coverage",
    "interval length",
    "RMSE",
    "held-out cover"
  ))
  for (which in meanFunctions) {
    for (armName in prefixArms) {
      keep <- results$meanFunction == which & results$arm == armName
      if (!any(keep)) {
        next
      }
      spec <- arms[[armName]]
      cat(sprintf(
        "%-12s %-22s %-9d %-20s %-20s %-20s %s\n",
        which,
        armName,
        spec$prefix * spec$nChains,
        surfacesRange(results$coveragePrefixTrain[keep]),
        surfacesRange(results$lengthPrefixTrain[keep], digits = 2L),
        surfacesRange(results$rmsePrefixTrain[keep], digits = 2L),
        surfacesRange(results$coveragePrefixTest[keep])
      ))
      cat(sprintf(
        "%-12s %-22s %-9d %-20s %-20s %-20s %s\n",
        which,
        armName,
        spec$nSamples * spec$nChains,
        surfacesRange(results$coverageTrain[keep]),
        surfacesRange(results$lengthTrain[keep], digits = 2L),
        surfacesRange(results$rmseTrain[keep], digits = 2L),
        surfacesRange(results$coverageTest[keep])
      ))
    }
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
