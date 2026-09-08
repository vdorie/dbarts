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
# Eight further arms sit on that same shipped four-chain configuration and the
# same twenty seeds, varying only the proposal mixture or the level-fibre
# flag, so the kernel contrast is read where the shipped chain default reads
# it rather than at the one-chain configuration the earlier move-set grid
# used:
#
#   independent75pool4bd           birth_death 1, swap 0, change 0
#   independent75pool4swap         birth_death 0.5, swap 0.1, change 0.4, the
#                                  former default
#   independent75pool4perturbB     birth_death 0.6, swap 0, change 0.24,
#                                  perturb 0.16, the shipped mixture with the
#                                  0.16 taken out of change alone
#   independent75pool4perturbMixed birth_death 0.5, swap 0, change 0.34,
#                                  perturb 0.16, which takes the 0.16 out of
#                                  change and birth/death both. It is NOT the
#                                  arm above and its readouts do not stand in
#                                  for it.
#   independent75pool4ruleGibbsB   birth_death 0.6, swap 0, change 0.24,
#                                  rule_gibbs 0.16, the same share taken out
#                                  of change alone. This is the dosage the
#                                  rule_gibbs benefit study's kill criterion
#                                  turns on.
#   independent75pool4ruleGibbs32  birth_death 0.6, swap 0, change 0.08,
#                                  rule_gibbs 0.32, the second dosage of that
#                                  study's grid. It is reported, not gated.
#   independent75pool4level        the shipped mixture unchanged, with the
#                                  level-fibre Gibbs step on. It varies no
#                                  proposal probability at all: the step is an
#                                  exact draw on the leaf fibre, taken once a
#                                  sweep ahead of the tree loop, and it is the
#                                  dosage the level-fibre benefit study's kill
#                                  criterion turns on.
#   independent75pool4ruleGibbsBlevel
#                                  the rule_gibbs arm above and the level step
#                                  together, one kernel structural and one on
#                                  the leaf fibre. It is reported, not gated,
#                                  and it is what says whether the two gains
#                                  stack.
#
# independent75pool4 is their control and is re-run beside them so the
# contrast is paired within one session. Both perturb arms are PILOT levels,
# not a confirmatory run of the perturb design; the two rule_gibbs arms are
# that design's own confirmatory run, its arm B being the first of them, and
# independent75pool4level is the level-fibre design's.
#
# A further arm changes no setting at all, only the seed:
#
#   independent75pool4sham         the shipped mixture again, on the same data
#                                  seeds but at sampler seeds offset by 1000.
#                                  It is the same kernel as its control, so
#                                  its paired difference measures the
#                                  harness's own seed-to-seed spread, which is
#                                  what a bar on that difference has to clear.
#
# Two run options steer the seeds rather than the arms. An arm's samplerOffset
# shifts that arm's sampler seed alone, leaving its data seed where the
# control's is, which is what makes the sham arm a second draw on the same
# data. `seedBlock=k` shifts the data-seed and sampler-seed index together by
# whole blocks of the replicate count, so block 2 runs seeds 21 to 40 and a
# flagged contrast can be re-run on seeds no arm has seen. A bare mean
# function name restricts the run to that mean function.
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
# Usage: Rscript C1-he-hahn.R [outputDir] [quick] [seedBlock=k] [meanFn ...]
#                             [arm ...]

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

# `seedBlock=k` moves the replicate index on by whole blocks, so block 1 is
# replicates 1 to nReplicates and block 2 the next nReplicates. Both the data
# seed and the sampler seed are indexed by it, so a block is a set of seeds no
# other block has drawn: this is the fresh-seed re-run the battery's rule
# requires of a flagged cell.
seedBlockArg <- grep("^seedBlock=", args, value = TRUE)
seedBlock <- if (length(seedBlockArg) > 0L) {
  as.integer(sub("^seedBlock=", "", seedBlockArg[1L]))
} else {
  1L
}
if (is.na(seedBlock) || seedBlock < 1L) {
  stop("seedBlock must be a positive integer")
}

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
# set of readouts from the first that many kept draws of each chain. `probs`
# is a proposal mixture; NULL leaves it unset, which is the shipped one.
# `levelGibbs` switches on the level-fibre Gibbs step and is named in every
# call: FALSE is bart2's own default, so naming it leaves every other arm's
# draws where they were. `samplerOffset` adds to the sampler seed and to
# nothing else, so an arm carrying one sees its control's data at a different
# MCMC stream.
armSpec <- function(
  design,
  nTrees,
  growSweeps = 0L,
  nChains = 1L,
  armBurn = nBurn,
  armSamples = nSamples,
  prefix = NA_integer_,
  probs = NULL,
  levelGibbs = FALSE,
  samplerOffset = 0L
) {
  list(
    design = design,
    nTrees = nTrees,
    growSweeps = growSweeps,
    nChains = nChains,
    nBurn = armBurn,
    nSamples = armSamples,
    prefix = prefix,
    probs = probs,
    levelGibbs = levelGibbs,
    samplerOffset = samplerOffset
  )
}

pool4Samples <- if (quick) 100L else 500L

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
    armSamples = pool4Samples
  ),
  independent75pool4long = armSpec("independent", 75L, nChains = 4L),
  independent75long = armSpec(
    "independent",
    75L,
    armSamples = if (quick) 2000L else 25000L,
    prefix = if (quick) 500L else 2500L
  ),
  independent75pool4bd = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(birth_death = 1, swap = 0, change = 0, perturb = 0, birth = 0.5)
  ),
  independent75pool4swap = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.5,
      swap = 0.1,
      change = 0.4,
      perturb = 0,
      birth = 0.5
    )
  ),
  independent75pool4perturbB = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.24,
      perturb = 0.16,
      birth = 0.5
    )
  ),
  independent75pool4perturbMixed = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.5,
      swap = 0,
      change = 0.34,
      perturb = 0.16,
      birth = 0.5
    )
  ),
  independent75pool4ruleGibbsB = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.24,
      perturb = 0,
      rule_gibbs = 0.16,
      birth = 0.5
    )
  ),
  independent75pool4ruleGibbs32 = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.08,
      perturb = 0,
      rule_gibbs = 0.32,
      birth = 0.5
    )
  ),
  independent75pool4level = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    levelGibbs = TRUE
  ),
  independent75pool4ruleGibbsBlevel = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.24,
      perturb = 0,
      rule_gibbs = 0.16,
      birth = 0.5
    ),
    levelGibbs = TRUE
  ),
  independent75pool4sham = armSpec(
    "independent",
    75L,
    nChains = 4L,
    armBurn = 500L,
    armSamples = pool4Samples,
    samplerOffset = 1000L
  )
)

# The move-set arms above are read against the shipped mixture at the same
# chain configuration, so that arm is their paired control. The sham arm is
# read the same way and is that reading's calibration: it differs from the
# control in nothing but its sampler seed.
movesetControl <- "independent75pool4"
movesetArms <- c(
  "independent75pool4bd",
  "independent75pool4swap",
  "independent75pool4perturbB",
  "independent75pool4perturbMixed",
  "independent75pool4ruleGibbsB",
  "independent75pool4ruleGibbs32",
  "independent75pool4level",
  "independent75pool4ruleGibbsBlevel",
  "independent75pool4sham"
)
armNames <- names(arms)
selectedArms <- intersect(armNames, args)
if (length(selectedArms) > 0L) {
  arms <- arms[selectedArms]
}

selectedMeanFunctions <- intersect(meanFunctions, args)
if (length(selectedMeanFunctions) > 0L) {
  meanFunctions <- selectedMeanFunctions
}

outputDir <- surfacesOutputDir(
  args,
  flags = c("quick", seedBlockArg, armNames, meanFunctions)
)

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

# The pooled interval taken apart chain by chain: each chain's own posterior
# mean of f and its own 95% interval at the ESS points, as points x chains.
# The between-chain ratio above says how far the chain means sit apart in
# units of the pooled spread; this says whether the intervals themselves
# still meet.
armChainSummary <- function(draws, nChains, nSamples, points) {
  means <- matrix(NA_real_, length(points), nChains)
  lower <- means
  upper <- means
  for (chain in seq_len(nChains)) {
    chainDraws <- armChainDraws(draws, chain, nSamples)[, points, drop = FALSE]
    interval <- armInterval(chainDraws)
    means[, chain] <- colMeans(chainDraws)
    lower[, chain] <- interval[1L, ]
    upper[, chain] <- interval[2L, ]
  }
  list(mean = means, lower = lower, upper = upper)
}

# How much two chains' intervals share at one point: the length they have in
# common over the length they span together, 0 when they are disjoint and 1
# when they coincide. Reported per point as the mean over the chain pairs and
# as the share of pairs that meet at all.
armIntervalOverlap <- function(summary) {
  pairs <- combn(ncol(summary$mean), 2L)
  shared <- vapply(
    seq_len(ncol(pairs)),
    function(k) {
      a <- pairs[1L, k]
      b <- pairs[2L, k]
      common <- pmin(summary$upper[, a], summary$upper[, b]) -
        pmax(summary$lower[, a], summary$lower[, b])
      spanned <- pmax(summary$upper[, a], summary$upper[, b]) -
        pmin(summary$lower[, a], summary$lower[, b])
      pmax(common, 0) / spanned
    },
    numeric(nrow(summary$mean))
  )
  list(shared = rowMeans(shared), meeting = rowMeans(shared > 0))
}

# A paired difference against the control arm, in the shape the move-set
# tables quote raw signals in: mean, standard deviation, and the count of
# seeds on which the difference is positive.
armPairedDifference <- function(x, digits = 2L) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return("      -")
  }
  fmt <- paste0("%+.", digits, "f +/- %.", digits, "f (%d/%d)")
  sprintf(fmt, mean(x), sd(x), sum(x > 0), length(x))
}

surfacesUptime("uptime before")

rows <- list()
chainSummaries <- list()
for (which in meanFunctions) {
  for (armName in names(arms)) {
    arm <- arms[[armName]]
    for (replicate in seq_len(nReplicates)) {
      seedIndex <- (seedBlock - 1L) * nReplicates + replicate
      set.seed(surfacesDataSeed("C1", paste0(which, arm$design), seedIndex))
      data <- surfacesHeHahn(n, nTest, p, which, kappa, design = arm$design)
      call <- list(
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
        levelGibbs = arm$levelGibbs,
        verbose = FALSE,
        seed = surfacesSamplerSeed(seedIndex) + arm$samplerOffset
      )
      if (!is.null(arm$probs)) {
        call$proposal.probs <- arm$probs
      }
      startedAt <- proc.time()
      fit <- do.call(bart2, call)
      elapsed <- (proc.time() - startedAt)[["elapsed"]]
      # Both extractions pool; the chains are recovered by row block below.
      trainDraws <- extract(fit, type = "ev", sample = "train")
      testDraws <- extract(fit, type = "ev", sample = "test")
      trainInterval <- armInterval(trainDraws)
      ess <- armChainEss(testDraws, arm$nChains, arm$nSamples, essPoints)
      # The chain-by-chain intervals are kept for the control arm alone, which
      # is the one the move-set arms are read against.
      chainSummary <- if (armName == movesetControl && arm$nChains > 1L) {
        armChainSummary(testDraws, arm$nChains, arm$nSamples, essPoints)
      } else {
        NULL
      }
      overlap <- if (is.null(chainSummary)) {
        NULL
      } else {
        armIntervalOverlap(chainSummary)
      }
      if (!is.null(chainSummary)) {
        chainSummaries[[length(chainSummaries) + 1L]] <- list(
          meanFunction = which,
          arm = armName,
          replicate = seedIndex,
          points = essPoints,
          mean = chainSummary$mean,
          lower = chainSummary$lower,
          upper = chainSummary$upper
        )
      }
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
        replicate = seedIndex,
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
        chainOverlap = if (is.null(overlap)) {
          NA_real_
        } else {
          median(overlap$shared)
        },
        chainOverlapMeeting = if (is.null(overlap)) {
          NA_real_
        } else {
          median(overlap$meeting)
        },
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
        seedIndex,
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
  "%-12s %-31s %-20s %-20s %-20s %-16s %s\n",
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
      "%-12s %-31s %-20s %-20s %-20s %-16s %s\n",
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
  "%-12s %-31s %-9s %-16s %-16s %-16s %s\n",
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
      "%-12s %-31s %-9s %-16s %-16s %-16s %.1f\n",
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
    "%-12s %-31s %-9s %-20s %-20s %-20s %s\n",
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
        "%-12s %-31s %-9d %-20s %-20s %-20s %s\n",
        which,
        armName,
        spec$prefix * spec$nChains,
        surfacesRange(results$coveragePrefixTrain[keep]),
        surfacesRange(results$lengthPrefixTrain[keep], digits = 2L),
        surfacesRange(results$rmsePrefixTrain[keep], digits = 2L),
        surfacesRange(results$coveragePrefixTest[keep])
      ))
      cat(sprintf(
        "%-12s %-31s %-9d %-20s %-20s %-20s %s\n",
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

contrastArms <- intersect(movesetArms, names(arms))
if (movesetControl %in% names(arms) && length(contrastArms) > 0L) {
  surfacesHeader(sprintf(
    "C1 move sets: paired difference against %s, seed by seed",
    movesetControl
  ))
  cat(sprintf(
    "%-12s %-31s %-31s %-31s %s\n",
    "mean fn",
    "arm",
    "summed min ESS",
    "95% coverage",
    "RMSE"
  ))
  for (which in meanFunctions) {
    control <- results[
      results$meanFunction == which & results$arm == movesetControl,
    ]
    for (armName in contrastArms) {
      contrast <- results[
        results$meanFunction == which & results$arm == armName,
      ]
      paired <- match(control$replicate, contrast$replicate)
      contrast <- contrast[paired[!is.na(paired)], ]
      base <- control[!is.na(paired), ]
      if (nrow(contrast) == 0L) {
        next
      }
      cat(sprintf(
        "%-12s %-31s %-31s %-31s %s\n",
        which,
        armName,
        armPairedDifference(contrast$minEss - base$minEss, digits = 1L),
        armPairedDifference(
          contrast$coverageTrain - base$coverageTrain,
          digits = 3L
        ),
        armPairedDifference(contrast$rmseTrain - base$rmseTrain, digits = 3L)
      ))
    }
  }
}

if (any(is.finite(results$chainOverlap))) {
  surfacesHeader(sprintf(
    "C1 chain agreement at the 25 ESS points: %s",
    movesetControl
  ))
  cat(sprintf(
    "%-12s %-31s %-24s %s\n",
    "mean fn",
    "arm",
    "interval overlap",
    "pairs that meet"
  ))
  for (which in meanFunctions) {
    keep <- results$meanFunction == which & is.finite(results$chainOverlap)
    if (!any(keep)) {
      next
    }
    cat(sprintf(
      "%-12s %-31s %-24s %s\n",
      which,
      movesetControl,
      surfacesRange(results$chainOverlap[keep], digits = 2L),
      surfacesRange(results$chainOverlapMeeting[keep], digits = 2L)
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
    chainSummaries = chainSummaries,
    settings = list(
      nReplicates = nReplicates,
      seedBlock = seedBlock,
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
