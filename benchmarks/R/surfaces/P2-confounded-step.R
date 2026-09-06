#!/usr/bin/env Rscript

# P2, the confounded step function: Pratola (2016) section 2.3, taking the
# problem from Wu, Tjelmeland and West (2007). Three columns, three hundred
# rows, one tree, and a design in which x1 and x3 are confounded, so
# root-on-x1 and root-on-x3 are two exactly equiprobable representations of
# the same fitted function.
#
# Published number to reproduce: "We fit BART to this dataset using only
# m = 1 trees and found that the acceptance rate of tree moves (after the
# initial few steps of the sampler) was 0" (arXiv 1312.1895 section 2.3).
#
# Primary statistic: the between-chain standard deviation of the fraction of
# draws whose root splits on x1, read against the exact 0.5 symmetry the
# design supplies. A mixing sampler puts every chain near 0.5 and the
# between-chain sd near its Monte Carlo floor; a sampler that locks each
# chain into one representation puts the fraction at 0 or 1 and the sd near
# 0.5.
#
# Secondary: an acceptance-rate proxy that needs no engine hook. At one tree
# with no thinning one kept draw is one sweep, so the tree structure differs
# from the previous draw's exactly when a structural move was accepted.
#
# Two arms on matched seeds, both reachable at runtime through
# proposal.probs, so the contrast costs a grid of fits and nothing else:
#
#   default        birth_death 0.5, swap 0.1, change 0.4, birth 0.5
#   birth/death    birth_death 1.0, swap 0.0, change 0.0, birth 0.5
#
# The duplicate-column design is the null control the battery requires
# alongside this cell: two exactly identical predictor columns, where the
# likelihood and prior ratios are both exactly 1, so the pooled fraction on
# the first of the pair must be 1/2 with a non-zero switch count in every
# chain. If it is not, the harness is measuring itself.
#
# Usage: Rscript P2-confounded-step.R [outputDir] [quick]

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

nReplicates <- if (quick) 2L else 5L
nChains <- 8L
nBurn <- 1000L
nSamples <- if (quick) 500L else 2000L

arms <- list(
  default = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  birthdeath = c(birth_death = 1, swap = 0, change = 0, birth = 0.5)
)

designs <- list(
  confounded = surfacesConfoundedStep,
  duplicate = surfacesDuplicateColumnNull
)

surfacesUptime("uptime before")

rows <- list()
for (designName in names(designs)) {
  for (replicate in seq_len(nReplicates)) {
    set.seed(surfacesDataSeed("P2", designName, replicate))
    data <- designs[[designName]]()
    for (armName in names(arms)) {
      startedAt <- proc.time()
      fit <- bart2(
        data$x,
        data$y,
        n.trees = 1L,
        n.chains = nChains,
        n.burn = nBurn,
        n.samples = nSamples,
        n.thin = 1L,
        n.threads = 1L,
        keepTrees = TRUE,
        verbose = FALSE,
        seed = surfacesSamplerSeed(replicate),
        proposal.probs = arms[[armName]]
      )
      elapsed <- (proc.time() - startedAt)[["elapsed"]]
      trees <- extract(fit, type = "trees")
      roots <- surfacesRootVariable(trees)

      # Per chain: the share of draws rooted on each column, the share
      # rooted on the first of the confounded (or duplicated) pair among
      # draws rooted on either, and the number of times the root variable
      # moves. The pair is (x1, x3) in the confounded design and (x1, x2) in
      # the duplicate-column null.
      pair <- if (designName == "confounded") c(1L, 3L) else c(1L, 2L)
      perChain <- lapply(sort(unique(roots$chain)), function(ch) {
        v <- roots$var[roots$chain == ch]
        onPair <- v %in% pair
        data.frame(
          chain = ch,
          fractionFirst = mean(v == pair[1L]),
          fractionSecond = mean(v == pair[2L]),
          fractionStump = mean(v == -1L),
          shareFirstGivenPair = if (any(onPair)) {
            mean(v[onPair] == pair[1L])
          } else {
            NA_real_
          },
          switches = sum(v[-1L] != v[-length(v)]),
          distinctRoots = length(unique(v))
        )
      })
      perChain <- do.call(rbind, perChain)

      rows[[length(rows) + 1L]] <- data.frame(
        design = designName,
        arm = armName,
        replicate = replicate,
        betweenChainSd = sd(perChain$fractionFirst),
        pooledFirst = mean(perChain$fractionFirst),
        pooledShareGivenPair = mean(perChain$shareFirstGivenPair, na.rm = TRUE),
        minChainFirst = min(perChain$fractionFirst),
        maxChainFirst = max(perChain$fractionFirst),
        minSwitches = min(perChain$switches),
        meanSwitches = mean(perChain$switches),
        meanDistinctRoots = mean(perChain$distinctRoots),
        acceptanceProxy = surfacesStructureChangeRate(trees),
        wall = elapsed,
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "%-11s %-11s rep %d  sd(p1) %.3f  pooled p1 %.3f  accept %.4f\n",
        designName,
        armName,
        replicate,
        rows[[length(rows)]]$betweenChainSd,
        rows[[length(rows)]]$pooledFirst,
        rows[[length(rows)]]$acceptanceProxy
      ))
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader("P2 confounded step: mean over seeds (min-max)")
cat(sprintf(
  "%-11s %-11s %-20s %-20s %-22s %-16s %s\n",
  "design",
  "arm",
  "between-chain sd",
  "pooled p(root x1)",
  "acceptance proxy",
  "root switches",
  "min switches"
))
for (designName in names(designs)) {
  for (armName in names(arms)) {
    keep <- results$design == designName & results$arm == armName
    cat(sprintf(
      "%-11s %-11s %-20s %-20s %-22s %-16s %d\n",
      designName,
      armName,
      surfacesRange(results$betweenChainSd[keep]),
      surfacesRange(results$pooledFirst[keep]),
      surfacesRange(results$acceptanceProxy[keep], digits = 4L),
      surfacesRange(results$meanSwitches[keep], digits = 1L),
      min(results$minSwitches[keep])
    ))
  }
}

surfacesHeader("published reference")
cat("Pratola 2016 sec 2.3: acceptance rate of tree moves 0 at m = 1\n")
cat("symmetry oracle: pooled share on x1 among {x1, x3} draws = 0.5\n")
cat(sprintf(
  "measured share on x1 given the pair, confounded design: %s (default), %s (birth/death)\n",
  surfacesRange(
    results$pooledShareGivenPair[
      results$design == "confounded" & results$arm == "default"
    ]
  ),
  surfacesRange(
    results$pooledShareGivenPair[
      results$design == "confounded" & results$arm == "birthdeath"
    ]
  )
))

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    settings = list(
      nReplicates = nReplicates,
      nChains = nChains,
      nBurn = nBurn,
      nSamples = nSamples,
      nTrees = 1L,
      arms = arms
    )
  ),
  outputDir,
  "P2-confounded-step"
)
