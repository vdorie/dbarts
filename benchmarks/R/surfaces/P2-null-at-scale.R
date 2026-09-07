#!/usr/bin/env Rscript

# P2's null control, re-run at production tree counts. P2-confounded-step.R
# runs the duplicate-column null (surfacesDuplicateColumnNull, two exactly
# identical columns x1/x2 and one irrelevant column x3) at m = 1, where the
# swap-carrying mixture switches representation freely (78.5 switches per
# chain) and the shipped no-swap move set leaves a minority of chains stuck
# at an x3 root for the whole chain. A single tree is the setting where a
# stuck root is visible at all; an ensemble can self-average a stuck tree's
# label away without any chain ever "mixing" it. This cell asks whether the
# x3-root representation survives at production tree counts or is washed out
# by the ensemble.
#
# Same duplicate-column design, same matched-seed idiom, two of P2's three
# move-set arms (swap and no-swap only - birth/death is P2's own contrast,
# not part of the tree-count question), n.trees in {50, 200}.
#
# Readout is tree-level, not chain-level: for every (arm, n.trees, replicate,
# chain), the share of (draw, tree) roots landing on x3 (the irrelevant
# column), on x1 and on x2, and the count of trees (out of n.trees) whose
# root sits on x3 for every kept draw of that chain - the direct multi-tree
# analogue of P2's single-tree "chain never leaves x3" finding.
#
# Usage: Rscript P2-null-at-scale.R [outputDir] [quick]

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
treeCounts <- c(50L, 200L)

arms <- list(
  swap = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  noswap = c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5)
)

surfacesUptime("uptime before")

rows <- list()
for (armName in names(arms)) {
  for (nTrees in treeCounts) {
    for (replicate in seq_len(nReplicates)) {
      set.seed(surfacesDataSeed("P2", "duplicate", replicate))
      data <- surfacesDuplicateColumnNull()
      startedAt <- proc.time()
      fit <- bart2(
        data$x,
        data$y,
        n.trees = nTrees,
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

      # Per chain: the share of (draw, tree) roots on x3, x1 and x2, and the
      # number of trees (of n.trees) whose root is on x3 in EVERY kept draw
      # of this chain - stuck for the chain's entire run, not just typical.
      perChain <- lapply(sort(unique(roots$chain)), function(ch) {
        v <- roots[roots$chain == ch, ]
        stuckOnX3 <- sum(tapply(v$var, v$tree, function(x) all(x == 3L)))
        data.frame(
          chain = ch,
          x3Share = mean(v$var == 3L),
          x1Share = mean(v$var == 1L),
          x2Share = mean(v$var == 2L),
          stuckOnX3 = stuckOnX3
        )
      })
      perChain <- do.call(rbind, perChain)

      rows[[length(rows) + 1L]] <- data.frame(
        arm = armName,
        nTrees = nTrees,
        replicate = replicate,
        chain = perChain$chain,
        x3Share = perChain$x3Share,
        x1Share = perChain$x1Share,
        x2Share = perChain$x2Share,
        stuckOnX3 = perChain$stuckOnX3,
        wall = elapsed,
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "%-8s trees %3d  rep %d  x3 share %.3f (%.3f-%.3f)  stuck %d  %.0fs\n",
        armName,
        nTrees,
        replicate,
        mean(perChain$x3Share),
        min(perChain$x3Share),
        max(perChain$x3Share),
        sum(perChain$stuckOnX3),
        elapsed
      ))
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader("P2 null control at production tree counts: per (arm, n.trees)")
cat(sprintf(
  "%-8s %-7s %-22s %-22s %-16s\n",
  "arm",
  "trees",
  "mean x3 share",
  "between-chain sd",
  "stuck-on-x3 trees"
))
for (armName in names(arms)) {
  for (nTrees in treeCounts) {
    keep <- results$arm == armName & results$nTrees == nTrees
    sdPerRep <- tapply(results$x3Share[keep], results$replicate[keep], sd)
    cat(sprintf(
      "%-8s %-7d %-22s %-22s %-16d\n",
      armName,
      nTrees,
      surfacesRange(results$x3Share[keep]),
      surfacesRange(sdPerRep),
      sum(results$stuckOnX3[keep])
    ))
  }
}

surfacesHeader(
  "context: mean x1 / x2 share (the duplicated, signal-carrying pair)"
)
cat(sprintf(
  "%-8s %-7s %-22s %-22s\n",
  "arm",
  "trees",
  "mean x1 share",
  "mean x2 share"
))
for (armName in names(arms)) {
  for (nTrees in treeCounts) {
    keep <- results$arm == armName & results$nTrees == nTrees
    cat(sprintf(
      "%-8s %-7d %-22s %-22s\n",
      armName,
      nTrees,
      surfacesRange(results$x1Share[keep]),
      surfacesRange(results$x2Share[keep])
    ))
  }
}

surfacesHeader("published reference")
cat("P2-confounded-step.R, m = 1: swap 0 stuck-on-x3 chains of 40, ")
cat("no-swap 5 of 40 (docs/design/benchmark-surfaces.md sec 10.1)\n")

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    settings = list(
      nReplicates = nReplicates,
      nChains = nChains,
      nBurn = nBurn,
      nSamples = nSamples,
      treeCounts = treeCounts,
      arms = arms
    )
  ),
  outputDir,
  "P2-null-at-scale"
)
