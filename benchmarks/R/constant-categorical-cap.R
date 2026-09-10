#!/usr/bin/env Rscript

# Measures categoricalExhaustiveCap (src/bartcore/scan.hpp): the present-level
# count at or below which the categorical rule scan enumerates every one of
# the 2^(P-1) - 1 balanced partitions, and above which it emits the P - 1
# sorted prefixes instead. The cap is 10, so 511 candidates is the widest
# enumeration a fit ever pays for and P = 11 drops to ten prefixes.
#
# Two things are reported per present-level count P. TIME is the sampler's
# msec per iteration. The grid is deliberately many trees over few rows: a
# proposal's histogram pass is O(node members) and its enumeration is O(2^P)
# below the cap, so small nodes are where the candidate count shows. A
# throwaway fit runs first, because the first fit of a session carries page
# faults that swamp the difference being measured.
# ACCEPTANCE is the share of (tree, adjacent kept draw) pairs whose structure
# signature moved, the coarse proposal-acceptance proxy available without the
# move-census build, alongside the share of interior nodes splitting on the
# factor: the prefix family is a subset of the partition family, so if the cap
# bound, the columns past it would show fewer accepted moves and less use of
# the factor, not merely a cheaper sweep.
#
# Usage: Rscript benchmarks/R/constant-categorical-cap.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cap is revisited.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

n <- if (quick) 1000L else 2000L
nTrees <- if (quick) 50L else 200L
nSamples <- if (quick) 40L else 200L
nBurn <- if (quick) 20L else 100L
levelCounts <- if (quick) c(8L, 12L) else c(8L, 10L, 12L, 14L)

structureChangeRate <- function(trees) {
  signatures <- tapply(
    trees$var,
    list(trees$sample, trees$tree),
    function(v) paste0(v, collapse = ",")
  )
  moved <- signatures[-1L, , drop = FALSE] !=
    signatures[-nrow(signatures), , drop = FALSE]
  mean(moved)
}

warmup <- bart(
  y ~ g + x1 + x2,
  data.frame(
    g = factor(sample.int(4L, 200L, replace = TRUE)),
    x1 = runif(200L),
    x2 = runif(200L),
    y = rnorm(200L)
  ),
  n.trees = nTrees,
  n.samples = 20L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE
)
rm(warmup)

cat(sprintf(
  "n = %d, %d trees, %d kept draws after %d burn\n",
  n,
  nTrees,
  nSamples,
  nBurn
))
cat(sprintf(
  "%8s %10s %12s %12s %12s\n",
  "levels",
  "candidates",
  "msec/iter",
  "accept",
  "factor_split"
))
for (P in levelCounts) {
  set.seed(20L + P)
  g <- factor(sample.int(P, n, replace = TRUE))
  x1 <- runif(n)
  x2 <- runif(n)
  effect <- rnorm(P, 0, 1)
  y <- effect[as.integer(g)] + 0.25 * x1 + rnorm(n, 0, 0.5)
  df <- data.frame(g = g, x1 = x1, x2 = x2, y = y)

  elapsed <- system.time(
    fit <- bart(
      y ~ g + x1 + x2,
      df,
      n.trees = nTrees,
      n.samples = nSamples,
      n.burn = nBurn,
      n.chains = 1L,
      n.threads = 1L,
      keepTrees = TRUE,
      seed = 7L,
      verbose = FALSE
    )
  )[["elapsed"]]

  trees <- extract(fit, type = "trees")
  interior <- trees$var > 0L
  candidates <- if (P <= 10L) 2^(P - 1L) - 1L else P - 1L

  cat(sprintf(
    "%8d %10d %12.3f %12.4f %12.4f\n",
    P,
    candidates,
    elapsed / (nSamples + nBurn) * 1000,
    structureChangeRate(trees),
    mean(trees$var[interior] == 1L)
  ))
}
