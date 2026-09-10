#!/usr/bin/env Rscript

# Measures testFitParallelCutoff (src/bartcore/chain.hpp): the test-row count
# below which a chain routes its test matrix through the trees on its own
# thread instead of borrowing its share of the thread budget. The value is
# 65536, uncalibrated.
#
# The comparison is a whole sampler sweep with a fixed training set and a
# growing test matrix, at one thread and at four, so the difference between
# the two columns is the test-fit routing and nothing else (a single chain
# has no other use for the extra workers). Below the cutoff both columns run
# the same serial code and must agree to noise; above it the threaded column
# is the only one that fans out. The cutoff BINDS if the threaded column wins
# materially at a test size below 65536, since a caller there is paying the
# serial path for no reason.
#
# Usage: Rscript benchmarks/R/constant-testfit-parallel-cutoff.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cutoff is revisited. Run it on a quiet machine.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nTrain <- 2000L
nTrees <- if (quick) 25L else 75L
nIter <- if (quick) 20L else 100L
nRounds <- if (quick) 2L else 5L
testSizes <- if (quick) {
  c(16384L, 65536L, 262144L)
} else {
  c(8192L, 16384L, 32768L, 65536L, 131072L, 262144L)
}
threadCounts <- c(1L, 4L)

set.seed(13)
p <- 10L
x <- matrix(runif(nTrain * p), nTrain, p)
colnames(x) <- paste0("x", seq_len(p))
y <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 10 * x[, 4L] + rnorm(nTrain)

runOnce <- function(nTest, nThreads) {
  set.seed(29)
  xTest <- matrix(runif(nTest * p), nTest, p)
  colnames(xTest) <- colnames(x)
  control <- dbartsControl(
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = nThreads,
    n.samples = nIter,
    n.burn = 0L,
    updateState = FALSE,
    verbose = FALSE
  )
  sampler <- dbarts(x, y, xTest, control = control)
  invisible(sampler$run(nIter, 0L)) # warm the trees and the pool
  best <- Inf
  for (round in seq_len(nRounds)) {
    start <- proc.time()[["elapsed"]]
    invisible(sampler$run(0L, nIter))
    best <- min(best, proc.time()[["elapsed"]] - start)
  }
  best / nIter * 1000
}

cat(sprintf(
  "n.train = %d, %d trees, %d iterations per round, %d rounds, min taken\n",
  nTrain,
  nTrees,
  nIter,
  nRounds
))
cat(sprintf(
  "%10s %14s %14s %10s %s\n",
  "n.test",
  "serial_ms",
  "threaded_ms",
  "speedup",
  "path"
))
for (nTest in testSizes) {
  serial <- runOnce(nTest, threadCounts[[1L]])
  threaded <- runOnce(nTest, threadCounts[[2L]])
  cat(sprintf(
    "%10d %14.3f %14.3f %10.3f %s\n",
    nTest,
    serial,
    threaded,
    serial / threaded,
    if (nTest >= 65536L) "threaded" else "serial (both)"
  ))
}
