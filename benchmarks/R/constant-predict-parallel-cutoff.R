#!/usr/bin/env Rscript

# Calibrates predictParallelCutoff (src/bartcore/sampler.hpp), the traversal
# count below which an out-of-sample replay runs inline on the caller's
# thread. A traversal is one (row, tree, slab) descent; the header derives
# the 1e7 default from a ~2.6 ns cost per traversal, so the number is an
# arithmetic estimate rather than a measurement. dec-B93 asks for it
# calibrated.
#
# The engine exposes the seam the tests use - bartcore_setPredictParallelCutoff,
# which replaces the constant and returns the value it replaced - so both arms
# run the SAME build: cutoff 1 forces the fan-out at every size, and a cutoff
# above every traversal count in the sweep forces the inline path. The
# crossover is where the threaded column first beats the serial one by more
# than measurement noise; that traversal count, not 1e7, is what the constant
# should read.
#
# Usage: Rscript benchmarks/R/constant-predict-parallel-cutoff.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cutoff is revisited. Run it on a quiet machine.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nTrain <- 1000L
nTrees <- 75L
nSamples <- if (quick) 50L else 200L
nRounds <- if (quick) 3L else 7L
nThreads <- 4L
rowCounts <- if (quick) {
  c(30L, 300L, 3000L)
} else {
  c(10L, 100L, 1000L, 10000L, 100000L)
}
lowRowCounts <- c(10L, 30L, 100L, 300L, 1000L, 3000L, 10000L)

set.seed(17)
p <- 10L
x <- matrix(runif(nTrain * p), nTrain, p)
colnames(x) <- paste0("x", seq_len(p))
y <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 10 * x[, 4L] + rnorm(nTrain)

makeFit <- function(nTrees, nSamples) {
  bart(
    x,
    y,
    n.trees = nTrees,
    n.samples = nSamples,
    n.burn = 100L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    seed = 7L,
    verbose = FALSE
  )
}

timePredict <- function(fit, newdata) {
  invisible(predict(fit, newdata, n.threads = nThreads))
  # Sys.time, not proc.time: the small end of the sweep runs in tens of
  # microseconds and proc.time's elapsed field quantizes to a millisecond
  best <- Inf
  for (round in seq_len(nRounds)) {
    start <- Sys.time()
    invisible(predict(fit, newdata, n.threads = nThreads))
    best <- min(best, as.numeric(Sys.time() - start, units = "secs"))
  }
  best * 1000
}

serialCutoff <- .Machine$integer.max
previous <- .Call(dbarts:::C_dbarts_bartcore_setPredictParallelCutoff, 1L)

cat(sprintf(
  "%d threads, %d rounds, min taken; the crossover is where speedup passes 1\n",
  nThreads,
  nRounds
))
cat(sprintf(
  "%7s %8s %9s %14s %14s %14s %10s\n",
  "trees",
  "draws",
  "rows",
  "traversals",
  "serial_ms",
  "threaded_ms",
  "speedup"
))
sweep <- function(fitTrees, fitDraws, rowCounts) {
  fit <- makeFit(fitTrees, fitDraws)
  for (rows in rowCounts) {
    set.seed(31)
    newdata <- matrix(runif(rows * p), rows, p)
    colnames(newdata) <- colnames(x)
    traversals <- as.numeric(rows) * fitTrees * fitDraws

    invisible(.Call(
      dbarts:::C_dbarts_bartcore_setPredictParallelCutoff,
      serialCutoff
    ))
    serial <- timePredict(fit, newdata)
    invisible(.Call(dbarts:::C_dbarts_bartcore_setPredictParallelCutoff, 1L))
    threaded <- timePredict(fit, newdata)

    cat(sprintf(
      "%7d %8d %9d %14.3e %14.3f %14.3f %10.3f\n",
      fitTrees,
      fitDraws,
      rows,
      traversals,
      serial,
      threaded,
      serial / threaded
    ))
  }
}
# the low arm: a small fit is the only way to reach traversal counts under
# 1e5, where the spawn and join are a real share of the replay
sweep(10L, if (quick) 10L else 10L, if (quick) c(10L, 100L) else lowRowCounts)
sweep(nTrees, nSamples, rowCounts)

invisible(.Call(dbarts:::C_dbarts_bartcore_setPredictParallelCutoff, previous))
