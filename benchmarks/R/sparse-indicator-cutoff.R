#!/usr/bin/env Rscript

# Measures the memory (and, incidentally, construction time) a factor's
# "indicators" dummy expansion costs dense against sparse, across a sweep of
# level counts, to justify sparseIndicatorLevelCutoff (R/utility.R): the
# level count past which makeModelMatrixFromDataFrame builds a wide factor's
# indicator block as a dgCMatrix automatically (dec-B100).
#
# One-hot storage costs O(n) regardless of level count (one stored entry per
# row, at most), against the dense block's O(n * K); this script reports the
# ratio at each K in the sweep and the level count where it first exceeds a
# few tens. Construction TIME is reported too, but is not the deciding
# factor: a real C builder (makeModelMatrixFromDataFrame's dense path) stays
# faster to build than this script's own from-scratch sparse assembly
# throughout the range swept, so the cutoff is a memory choice, not a speed
# one - the comment in R/utility.R states the same conclusion this script
# reproduces.
#
# Usage: Rscript benchmarks/R/sparse-indicator-cutoff.R [n]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cutoff is revisited.

suppressPackageStartupMessages(library(dbarts))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1L) as.integer(args[[1L]]) else 5000L

levelsToSweep <- c(5L, 10L, 25L, 50L, 75L, 100L, 150L, 200L, 300L, 500L)

cat(sprintf(
  "%8s %12s %12s %12s %10s\n",
  "levels",
  "dense_MB",
  "sparse_MB",
  "ratio",
  "dense/sparse_ms"
))
for (K in levelsToSweep) {
  set.seed(1)
  f <- factor(sample.int(K, n, replace = TRUE))
  df <- data.frame(f = f)

  denseTime <- system.time(
    for (r in 1:5) {
      mmDense <- .Call(
        dbarts:::C_dbarts_makeModelMatrixFromDataFrame,
        df,
        TRUE
      )
    }
  )[["elapsed"]] /
    5 *
    1000

  sparseTime <- system.time(
    for (r in 1:5) {
      mmSparse <- dbarts:::sparseFactorIndicatorSlices(f, "f", TRUE)
    }
  )[["elapsed"]] /
    5 *
    1000

  denseMB <- as.numeric(object.size(mmDense)) / 1e6
  # the sparse block's own footprint: one stored entry (index + value) per
  # present-level row, the shape assembleMixedMatrix's dgCMatrix ends up with
  nnz <- sum(lengths(mmSparse$i))
  sparseMB <- (nnz * 12) / 1e6 # 4-byte index + 8-byte double per entry

  cat(sprintf(
    "%8d %12.4f %12.4f %12.3f %10.4f/%0.4f\n",
    K,
    denseMB,
    sparseMB,
    denseMB / sparseMB,
    denseTime,
    sparseTime
  ))
}
