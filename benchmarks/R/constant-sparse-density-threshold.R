#!/usr/bin/env Rscript

# Measures sparseDensityThreshold (src/bartcore/data.hpp): the nonzero
# fraction at or below which a CSC-built column keeps rank-bitmap hot storage
# (a bitmap, a word-rank index and the packed nonzero codes) instead of
# densifying to one xint_t per row. The value is 0.2, stated without
# attribution.
#
# Two things are reported per density, on a design whose columns all sit on
# the same side of the threshold. MEMORY is the two layouts' analytic size
# (dense 2 bytes per row per column; sparse 1 bit plus a 4-byte word rank per
# 64 rows, plus 2 bytes per stored nonzero) beside the resident-set growth
# across sampler construction, each density in a FRESH R process because an
# allocator that has already grown the heap reports nothing on a second
# measurement. The resident figure is an UPPER BOUND on the store and does not
# isolate it - everything else the sampler allocates costs the same in either
# layout and swamps it - so the memory verdict rests on the model. TIME is
# msec per sweep: the gather cost of the rank decode against a dense read, and
# it is where the threshold shows, the two rows either side of 0.2 differing
# by the whole cost of the rank machinery. The threshold BINDS if that
# crossover sits far from 0.2.
#
# Usage: Rscript benchmarks/R/constant-sparse-density-threshold.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the threshold is revisited. Run it on a quiet machine.

suppressPackageStartupMessages(library(dbarts))
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("this sweep needs the Matrix package")
}

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
childArg <- grep("^child=", args, value = TRUE)

n <- if (quick) 50000L else 100000L
p <- if (quick) 40L else 100L
nTrees <- if (quick) 25L else 50L
nIter <- if (quick) 10L else 25L
nRounds <- if (quick) 2L else 4L
densities <- if (quick) {
  c(0.10, 0.19, 0.21, 0.40)
} else {
  c(0.05, 0.10, 0.15, 0.19, 0.21, 0.25, 0.35, 0.50)
}

residentKB <- function() {
  as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE))
}

measure <- function(density) {
  set.seed(23)
  # exactly round(n * density) DISTINCT rows per column: sampling (i, j) pairs
  # with replacement would let sparseMatrix sum duplicates, putting the
  # realized density below the requested one and moving the column across the
  # threshold without the label saying so
  perColumn <- round(n * density)
  rows <- unlist(lapply(seq_len(p), function(j) sample.int(n, perColumn)))
  cols <- rep(seq_len(p), each = perColumn)
  x <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 0.5 + runif(length(rows)),
    dims = c(n, p)
  )
  colnames(x) <- paste0("x", seq_len(p))
  y <- as.vector(2 * x[, 1L] - 1.5 * x[, 2L]) + rnorm(n, 0, 0.3)
  storedNonzeros <- length(x@x)

  control <- dbartsControl(
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = nIter,
    n.burn = 0L,
    updateState = FALSE,
    verbose = FALSE
  )
  invisible(gc(FALSE))
  before <- residentKB()
  # read the resident set BEFORE the first sweep: run() allocates the index
  # buffers and the per-tree fits, which cost the same in either layout and
  # would swamp the difference the column store makes
  sampler <- suppressWarnings(dbarts(x, y, control = control))
  after <- residentKB()
  invisible(sampler$run(nIter, 0L))

  best <- Inf
  for (round in seq_len(nRounds)) {
    start <- Sys.time()
    invisible(sampler$run(0L, nIter))
    best <- min(best, as.numeric(Sys.time() - start, units = "secs"))
  }
  c(
    density = length(x@x) / (as.numeric(n) * p),
    rssMB = (after - before) / 1024,
    denseMB = (2 * as.numeric(n) * p) / 1e6,
    sparseMB = (as.numeric(p) *
      (n / 8 + ceiling(n / 64) * 4) +
      2 * storedNonzeros) /
      1e6,
    msec = best / nIter * 1000
  )
}

if (length(childArg) > 0L) {
  values <- measure(as.numeric(sub("^child=", "", childArg[[1L]])))
  cat("CHILD", paste0(values, collapse = " "), "\n")
} else {
  scriptPath <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
  )
  cat(sprintf(
    "n = %d, %d columns, %d trees, %d iterations, %d rounds, min taken\n",
    n,
    p,
    nTrees,
    nIter,
    nRounds
  ))
  cat(sprintf(
    "%9s %10s %12s %12s %12s %12s\n",
    "density",
    "storage",
    "rss_MB",
    "dense_MB",
    "sparse_MB",
    "msec/iter"
  ))
  for (density in densities) {
    childArgs <- c(scriptPath, paste0("child=", density))
    if (quick) {
      childArgs <- c(childArgs, "quick")
    }
    output <- system2(
      file.path(R.home("bin"), "Rscript"),
      childArgs,
      stdout = TRUE,
      stderr = FALSE
    )
    values <- as.numeric(strsplit(
      grep("^CHILD ", output, value = TRUE)[[1L]],
      " +"
    )[[1L]][-1L])
    cat(sprintf(
      "%9.3f %10s %12.1f %12.1f %12.1f %12.3f\n",
      values[[1L]],
      if (values[[1L]] <= 0.2) "sparse" else "dense",
      values[[2L]],
      values[[3L]],
      values[[4L]],
      values[[5L]]
    ))
  }
}
