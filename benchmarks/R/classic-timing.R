#!/usr/bin/env Rscript

# Wall-time comparison of this release against dbarts 0.9-34, the timing
# half of benchmarks/R/classic-compare.R (which compares posteriors and
# times nothing). Measurement, not a gate: bench-sampler.R is the
# zero-regression gate and compares this tree against itself, while this
# script compares two INSTALLED releases and so cannot run in CI.
#
# One invocation fits ONE design cell under whichever dbarts R_LIBS points
# at, and prints the wall-clock seconds of the fit call alone. A fresh
# process per fit is the point: only one dbarts ever loads, and the two
# releases cannot be linked into the same session anyway.
#
#   R_LIBS=<lib-0.9-34> Rscript benchmarks/R/classic-timing.R a
#   R_LIBS=<lib-1.0-0>  Rscript benchmarks/R/classic-timing.R a
#
# The four cells, Friedman data at ten predictors, 200 trees, 1000 draws
# after 500 discarded:
#
#   a  n = 1000,  one chain              c  n = 1000, four chains, four threads
#   b  n = 10000, one chain              d  n = 1000, probit, one chain
#
# Alternate which library runs first on each repetition so drift over the
# run falls on both alike, and take medians over the repetitions; the
# driver loop and the recorded numbers are in benchmarks/README.md and
# docs/plans/classic-compare.md, "Wall time".

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: classic-timing.R <a|b|c|d>")
}
cell <- args[[1L]]

# 0.9-34 spells the BayesTree-style door 'bart'; 1.0-0 spells it 'bartBT'
# and forwards a BayesTree-shaped 'bart' call to it with a once-per-session
# warning. Naming the successor where it exists keeps that forwarding path,
# and its warning, out of the timed call.
bartFn <- if (
  exists("bartBT", where = asNamespace("dbarts"), inherits = FALSE)
) {
  dbarts::bartBT
} else {
  dbarts::bart
}

# 0.9-x's tree-move mixture, pinned on both sides (1.0-0 defaults to
# birth_death 0.6 / swap 0 / change 0.4 / birth 0.5 instead). Every other
# default that moved between the releases is pinned in the call below, the
# same pinning classic-compare.R uses.
proposalProbs <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)

friedman <- function(x) {
  10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
}

p <- 10L
specs <- list(
  a = list(n = 1000L, nchain = 1L, nthread = 1L, binary = FALSE, seed = 71001L),
  b = list(n = 1e4L, nchain = 1L, nthread = 1L, binary = FALSE, seed = 71002L),
  c = list(n = 1000L, nchain = 4L, nthread = 4L, binary = FALSE, seed = 71003L),
  d = list(n = 1000L, nchain = 1L, nthread = 1L, binary = TRUE, seed = 71004L)
)
spec <- specs[[cell]]
if (is.null(spec)) {
  stop("unknown cell '", cell, "'; expected one of a, b, c, d")
}

# data are fixed per cell - same seed every repetition and every library -
# so only the fit itself varies across timings. No cell carries a factor
# predictor, so this door's indicator expansion is inert on both sides.
set.seed(spec$seed)
x <- matrix(runif(spec$n * p), spec$n)
y <- if (spec$binary) {
  rbinom(spec$n, 1L, pnorm(scale(friedman(x))))
} else {
  friedman(x) + rnorm(spec$n)
}

# the forwarding and argument-rename warnings are not part of what is being
# timed and would otherwise land mid-call
muffleBenign <- function(w) {
  msg <- conditionMessage(w)
  if (
    grepl("deprecated", msg) ||
      grepl("bartBT", msg) ||
      grepl("proposalprobs", msg)
  ) {
    invokeRestart("muffleWarning")
  }
}

t0 <- proc.time()[["elapsed"]]
invisible(withCallingHandlers(
  bartFn(
    x.train = x,
    y.train = y,
    sigest = NA_real_,
    sigdf = 3.0,
    sigquant = 0.90,
    k = 2.0,
    power = 2.0,
    base = 0.95,
    binaryOffset = 0.0,
    ntree = 200L,
    ndpost = 1000L,
    nskip = 500L,
    keepevery = 1L,
    keeptrainfits = TRUE,
    usequants = FALSE,
    numcut = 100L,
    verbose = FALSE,
    nchain = spec$nchain,
    nthread = spec$nthread,
    combinechains = TRUE,
    keeptrees = FALSE,
    keepcall = FALSE,
    proposalprobs = proposalProbs
  ),
  warning = muffleBenign
))
elapsed <- proc.time()[["elapsed"]] - t0

cat(sprintf(
  "ELAPSED %.4f VERSION %s CELL %s\n",
  elapsed,
  as.character(utils::packageVersion("dbarts")),
  cell
))
