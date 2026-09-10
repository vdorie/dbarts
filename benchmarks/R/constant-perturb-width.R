#!/usr/bin/env Rscript

# Measures perturbWidth (src/bartcore/moves.hpp): the half-width, in grid
# positions, of the window the perturb move draws a cut displacement from.
# It is a COMPILE-TIME constant with no knob, so a width arm needs its own
# build. This script measures ONE width - the width of the build it is run
# against - and takes the number only as a label for the table; run it once
# per build and read the rows together:
#
#   cp -r <worktree> /tmp/perturb-w2 && edit perturbWidth in
#     /tmp/perturb-w2/src/bartcore/moves.hpp
#   R CMD INSTALL --preclean -l /tmp/lib-w2 /tmp/perturb-w2
#   R_LIBS=/tmp/lib-w2 Rscript benchmarks/R/constant-perturb-width.R width=2
#
# Never edit the width in the worktree itself. Two things are reported, both
# on a perturb-dominant mixture (birth_death 0.1, change 0.1, perturb 0.8) so
# that most accepted structure moves are perturbs: ACCEPTANCE, the share of
# (tree, adjacent kept draw) pairs whose cut vector moved, and ESS per kept
# draw for sigma and for the fitted value at individual training points (the
# median over a sample of them - the fit MEAN mixes too fast to separate the
# arms). A wider window proposes bigger displacements, so acceptance must
# fall; the width binds only if the mixing it buys rises faster than
# acceptance falls.
#
# Usage: Rscript benchmarks/R/constant-perturb-width.R [quick] [width=<w>]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the width is revisited.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
widthArg <- grep("^width=", args, value = TRUE)
width <- if (length(widthArg) > 0L) {
  as.integer(sub("^width=", "", widthArg[[1L]]))
} else {
  NA_integer_
}

nSamples <- if (quick) 200L else 2000L
nBurn <- if (quick) 100L else 500L
nTrees <- if (quick) 20L else 50L
sizes <- if (quick) c(500L) else c(500L, 2000L)

proposalProbs <- c(
  birth_death = 0.10,
  swap = 0,
  change = 0.10,
  perturb = 0.80,
  birth = 0.5
)

# initial-positive-sequence effective sample size, the standard
# autocorrelation-time estimator; returns draws when the chain is white
effectiveSize <- function(z) {
  m <- length(z)
  if (stats::var(z) <= 0) {
    return(NA_real_)
  }
  rho <- stats::acf(z, lag.max = min(m - 1L, 500L), plot = FALSE)$acf[-1L]
  total <- 0
  for (i in seq(1L, length(rho) - 1L, by = 2L)) {
    pair <- rho[[i]] + rho[[i + 1L]]
    if (pair <= 0) {
      break
    }
    total <- total + pair
  }
  m / (1 + 2 * total)
}

cutChangeRate <- function(trees) {
  signatures <- tapply(
    seq_len(nrow(trees)),
    list(trees$sample, trees$tree),
    function(rows) {
      interior <- trees$var[rows] > 0L
      paste0(
        trees$var[rows],
        ":",
        ifelse(interior, format(trees$value[rows], digits = 12), ""),
        collapse = ","
      )
    }
  )
  moved <- signatures[-1L, , drop = FALSE] !=
    signatures[-nrow(signatures), , drop = FALSE]
  mean(moved)
}

cat(sprintf(
  "perturbWidth = %s (build label), %d trees, %d kept draws after %d burn\n",
  ifelse(is.na(width), "unlabelled", width),
  nTrees,
  nSamples,
  nBurn
))
cat(sprintf(
  "%6s %8s %12s %12s %14s %14s\n",
  "width",
  "n",
  "msec/iter",
  "cut_accept",
  "ESS_sigma",
  "ESS_point"
))
for (n in sizes) {
  set.seed(5)
  x <- matrix(runif(n * 5L), n, 5L)
  colnames(x) <- paste0("x", 1:5)
  y <- 10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L] +
    rnorm(n)

  start <- proc.time()[["elapsed"]]
  fit <- bart(
    x,
    y,
    n.trees = nTrees,
    n.samples = nSamples,
    n.burn = nBurn,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    proposal.probs = proposalProbs,
    seed = 7L,
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - start

  trees <- extract(fit, type = "trees")
  sigma <- as.vector(fit$sigma)
  sigma <- sigma[(length(sigma) - nSamples + 1L):length(sigma)]
  # individual fitted values, not their mean: the mean over training points
  # averages the tree-structure noise away and saturates the estimator
  ev <- extract(fit, type = "ev")
  set.seed(53)
  probes <- sample.int(ncol(ev), min(20L, ncol(ev)))
  pointEss <- vapply(probes, function(j) effectiveSize(ev[, j]), numeric(1L))

  cat(sprintf(
    "%6s %8d %12.3f %12.4f %14.1f %14.1f\n",
    ifelse(is.na(width), "?", width),
    n,
    elapsed / (nSamples + nBurn) * 1000,
    cutChangeRate(trees),
    effectiveSize(sigma),
    median(pointEss, na.rm = TRUE)
  ))
}
