#!/usr/bin/env Rscript

# Measures maxLeafSize_ (src/bartcore/model.hpp), the GP leaf's member-count
# ceiling, above which a leaf evaluation abandons the Gaussian process and
# scores as a constant leaf. Its default is 256, set from the gp() argument
# max.leaf.size (R/model.R), so unlike the other engine constants in the
# audit this one is already settable and the sweep is a plain argument sweep.
# The GP leaf draw is O(m^3) in the leaf's member count, so the ceiling is a
# time bound, and dec-B110 asks how much of a fit it silently disables.
#
# FIT TIME is msec per iteration. FALLBACK SHARE is the share of terminal
# nodes across the kept draws holding more than max.leaf.size members - the
# leaf evaluations that took the constant-leaf path - with the share of
# TRAINING ROWS sitting in such a leaf beside it, which is the quantity a user
# cares about: a small share of oversized leaves can still cover most of the
# data. The value BINDS if a fit at the recommended 10 to 25 trees falls back
# on a large share of its rows at 256 and stops doing so at 512 or 1024
# without an unacceptable time cost.
#
# The full grid is deliberately small (n = 1200, 20 trees, 60 iterations):
# the 1024 arm is cubic in a leaf that can hold most of the design, so a
# training set of the size the other sweeps use would not finish.
#
# Usage: Rscript benchmarks/R/constant-gp-max-leaf-size.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the ceiling is revisited. Run it on a quiet machine.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

n <- if (quick) 1000L else 1200L
nTrees <- if (quick) 10L else 20L
nSamples <- if (quick) 25L else 40L
nBurn <- if (quick) 25L else 20L
sizes <- if (quick) c(128L, 512L) else c(128L, 256L, 512L, 1024L)

set.seed(19)
x1 <- runif(n)
x2 <- runif(n, -1, 1)
x3 <- runif(n)
y <- sin(4 * pi * x1) + x2^2 + rnorm(n, 0, 0.2)
df <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)

cat(sprintf(
  "n = %d, %d trees, %d kept draws after %d burn\n",
  n,
  nTrees,
  nSamples,
  nBurn
))
cat(sprintf(
  "%14s %12s %14s %14s %10s\n",
  "max.leaf.size",
  "msec/iter",
  "leaf_fallback",
  "row_fallback",
  "rmse"
))
for (maxLeafSize in sizes) {
  start <- Sys.time()
  fit <- bart(
    y ~ x1 + x2 + x3,
    df,
    node.prior = dbartsPriors$gp(c("x1", "x2"), max.leaf.size = maxLeafSize),
    n.trees = nTrees,
    n.samples = nSamples,
    n.burn = nBurn,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    seed = 7L,
    verbose = FALSE
  )
  elapsed <- as.numeric(Sys.time() - start, units = "secs")

  trees <- extract(fit, type = "trees")
  leafSizes <- trees$n[trees$var < 0L]
  oversized <- leafSizes > maxLeafSize
  fitted <- apply(extract(fit, type = "ev"), 2L, mean)

  cat(sprintf(
    "%14d %12.3f %14.4f %14.4f %10.4f\n",
    maxLeafSize,
    elapsed / (nSamples + nBurn) * 1000,
    mean(oversized),
    sum(leafSizes[oversized]) / sum(leafSizes),
    sqrt(mean((fitted - y)^2))
  ))
}
