#!/usr/bin/env Rscript

# Mixing gate for the negative-binomial dispersion r
# (docs/plans/nbinom-log-mean.md, Gates). negbin-exact.R checks the stationary
# law on one tree and n = 50; it cannot see a chain that never leaves its cold
# start. This gate fits the default forest at realistic n and checks that r
# moves, that two chains agree on it, and that the predictive law covers fresh
# counts.
#
# Design: mu = 8 exp(x1), five uniform predictors, default bart() settings
# except n.chains = 2, on two cells where r is identified at these means:
# r0 = 5 at n = 2000 and r0 = 2 at n = 500. r0 = 30 is left out on purpose: at
# these means 30 and 50 are hard to tell apart, the forest absorbs the variance
# difference, and a correct sampler can put most mass on 50.
#
# Pass, per cell and seed:
#   (i)   each chain left the cold start: under half its draws at r = 8;
#   (ii)  the chains agree: split-Rhat on r below 1.05. When every split half
#         is constant (zero within-half variance) Rhat is undefined: all halves
#         at one shared value passes (Rhat taken as 1), anything else fails;
#   (iii) the pooled central 95% set of r contains r0;
#   (iv)  90% predictive coverage of 1000 fresh counts by randomized PIT,
#         u = F(y - 1) + V (F(y) - F(y - 1)) with F the ppd draws' ecdf and V
#         uniform, lies in 0.90 +- 0.04. The randomization removes the
#         over-coverage of equal-tailed integer quantiles at small means.
#
# A failure is a finding; never widen a band to pass. quick and full differ in
# seeds only (full adds two).
#
# Usage: Rscript negbin-mixing.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

seeds <- if (quick) 1L else 1:3
cells <- list(
  list(r0 = 5, n = 2000L),
  list(r0 = 2, n = 500L)
)
numPredictors <- 5L
numTest <- 1000L
coldStart <- 8
rhatMax <- 1.05
coverageTarget <- 0.90
coverageBand <- 0.04

simulate <- function(n, r0) {
  x <- matrix(runif(n * numPredictors), n, numPredictors)
  colnames(x) <- paste0("x", seq_len(numPredictors))
  mu <- 8 * exp(x[, 1L])
  list(x = x, y = rnbinom(n, size = r0, mu = mu))
}

# split-Rhat (BDA3) over the chains' halves; draws is samples x chains
splitRhat <- function(draws) {
  half <- nrow(draws) %/% 2L
  halves <- cbind(draws[seq_len(half), ], draws[half + seq_len(half), ])
  within <- mean(apply(halves, 2L, var))
  if (within == 0) {
    return(if (length(unique(as.vector(halves))) == 1L) 1 else Inf)
  }
  between <- half * var(colMeans(halves))
  sqrt(((half - 1) / half * within + between / half) / within)
}

# randomized PIT of each fresh count against its ppd draws (draws x rows)
randomizedPit <- function(ppd, y) {
  below <- colMeans(ppd < rep(y, each = nrow(ppd)))
  atOrBelow <- colMeans(ppd <= rep(y, each = nrow(ppd)))
  below + runif(length(y)) * (atOrBelow - below)
}

runCell <- function(cell, seed) {
  set.seed(seed)
  train <- simulate(cell$n, cell$r0)
  test <- simulate(numTest, cell$r0)
  fit <- bart(
    train$x,
    train$y,
    test = test$x,
    family = "nbinom",
    n.chains = 2L,
    seed = seed,
    verbose = FALSE
  )
  r <- extract(fit, type = "dispersion", combineChains = FALSE)
  if (is.null(dim(r))) {
    r <- matrix(r, ncol = 2L)
  } else if (nrow(r) == 2L) {
    r <- t(r)
  }
  atCold <- colMeans(r == coldStart)
  rhat <- splitRhat(r)
  band <- quantile(as.vector(r), c(0.025, 0.975), names = FALSE, type = 1L)
  ppd <- extract(fit, type = "ppd", sample = "test", combineChains = TRUE)
  u <- randomizedPit(ppd, test$y)
  coverage <- mean(u >= 0.05 & u <= 0.95)
  pass <- c(
    cold = all(atCold < 0.5),
    rhat = rhat < rhatMax,
    band = band[1L] <= cell$r0 && cell$r0 <= band[2L],
    coverage = abs(coverage - coverageTarget) <= coverageBand
  )
  tab <- vapply(
    seq_len(ncol(r)),
    function(chain) {
      counts <- table(r[, chain])
      paste(names(counts), counts, sep = ":", collapse = " ")
    },
    ""
  )
  cat(sprintf(
    "r0 = %g, n = %d, seed %d\n  r by chain: %s\n",
    cell$r0,
    cell$n,
    seed,
    paste(tab, collapse = " | ")
  ))
  cat(sprintf(
    "  at r = 8 %s%s; split-Rhat %.3f%s; 95%% set [%g, %g]%s; coverage %.3f%s\n",
    paste(sprintf("%.2f", atCold), collapse = "/"),
    if (pass[["cold"]]) "" else " <- FAIL",
    rhat,
    if (pass[["rhat"]]) "" else " <- FAIL",
    band[1L],
    band[2L],
    if (pass[["band"]]) "" else " <- FAIL",
    coverage,
    if (pass[["coverage"]]) "" else " <- FAIL"
  ))
  all(pass)
}

cat("Negative-binomial mixing gate (default forest, two chains):\n")
results <- unlist(lapply(cells, function(cell) {
  vapply(seeds, function(seed) runCell(cell, seed), TRUE)
}))

if (!all(results)) {
  quit(status = 1L)
}
cat("\nOK: the dispersion mixes and the predictive law covers fresh counts\n")
