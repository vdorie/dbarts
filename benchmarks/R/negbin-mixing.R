#!/usr/bin/env Rscript

# Mixing gate for the negative-binomial shape r
# (docs/plans/nbinom-log-mean.md, Gates; docs/plans/check-shape-fixes.md).
# negbin-exact.R checks the stationary law on one tree and n = 50; it cannot
# see a chain that never leaves where it began. This gate fits the default
# forest at realistic n, with the two chains' shapes SET APART at the start,
# and checks that r moves, that the two chains agree, and that the predictive
# law covers fresh counts.
#
# Design: mu = 8 exp(x1), five uniform predictors, default bart() settings
# except n.chains = 2, on two cells whose posterior on r spreads over several
# grid values, so agreement between chains can be measured: r0 = 8 at n = 400
# (mass on 5, 6, 8 and 10) and r0 = 10 at n = 500 (mass on 8 to 20). A cell whose
# posterior is a single value is not a cell: split-Rhat is undefined there and
# the gate reports it as a failure of the cell, not a pass. r0 = 30 is left
# out on purpose: at these means 30 and 50 are hard to tell apart, the forest
# absorbs the variance difference, and a correct sampler can put most mass on
# 50.
#
# Each chain keeps 2000 draws, not bart()'s 500: r is a grid draw that drifts
# with the forest, its split halves of 500 draws differ by more than 1.05 in
# Rhat on a correct sampler (seed 3 of the r0 = 10 cell does), and 2000 is
# where the chains' block means settle together.
#
# The sampler has no setter for the shape; its stored state carries one per
# chain and setState installs it, so the fit is built with samplerOnly, chain 1
# is set to the grid's smallest value (1) and chain 2 to its largest (50), and
# the sampler is run. A sampler that does not move r leaves chain 1 at 1 and
# chain 2 at 50, which fails (i) to (iii).
#
# Pass, per cell and seed:
#   (i)   each chain left its start: under half its draws at its start value;
#   (ii)  the chains agree: split-Rhat on r below 1.05 (undefined, from no
#         spread in any split half, fails);
#   (iii) the pooled central 95% set of r contains r0;
#   (iv)  90% predictive coverage of 1000 fresh counts by randomized PIT,
#         u = F(y - 1) + V (F(y) - F(y - 1)) with F the ppd draws' ecdf and V
#         uniform, lies in 0.90 +- 0.04. The randomization removes the
#         over-coverage of equal-tailed integer quantiles at small means. The
#         ppd is drawn here from the run's own test log means and shapes, as
#         predict(type = "ppd") draws it.
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
  list(r0 = 8, n = 400L),
  list(r0 = 10, n = 500L)
)
# the shape each chain is set to before the run: the grid's ends
startShapes <- c(1, 50)
numPredictors <- 5L
numTest <- 1000L
numSamples <- 2000L
rhatMax <- 1.05
coverageTarget <- 0.90
coverageBand <- 0.04

simulate <- function(n, r0) {
  x <- matrix(runif(n * numPredictors), n, numPredictors)
  colnames(x) <- paste0("x", seq_len(numPredictors))
  mu <- 8 * exp(x[, 1L])
  list(x = x, y = rnbinom(n, size = r0, mu = mu))
}

# split-Rhat (BDA3) over the chains' halves; draws is samples x chains. NA when
# no half varies: the agreement of chains with no spread cannot be measured.
splitRhat <- function(draws) {
  half <- nrow(draws) %/% 2L
  halves <- cbind(draws[seq_len(half), ], draws[half + seq_len(half), ])
  within <- mean(apply(halves, 2L, var))
  if (within == 0) {
    return(NA_real_)
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

# one posterior-predictive count per kept draw and fresh row, drawn at the
# draw's own shape; test is rows x draws x chains of log means, shape is
# draws x chains. Returns (draws * chains) x rows.
posteriorPredictive <- function(test, shape) {
  do.call(
    rbind,
    lapply(seq_len(ncol(shape)), function(chain) {
      mu <- exp(t(test[,, chain]))
      matrix(
        rnbinom(
          length(mu),
          size = rep(shape[, chain], times = ncol(mu)),
          mu = mu
        ),
        nrow(mu)
      )
    })
  )
}

runCell <- function(cell, seed) {
  set.seed(seed)
  train <- simulate(cell$n, cell$r0)
  test <- simulate(numTest, cell$r0)
  sampler <- suppressMessages(bart(
    train$x,
    train$y,
    test = test$x,
    family = "nbinom",
    n.chains = 2L,
    seed = seed,
    verbose = FALSE,
    samplerOnly = TRUE
  ))
  state <- sampler$state
  for (chain in seq_along(startShapes)) {
    state[[chain]]$shape <- startShapes[chain]
  }
  sampler$setState(state)
  if (!identical(as.numeric(sampler$getShape()), startShapes)) {
    stop("the chains' shapes were not set apart")
  }
  run <- sampler$run(sampler$control@n.burn, numSamples)
  r <- matrix(run$shape, ncol = 2L)
  atStart <- vapply(
    seq_along(startShapes),
    function(chain) mean(r[, chain] == startShapes[chain]),
    0
  )
  rhat <- splitRhat(r)
  band <- quantile(as.vector(r), c(0.025, 0.975), names = FALSE, type = 1L)
  ppd <- posteriorPredictive(run$test, r)
  u <- randomizedPit(ppd, test$y)
  coverage <- mean(u >= 0.05 & u <= 0.95)
  pass <- c(
    moved = all(atStart < 0.5),
    rhat = !is.na(rhat) && rhat < rhatMax,
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
    "r0 = %g, n = %d, seed %d, chains set to %s\n  r by chain: %s\n",
    cell$r0,
    cell$n,
    seed,
    paste(startShapes, collapse = "/"),
    paste(tab, collapse = " | ")
  ))
  cat(sprintf(
    "  at start %s%s; split-Rhat %s%s; 95%% set [%g, %g]%s; coverage %.3f%s\n",
    paste(sprintf("%.2f", atStart), collapse = "/"),
    if (pass[["moved"]]) "" else " <- FAIL",
    if (is.na(rhat)) "undefined (no spread)" else sprintf("%.3f", rhat),
    if (pass[["rhat"]]) "" else " <- FAIL",
    band[1L],
    band[2L],
    if (pass[["band"]]) "" else " <- FAIL",
    coverage,
    if (pass[["coverage"]]) "" else " <- FAIL"
  ))
  all(pass)
}

cat("Negative-binomial mixing gate (default forest, two chains set apart):\n")
results <- unlist(lapply(cells, function(cell) {
  vapply(seeds, function(seed) runCell(cell, seed), TRUE)
}))

if (!all(results)) {
  quit(status = 1L)
}
cat(
  "\nOK: the shape moves and the chains agree on it, and the predictive law covers fresh counts\n"
)
