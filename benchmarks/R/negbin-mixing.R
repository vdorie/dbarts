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
# with the forest, and 500 split halves can differ by more than 1.05 in Rhat on
# a correct sampler.
#
# The sampler has no setter for the shape; its stored state carries one per
# chain and setState installs it, so the fit is built with samplerOnly, chain 1
# is set to the grid's smallest value (1) and chain 2 to its largest (50), and
# the sampler is run. A sampler that does not move r leaves chain 1 at 1 and
# chain 2 at 50, which fails (i), (ii) and (iv), the set [1, 50] holding r0.
#
# Pass, per cell and seed:
#   (i)   each chain left its start: under half its draws at its start value;
#   (ii)  the chains agree: split-Rhat on r below 1.05 (undefined, from no
#         spread in any split half, fails);
#   (iii) the pooled central 99.9% set of r contains r0;
#   (iv)  90% predictive coverage of 1000 fresh counts by randomized PIT,
#         u = F(y - 1) + V (F(y) - F(y - 1)) with F the ppd draws' ecdf and V
#         uniform, lies in [0.84, 0.95]. The randomization removes the
#         over-coverage of equal-tailed integer quantiles at small means. The
#         ppd is drawn here from the run's own test log means and shapes, as
#         predict(type = "ppd") draws it;
#   (v)   full only: the pooled central 95% set contains r0 in all but at most
#         3 of the 6 runs. A 95% set misses the truth in about one dataset in
#         twenty by construction, so it is read over the runs, not in each.
# quick runs seed 1 of each cell and checks (i) to (iv); full runs seeds 1 to 3
# and adds (v).
#
# False failures, measured on the correct sampler. The coverage of a run is set
# by its data far more than by the engine's stream (sd 0.016 over datasets,
# 0.002 over streams), so its limits are the mean 0.897 of 90 fresh datasets
# plus and minus 3.5 sd, rounded. Over 45 fresh datasets a cell (seeds 1 to
# 45, engine seed the data seed) no run failed (i) to (iv): 0 of 90, the
# largest Rhat 1.018, coverage 0.858 to 0.929, and the 95% set missed r0 in
# 5 of 90. Over 12 engine streams on each pinned dataset (seeds 1 to 3, both
# cells) no run failed: 0 of 72, the largest Rhat 1.035. Resampling the
# measured runs into full-mode gates gives a false-failure rate under 1e-4;
# the per-run counts alone bound it at 3% (95% upper limit for 0 of 90).
#
# What it catches (the shape draw mutated in a build of its own; quick, full):
# a shape that never moves (caught, caught); a move applied one sweep in 100
# (caught, caught); one sweep in 20 (missed, caught by Rhat in one run); every
# draw one grid step up (missed, caught by (iii), (iv) and (v)). A shape draw
# off by less than that, or a move applied one sweep in 20 or more often, is
# not caught here; negbin-exact.R holds the stationary law at n = 50.
#
# A failure is a finding.
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
coverageLow <- 0.84
coverageHigh <- 0.95
set95Misses <- 3L

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

# one fit; returns the numbers the checks read. engineSeed seeds the sampler,
# seed the data.
runCell <- function(cell, seed, engineSeed = seed) {
  set.seed(seed)
  train <- simulate(cell$n, cell$r0)
  test <- simulate(numTest, cell$r0)
  sampler <- suppressMessages(bart(
    train$x,
    train$y,
    test = test$x,
    family = "nbinom",
    n.chains = 2L,
    seed = engineSeed,
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
  u <- randomizedPit(posteriorPredictive(run$test, r), test$y)
  tab <- vapply(
    seq_len(ncol(r)),
    function(chain) {
      counts <- table(r[, chain])
      paste(names(counts), counts, sep = ":", collapse = " ")
    },
    ""
  )
  list(
    cell = cell,
    seed = seed,
    atStart = vapply(
      seq_along(startShapes),
      function(chain) mean(r[, chain] == startShapes[chain]),
      0
    ),
    rhat = splitRhat(r),
    wide = quantile(as.vector(r), c(0.0005, 0.9995), names = FALSE, type = 1L),
    set95 = quantile(as.vector(r), c(0.025, 0.975), names = FALSE, type = 1L),
    coverage = mean(u >= 0.05 & u <= 0.95),
    table = tab
  )
}

# the per-run checks; set95 is read in aggregate, not here
judge <- function(res) {
  c(
    moved = all(res$atStart < 0.5),
    rhat = !is.na(res$rhat) && res$rhat < rhatMax,
    wide = res$wide[1L] <= res$cell$r0 && res$cell$r0 <= res$wide[2L],
    coverage = res$coverage >= coverageLow && res$coverage <= coverageHigh
  )
}

report <- function(res, pass) {
  cat(sprintf(
    "r0 = %g, n = %d, seed %d, chains set to %s\n  r by chain: %s\n",
    res$cell$r0,
    res$cell$n,
    res$seed,
    paste(startShapes, collapse = "/"),
    paste(res$table, collapse = " | ")
  ))
  flag <- function(name) if (pass[[name]]) "" else " <- FAIL"
  cat(sprintf(
    "  at start %s%s; split-Rhat %s%s; 99.9%% set [%g, %g]%s; 95%% set [%g, %g]; coverage %.3f%s\n",
    paste(sprintf("%.2f", res$atStart), collapse = "/"),
    flag("moved"),
    if (is.na(res$rhat)) "undefined (no spread)" else sprintf("%.3f", res$rhat),
    flag("rhat"),
    res$wide[1L],
    res$wide[2L],
    flag("wide"),
    res$set95[1L],
    res$set95[2L],
    res$coverage,
    flag("coverage")
  ))
}

if (sys.nframe() == 0L) {
  cat("Negative-binomial mixing gate (default forest, two chains set apart):\n")
  runs <- unlist(
    lapply(cells, function(cell) {
      lapply(seeds, function(seed) runCell(cell, seed))
    }),
    recursive = FALSE
  )
  passes <- lapply(runs, judge)
  for (i in seq_along(runs)) {
    report(runs[[i]], passes[[i]])
  }
  ok <- all(vapply(passes, all, TRUE))
  if (!quick) {
    held <- vapply(
      runs,
      function(res) {
        res$set95[1L] <= res$cell$r0 && res$cell$r0 <= res$set95[2L]
      },
      TRUE
    )
    cat(sprintf(
      "95%% set holds r0 in %d of %d runs (at least %d needed)\n",
      sum(held),
      length(held),
      length(held) - set95Misses
    ))
    ok <- ok && sum(!held) <= set95Misses
  }
  if (!ok) {
    quit(status = 1L)
  }
  cat(
    "\nOK: the shape moves and the chains agree on it, and the predictive law covers fresh counts\n"
  )
}
