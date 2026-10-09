#!/usr/bin/env Rscript

# How fast k mixes under probit with the probit rescaling step
# (docs/design/probit-k-scale-move.md). A measurement, not a gate: no baseline
# and no pass/fail exit. Its runs are that design's kill-criteria evidence.
#
#   census     k's and the average fit's integrated autocorrelation time over
#              the SBC probit-k arm's prior-drawn datasets: per dataset, the
#              truth and the data drawn as the arm draws them (seeded by the
#              dataset's index, so every run sees the same data), the arm's
#              start, 30,000 sweeps of burn-in and 200,000 recorded at thin 10.
#              Geyer's initial monotone sequence, in sweeps.
#   truth      the start-from-truth invariance statistic: the state is an exact
#              posterior draw (trees, leaves and k0 from the prior, y given the
#              fit, latents given y), the chain runs on, and log k at lags 300
#              to 3000 averaged within a replication must not drift from log k0.
#              Reports its z; a correct kernel sits within |z| 3 at R = 2000.
#   summarize  the census tables and the truth z from saved runs.
#
# Usage:
#   Rscript probit-k-mixing.R census <first> <last> <out.rds>
#   Rscript probit-k-mixing.R truth <seed> <reps> <out.rds>
#   Rscript probit-k-mixing.R summarize <file.rds>...
# Run from the repository root (it sources benchmarks/R/sbc.R for the arm's
# configuration).

suppressPackageStartupMessages(library(dbarts))
source(file.path("benchmarks", "R", "sbc.R"))

args <- commandArgs(trailingOnly = TRUE)
action <- args[1L]

# Geyer's initial monotone sequence estimate of the integrated autocorrelation
# time, in units of the series' own spacing
iatGeyer <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  if (!(var(x) > 0)) {
    return(NA_real_)
  }
  m <- 2^ceiling(log2(2 * n))
  fx <- fft(c(x, numeric(m - n)))
  ac <- Re(fft(Mod(fx)^2, inverse = TRUE))[seq_len(n)]
  ac <- ac / ac[1L]
  numPairs <- floor((n - 1) / 2)
  pairSums <- ac[2 * (0:(numPairs - 1)) + 1] + ac[2 * (0:(numPairs - 1)) + 2]
  firstNonPositive <- which(pairSums <= 0)[1L]
  if (!is.na(firstNonPositive)) {
    pairSums <- pairSums[seq_len(firstNonPositive - 1L)]
  }
  if (length(pairSums) == 0L) {
    return(1)
  }
  max(-1 + 2 * sum(cummin(pairSums)), 1 / n)
}

# theta0 drawn from the prior: trees, k0, leaves given k0, and the fit
drawTruth <- function(sampler, config) {
  sampler$sampleTreesFromPrior()
  k0 <- sbcKDraw(config)
  sampler$setLeafPrior(dbartsPriors$normal(k = k0))
  sampler$sampleLeafParametersFromPrior()
  list(k0 = k0, f0 = as.numeric(sampler$predict(config$x)))
}

runCensus <- function(first, last, out, burn = 30000L, numRecorded = 200000L) {
  thin <- 10L
  config <- sbcConfigProbitK()
  set.seed(1L)
  sampler <- sbcMakeSampler(config, numRecorded %/% thin, thin, 1L)
  results <- list()
  for (r in first:last) {
    set.seed(100000L + r)
    truth <- drawTruth(sampler, config)
    y0 <- sbcSimulate(config, truth$f0, 1.0)
    # the arm's start: fresh prior trees and a k from the prior
    set.seed(200000L + r)
    start <- drawTruth(sampler, config)
    sampler$setLeafPrior(config$nodePrior)
    sampler$setResponse(y0)
    started <- proc.time()[[3L]]
    sampler$run(burn - 1L, 1L)
    burned <- proc.time()[[3L]]
    fit <- sampler$run(0L, numRecorded %/% thin)
    ended <- proc.time()[[3L]]
    k <- as.numeric(fit$k)
    result <- list(
      r = r,
      k0 = truth$k0,
      kStart = start$k0,
      tauK = thin * iatGeyer(log(k)),
      tauF = thin * iatGeyer(colMeans(fit$train)),
      pSmall = mean(k < 0.135),
      secondsBurn = burned - started,
      secondsRecorded = ended - burned,
      kLagThin = k[seq(1L, length(k), by = 100L)]
    )
    results[[length(results) + 1L]] <- result
    cat(sprintf(
      "r %d k0 %.3f tau(k) %.0f tau(avg f) %.0f %.1f s\n",
      r,
      truth$k0,
      result$tauK,
      result$tauF,
      ended - started
    ))
  }
  saveRDS(
    list(
      action = "census",
      burn = burn,
      numRecorded = numRecorded,
      thin = thin,
      results = results
    ),
    out
  )
}

runTruth <- function(seed, numReps, out) {
  config <- sbcConfigProbitK()
  set.seed(seed)
  sampler <- sbcMakeSampler(config, 1L, 1L, seed)
  lags <- c(1, 3, 10, 30, 100, 300, 1000, 2000, 3000)
  results <- vector("list", numReps)
  started <- proc.time()[[3L]]
  for (r in seq_len(numReps)) {
    truth <- drawTruth(sampler, config)
    sampler$setLeafPrior(config$nodePrior)
    sampler$setResponse(sbcSimulate(config, truth$f0, 1.0))
    k <- numeric(length(lags))
    previous <- 0
    for (i in seq_along(lags)) {
      fit <- sampler$run(as.integer(lags[i] - previous) - 1L, 1L)
      k[i] <- as.numeric(fit$k)
      previous <- lags[i]
    }
    results[[r]] <- list(k0 = truth$k0, k = k)
  }
  cat(sprintf(
    "seed %d: %d replications in %.0f s\n",
    seed,
    numReps,
    proc.time()[[3L]] - started
  ))
  saveRDS(
    list(action = "truth", lags = lags, results = results),
    out
  )
}

summarize <- function(files) {
  # a run saved while the step had a switch and recorded with it off is not
  # this script's measurement
  runs <- Filter(
    function(run) !identical(run$rescale, FALSE),
    lapply(files, readRDS)
  )
  quantileAt <- function(v, p) as.numeric(quantile(v, p, na.rm = TRUE))
  census <- Filter(function(run) run$action == "census", runs)
  if (length(census) > 0L) {
    results <- do.call(c, lapply(census, `[[`, "results"))
    tauK <- vapply(results, `[[`, 0, "tauK")
    tauF <- vapply(results, `[[`, 0, "tauF")
    seconds <- vapply(results, `[[`, 0, "secondsRecorded")
    essPerSecond <- census[[1L]]$numRecorded / tauK / seconds
    cat(sprintf(
      "census: %d datasets\n  tau(k) median %.0f, 90th %.0f, 99th %.0f, max %.0f\n  tau(avg f) median %.0f, 90th %.0f\n  ESS(k)/s median %.0f, 10th %.0f; %.3f s a 10,000 recorded sweeps\n",
      length(results),
      median(tauK),
      quantileAt(tauK, 0.9),
      quantileAt(tauK, 0.99),
      max(tauK),
      median(tauF),
      quantileAt(tauF, 0.9),
      median(essPerSecond),
      quantileAt(essPerSecond, 0.1),
      1e4 * median(seconds) / census[[1L]]$numRecorded
    ))
  }
  truth <- Filter(function(run) run$action == "truth", runs)
  if (length(truth) > 0L) {
    results <- do.call(c, lapply(truth, `[[`, "results"))
    lags <- truth[[1L]]$lags
    k <- t(vapply(results, `[[`, numeric(length(lags)), "k"))
    k0 <- vapply(results, `[[`, 0, "k0")
    drift <- rowMeans(log(k[, lags >= 300, drop = FALSE])) - log(k0)
    cat(sprintf(
      "truth: R = %d, window lags 300 to 3000: drift of log k %+.4f (se %.4f), z %+.2f\n",
      length(drift),
      mean(drift),
      sd(drift) / sqrt(length(drift)),
      mean(drift) / (sd(drift) / sqrt(length(drift)))
    ))
  }
}

switch(
  action,
  census = runCensus(as.integer(args[2L]), as.integer(args[3L]), args[4L]),
  truth = runTruth(as.integer(args[2L]), as.integer(args[3L]), args[4L]),
  summarize = summarize(args[-1L]),
  stop("the first argument is census, truth or summarize")
)
