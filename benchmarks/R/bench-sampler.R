#!/usr/bin/env Rscript

# End-to-end sampler timing benchmarks with baseline record/compare.
# The zero-regression gate for the core-generalization work
# (docs/design/core-generalization.md).
#
# Usage:
#   Rscript bench-sampler.R                    run and print
#   Rscript bench-sampler.R record [out.csv]   run and write a baseline
#   Rscript bench-sampler.R compare base.csv   run and compare; exits 1 on
#                                              any metric > 5% slower
# Append 'quick' for a fast smoke test (not comparable to full runs).

suppressPackageStartupMessages(library(dbarts))

args  <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args  <- setdiff(args, "quick")
# engine=new times the bartcore engine, so recording under classic and
# comparing under engine=new measures the cross-engine gap on one build
useNewEngine <- "engine=new" %in% args
args  <- setdiff(args, "engine=new")
engine <- if (useNewEngine) "bartcore" else "classic"
mode  <- if (length(args) >= 1L) args[[1L]] else "print"

genFriedman <- function(n, p = 10L) {
  x <- matrix(runif(n * p), n, p)
  f <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] + 5 * x[, 5L]
  list(x = x, f = f, y = f + rnorm(n))
}

newSampler <- function(x, y, n.trees) {
  control <- dbartsControl(
    verbose = FALSE, n.trees = n.trees, n.chains = 1L, n.threads = 1L,
    updateState = FALSE, engine = engine
  )
  dbarts(x, y, control = control)
}

timeMedian <- function(fn, reps) {
  median(vapply(
    seq_len(reps),
    function(i) system.time(fn())[["elapsed"]],
    numeric(1L)
  ))
}

runBenchmarks <- function(quick) {
  reps    <- if (quick) 1L else 3L
  n.samps <- if (quick) 50L else 500L
  rows    <- data.frame()
  addRow <- function(scenario, metric, value)
    rows <<- rbind(rows, data.frame(scenario = scenario, metric = metric, value = value))

  # Plain run throughput.
  runScenarios <- list(
    list(name = "run-n1000-p10-t75",   n = 1000L,  n.trees = 75L),
    list(name = "run-n1000-p10-t200",  n = 1000L,  n.trees = 200L),
    list(name = "run-n10000-p10-t75",  n = 10000L, n.trees = 75L)
  )
  if (quick) runScenarios <- runScenarios[1L]

  for (scenario in runScenarios) {
    set.seed(4001L)
    data <- genFriedman(scenario$n)
    sampler <- newSampler(data$x, data$y, scenario$n.trees)
    invisible(sampler$run(200L, 1L))
    elapsed <- timeMedian(function() invisible(sampler$run(0L, n.samps)), reps)
    addRow(scenario$name, "ms_per_iteration", 1000 * elapsed / n.samps)
  }

  # Binary probit throughput.
  set.seed(4002L)
  data <- genFriedman(1000L)
  y.binary <- rbinom(1000L, 1L, pnorm(scale(data$f)))
  sampler <- newSampler(data$x, y.binary, 75L)
  invisible(sampler$run(200L, 1L))
  elapsed <- timeMedian(function() invisible(sampler$run(0L, n.samps)), reps)
  addRow("run-binary-n1000-p10-t75", "ms_per_iteration", 1000 * elapsed / n.samps)

  # Embedded-Gibbs pattern: mutate offset, draw a single sample, repeat.
  set.seed(4003L)
  data <- genFriedman(1000L)
  sampler <- newSampler(data$x, data$y, 75L)
  invisible(sampler$run(200L, 1L))
  offsets <- matrix(rnorm(1000L * 20L, sd = 0.1), 1000L)
  # long enough that system.time's millisecond granularity stays well under
  # the 5% regression threshold
  n.gibbs <- if (quick) 20L else 250L
  elapsed <- timeMedian(function() {
    for (i in seq_len(n.gibbs)) {
      sampler$setOffset(offsets[, 1L + i %% 20L])
      invisible(sampler$run(0L, 1L))
    }
  }, reps)
  addRow("embedded-offset-run1-n1000-t75", "ms_per_gibbs_step", 1000 * elapsed / n.gibbs)

  # Single-column predictor replacement with tree validation/rollback. The
  # accept/reject mix of random replacements depends on the chain state (one
  # tiny leaf rejects most candidates), so time the two paths separately with
  # deterministic workloads: an identity swap always accepts (full
  # revalidation + fits rebuild) and a degenerate column always rejects
  # (early exit + rollback).
  set.seed(4004L)
  data <- genFriedman(1000L)
  sampler <- newSampler(data$x, data$y, 75L)
  invisible(sampler$run(200L, 1L))
  n.updates <- if (quick) 40L else 500L
  x2 <- data$x[, 2L]
  elapsed <- timeMedian(function() {
    for (i in seq_len(n.updates))
      invisible(sampler$setPredictor(x2, column = 2L, forceUpdate = FALSE))
  }, reps)
  addRow("setPredictor-accept-n1000-t75", "ms_per_update", 1000 * elapsed / n.updates)

  x2.degenerate <- rep(0.5, 1000L)
  elapsed <- timeMedian(function() {
    for (i in seq_len(n.updates))
      invisible(sampler$setPredictor(x2.degenerate, column = 2L, forceUpdate = FALSE))
  }, reps)
  addRow("setPredictor-reject-n1000-t75", "ms_per_update", 1000 * elapsed / n.updates)

  rows$value <- round(rows$value, 4L)
  rows$rev   <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
  rows$date  <- format(Sys.Date())
  rows$quick <- quick
  rows
}

results <- runBenchmarks(quick)

if (mode == "record") {
  out.file <- if (length(args) >= 2L) args[[2L]] else "sampler-baseline.csv"
  write.csv(results, out.file, row.names = FALSE)
  cat("wrote", nrow(results), "measurements to", out.file, "\n")
  print(results[c("scenario", "metric", "value")], row.names = FALSE)
} else if (mode == "compare") {
  if (length(args) < 2L) stop("usage: bench-sampler.R compare baseline.csv")
  baseline <- read.csv(args[[2L]])
  if (!identical(unique(baseline$quick), quick))
    warning("comparing against a baseline recorded at a different quick setting")
  merged <- merge(
    baseline[c("scenario", "metric", "value")], results[c("scenario", "metric", "value")],
    by = c("scenario", "metric"), suffixes = c(".base", ".curr")
  )
  merged$ratio <- round(merged$value.curr / merged$value.base, 3L)
  merged$flag  <- ifelse(merged$ratio > 1.05, "REGRESSION", "")
  print(merged, row.names = FALSE)
  if (any(merged$flag != "")) {
    cat("\nFAIL:", sum(merged$flag != ""), "metric(s) regressed more than 5%\n")
    quit(status = 1L)
  }
  cat("\nOK: no metric regressed more than 5%\n")
} else {
  print(results[c("scenario", "metric", "value")], row.names = FALSE)
}
