#!/usr/bin/env Rscript

# Statistical-equivalence harness: verifies that two engine builds target the
# same posterior (docs/design/core-generalization.md). Data are fixed per
# scenario; only the MCMC seed varies, so across-seed spread is pure Monte
# Carlo variability. Record a baseline with build A, then compare from a tree
# with build B installed; per-summary Welch z-statistics should look standard
# normal if the posteriors agree.
#
# If both builds produce identical RNG streams the comparison reports an
# exact match instead, which doubles as a no-unintended-RNG-shift check for
# refactors of the current engine.
#
# Usage:
#   Rscript equivalence.R record [out.rds]
#   Rscript equivalence.R compare baseline.rds
# Append 'quick' for a fast smoke test (not comparable to full runs).

suppressPackageStartupMessages(library(dbarts))

args  <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args  <- setdiff(args, "quick")
useNewEngine <- "engine=new" %in% args
args  <- setdiff(args, "engine=new")
mode  <- if (length(args) >= 1L) args[[1L]] else "record"

if (useNewEngine) {
  source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                   "bartcore-shim.R"))
  loadBartcoreShim()
}

n.seeds <- if (quick) 3L else 20L
ndpost  <- if (quick) 250L else 1000L
nskip   <- if (quick) 100L else 500L
ntree   <- if (quick) 50L else 200L
n.test  <- 25L

friedman <- function(x)
  10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] + 5 * x[, 5L]

makeScenarios <- function() {
  result <- list()

  set.seed(5101L)
  x <- matrix(runif(500L * 10L), 500L)
  result$friedman <- list(
    x = x, y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test), binary = FALSE
  )

  set.seed(5102L)
  x <- matrix(runif(1000L * 10L), 1000L)
  result$probit <- list(
    x = x, y = rbinom(1000L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test), binary = TRUE
  )

  set.seed(5103L)
  x <- matrix(runif(500L * 10L), 500L)
  weights <- sample(c(1, 2), 500L, replace = TRUE)
  result$weighted <- list(
    x = x, y = friedman(x) + rnorm(500L) / sqrt(weights), weights = weights,
    x.test = matrix(runif(n.test * 10L), n.test), binary = FALSE
  )

  # exercises the non-uniform split-variable selection path (DART's seam)
  set.seed(5104L)
  x <- matrix(runif(500L * 10L), 500L)
  result$splitprobs <- list(
    x = x, y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test), binary = FALSE,
    splitprobs = c(rep(0.15, 5L), rep(0.05, 5L))
  )

  # binary with the default chi hyperprior on k; driven through the sampler
  # API (old engine) and the package's bartcore bridge (new engine), since
  # bart() cannot express an adaptive k
  set.seed(5105L)
  x <- matrix(runif(600L * 10L), 600L)
  result$chik <- list(
    x = x, y = rbinom(600L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test), binary = TRUE,
    samplerApi = TRUE
  )

  # quantile cut points over a mix of continuous columns (thinned to numcut)
  # and discrete ones (columns 4-5 and 8-10 on an 11-level grid, inducing
  # 10 unique-midpoint cuts each)
  set.seed(5106L)
  x <- matrix(runif(500L * 10L), 500L)
  x.test <- matrix(runif(n.test * 10L), n.test)
  discrete <- c(4L, 5L, 8L, 9L, 10L)
  x[, discrete] <- round(x[, discrete], 1L)
  x.test[, discrete] <- round(x.test[, discrete], 1L)
  result$quants <- list(
    x = x, y = friedman(x) + rnorm(500L),
    x.test = x.test, binary = FALSE, usequants = TRUE
  )

  result
}

fitViaSamplerApi <- function(scenario, engineIsNew) {
  control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = ntree,
                           updateState = FALSE)
  sampler <- dbarts(scenario$x, scenario$y, test = scenario$x.test,
                    control = control)
  if (engineIsNew) {
    bcSampler <- dbarts:::bartcoreSampler(sampler)
    r <- dbarts:::bartcoreRun(bcSampler, nskip, ndpost)
    list(yhat.test = t(r$yhat.test), varcount = t(r$varcount),
         sigma = r$sigma, k = r$k)
  } else {
    r <- sampler$run(nskip, ndpost)
    list(yhat.test = t(r$test), varcount = t(r$varcount),
         sigma = r$sigma, k = r$k)
  }
}

fitSummaries <- function(scenario, seed) {
  set.seed(seed)
  # Test-data weights are irrelevant here (no posterior-predictive use);
  # muffle only that warning so real ones stay visible.
  splitprobs <- scenario$splitprobs
  fit <- if (!is.null(scenario$samplerApi)) {
    fitViaSamplerApi(scenario, useNewEngine)
  } else if (useNewEngine) {
    bartcore_bart(
      scenario$x, scenario$y, x.test = scenario$x.test,
      weights = scenario$weights, splitprobs = splitprobs,
      usequants = isTRUE(scenario$usequants),
      ntree = ntree, ndpost = ndpost, nskip = nskip
    )
  } else if (!is.null(splitprobs)) withCallingHandlers(
    bart(
      scenario$x, scenario$y, x.test = scenario$x.test,
      weights = scenario$weights, splitprobs = splitprobs,
      usequants = isTRUE(scenario$usequants),
      ntree = ntree, ndpost = ndpost, nskip = nskip,
      nchain = 1L, nthread = 1L, verbose = FALSE
    ),
    warning = function(w) {
      if (grepl("'weights' are ignored for test data", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  ) else withCallingHandlers(
    bart(
      scenario$x, scenario$y, x.test = scenario$x.test,
      weights = scenario$weights,
      usequants = isTRUE(scenario$usequants),
      ntree = ntree, ndpost = ndpost, nskip = nskip,
      nchain = 1L, nthread = 1L, verbose = FALSE
    ),
    warning = function(w) {
      if (grepl("'weights' are ignored for test data", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
  p <- ncol(scenario$x)
  result <- c(
    setNames(colMeans(fit$yhat.test), paste0("fhat.test.", seq_len(n.test))),
    setNames(colMeans(fit$varcount / rowSums(fit$varcount)), paste0("vprop.", seq_len(p)))
  )
  if (!scenario$binary)
    result <- c(result, sigma.mean = mean(fit$sigma), sigma.sd = sd(fit$sigma))
  if (!is.null(fit$k))
    result <- c(result, k.mean = mean(fit$k), k.sd = sd(fit$k))
  result
}

runAll <- function(scenarios) {
  lapply(scenarios, function(scenario)
    do.call(rbind, lapply(seq_len(n.seeds), function(seed) fitSummaries(scenario, seed)))
  )
}

scenarios <- makeScenarios()

if (mode == "record") {
  out.file <- if (length(args) >= 2L) args[[2L]] else "equivalence-baseline.rds"
  results <- runAll(scenarios)
  meta <- list(
    rev = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    date = format(Sys.Date()), quick = quick,
    n.seeds = n.seeds, ndpost = ndpost, nskip = nskip, ntree = ntree
  )
  saveRDS(list(meta = meta, results = results), out.file)
  cat("wrote baseline for", length(results), "scenarios x", n.seeds, "seeds to", out.file, "\n")
} else if (mode == "compare") {
  if (length(args) < 2L) stop("usage: equivalence.R compare baseline.rds")
  baseline <- readRDS(args[[2L]])
  settings <- baseline$meta[c("quick", "n.seeds", "ndpost", "nskip", "ntree")]
  if (!identical(settings, list(quick = quick, n.seeds = n.seeds, ndpost = ndpost,
                                nskip = nskip, ntree = ntree)))
    stop("baseline was recorded with different settings: ",
         paste(names(settings), unlist(settings), sep = "=", collapse = ", "))

  results <- runAll(scenarios)
  anyFailure <- FALSE
  for (name in names(baseline$results)) {
    a <- baseline$results[[name]]
    b <- results[[name]]
    if (identical(a, b)) {
      cat(sprintf("%-10s identical draws (same RNG stream)\n", name))
      next
    }
    z <- (colMeans(a) - colMeans(b)) /
      sqrt(apply(a, 2L, var) / nrow(a) + apply(b, 2L, var) / nrow(b))
    n.warn <- sum(abs(z) > 3, na.rm = TRUE)
    n.fail <- sum(abs(z) > 4, na.rm = TRUE)
    cat(sprintf(
      "%-10s %d summaries, max |z| = %.2f%s%s\n",
      name, length(z), max(abs(z), na.rm = TRUE),
      if (n.warn > 0L) sprintf(", %d with |z| > 3", n.warn) else "",
      if (n.fail > 0L) sprintf(", %d with |z| > 4 <- FAIL", n.fail) else ""
    ))
    if (n.fail > 0L) {
      anyFailure <- TRUE
      failed <- which(abs(z) > 4)
      cat("  worst offenders:", paste0(names(z)[failed], " (z=", round(z[failed], 2L), ")",
                                       collapse = ", "), "\n")
    }
  }
  if (anyFailure) quit(status = 1L)
  cat("\nOK: posteriors statistically indistinguishable at |z| > 4\n")
} else {
  stop("unknown mode '", mode, "'; use record or compare")
}
