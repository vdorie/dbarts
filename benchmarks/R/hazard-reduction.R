#!/usr/bin/env Rscript

# The discrete-time hazard reduction gate (docs/design/survival.md section 5) -
# the single load-bearing correctness gate, and it is free. A discrete-time
# hazard fit IS an ordinary probit/logistic fit on a person-period-expanded
# design (the family adds no engine code and remaps to the binary token before
# the sampler is built), so a hazard fit must equal, DRAW FOR DRAW, the binary
# fit a user would get by expanding the data by hand and calling
# bart2(family = "probit"/"logistic") on the same design with the same seed and
# control. This makes the "hazard is sugar" thesis testable: the two fits
# consume the identical RNG stream, so the trees, latent fits, and varcount
# draws are BITWISE identical. The comparison is over the draw components (not
# the packaged objects, which differ by construction in the hazard-only
# $periods marker). Both links are covered.
#
# Usage: Rscript benchmarks/R/hazard-reduction.R

suppressPackageStartupMessages(library(dbarts))

set.seed(919L)
n <- 300L
p <- 4L
x <- matrix(runif(n * p), n, p)
f <- 0.9 * sin(pi * x[, 1L]) + 0.6 * x[, 2L] - 0.4 * x[, 3L]

# a genuinely discrete grid: a small integer number of periods, with a
# baseline hazard that rises over time and a covariate effect
K <- 8L
baseline <- seq(-2.2, -0.6, length.out = K)
hazard <- plogis(outer(baseline, rep(0, n)) + matrix(f, K, n, byrow = TRUE))
event <- rep(K + 1L, n) # K + 1 == survived the whole grid (censored at K)
for (i in seq_len(n)) {
  for (k in seq_len(K)) {
    if (runif(1L) < hazard[k, i]) {
      event[i] <- k
      break
    }
  }
}
# right-censor at the horizon, plus some random early censoring
censorPeriod <- pmin(event, sample.int(K, n, replace = TRUE))
status <- as.double(event <= censorPeriod & event <= K)
time <- as.double(pmin(event, censorPeriod, K))

args <- list(
  n.trees = 40L,
  n.burn = 60L,
  n.samples = 120L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 271828L
)

# the by-hand expanded design, built by the same exported expander the hazard
# ingestion uses - same column order (period appended last), same values, so
# the two fits build the identical cut grid
expansion <- dbarts:::expandDiscreteTimeHazard(x, time, status)

compareLink <- function(hazardToken, binaryFamily) {
  fitHazard <- do.call(
    bart2,
    c(list(x, cbind(time, status), family = hazardToken), args)
  )
  fitBinary <- do.call(
    bart2,
    c(list(expansion$x, expansion$y, family = binaryFamily), args)
  )

  results <- c(
    yhat.train = identical(
      extract(fitHazard, type = "bart", sample = "train"),
      extract(fitBinary, type = "bart", sample = "train")
    ),
    varcount = identical(fitHazard$varcount, fitBinary$varcount),
    trees = identical(
      extract(fitHazard, type = "trees"),
      extract(fitBinary, type = "trees")
    )
  )
  # the packaged objects differ ONLY in the hazard marker: $periods present on
  # the hazard fit, absent on the binary one, and $family reads the same binary
  # token on both
  markerOnly <-
    identical(fitHazard$family, binaryFamily) &&
    identical(fitBinary$family, binaryFamily) &&
    !is.null(fitHazard[["periods"]]) &&
    is.null(fitBinary[["periods"]])

  bitwise <- all(results) && markerOnly
  cat(sprintf(
    "%-16s reduction to %-9s %s (%s)\n",
    hazardToken,
    binaryFamily,
    if (bitwise) "BITWISE IDENTICAL" else "DIVERGED",
    paste(sprintf("%s=%s", names(results), results), collapse = ", ")
  ))
  if (!markerOnly) {
    cat("  marker check FAILED: the objects differ beyond $periods\n")
  }
  bitwise
}

okProbit <- compareLink("hazard", "probit")
okLogistic <- compareLink("hazard.logistic", "logistic")

if (okProbit && okLogistic) {
  cat(
    "\nOK: discrete-time hazard reduces bitwise to the binary fit (both links)\n"
  )
} else {
  cat("\nFAIL: the hazard fit diverged from its by-hand binary reduction\n")
  quit(status = 1L)
}
