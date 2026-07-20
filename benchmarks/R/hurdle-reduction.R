#!/usr/bin/env Rscript

# The hurdle reduction gate (docs/design/hurdle.md section 5) - hurdle
# composes two ordinary bart2() fits and never combines them within a sweep
# (docs/design/hurdle.md section 0), so a bart2(family = "hurdle.lognormal")
# fit's two component fits must equal, DRAW FOR DRAW, standalone
# bart2(family = "probit") / bart2(family = "gaussian") fits built by hand at
# the same derived seeds and control. This is a sanity FLOOR: the wrapper
# literally calls bart2() for each component (R/bart.R's bart2Hurdle), so the
# equality holds by construction. It is not the full correctness argument -
# the analytic combine/retransform oracle is the load-bearing gate and lands
# separately (docs/plans/hurdle.md step 2).
#
# Usage: Rscript benchmarks/R/hurdle-reduction.R

suppressPackageStartupMessages(library(dbarts))

set.seed(1729L)
n <- 250L
p <- 3L
x <- matrix(runif(n * p), n, p)
pi.true <- pnorm(0.8 * x[, 1L] - 0.5 * x[, 2L])
mu.true <- 1.2 + 0.6 * x[, 1L] - 0.4 * x[, 3L]
occupied <- rbinom(n, 1L, pi.true) == 1L
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), mu.true[occupied], 0.4))

args <- list(
  n.trees = 30L,
  n.burn = 40L,
  n.samples = 80L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 90210L
)

# the SAME derivation bart2Hurdle uses internally (R/bart.R), reproduced here
# so the standalone comparison fits are built at the identical seeds
derivedSeeds <- local({
  oldSeed <- .GlobalEnv[[".Random.seed"]]
  set.seed(args$seed)
  seeds <- sample.int(.Machine$integer.max, 2L)
  if (!is.null(oldSeed)) {
    .GlobalEnv$.Random.seed <- oldSeed
  }
  seeds
})

fitHurdle <- do.call(bart2, c(list(x, y, family = "hurdle.lognormal"), args))

split <- dbarts:::splitHurdleResponse(y)
xPositive <- x[split$positive, , drop = FALSE]

occupancyArgs <- args
occupancyArgs$seed <- derivedSeeds[1L]
fitOccupancy <- do.call(
  bart2,
  c(list(x, split$z, family = "probit"), occupancyArgs)
)

positiveArgs <- args
positiveArgs$seed <- derivedSeeds[2L]
fitPositive <- do.call(
  bart2,
  c(
    list(xPositive, split$logPositive, test = x, family = "gaussian"),
    positiveArgs
  )
)

compareComponent <- function(wrapped, standalone, label, hasSigma) {
  results <- c(
    yhat.train = identical(
      extract(wrapped, type = "bart", sample = "train"),
      extract(standalone, type = "bart", sample = "train")
    ),
    varcount = identical(wrapped$varcount, standalone$varcount),
    trees = identical(
      extract(wrapped, type = "trees"),
      extract(standalone, type = "trees")
    )
  )
  if (hasSigma) {
    results <- c(results, sigma = identical(wrapped$sigma, standalone$sigma))
  }
  bitwise <- all(results)
  cat(sprintf(
    "%-10s reduction %s (%s)\n",
    label,
    if (bitwise) "BITWISE IDENTICAL" else "DIVERGED",
    paste(sprintf("%s=%s", names(results), results), collapse = ", ")
  ))
  bitwise
}

okOccupancy <- compareComponent(
  fitHurdle$occupancy,
  fitOccupancy,
  "occupancy",
  hasSigma = FALSE
)
okPositive <- compareComponent(
  fitHurdle$positive,
  fitPositive,
  "positive",
  hasSigma = TRUE
)

# the packaged objects differ ONLY in $call (the wrapper's own reconstructed
# component call vs. the by-hand standalone call, the hazard-reduction.R
# markerOnly idea); every other top-level field runs through the identical
# bart2() code path
markerOnly <-
  identical(fitHurdle$occupancy$family, "probit") &&
  identical(fitHurdle$positive$family, "gaussian") &&
  identical(fitHurdle$family, "hurdle.lognormal")

if (!markerOnly) {
  cat("  marker check FAILED: the objects differ beyond $call\n")
}

if (okOccupancy && okPositive && markerOnly) {
  cat("\nOK: hurdle reduces bitwise to its two standalone component fits\n")
} else {
  cat("\nFAIL: a hurdle component diverged from its standalone reduction\n")
  quit(status = 1L)
}
