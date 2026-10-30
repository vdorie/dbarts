#!/usr/bin/env Rscript

# Statistical equivalence of forests that share one total and blocks(): with
# the same trees per group, y ~ forest(x1) + forest(x2) + forest(x3) and
# y ~ forest(x1 + x2 + x3, blocks = blocks(list("x1", "x2", "x3"))) are one
# model, each group's leaves over its own trees with the total's variance
# divided evenly among the groups. Four scenarios: three forests at 75 trees
# (25 each) under gaussian at k 2 and under probit with k drawn; two forests
# at 75 trees (38 and 37); and a stated count, forest(x1, n.trees = 10) beside
# two forests that split the rest (38 and 37), against trees.per.group =
# c(10, 38, 37) on a forest of 85. The functionals are the posterior mean of
# the fit at each row and of sigma (gaussian) or k (probit), each spelling run
# on its own seeds; a Welch z over seeds past 4 fails. With about a hundred
# functionals a z between 4 and 5 is seen now and then (a spelling against
# itself reaches 3.3), so such a scenario is run once more on a second set of
# seeds and fails only if the same functional is out again.
#
# Usage: Rscript forests-blocks-equivalence.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
numSeeds <- if (quick) 40L else 80L
numBurn <- 500L
numKept <- if (quick) 1000L else 2000L
numCores <- as.integer(Sys.getenv("EQUIVALENCE_CORES", "2"))
zBound <- 4
zRerun <- 5

set.seed(20261030L)
n <- 100L
data <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
signal <- with(data, sin(2 * pi * x1) + 2 * (x2 - 0.5)^2 + x3)
data$y <- signal + rnorm(n, sd = 0.5)
data$z <- as.double(rbinom(n, 1L, pnorm(signal - mean(signal))))

scenarios <- list(
  "three forests, gaussian" = list(
    forests = y ~ forest(x1) + forest(x2) + forest(x3),
    blocks = y ~ forest(x1 + x2 + x3, blocks = blocks(list("x1", "x2", "x3"))),
    family = "gaussian"
  ),
  "three forests, probit" = list(
    forests = z ~ forest(x1) + forest(x2) + forest(x3),
    blocks = z ~ forest(x1 + x2 + x3, blocks = blocks(list("x1", "x2", "x3"))),
    family = "probit"
  ),
  "two forests, 38 and 37" = list(
    forests = y ~ forest(x1) + forest(x2 + x3),
    blocks = y ~
      forest(x1 + x2 + x3, blocks = blocks(list("x1", c("x2", "x3")))),
    family = "gaussian"
  ),
  "a stated count, 10, 38 and 37" = list(
    forests = y ~ forest(x1, n.trees = 10L) + forest(x2) + forest(x3),
    blocks = y ~
      forest(
        x1 + x2 + x3,
        blocks = blocks(
          list("x1", "x2", "x3"),
          trees.per.group = c(10L, 38L, 37L)
        )
      ),
    family = "gaussian",
    blocksTrees = 85L
  )
)

posteriorMeans <- function(seed, formula, family, numTrees) {
  sampler <- dbarts(
    formula,
    data,
    family = family,
    control = dbartsControl(
      n.trees = numTrees,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE,
      seed = seed
    )
  )
  run <- sampler$run(numBurn, numKept)
  scalar <- if (family == "gaussian") {
    c(sigma = mean(run$sigma))
  } else {
    c(k = mean(run$k))
  }
  fit <- rowMeans(run$train)
  names(fit) <- paste0("fit[", seq_len(n), "]")
  c(scalar, fit)
}

draws <- function(formula, family, offset, numTrees = 75L) {
  do.call(
    rbind,
    parallel::mclapply(
      offset + seq_len(numSeeds),
      posteriorMeans,
      formula = formula,
      family = family,
      numTrees = numTrees,
      mc.cores = numCores
    )
  )
}

cat(sprintf(
  "forests against blocks() (%s): %d seeds a spelling, %d kept after %d\n",
  if (quick) "quick" else "full",
  numSeeds,
  numKept,
  numBurn
))
anyFailure <- FALSE
compare <- function(scenario, offset) {
  forests <- draws(scenario$forests, scenario$family, offset)
  blocked <- draws(
    scenario$blocks,
    scenario$family,
    offset + 1000L,
    if (is.null(scenario$blocksTrees)) 75L else scenario$blocksTrees
  )
  z <- (colMeans(forests) - colMeans(blocked)) /
    sqrt((apply(forests, 2L, var) + apply(blocked, 2L, var)) / numSeeds)
  list(forests = forests[, 1L], blocked = blocked[, 1L], z = z)
}
for (name in names(scenarios)) {
  first <- compare(scenarios[[name]], 0L)
  out <- abs(first$z) > zBound
  rerun <- any(out) && all(abs(first$z[out]) < zRerun)
  if (rerun) {
    out <- out & abs(compare(scenarios[[name]], 5000L)$z) > zBound
  }
  z <- first$z
  cat(sprintf(
    "%-30s %s %.4f vs %.4f (z %+.2f) | max |z| over %d fits %.2f%s%s\n",
    name,
    names(z)[1L],
    mean(first$forests),
    mean(first$blocked),
    z[[1L]],
    n,
    max(abs(z[-1L])),
    if (rerun) ", run again" else "",
    if (any(out)) "  <- FAIL" else ""
  ))
  anyFailure <- anyFailure || any(out)
}
if (anyFailure) {
  quit(status = 1L)
}
cat("OK: the forests and blocks() spellings agree\n")
