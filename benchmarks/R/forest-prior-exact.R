#!/usr/bin/env Rscript

# Exact-posterior gate for the rule a forest's coefficient and size follow in
# a model of several forests (docs/plans/forest-multiplier.md, "The rule"),
# its held arms. One binary predictor and sigma held, so a single-tree forest
# is a stump or one split and its leaves integrate out in closed form: given
# how many of each forest's trees split, the forests' cell values are jointly
# normal, the response is their sum through the held multipliers plus noise,
# and the posterior over the configurations is a finite sum. No quadrature.
#
# The oracle is written from the rule and reads nothing off the sampler but
# the unit and the shift: a held coefficient is 1 with no basis and (0, 1)
# for the two levels of a factor, so the forest enters with multiplier 1 or
# z; the forest's own fit then has prior sd its sd, exactly, each of T trees'
# leaves N(0, sd^2 / T); a held forest with no basis that states no sd has a
# quarter of the response's range, taken here from the response drawn, and a
# held factor one unit. Each arm compares the posterior mean of each forest's
# term in each cell, the second level's for a factor, and the probability
# that each forest's tree splits.
#
# Verdict: an arm fails when a gap exceeds 0.0005 plus 4 standard errors over
# seeds. With about forty functionals a run a lone z between 4 and 5 is
# expected now and then, so such an arm is run once more on a second set of
# seeds and fails only if the same functional is out again.
#
# Not seen here: a drawn coefficient, a numeric multiplier (both push 2's),
# probit and logistic, and weights.
#
# Usage: Rscript forest-prior-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
numSeeds <- if (quick) 32L else 64L
numKept <- if (quick) 25000L else 50000L
numBurn <- 500L
numCores <- as.integer(Sys.getenv("EQUIVALENCE_CORES", "2"))
gapFloor <- 0.0005
zBound <- 4
zRerun <- 5

# ---- the design ----

set.seed(20261017L)
perGroup <- 4L
cell <- rep(0:1, each = 4L * perGroup)
z1 <- rep(rep(c(0, 1, 0, 1), each = perGroup), 2L)
z2 <- rep(rep(c(0, 0, 1, 1), each = perGroup), 2L)
n <- length(cell)
x <- matrix(as.double(cell), ncol = 1L, dimnames = list(NULL, "x"))
sigma <- 0.5 # the residual sd, held; fixed() takes its square
y <- 1 + 0.8 * cell + z1 * (0.5 + 0.4 * cell) - 0.3 * z2 + rnorm(n, sd = sigma)
base <- 0.5
power <- 2

# ---- the oracle ----

# forests: a list of list(multiplier = <n-vector>, trees = T, sd = <response
# units>). The forest's two cell values G_f = (G_f(0), G_f(1)) have, given k
# of its T trees split, covariance (sd^2 / T) [[T, T - k], [T - k, T]], and k
# is binomial(T, base). y - shift = sum_f m_f G_f(cell) + N(0, sigma^2).
exactPosterior <- function(forests, shift) {
  K <- length(forests)
  design <- matrix(0, n, 2L * K)
  for (f in seq_len(K)) {
    design[cbind(seq_len(n), 2L * (f - 1L) + cell + 1L)] <- forests[[
      f
    ]]$multiplier
  }
  counts <- as.matrix(expand.grid(lapply(forests, function(forest) {
    0:forest$trees
  })))
  residual <- y - shift
  logWeights <- numeric(nrow(counts))
  means <- matrix(0, nrow(counts), 2L * K)
  for (r in seq_len(nrow(counts))) {
    prior <- matrix(0, 2L * K, 2L * K)
    logPrior <- 0
    for (f in seq_len(K)) {
      trees <- forests[[f]]$trees
      k <- counts[r, f]
      leaf <- forests[[f]]$sd^2 / trees
      block <- leaf * matrix(c(trees, trees - k, trees - k, trees), 2L)
      at <- 2L * (f - 1L) + 1:2
      prior[at, at] <- block
      logPrior <- logPrior + dbinom(k, trees, base, log = TRUE)
    }
    marginal <- design %*% prior %*% t(design) + diag(sigma^2, n)
    chol <- chol(marginal)
    solved <- backsolve(chol, forwardsolve(t(chol), residual))
    logWeights[r] <- logPrior -
      sum(log(diag(chol))) -
      0.5 * sum(residual * solved)
    means[r, ] <- prior %*% t(design) %*% solved
  }
  weights <- exp(logWeights - max(logWeights))
  weights <- weights / sum(weights)
  term <- colSums(weights * means)
  splits <- vapply(
    seq_len(K),
    function(f) sum(weights * counts[, f]) / forests[[f]]$trees,
    0
  )
  names(term) <- paste0(
    "term",
    rep(seq_len(K), each = 2L),
    "[",
    rep(0:1, K),
    "]"
  )
  names(splits) <- paste0("split", seq_len(K))
  c(term, splits)
}

# ---- the sampler ----

# `arm`: a list of forests, each list(basis = NULL, z1 or z2, sd = NULL or a
# number of units), and n.trees
heldForest <- function(spec, trees) {
  arguments <- list(
    n.trees = trees,
    base = base,
    power = power,
    amplitude = dbarts::dbartsPriors$fixed()
  )
  if (!is.null(spec$basis)) {
    arguments$basis <- factor(spec$basis, levels = c(0, 1))
  }
  if (!is.null(spec$sd)) {
    arguments$sd <- spec$sd
  }
  do.call(dbarts::dbartsForests$forest, arguments)
}

# the forests' terms in each cell from the combined train draws: a held
# forest's multiplier is 1 or its z, so the rows at z = 0 of every factor
# carry the forest with no basis and each factor adds its rows at z = 1
termDraws <- function(train, shift, arm) {
  K <- length(arm$forests)
  plain <- matrix(0, 2L, ncol(train))
  out <- list()
  baseRows <- function(c) which(cell == c & z1 == 0 & z2 == 0)
  for (c in 0:1) {
    plain[c + 1L, ] <- colMeans(train[baseRows(c), , drop = FALSE]) - shift
  }
  for (f in seq_len(K)) {
    spec <- arm$forests[[f]]
    for (c in 0:1) {
      name <- sprintf("term%d[%d]", f, c)
      if (is.null(spec$basis)) {
        out[[name]] <- plain[c + 1L, ]
      } else {
        on <- if (identical(spec$basis, z1)) {
          z1 == 1 & z2 == 0
        } else {
          z2 == 1 & z1 == 0
        }
        rows <- which(cell == c & on)
        out[[name]] <- colMeans(train[rows, , drop = FALSE]) -
          shift -
          plain[c + 1L, ]
      }
    }
  }
  out
}

samplerFit <- function(seed, arm, kept) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.cuts = 1L,
    updateState = FALSE,
    keepTrainingFits = TRUE
  )
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(sigma^2)),
    forests = lapply(arm$forests, heldForest, trees = arm$trees),
    control = control
  )
  result <- sampler$run(numBurn, kept)
  shift <- sampler$getLeafPrior(1L)$response.shift
  terms <- termDraws(result$train, shift, arm)
  counts <- result$varcount
  K <- length(arm$forests)
  splits <- lapply(seq_len(K), function(f) counts[1L, f, ] / arm$trees)
  names(splits) <- paste0("split", seq_len(K))
  vapply(c(terms, splits), mean, 0)
}

# the unit and the shift, the two numbers the oracle takes off a sampler
samplerConstants <- function(arm) {
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(sigma^2)),
    forests = lapply(arm$forests, heldForest, trees = arm$trees),
    control = dbartsControl(n.chains = 1L, n.threads = 1L, n.cuts = 1L)
  )
  list(
    unit = attr(sampler$control, "bartcore.forests", exact = TRUE)$unit,
    shift = sampler$getLeafPrior(1L)$response.shift
  )
}

# the oracle's forests for an arm: a stated sd as written, a default by the
# rule, from the response drawn
oracleForests <- function(arm, unit) {
  lapply(arm$forests, function(spec) {
    sd <- if (!is.null(spec$sd)) {
      spec$sd
    } else if (is.null(spec$basis)) {
      diff(range(y)) / 4
    } else {
      unit
    }
    list(
      multiplier = if (is.null(spec$basis)) rep(1, n) else spec$basis,
      trees = arm$trees,
      sd = sd
    )
  })
}

runArm <- function(name, arm) {
  constants <- samplerConstants(arm)
  exact <- exactPosterior(oracleForests(arm, constants$unit), constants$shift)
  verdict <- function(seeds) {
    fits <- parallel::mclapply(
      seeds,
      samplerFit,
      arm = arm,
      kept = numKept,
      mc.cores = numCores
    )
    fits <- do.call(rbind, fits)
    mean <- colMeans(fits)[names(exact)]
    se <- apply(fits, 2L, sd)[names(exact)] / sqrt(nrow(fits))
    gap <- mean - exact
    list(
      mean = mean,
      se = se,
      gap = gap,
      z = gap / se,
      out = abs(gap) > gapFloor + zBound * se
    )
  }
  started <- proc.time()[["elapsed"]]
  first <- verdict(seq_len(numSeeds))
  failed <- first$out
  rerun <- any(failed) && all(abs(first$z[failed]) < zRerun)
  if (rerun) {
    second <- verdict(numSeeds + seq_len(numSeeds))
    failed <- failed & second$out
  }
  cat(sprintf("%s (%.1f s)\n", name, proc.time()[["elapsed"]] - started))
  cat(sprintf(
    "  %-10s %9s %9s %9s %7s\n",
    "quantity",
    "sampler",
    "exact",
    "se",
    "z"
  ))
  for (label in names(exact)) {
    cat(sprintf(
      "  %-10s %+9.5f %+9.5f %9.2e %+7.2f%s\n",
      label,
      first$mean[[label]],
      exact[[label]],
      first$se[[label]],
      first$z[[label]],
      if (failed[[label]]) "  <- FAIL" else ""
    ))
  }
  if (rerun) {
    cat("  a lone z under", zRerun, "was run again on a second set of seeds\n")
  }
  any(failed)
}

# ---- the arms ----

constants <- samplerConstants(list(
  forests = list(list(), list(basis = z1)),
  trees = 1L
))
unit <- constants$unit
cat(sprintf(
  "forest prior exact gate (%s): %d seeds of %d kept sweeps, unit %.4f\n",
  if (quick) "quick" else "full",
  numSeeds,
  numKept,
  unit
))
arms <- list(
  "no basis held at 0.3 of the unit" = list(
    forests = list(list(sd = 0.3 * unit), list(basis = z1)),
    trees = 1L
  ),
  "no basis held at 2 of the unit" = list(
    forests = list(list(sd = 2 * unit), list(basis = z1)),
    trees = 1L
  ),
  "no basis held at its default" = list(
    forests = list(list(), list(basis = z1, sd = 0.7 * unit)),
    trees = 1L
  ),
  "factor held at 0.6 of the unit" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1, sd = 0.6 * unit)),
    trees = 1L
  ),
  "factor held at 0.2 of the unit" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1, sd = 0.2 * unit)),
    trees = 1L
  ),
  "factor held at its default" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1)),
    trees = 1L
  ),
  "both held at their defaults" = list(
    forests = list(list(), list(basis = z1)),
    trees = 1L
  ),
  "three forests of three trees" = list(
    forests = list(list(), list(basis = z1), list(basis = z2, sd = 0.4 * unit)),
    trees = 3L
  )
)
anyFailure <- FALSE
for (name in names(arms)) {
  anyFailure <- runArm(name, arms[[name]]) || anyFailure
}
if (anyFailure) {
  quit(status = 1L)
}
cat("OK: every held arm matches the exact posterior\n")
