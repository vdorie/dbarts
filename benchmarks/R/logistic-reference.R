#!/usr/bin/env Rscript

# Validation of the bartcore binary samplers against the exact posterior,
# plus an informational comparison to the BART package.
#
# Part 1 (the gate): single-tree BART on one predictor with 3 cut points
# admits brute-force enumeration of the tree space and 1-D quadrature for
# the leaf marginals, so the posterior predictive is computable exactly.
# Both binary families are checked; failure means the sampler's stationary
# distribution is wrong. Data are a fixed grid with 25 observations per
# inter-cut cell so every realizable leaf is well occupied.
#
# Part 2 (context, requires the BART package): posterior mean test
# probabilities against BART::lbart / BART::pbart under exactly matched
# priors (leaf sd qlogis(0.975) / (k sqrt(ntree)) for logistic - lbart's
# tau.interval = 0.95 convention - and 3 / (k sqrt(ntree)) for probit,
# k fixed at 2, tree prior 0.95 / (1 + depth)^2, 100 uniform cuts, fixed
# centering offsets). This part does not gate: on the part-1 problem lbart
# itself deviates from the exact posterior by ~0.03 in the anti-shrinkage
# direction (bartcore matches to MC error), and that bias reappears here
# as a ~0.02 mean probability difference; pbart agrees to ~0.005.
#
# Usage: Rscript logistic-reference.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

# ---- part 1: exact single-tree posteriors ----

exactNdpost <- if (quick) 10000L else 50000L
exactSeeds <- if (quick) 1L else 3L
exactTolerance <- if (quick) 0.01 else 0.005  # vs single-chain MC error

set.seed(999L)
nPerCell <- 25L
cellCenters <- c(0.125, 0.375, 0.625, 0.875)
x1 <- matrix(rep(cellCenters, each = nPerCell) +
               runif(4L * nPerCell, -0.1, 0.1), ncol = 1L)
cell <- rep(1:4, each = nPerCell)
y1 <- rbinom(length(cell), 1L, c(0.2, 0.4, 0.6, 0.8)[cell])
x1.test <- matrix(cellCenters, ncol = 1L)

# dbarts' uniform cut rule with n.cuts = 3; each cut separates two cells
cuts <- min(x1) + (1:3) * (max(x1) - min(x1)) / 4
stopifnot(identical(cell, findInterval(x1[, 1L], cuts) + 1L))

base <- 0.5
power <- 2
k <- 2

successesByCell <- as.vector(tapply(y1, cell, sum))
countsByCell <- as.vector(table(cell))

# CGM prior over trees on cut interval [loCut, hiCut]: leaf w.p.
# 1 - base / (1 + depth)^power (forced once the interval empties), else a
# uniformly drawn cut splits the cell range. Leaves are contiguous cell
# ranges.
enumerate <- function(loCell, hiCell, loCut, hiCut, depth) {
  growth <- if (hiCut >= loCut) base / (1 + depth)^power else 0
  result <- list(list(leaves = list(c(loCell, hiCell)),
                      logPrior = log(1 - growth)))
  if (hiCut < loCut) return(result)
  for (j in loCut:hiCut) {
    lefts <- enumerate(loCell, j, loCut, j - 1L, depth + 1L)
    rights <- enumerate(j + 1L, hiCell, j + 1L, hiCut, depth + 1L)
    for (left in lefts) for (right in rights) {
      result[[length(result) + 1L]] <- list(
        leaves = c(left$leaves, right$leaves),
        logPrior = log(growth) - log(hiCut - loCut + 1) +
          left$logPrior + right$logPrior
      )
    }
  }
  result
}
trees <- enumerate(1L, 4L, 1L, 3L, 0L)

exactPredictive <- function(linkLogProbability, linkProbability, offset, tau) {
  muGrid <- seq(-10, 10, by = 0.005)
  muDensity <- dnorm(muGrid, 0, tau)
  leafQuantities <- function(s, m) {
    w <- muDensity * exp(s * linkLogProbability(muGrid + offset) +
                           (m - s) * linkLogProbability(-(muGrid + offset)))
    list(logMarginal = log(sum(w) * 0.005),
         meanProbability = sum(linkProbability(muGrid + offset) * w) / sum(w))
  }
  logWeights <- numeric(length(trees))
  predictions <- matrix(0, length(trees), 4L)
  for (t in seq_along(trees)) {
    logWeight <- trees[[t]]$logPrior
    for (leaf in trees[[t]]$leaves) {
      cells <- leaf[1L]:leaf[2L]
      q <- leafQuantities(sum(successesByCell[cells]), sum(countsByCell[cells]))
      logWeight <- logWeight + q$logMarginal
      predictions[t, cells] <- q$meanProbability
    }
    logWeights[t] <- logWeight
  }
  weights <- exp(logWeights - max(logWeights))
  colSums(weights * predictions) / sum(weights)
}

fitSingleTree <- function(seed, family, nodeScale, offset, linkinv) {
  set.seed(seed)
  control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 1L,
                           n.cuts = 3L, updateState = FALSE)
  sampler <- dbarts(x1, y1, test = x1.test, offset = offset,
                    control = control, node.prior = normal(k),
                    tree.prior = cgm(power, base))
  sampler$model@node.scale <- nodeScale
  sampler$data@offset.test <- NULL  # engine fits exclude the offset
  bc <- dbarts:::bartcoreSampler(sampler, family = family)
  r <- dbarts:::bartcoreRun(bc, 5000L, exactNdpost)
  colMeans(linkinv(t(r$test) + offset))
}

checkExact <- function(name, family, nodeScale, offset,
                       linkLogProbability, linkinv) {
  exact <- exactPredictive(linkLogProbability, linkinv, offset,
                           nodeScale / k)
  fit <- colMeans(do.call(rbind, lapply(
    seq_len(exactSeeds), fitSingleTree,
    family = family, nodeScale = nodeScale, offset = offset,
    linkinv = linkinv
  )))
  gap <- max(abs(fit - exact))
  cat(sprintf("%-10s exact %s | sampler %s | max gap %.4f%s\n",
              name,
              paste(sprintf("%.4f", exact), collapse = " "),
              paste(sprintf("%.4f", fit), collapse = " "),
              gap, if (gap > exactTolerance) " <- FAIL" else ""))
  gap > exactTolerance
}

cat("part 1: single-tree posteriors vs exact enumeration\n")
anyFailure <- checkExact(
  "logistic", "logistic", qlogis(0.975), qlogis(mean(y1)),
  function(q) plogis(q, log.p = TRUE), plogis
)
anyFailure <- checkExact(
  "probit", "probit", 3.0, qnorm(mean(y1)),
  function(q) pnorm(q, log.p = TRUE), pnorm
) || anyFailure

# ---- part 2: BART package comparison (informational) ----

if (requireNamespace("BART", quietly = TRUE)) {
  n.seeds <- if (quick) 3L else 20L
  ndpost  <- if (quick) 250L else 1000L
  nskip   <- if (quick) 100L else 500L
  ntree   <- if (quick) 50L else 200L
  n.test  <- 25L

  set.seed(5109L)
  n <- 500L
  p <- 5L
  x <- matrix(runif(n * p), n)
  x.test <- matrix(runif(n.test * p), n.test)
  linearPredictor <-
    (10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
       10 * x[, 4L] + 5 * x[, 5L] - 14) / 4
  y <- rbinom(n, 1L, plogis(linearPredictor))

  fitBartcore <- function(seed, family, nodeScale, offset, linkinv) {
    set.seed(seed)
    control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = ntree,
                             updateState = FALSE)
    sampler <- dbarts(x, y, test = x.test, offset = offset,
                      control = control, node.prior = normal(k))
    sampler$model@node.scale <- nodeScale
    sampler$data@offset.test <- NULL
    bc <- dbarts:::bartcoreSampler(sampler, family = family)
    r <- dbarts:::bartcoreRun(bc, nskip, ndpost)
    colMeans(linkinv(t(r$test) + offset))
  }

  fitBART <- function(seed, fitter, offset) {
    set.seed(seed)
    log <- capture.output(
      r <- fitter(x, y, x.test = x.test, k = k, ntree = ntree, numcut = 100L,
                  ndpost = ndpost, nskip = nskip, binaryOffset = offset)
    )
    colMeans(r$prob.test)
  }

  report <- function(name, a, b) {
    cat(sprintf(
      "%-24s mean |prob difference| = %.4f, max = %.4f (base rate %.2f)\n",
      name, mean(abs(colMeans(a) - colMeans(b))),
      max(abs(colMeans(a) - colMeans(b))), mean(y)
    ))
  }

  cat("\npart 2: BART package comparison (informational; see header)\n")
  offset <- qlogis(mean(y))
  a <- do.call(rbind, lapply(seq_len(n.seeds), fitBartcore,
                             family = "logistic", nodeScale = qlogis(0.975),
                             offset = offset, linkinv = plogis))
  b <- do.call(rbind, lapply(seq_len(n.seeds), fitBART,
                             fitter = BART::lbart, offset = offset))
  report("logistic vs BART::lbart", a, b)

  offset <- qnorm(mean(y))
  a <- do.call(rbind, lapply(seq_len(n.seeds), fitBartcore,
                             family = "probit", nodeScale = 3.0,
                             offset = offset, linkinv = pnorm))
  b <- do.call(rbind, lapply(seq_len(n.seeds), fitBART,
                             fitter = BART::pbart, offset = offset))
  report("probit vs BART::pbart", a, b)
} else {
  cat("\npart 2 skipped: BART package not installed\n")
}

if (anyFailure) quit(status = 1L)
cat("\nOK: samplers match the exact posterior\n")
