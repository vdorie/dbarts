#!/usr/bin/env Rscript

# Exact-posterior gate for where a sweep draws its augmentation variables
# relative to an outer step that reads the fit. The BART sampler sits inside a
# larger one: the offset is gamma * v, gamma has a uniform prior on a grid, and
# after every sweep the outer step draws gamma given the fit and installs the
# new offset with setOffset. Two outer steps are valid for a latent-variable
# family. One conditions on the latent (z for probit, omega for logistic); it
# is exact whenever the latent it reads is a draw given the fit. The other
# integrates the latent out (pnorm, plogis of the fit); it is exact only if the
# state the next sweep starts from carries a latent drawn given the fit and the
# new offset, which holds when the sweep draws the latents before anything
# reads them.
#
# Thirty-two rows in four cells of one ordinal predictor, one tree,
# n.cuts = 3, so the tree space is enumerated under the CGM prior: every tree
# times every gamma on the grid -1, -0.5, ..., 4, the leaves integrated by
# quadrature and k by quadrature under its chi(1.5, 2) prior. Four arms, all
# with k drawn, so the probit rescaling step runs:
#   probit-integrated    probit, outer step integrating z out;
#   probit-conditional   probit, outer step conditioning on z;
#   logistic-integrated  logistic, outer step integrating omega out;
#   logistic-conditional logistic, outer step conditioning on omega.
# Compared, per arm: the long-run P(gamma = g) at each of the 11 grid points
# (Rao-Blackwellized: the mean of the conditional gamma is drawn from) and the
# posterior mean of the fit without the offset at each cell, by batch-means z.
#
# What it checks. Draw placement and nothing else: where in the sweep the
# latents are drawn, against the outer step's read of them. A stale working
# pointer, a missed cache invalidation or a latent split from the parameter
# that shares its draw leaves these arms exact or out of their reach, and is
# for the unit tests.
#
# Sizes. The standard errors are batch-means errors over the pooled batches.
# Quick mode keeps 4 chains of 250,000 sweeps after 2,000 of burn-in; full mode
# 4 chains of 2,000,000. The bound is |z| <= 4.5 on each of an arm's 12
# statistics: the P(gamma) grid points holding at least 1e-4 of the exact mass
# (eight of the eleven; the rest are mostly zero with rare visits, which no
# normal error describes) and the four cell fits. By a union bound with
# Student-t tails a correct sampler fails a fresh seed about once in 2000 quick
# runs; the seeds are fixed, so a
# given build either passes or does not. The sizes were set from the size of
# the effects the arms exist to find (the integrated probit arm's P(gamma = 4)
# is off by 0.03, about 15 standard errors at 400,000 sweeps).
#
# Usage: Rscript outer-step-exact.R [quick] [probit-integrated]
#   [probit-conditional] [logistic-integrated] [logistic-conditional]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
allArms <- c(
  "probit-integrated",
  "probit-conditional",
  "logistic-integrated",
  "logistic-conditional"
)
armNames <- intersect(allArms, args)
if (length(armNames) == 0L) {
  armNames <- allArms
}

numChains <- 4L
batchSize <- if (quick) 5000L else 20000L
numBatches <- if (quick) 50L else 100L
burnIn <- 2000L
engineSeed <- 20261105L
zBound <- 4.5
# a grid point with less exact mass than this is printed and not tested: its
# batch means are mostly zero with rare visits, and no normal error describes
# them (one such point read z = -22 on a correct sampler)
minMass <- 1e-4

# ---- the design ----

set.seed(20261103L)
cells <- rep(1:4, each = 8L)
n <- length(cells)
K <- 4L
v <- rnorm(n)
y <- rbinom(n, 1L, pnorm(c(-0.6, -0.2, 0.3, 0.8)[cells] + v))
sgn <- 2 * y - 1
grid <- seq(-1, 4, by = 0.5)
numGrid <- length(grid)
offsetMatrix <- outer(v, grid)
base <- 0.95
power <- 2
kPrior <- c(df = 1.5, scale = 2)

# ---- tree enumeration ----
#
# Every tree as its leaves, each a contiguous range of cells, and its log
# prior: bd-balance.R's recursion over runs of adjacent cells.

enumerateTrees <- function(loCell, hiCell, depth) {
  growth <- if (hiCell > loCell) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(loCell, hiCell)),
    logPrior = log(1 - growth)
  ))
  for (j in seq_len(hiCell - loCell) + loCell - 1L) {
    for (left in enumerateTrees(loCell, j, depth + 1L)) {
      for (right in enumerateTrees(j + 1L, hiCell, depth + 1L)) {
        result[[length(result) + 1L]] <- list(
          leaves = c(left$leaves, right$leaves),
          logPrior = log(growth) -
            log(hiCell - loCell) +
            left$logPrior +
            right$logPrior
        )
      }
    }
  }
  result
}

# ---- the exact posterior ----
#
# p(gamma, tree, k, leaves | y) is proportional to the tree prior, the chi
# prior of k, and for each leaf its N(0, (scale / k)^2) prior times the
# likelihood of its rows at the leaf value plus gamma * v. Leaves are
# independent given the tree, k and gamma, so each is a one-dimensional
# integral, taken by the trapezoid rule on a fine grid in the leaf value (the
# integrand is smooth and negligible at the ends, every leaf holding both
# classes), and k by the trapezoid rule in log k.

exactPosterior <- function(linkLogProbability, kScale) {
  muStep <- 0.01
  muGrid <- seq(-40, 40, by = muStep)
  logK <- seq(log(0.01), log(60), length.out = 600L)
  kValues <- exp(logK)
  u <- kValues / kPrior[["scale"]]
  logKPrior <- (kPrior[["df"]] - 1) * log(u) - 0.5 * u^2 + log(kValues)
  # density of the leaf value under each k, one row a k
  prior <- t(vapply(
    kValues,
    function(k) dnorm(muGrid, 0, kScale / k) * muStep,
    muGrid
  ))
  trees <- enumerateTrees(1L, K, 0L)
  ranges <- unique(unlist(lapply(trees, `[[`, "leaves"), recursive = FALSE))
  rangeKey <- vapply(ranges, paste, "", collapse = "-")

  # per leaf range and gamma: the log marginal and the posterior mean of the
  # leaf value, over k
  logMarginal <- array(0, c(length(ranges), length(logK), numGrid))
  leafMean <- array(0, c(length(ranges), length(logK), numGrid))
  for (r in seq_along(ranges)) {
    rows <- cells >= ranges[[r]][1L] & cells <= ranges[[r]][2L]
    for (g in seq_len(numGrid)) {
      logLik <- numeric(length(muGrid))
      for (i in which(rows)) {
        logLik <- logLik +
          linkLogProbability(sgn[i] * (muGrid + grid[g] * v[i]))
      }
      top <- max(logLik)
      lik <- exp(logLik - top)
      mass <- as.vector(prior %*% lik)
      logMarginal[r, , g] <- top + log(mass)
      leafMean[r, , g] <- as.vector(prior %*% (lik * muGrid)) / mass
    }
  }

  # sum over trees; fit[c, g] and gammaMass[g] accumulate unnormalized
  logWeights <- vector("list", length(trees))
  for (t in seq_along(trees)) {
    kinds <- match(
      vapply(trees[[t]]$leaves, paste, "", collapse = "-"),
      rangeKey
    )
    total <- matrix(trees[[t]]$logPrior + logKPrior, length(logK), numGrid)
    for (kind in kinds) {
      total <- total + logMarginal[kind, , ]
    }
    logWeights[[t]] <- list(kinds = kinds, total = total)
  }
  top <- max(vapply(logWeights, function(w) max(w$total), 0))
  gammaMass <- numeric(numGrid)
  fit <- matrix(0, K, numGrid)
  for (t in seq_along(trees)) {
    w <- exp(logWeights[[t]]$total - top)
    gammaMass <- gammaMass + colSums(w)
    for (j in seq_along(trees[[t]]$leaves)) {
      range <- trees[[t]]$leaves[[j]]
      fitHere <- colSums(w * leafMean[logWeights[[t]]$kinds[j], , ])
      for (cell in range[1L]:range[2L]) {
        fit[cell, ] <- fit[cell, ] + fitHere
      }
    }
  }
  list(
    numTrees = length(trees),
    gamma = gammaMass / sum(gammaMass),
    fit = rowSums(fit) / sum(gammaMass)
  )
}

# ---- the combined sampler ----

makeSampler <- function(family, seed) {
  dbarts(
    matrix(as.double(cells), n, 1L),
    y,
    family = family,
    n.trees = 1L,
    n.cuts = 3L,
    power = power,
    base = base,
    offset = grid[1L] * v,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE,
      seed = seed
    )
  )
}

# log conditional of gamma over the grid, given the fit f (without the
# offset) and the latent latent
logConditional <- function(arm, f, latent) {
  switch(
    arm,
    "probit-integrated" = colSums(pnorm(
      sgn * (f + offsetMatrix),
      log.p = TRUE
    )),
    "logistic-integrated" = colSums(plogis(
      sgn * (f + offsetMatrix),
      log.p = TRUE
    )),
    "probit-conditional" = -0.5 * colSums((latent - f - offsetMatrix)^2),
    "logistic-conditional" = -0.5 *
      colSums(latent * ((y - 0.5) / latent - f - offsetMatrix)^2)
  )
}

# batch means of the conditional P(gamma = g) and of the fit at each cell, one
# row per batch; a chain's burn-in sweeps are dropped
runChain <- function(arm, seed) {
  family <- sub("-.*", "", arm)
  needsLatent <- grepl("conditional", arm)
  sampler <- makeSampler(family, seed)
  set.seed(seed + 1000000L)
  # the bridge entries $run, $setOffset and the reads call, without the
  # per-call argument handling, which at one sweep a call is most of the run
  # time
  pointer <- sampler$getPointer()
  keepFits <- sampler$control@keepFits
  runEntry <- dbarts:::C_dbarts_bartcore_run
  offsetEntry <- dbarts:::C_dbarts_bartcore_setOffset
  fitEntry <- dbarts:::C_dbarts_bartcore_getFitsWithoutOffset
  latentEntry <- dbarts:::C_dbarts_bartcore_getLatents
  firstRow <- match(seq_len(K), cells)
  latent <- NULL

  g <- 1L
  step <- function() {
    .Call(offsetEntry, pointer, grid[g] * v, FALSE)
    .Call(runEntry, pointer, 0L, 1L, NULL, NULL, keepFits)
    f <- as.numeric(.Call(fitEntry, pointer))
    if (needsLatent) {
      latent <<- as.numeric(.Call(latentEntry, pointer, NULL))
    }
    lp <- logConditional(arm, f, latent)
    p <- exp(lp - max(lp))
    p <- p / sum(p)
    g <<- sample.int(numGrid, 1L, prob = p)
    list(p = p, fit = f[firstRow])
  }

  for (i in seq_len(burnIn)) {
    step()
  }
  pBatch <- matrix(0, numBatches, numGrid)
  fitBatch <- matrix(0, numBatches, K)
  for (b in seq_len(numBatches)) {
    pSum <- numeric(numGrid)
    fitSum <- numeric(K)
    for (i in seq_len(batchSize)) {
      s <- step()
      pSum <- pSum + s$p
      fitSum <- fitSum + s$fit
    }
    pBatch[b, ] <- pSum / batchSize
    fitBatch[b, ] <- fitSum / batchSize
  }
  list(gamma = pBatch, fit = fitBatch)
}

# ---- verdict ----

compare <- function(label, batches, exact, tested = rep(TRUE, length(exact))) {
  estimate <- colMeans(batches)
  se <- apply(batches, 2L, sd) / sqrt(nrow(batches))
  z <- (estimate - exact) / se
  cat(sprintf("  %s\n", label))
  cat(sprintf("    %9s %9s %9s %7s\n", "exact", "sampler", "MCse", "z"))
  for (j in seq_along(exact)) {
    cat(sprintf(
      "    %9.5f %9.5f %9.5f %+7.1f%s\n",
      exact[j],
      estimate[j],
      se[j],
      z[j],
      if (!tested[j]) {
        "  (not tested)"
      } else if (is.na(z[j]) || abs(z[j]) > zBound) {
        " <- FAIL"
      } else {
        ""
      }
    ))
  }
  z[!tested] <- 0
  cat(sprintf("    largest |z| %.1f\n", max(abs(z))))
  list(failed = anyNA(z) || any(abs(z) > zBound), z = max(abs(z)))
}

exactCache <- list()
anyFailure <- FALSE
for (armIndex in seq_along(armNames)) {
  arm <- armNames[armIndex]
  family <- sub("-.*", "", arm)
  if (is.null(exactCache[[family]])) {
    kScale <- makeSampler(family, 1L)$getLeafPrior()$k.scale
    exactCache[[family]] <- exactPosterior(
      if (family == "probit") {
        function(q) pnorm(q, log.p = TRUE)
      } else {
        function(q) plogis(q, log.p = TRUE)
      },
      kScale
    )
  }
  exact <- exactCache[[family]]
  started <- proc.time()[3L]
  chains <- lapply(
    engineSeed + 100L * match(arm, allArms) + seq_len(numChains),
    function(seed) runChain(arm, seed)
  )
  cat(sprintf(
    "\n%s arm: %d trees enumerated; %d chains x %d batches x %d sweeps in %.0f s\n",
    arm,
    exact$numTrees,
    numChains,
    numBatches,
    batchSize,
    proc.time()[3L] - started
  ))
  gammaResult <- compare(
    "P(gamma = g | y), by grid point",
    do.call(rbind, lapply(chains, function(chain) chain$gamma)),
    exact$gamma,
    exact$gamma >= minMass
  )
  fitResult <- compare(
    "posterior mean of the fit without the offset, by cell",
    do.call(rbind, lapply(chains, function(chain) chain$fit)),
    exact$fit
  )
  armFailed <- gammaResult$failed || fitResult$failed
  cat(sprintf(
    "  %s: largest |z| %.1f%s\n",
    arm,
    max(gammaResult$z, fitResult$z),
    if (armFailed) " <- FAIL" else ""
  ))
  anyFailure <- anyFailure || armFailed
}

if (anyFailure) {
  cat("\nFAIL: the combined sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("\nOK: the combined sampler matches the exact posterior\n")
