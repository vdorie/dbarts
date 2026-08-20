#!/usr/bin/env Rscript

# Exact-posterior gate for the discrete-time hazard family
# (docs/design/survival.md, "Discrete-time hazard", sections 1 and 4). The
# companion reduction gate (hazard-reduction.R) proves the fit equals a binary
# fit on the expanded design DRAW FOR DRAW, but it builds its comparison design
# with the package's OWN expander, so it is blind by construction to a defect in
# the expansion itself: a wrong expander still reduces bitwise to a correct
# probit on the wrong rows. This gate closes that hole - the target is derived
# from the subject-level hazard likelihood, and no person-period expansion
# appears anywhere in the reference.
#
# THE DERIVATION. For K ordered periods the discrete hazard is h(k | x) =
# Phi(f(x, k)), and a subject observed to period t with status s contributes
#
#     L_i = h(t)^s prod_{k < t} (1 - h(k))   for an event (s = 1),
#         = prod_{k <= t} (1 - h(k))         for censoring (s = 0),
#
# so on a design where f is constant within each (covariate cell, period) the
# whole sample's likelihood collapses to ONE Bernoulli factor per cell, with
# sufficient statistics read straight off (x, time, status):
#
#     e[a, k] = #{i in cell a : s_i = 1, t_i = k}   terminal events,
#     R[a, k] = #{i in cell a : t_i >= k}           the period-k risk set.
#
# R is the at-risk count under the same-period convention the design settles on
# (section 1): a subject censored in period t IS at risk in t and contributes
# its (1 - h(t)) factor, so censoring never leaves the risk set early. That one
# count is what an off-by-one expansion gets wrong, and it is the reason this
# gate exists. e likewise pins the terminal-row indicator: only an EVENT's last
# period is a success. Cell membership pins the third piece of the expansion,
# that a subject's covariates are replicated onto all of its rows.
#
# THE CONFIGURATION. Two covariate cells (one predictor, one cut) and K = 3
# periods (two cuts on the appended period column), so both columns are live -
# a constant filler column would NOT be inert, since a non-quantile grid gives
# even a degenerate column its full complement of (repeated) cut points and the
# engine's availability test is index-based, which perturbs the tree prior. The
# single-tree space over the resulting 2 x 3 cell grid is enumerated exactly,
# each node's rule prior being uniform over its available variables and then
# over that variable's remaining cuts (src/bartcore/model.hpp's
# drawSplitVariable / drawRuleForVariable, whose density cancels against the
# birth proposal's). Every cell is occupied, so no enumerated split is empty.
#
# THE FUNCTIONALS. Leaves are conditionally independent given the structure, so
# both are exact:
#
#     E[h(k | x)]  per-cell hazard, a 1-D quadrature over the leaf parameter;
#     E[S(t | x)]  survival, S(t | x) = prod_{k <= t} (1 - h(k | x)), the
#                  original-data-scale quantity, computed as a product of
#                  per-leaf moments E[(1 - Phi(mu))^m] over the leaves the
#                  horizon crosses.
#
# What the survival channel actually DISCRIMINATES, stated honestly: the
# cumulative-product composition of section 4, which leaf each period maps to,
# and the risk sets underneath. It does NOT discriminate the moment refinement
# itself - on this design the leaf posteriors are tight enough that replacing
# E[(1 - h)^m] with the plug-in (E[1 - h])^m still reads a max survival gap of
# 0.0009 against the 0.004 tolerance, so the gate does not fire on it. The
# exact moment is used because it is the correct quantity, not because this
# configuration resolves it.
#
# Only the probit link is gated: both links are exactly gated in
# logistic-reference.R and bitwise in hazard-reduction.R, so the new content
# here is the expansion, which is link-free.
#
# Sensitivity, measured by poisoning the reference rather than the engine, at
# full size (against a clean gate reading 0.0008 / 0.0005): (a) dropping the k
# divisor from the leaf-prior sd (tau = 3.0 rather than node.scale / k = 1.5)
# reads max hazard gap 0.0121 and max survival gap 0.0077; (b) deriving R under
# the OTHER censoring convention - a censored subject at risk only through
# t_i - 1, the off-by-one an expander defect produces - reads max hazard gap
# 0.6909 and max survival gap 0.5375.
# Arm (b) is the one the reduction gate cannot see. Both sit far above the
# tolerances below, which bound sampler MC plus quadrature error only.
#
# Usage: Rscript hazard-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

ndpost <- if (quick) 30000L else 100000L
nburn <- 5000L
nSeeds <- if (quick) 2L else 3L
tolHazard <- if (quick) 0.008 else 0.004 # vs single-chain MC error
tolSurvival <- if (quick) 0.008 else 0.004

# ---- shipped constants ----

K <- 3L # periods
nCells <- 2L # covariate cells
base <- 0.5
power <- 2
k <- 2
nodeScale <- 3.0 # probit node.scale (R/model.R's defaultNodeScale)
tau <- nodeScale / k # one tree, so no sqrt(n.trees) divisor

# ---- fixed data: two covariate cells on a K-period grid ----

set.seed(3141L)
nPerCell <- 200L
n <- nCells * nPerCell
cellOf <- rep(seq_len(nCells), each = nPerCell)
cellValue <- c(0.25, 0.75) # one cut at the midpoint separates them
hazardTrue <- rbind(c(0.15, 0.22, 0.32), c(0.30, 0.38, 0.46))
eventPeriod <- rep(K + 1L, n)
for (i in seq_len(n)) {
  for (j in seq_len(K)) {
    if (runif(1L) < hazardTrue[cellOf[i], j]) {
      eventPeriod[i] <- j
      break
    }
  }
}
censorPeriod <- sample.int(K, n, replace = TRUE)
time <- as.double(pmin(eventPeriod, censorPeriod, K))
status <- as.double(eventPeriod <= censorPeriod & eventPeriod <= K)
x <- matrix(cellValue[cellOf], ncol = 1L)
# the grid is the distinct observed times (the family's default), so it must be
# exactly 1..K for the period column to carry K - 1 usable cuts
stopifnot(identical(sort(unique(time)), as.double(seq_len(K))))

# the two sufficient statistics, straight off (x, time, status)
events <- outer(
  seq_len(nCells),
  seq_len(K),
  Vectorize(function(a, j) {
    sum(cellOf == a & status == 1 & time == j)
  })
)
atRisk <- outer(
  seq_len(nCells),
  seq_len(K),
  Vectorize(function(a, j) {
    sum(cellOf == a & time >= j)
  })
)
stopifnot(all(atRisk > 0)) # every cell occupied: no enumerated split is empty
cat(sprintf(
  "events %s | risk sets %s | censoring rate %.2f\n",
  paste(as.vector(t(events)), collapse = " "),
  paste(as.vector(t(atRisk)), collapse = " "),
  mean(status == 0)
))

# ---- exact posterior over the enumerated single-tree space ----

# CGM prior over trees on the 2-D cut grid: a node holds a cut-index interval
# per variable (cut j of a variable separates its cells <= j from those above),
# grows with probability base / (1 + depth)^power whenever some variable still
# has a cut, and then draws uniformly over the available variables and that
# variable's remaining cuts. Leaves are cell rectangles.
enumerate <- function(cellLo, cellHi, cutLo, cutHi, depth) {
  nCuts <- pmax(cutHi - cutLo + 1L, 0L)
  available <- which(nCuts > 0L)
  growth <- if (length(available) > 0L) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(cellLo[1L], cellHi[1L], cellLo[2L], cellHi[2L])),
    logPrior = log(1 - growth)
  ))
  for (v in available) {
    logRule <- log(growth) - log(length(available)) - log(nCuts[v])
    for (j in cutLo[v]:cutHi[v]) {
      leftCellHi <- cellHi
      leftCellHi[v] <- j
      leftCutHi <- cutHi
      leftCutHi[v] <- j - 1L
      rightCellLo <- cellLo
      rightCellLo[v] <- j + 1L
      rightCutLo <- cutLo
      rightCutLo[v] <- j + 1L
      lefts <- enumerate(cellLo, leftCellHi, cutLo, leftCutHi, depth + 1L)
      rights <- enumerate(rightCellLo, cellHi, rightCutLo, cutHi, depth + 1L)
      for (left in lefts) {
        for (right in rights) {
          result[[length(result) + 1L]] <- list(
            leaves = c(left$leaves, right$leaves),
            logPrior = logRule + left$logPrior + right$logPrior
          )
        }
      }
    }
  }
  result
}
trees <- enumerate(
  c(1L, 1L),
  c(nCells, K),
  c(1L, 1L),
  c(nCells - 1L, K - 1L),
  0L
)
cat(sprintf("enumerated %d trees\n", length(trees)))

muGrid <- seq(-10, 10, by = 0.002)
dmu <- muGrid[2L] - muGrid[1L]
logMuPrior <- dnorm(muGrid, 0, tau, log = TRUE)
hazardOfMu <- pnorm(muGrid)
survivalOfMu <- pnorm(muGrid, lower.tail = FALSE)
logHazardOfMu <- pnorm(muGrid, log.p = TRUE)
logSurvivalOfMu <- pnorm(muGrid, lower.tail = FALSE, log.p = TRUE)

# Posterior quantities for one leaf rectangle, cached by its cell bounds for
# the duration of one exactPosterior call (the enumeration revisits each
# rectangle many times, and the statistics are fixed within a call). The leaf's
# Bernoulli data are the summed events and risk sets over its cells; the
# survival moments
# E[(1 - h)^m], m = 0..K, are what the cumulative product needs, since a
# horizon crosses only part of a leaf's periods.
leafCache <- new.env(parent = emptyenv())
leafQuantities <- function(events, atRisk, cells) {
  key <- paste(cells, collapse = ",")
  hit <- leafCache[[key]]
  if (!is.null(hit)) {
    return(hit)
  }
  rows <- cells[1L]:cells[2L]
  columns <- cells[3L]:cells[4L]
  e <- sum(events[rows, columns])
  r <- sum(atRisk[rows, columns])
  logIntegrand <- logMuPrior + e * logHazardOfMu + (r - e) * logSurvivalOfMu
  m <- max(logIntegrand)
  w <- exp(logIntegrand - m)
  sw <- sum(w)
  value <- list(
    logMarginal = log(sw * dmu) + m,
    meanHazard = sum(hazardOfMu * w) / sw,
    survivalMoment = vapply(0:K, function(j) sum(survivalOfMu^j * w) / sw, 0)
  )
  leafCache[[key]] <- value
  value
}

# Posterior-mean hazards and survival probabilities under the enumerated tree
# space, as a function of the sufficient statistics, which are the only inputs
# either poison arm in the header perturbs.
exactPosterior <- function(events, atRisk) {
  rm(list = ls(leafCache), envir = leafCache)
  logWeights <- numeric(length(trees))
  hazardByTree <- array(0, c(length(trees), nCells, K))
  survivalByTree <- array(1, c(length(trees), nCells, K))
  for (t in seq_along(trees)) {
    logWeight <- trees[[t]]$logPrior
    for (leaf in trees[[t]]$leaves) {
      q <- leafQuantities(events, atRisk, leaf)
      logWeight <- logWeight + q$logMarginal
      rows <- leaf[1L]:leaf[2L]
      periods <- leaf[3L]:leaf[4L]
      hazardByTree[t, rows, periods] <- q$meanHazard
      # each covariate cell in the leaf picks up E[(1 - h)^m], m the leaf's
      # periods at or below the horizon; leaves the horizon misses give 1
      for (j in seq_len(K)) {
        moment <- q$survivalMoment[sum(periods <= j) + 1L]
        survivalByTree[t, rows, j] <- survivalByTree[t, rows, j] * moment
      }
    }
    logWeights[t] <- logWeight
  }
  weights <- exp(logWeights - max(logWeights))
  weights <- weights / sum(weights)
  list(
    hazard = apply(hazardByTree * weights, c(2L, 3L), sum),
    survival = apply(survivalByTree * weights, c(2L, 3L), sum)
  )
}

# ---- sampler fit ----

# both newdata grids are built here, not by the package's expander: the hazard
# grid is one row per (cell, period) with the period column LAST, and the
# survival grid is one row per cell, which survivalProbabilities expands itself
hazardGrid <- cbind(
  rep(cellValue, each = K),
  rep(seq_len(K), times = nCells)
)
profiles <- matrix(cellValue, ncol = 1L)

fitSeed <- function(seed) {
  fit <- bart2(
    x,
    cbind(time, status),
    family = "hazard",
    n.trees = 1L,
    n.cuts = c(nCells - 1L, K - 1L),
    n.burn = nburn,
    n.samples = ndpost,
    n.chains = 1L,
    n.threads = 1L,
    k = k,
    power = power,
    base = base,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = seed
  )
  stopifnot(identical(fit$periods, as.double(seq_len(K))))
  hazard <- colMeans(predict(fit, hazardGrid, type = "ev"))
  # survivalProbabilities reports draws x times x observations, so the mean is
  # a (times, cells) matrix whose column-major order is the cell-major layout
  # the hazard grid above is already built in
  survival <- apply(
    survivalProbabilities(fit, seq_len(K), newdata = profiles),
    c(2L, 3L),
    mean
  )
  c(hazard, as.vector(survival))
}

exact <- exactPosterior(events, atRisk)
fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
fitHazard <- matrix(fit[seq_len(nCells * K)], nCells, K, byrow = TRUE)
fitSurvival <- matrix(fit[-seq_len(nCells * K)], nCells, K, byrow = TRUE)

gapHazard <- max(abs(fitHazard - exact$hazard))
gapSurvival <- max(abs(fitSurvival - exact$survival))

cat("Discrete-time hazard exact-posterior gate (single tree, 2 x 3 cells):\n")
for (a in seq_len(nCells)) {
  cat(sprintf(
    "  cell %d hazard    exact %s | sampler %s\n",
    a,
    paste(sprintf("%.4f", exact$hazard[a, ]), collapse = " "),
    paste(sprintf("%.4f", fitHazard[a, ]), collapse = " ")
  ))
  cat(sprintf(
    "  cell %d survival  exact %s | sampler %s\n",
    a,
    paste(sprintf("%.4f", exact$survival[a, ]), collapse = " "),
    paste(sprintf("%.4f", fitSurvival[a, ]), collapse = " ")
  ))
}
cat(sprintf(
  "  max hazard gap %.4f (tol %.3f)%s\n",
  gapHazard,
  tolHazard,
  if (gapHazard > tolHazard) "  <- FAIL" else ""
))
cat(sprintf(
  "  max survival gap %.4f (tol %.3f)%s\n",
  gapSurvival,
  tolSurvival,
  if (gapSurvival > tolSurvival) "  <- FAIL" else ""
))

if (gapHazard > tolHazard || gapSurvival > tolSurvival) {
  quit(status = 1L)
}
cat("\nOK: the hazard sampler matches the exact subject-level posterior\n")
