#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore ordinal (cumulative-probit) sampler
# (docs/design/ordinal.md section 7). Single tree, K = 3, small n over one
# binary categorical predictor (two cells): the tree space is exactly two
# structures - a shared root leaf, or a split into one leaf per cell - so the
# posterior is a closed-form quadrature.
#
# The key simplification: the cumulative-probit likelihood MARGINALIZES the
# latents analytically (section 1 is a Phi-difference), so the reference is
# z-FREE and integrates only over the leaf means mu and the ONE free cutpoint
# gamma_2 (gamma_1 = 0 pinned at K = 3). Agreement therefore validates BOTH the
# truncated-normal augmentation AND the cutpoint Metropolis update at once.
#
# Three details keep the target the sampler's ACTUAL posterior (section 7):
#   - the free cutpoint prior is the shipped log-gap prior, normal on
#     delta = log(gamma_2 - gamma_1) = log gamma_2 with sd 1.5; integrating in
#     delta space (gamma_2 = exp(delta), weight dnorm(delta, 0, 1.5)) IS the
#     1/gamma_2 Jacobian pushforward - the gate integrates the real prior, not a
#     flat stand-in;
#   - a leaf holding no observations contributes a marginal factor of 1 (its
#     empty likelihood product integrates its mu prior to unity) - not exercised
#     here (both cells are occupied) but the cellIntegral of a zero-count cell
#     returns 1 by construction;
#   - each structure's posterior weight is its tree prior TIMES its computed
#     marginal, renormalized over the two structures. The two cells share
#     gamma_2, so under the split they are conditionally independent given
#     gamma_2: an outer gamma_2 quadrature nesting an inner per-cell mu one.
#
# The sampler side needs per-draw cutpoints, which the engine exposes only
# through its state block (no run-output channel), so - as bart2's ordinal fit
# does - the run is driven one kept sample at a time and gamma_2 read from the
# state after each sweep, paired with that sweep's leaf latents. A single tree
# makes the per-draw state serialization cheap.
#
# Tolerances bound sampler MC plus quadrature error; never widen one to pass -
# a failure means the augmentation or the cutpoint update is wrong (perturbing
# the log-gap prior sd or the proposal scale in the engine moves the posterior
# and this gate catches it).
#
# Usage: Rscript ordinal-exact.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

ndpost <- if (quick) 20000L else 60000L
nburn <- 5000L
nSeeds <- if (quick) 2L else 3L
tolProb <- if (quick) 0.020 else 0.012 # category probabilities
tolGamma <- if (quick) 0.060 else 0.040 # the free cutpoint gamma_2

# ---- shipped constants (docs/design/ordinal.md section 3, "Shipped (C1)") ----

K <- 3L
k <- 2
numTrees <- 1L
nodeScale <- 3.0
tau <- nodeScale / (k * sqrt(numTrees)) # leaf-prior sd, 1.5
deltaSd <- 1.5 # log-gap prior N(0, 1.5^2)
power <- 2.0
base <- 0.95

# ---- fixed data: K = 3, two cells, deliberately imbalanced ----

set.seed(20260718L)
nPerCell <- 30L
cellProb <- list(c(0.60, 0.30, 0.10), c(0.15, 0.35, 0.50))
cell <- rep(0:1, each = nPerCell)
labels <- integer(length(cell))
for (c in 0:1) {
  idx <- which(cell == c)
  labels[idx] <- sample.int(
    K,
    length(idx),
    replace = TRUE,
    prob = cellProb[[c + 1L]]
  )
}
countsA <- tabulate(labels[cell == 0L], K)
countsB <- tabulate(labels[cell == 1L], K)
x <- matrix(as.double(cell), ncol = 1L)
y <- ordered(c("lo", "mid", "hi")[labels], levels = c("lo", "mid", "hi"))

# ---- exact posterior by structure enumeration + nested quadrature ----

muGrid <- seq(-8, 8, by = 0.02)
dmu <- muGrid[2L] - muGrid[1L]
wMu <- dnorm(muGrid, 0, tau) * dmu

deltaGrid <- seq(-9, 9, length.out = 601L)
dd <- deltaGrid[2L] - deltaGrid[1L]
wDelta <- dnorm(deltaGrid, 0, deltaSd) * dd
gamma2Grid <- exp(deltaGrid)

# For a cell with category counts cnt = (n1, n2, n3) and a fixed gamma_2,
# integrate over the leaf mean mu: g = marginal likelihood, h[k] = its k-th
# predictive moment E[P(y = k)]. P(y = 1) = Phi(-mu), P(y = 2) = Phi(gamma_2 -
# mu) - Phi(-mu), P(y = 3) = 1 - Phi(gamma_2 - mu). A zero-count cell returns
# g = 1 (the mu prior integrates to unity).
cellIntegral <- function(cnt, gamma2) {
  p1 <- pnorm(-muGrid)
  c2 <- pnorm(gamma2 - muGrid)
  p2 <- c2 - p1
  p3 <- 1 - c2
  loglik <- cnt[1L] *
    log(pmax(p1, 1e-300)) +
    cnt[2L] * log(pmax(p2, 1e-300)) +
    cnt[3L] * log(pmax(p3, 1e-300))
  wl <- wMu * exp(loglik)
  list(
    g = sum(wl),
    h = c(sum(wl * p1), sum(wl * p2), sum(wl * p3))
  )
}

# root structure: one leaf shared by both cells, so the combined counts drive a
# single mu integral and both cells share the predictive
rootDen <- 0
rootPred <- numeric(K)
rootGamma <- 0
# split structure: independent mu per cell given the shared gamma_2
splitDen <- 0
predA <- numeric(K)
predB <- numeric(K)
splitGamma <- 0
combined <- countsA + countsB
for (i in seq_along(gamma2Grid)) {
  g2 <- gamma2Grid[i]
  cRoot <- cellIntegral(combined, g2)
  rootDen <- rootDen + wDelta[i] * cRoot$g
  rootPred <- rootPred + wDelta[i] * cRoot$h
  rootGamma <- rootGamma + wDelta[i] * g2 * cRoot$g

  cA <- cellIntegral(countsA, g2)
  cB <- cellIntegral(countsB, g2)
  splitDen <- splitDen + wDelta[i] * cA$g * cB$g
  predA <- predA + wDelta[i] * cA$h * cB$g
  predB <- predB + wDelta[i] * cB$h * cA$g
  splitGamma <- splitGamma + wDelta[i] * g2 * cA$g * cB$g
}

# tree prior: a single binary predictor exhausts its one cut, so the children of
# a split are forced leaves - root has prior 1 - base, split has prior base
priorRoot <- 1 - base
priorSplit <- base
den <- priorRoot * rootDen + priorSplit * splitDen
exactA <- (priorRoot * rootPred + priorSplit * predA) / den
exactB <- (priorRoot * rootPred + priorSplit * predB) / den
exactGamma <- (priorRoot * rootGamma + priorSplit * splitGamma) / den

# ---- sampler fit: single tree, per-draw cutpoints from the state ----

fitSeed <- function(seed) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = numTrees,
    updateState = FALSE
  )
  host <- dbarts(
    x,
    y,
    family = "ordinal",
    control = control,
    tree.prior = cgm(power, base),
    node.prior = normal(k),
    verbose = FALSE
  )
  host$data@varTypes[1L] <- 1L # mark the predictor categorical
  bc <- dbarts:::bartcoreSampler(host, family = "ordinal")

  iA <- which(cell == 0L)[1L]
  iB <- which(cell == 1L)[1L]
  probsA <- numeric(K)
  probsB <- numeric(K)
  gammaSum <- 0
  for (s in seq_len(ndpost)) {
    r <- bartcoreRun(bc, if (s == 1L) nburn else 0L, 1L)
    st <- bartcoreStoreState(bc)
    # the K-1 finite cutpoints c(gamma_1 = 0, gamma_2); the probability transform
    # needs the whole vector, the gamma_2 tracking only its free entry
    gammaVec <- st[[1L]]$cutpoints
    gammaSum <- gammaSum + gammaVec[2L]
    etaA <- r$train[iA, 1L]
    etaB <- r$train[iB, 1L]
    probsA <- probsA +
      as.vector(dbarts:::ordinalCategoryProbabilities(etaA, gammaVec))
    probsB <- probsB +
      as.vector(dbarts:::ordinalCategoryProbabilities(etaB, gammaVec))
  }
  c(probsA / ndpost, probsB / ndpost, gammaSum / ndpost)
}

fit <- colMeans(do.call(rbind, lapply(seq_len(nSeeds), fitSeed)))
fitA <- fit[1:K]
fitB <- fit[(K + 1L):(2L * K)]
fitGamma <- fit[2L * K + 1L]

gapProb <- max(abs(c(fitA - exactA, fitB - exactB)))
gapGamma <- abs(fitGamma - exactGamma)

cat("Ordinal exact-posterior gate (single tree, K = 3, two cells):\n")
cat(sprintf(
  "  exact   A  %s\n",
  paste(sprintf("%.4f", exactA), collapse = " ")
))
cat(sprintf("  sampler A  %s\n", paste(sprintf("%.4f", fitA), collapse = " ")))
cat(sprintf(
  "  exact   B  %s\n",
  paste(sprintf("%.4f", exactB), collapse = " ")
))
cat(sprintf("  sampler B  %s\n", paste(sprintf("%.4f", fitB), collapse = " ")))
cat(sprintf(
  "  gamma_2  exact %.4f  sampler %.4f\n",
  exactGamma,
  fitGamma
))
cat(sprintf(
  "  max category-prob gap %.4f (tol %.3f)%s\n",
  gapProb,
  tolProb,
  if (gapProb > tolProb) "  <- FAIL" else ""
))
cat(sprintf(
  "  gamma_2 gap %.4f (tol %.3f)%s\n",
  gapGamma,
  tolGamma,
  if (gapGamma > tolGamma) "  <- FAIL" else ""
))

if (gapProb > tolProb || gapGamma > tolGamma) {
  quit(status = 1L)
}
cat("\nOK: ordinal sampler matches the exact posterior\n")
