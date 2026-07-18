#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore negative-binomial sampler
# (docs/design/negative-binomial.md section 6). Single tree, small n over one
# binary categorical predictor (two cells): the tree space is exactly two
# structures - a shared root leaf, or a split into one leaf per cell - so the
# posterior is a closed-form quadrature.
#
# The key simplification: the negative-binomial likelihood is CLOSED FORM in
# (leaf log-odds mu, dispersion r) - the Polya-Gamma augmentation omega
# integrates out - so the reference is omega-FREE. It integrates over the leaf
# mean mu and, in the estimated arm, sums over the r GRID (a finite sum, cleaner
# than a continuous integral). Agreement therefore validates the PG mean
# augmentation, the grid r update, AND their composition into the sweep (the
# section 5 r-first ordering; an invalid scan shifts the stationary law, which
# this gate can see).
#
# Two details keep the target the sampler's ACTUAL posterior:
#   - the r grid and its prior weights are the shipped ones (NBDispersionPrior):
#     grid {1,2,3,4,5,6,8,10,12,15,20,30,50}, weights the gamma(2, 0.1) kernel
#     r exp(-0.1 r) renormalized over the grid;
#   - each structure's posterior weight is its tree prior TIMES its computed
#     marginal, renormalized over the two structures. The two cells share r, so
#     under the split they are conditionally independent given r: a per-r sum
#     nesting an inner per-cell mu quadrature.
#
# BOTH modes are gated: the estimated arm (the grid r posterior AND the mean
# counts) and a fixed-r arm (r pinned, mean counts). r is read from the engine's
# state block (no run-output channel), so - as bart2's nbinom fit does - the run
# is driven one kept sample at a time and r read from the state after each sweep.
#
# STATED LIMITATION: fork (A) draws only integer-shape (exact) PG variates, so
# this gate exercises NO approximate path; the reference is omega-free.
#
# Tolerances bound sampler MC plus quadrature error; never widen one to pass.
#
# Usage: Rscript negbin-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

ndpost <- if (quick) 12000L else 30000L
nburn <- 4000L
nSeeds <- 2L
tolMean <- if (quick) 0.12 else 0.07 # mean counts mu = r exp(mu_leaf)
tolGrid <- if (quick) 0.045 else 0.025 # the grid r posterior distribution

# ---- shipped constants (docs/design/negative-binomial.md sections 1, 3) ----

k <- 2
numTrees <- 1L
nodeScale <- pi * sqrt(3) # nbinom reuses logistic's node.scale
tau <- nodeScale / (k * sqrt(numTrees)) # leaf-prior sd, pi sqrt(3) / 2
power <- 2.0
base <- 0.95
rFixed <- 5 # the fixed-r arm's pinned dispersion (a grid member)

grid <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
priorKernel <- grid * exp(-0.1 * grid) # gamma(2, 0.1) kernel
priorW <- priorKernel / sum(priorKernel) # renormalized grid prior

# ---- fixed data: two cells, deliberately different means, small counts ----

set.seed(20260718L)
nPerCell <- 25L
cell <- rep(0:1, each = nPerCell)
cntA <- rnbinom(nPerCell, size = 5L, mu = 1.5)
cntB <- rnbinom(nPerCell, size = 5L, mu = 4.0)
y <- as.double(c(cntA, cntB))
x <- matrix(as.double(cell), ncol = 1L)

# ---- exact posterior by structure enumeration + nested quadrature ----

muGrid <- seq(-10, 4, by = 0.01)
dmu <- muGrid[2L] - muGrid[1L]
wMu <- dnorm(muGrid, 0, tau) * dmu
logP <- -log1p(exp(-muGrid)) # log plogis(mu)
log1mP <- -log1p(exp(muGrid)) # log(1 - plogis(mu))
meanByMu <- exp(muGrid) # mu = r exp(mu_leaf), the r factor applied per r

# For a cell with integer counts cnt and dispersion r, integrate over the leaf
# mean mu: g = marginal likelihood, mc = its mean-count numerator E[r exp(mu)].
# An empty cell (length 0) returns g = 1 (the mu prior integrates to unity).
cellIntegral <- function(cnt, r) {
  nobs <- length(cnt)
  if (nobs == 0L) {
    return(list(g = 1, mc = 0))
  }
  # the mu-independent combinatorial term lgamma(y+r) - lgamma(r) - lgamma(y+1)
  combTerm <- sum(lgamma(cnt + r)) - nobs * lgamma(r) - sum(lgamma(cnt + 1))
  loglik <- combTerm + sum(cnt) * logP + nobs * r * log1mP
  wl <- wMu * exp(loglik)
  list(g = sum(wl), mc = r * sum(wl * meanByMu))
}

combined <- c(cntA, cntB)
nGrid <- length(grid)
rootG <- numeric(nGrid)
rootMC <- numeric(nGrid)
gA <- numeric(nGrid)
gB <- numeric(nGrid)
mcA <- numeric(nGrid)
mcB <- numeric(nGrid)
for (ki in seq_len(nGrid)) {
  r <- grid[ki]
  cr <- cellIntegral(combined, r)
  rootG[ki] <- cr$g
  rootMC[ki] <- cr$mc
  ca <- cellIntegral(cntA, r)
  cb <- cellIntegral(cntB, r)
  gA[ki] <- ca$g
  gB[ki] <- cb$g
  mcA[ki] <- ca$mc
  mcB[ki] <- cb$mc
}

# tree prior: a single binary predictor exhausts its one cut, so a split's
# children are forced leaves - root has prior 1 - base, split has prior base
priorRoot <- 1 - base
priorSplit <- base

# ---- estimated arm: sum over the r grid ----

rootMass <- priorRoot * priorW * rootG # per grid point
splitMass <- priorSplit * priorW * gA * gB
den <- sum(rootMass) + sum(splitMass)
gridPost <- (rootMass + splitMass) / den
# root's mean count is shared by both cells; split's is per-cell (the other
# leaf's marginal integrates out to its g factor)
mcRootNum <- priorRoot * sum(priorW * rootMC)
exactMeanA <- (mcRootNum + priorSplit * sum(priorW * mcA * gB)) / den
exactMeanB <- (mcRootNum + priorSplit * sum(priorW * mcB * gA)) / den

# ---- fixed arm: r pinned at rFixed, no grid sum ----

kf <- which(grid == rFixed)
denF <- priorRoot * rootG[kf] + priorSplit * gA[kf] * gB[kf]
mcRootNumF <- priorRoot * rootMC[kf]
exactFixedA <- (mcRootNumF + priorSplit * mcA[kf] * gB[kf]) / denF
exactFixedB <- (mcRootNumF + priorSplit * mcB[kf] * gA[kf]) / denF

# ---- sampler fits: single tree, per-draw r and mean counts from the state ----

fitSeed <- function(seed, dispersion) {
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
    family = "nbinom",
    dispersion = dispersion,
    control = control,
    tree.prior = cgm(power, base),
    node.prior = normal(k),
    verbose = FALSE
  )
  host$data@varTypes[1L] <- 1L # mark the predictor categorical
  bc <- dbarts:::bartcoreSampler(host, family = "nbinom")

  iA <- which(cell == 0L)[1L]
  iB <- which(cell == 1L)[1L]
  meanA <- 0
  meanB <- 0
  gridCounts <- numeric(nGrid)
  for (s in seq_len(ndpost)) {
    r <- dbarts:::bartcoreRun(bc, if (s == 1L) nburn else 0L, 1L)
    st <- dbarts:::bartcoreStoreState(bc)
    rDraw <- st[[1L]]$dispersion
    gridCounts[match(rDraw, grid)] <- gridCounts[match(rDraw, grid)] + 1
    meanA <- meanA + dbarts:::negbinMeanCounts(r$train[iA, 1L], rDraw)
    meanB <- meanB + dbarts:::negbinMeanCounts(r$train[iB, 1L], rDraw)
  }
  c(meanA / ndpost, meanB / ndpost, gridCounts / ndpost)
}

runArm <- function(dispersion) {
  rows <- do.call(
    rbind,
    lapply(seq_len(nSeeds), function(sd) fitSeed(sd, dispersion))
  )
  colMeans(rows)
}

# estimated arm
est <- runArm(NA_real_)
fitMeanA <- est[1L]
fitMeanB <- est[2L]
fitGrid <- est[-(1:2)]

# fixed arm
fx <- runArm(rFixed)
fitFixedA <- fx[1L]
fitFixedB <- fx[2L]

gapMeanEst <- max(abs(c(fitMeanA - exactMeanA, fitMeanB - exactMeanB)))
gapGrid <- max(abs(fitGrid - gridPost))
gapMeanFixed <- max(abs(c(fitFixedA - exactFixedA, fitFixedB - exactFixedB)))

cat("Negative-binomial exact-posterior gate (single tree, two cells):\n")
cat("--- estimated r (grid posterior) ---\n")
cat(sprintf(
  "  mean count A  exact %.4f  sampler %.4f\n",
  exactMeanA,
  fitMeanA
))
cat(sprintf(
  "  mean count B  exact %.4f  sampler %.4f\n",
  exactMeanB,
  fitMeanB
))
cat("  grid r posterior (r: exact / sampler):\n")
for (ki in seq_len(nGrid)) {
  if (gridPost[ki] > 0.005 || fitGrid[ki] > 0.005) {
    cat(sprintf(
      "    r = %-2g  %.4f / %.4f\n",
      grid[ki],
      gridPost[ki],
      fitGrid[ki]
    ))
  }
}
cat(sprintf(
  "  max mean-count gap %.4f (tol %.3f)%s\n",
  gapMeanEst,
  tolMean,
  if (gapMeanEst > tolMean) "  <- FAIL" else ""
))
cat(sprintf(
  "  max grid-posterior gap %.4f (tol %.3f)%s\n",
  gapGrid,
  tolGrid,
  if (gapGrid > tolGrid) "  <- FAIL" else ""
))
cat(sprintf("--- fixed r = %g ---\n", rFixed))
cat(sprintf(
  "  mean count A  exact %.4f  sampler %.4f\n",
  exactFixedA,
  fitFixedA
))
cat(sprintf(
  "  mean count B  exact %.4f  sampler %.4f\n",
  exactFixedB,
  fitFixedB
))
cat(sprintf(
  "  max mean-count gap %.4f (tol %.3f)%s\n",
  gapMeanFixed,
  tolMean,
  if (gapMeanFixed > tolMean) "  <- FAIL" else ""
))

if (gapMeanEst > tolMean || gapGrid > tolGrid || gapMeanFixed > tolMean) {
  quit(status = 1L)
}
cat("\nOK: negative-binomial sampler matches the exact posterior\n")
