#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore negative-binomial sampler
# (docs/design/negative-binomial.md section 6). Single tree, small n over one
# binary categorical predictor (two cells): the tree space is exactly two
# structures - a shared root leaf, or a split into one leaf per cell - so the
# posterior is a closed-form quadrature.
#
# The forest models the log mean: a cell's log mean (offset aside) is
# m = leaf + c, with c = log(max(sum y, 1/2) / sum exp(o)) the data transform,
# so m ~ N(c, tau^2) a priori, tau = A / (k sqrt(numTrees)). The
# negative-binomial likelihood is CLOSED FORM in (m, r) - the Polya-Gamma
# augmentation omega integrates out:
#   sum lgamma(y + r) - n lgamma(r) - sum lgamma(y + 1) + sum y (m + o)
#     + n r log r - sum (y + r) log(r + exp(m + o)),
# so the reference is omega-FREE. It integrates over m on a grid centred on c
# and, in the estimated arm, sums over the r GRID. Agreement therefore
# validates the PG mean augmentation, the grid r update given the means, AND
# their composition into the sweep (an invalid scan shifts the stationary law,
# which this gate can see).
#
# Two details keep the target the sampler's ACTUAL posterior:
#   - the r grid and its prior weights are the shipped ones (NBDispersionPrior):
#     grid {1,2,3,4,5,6,8,10,12,15,20,30,50}, weights the gamma(2, 0.1) kernel
#     r exp(-0.1 r) renormalized over the grid;
#   - each structure's posterior weight is its tree prior TIMES its computed
#     marginal, renormalized over the two structures. The two cells share r, so
#     under the split they are conditionally independent given r: a per-r sum
#     nesting an inner per-cell m quadrature.
#
# BOTH modes are gated: the estimated arm (the grid r posterior AND the mean
# counts exp(m)) under a two-level exposure offset (o in {0, log 2} within each
# cell, so the anchor's offset, c and log r terms are all exercised), and a
# fixed-r arm (r pinned, mean counts) without one. The engine's c is read from
# the sampler and checked against the formula first. r is read from the
# engine's state after each sweep, so the run is driven one kept sample at a
# time; the gated mean is exp(train - offset) at one row of each cell.
#
# STATED LIMITATION: fork (A) draws only integer-shape (exact) PG variates, so
# this gate exercises NO approximate path; the reference is omega-free.
#
# Tolerances bound sampler MC plus quadrature error; never widen one to pass.
#
# The single-tree fits below run swap at 0.1 rather than the shipped zero:
# swap is the only proposal that rotates a child's rule up the tree, and with
# one tree the change move cannot re-root once the splits beneath the root
# depend on the root's own variable, so nothing else crosses between
# rootings. The kernel's correctness at the shipped mixture is the balance
# gates' job, not this one's.
#
# Usage: Rscript negbin-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

ndpost <- if (quick) 12000L else 30000L
nburn <- 4000L
nSeeds <- 2L
tolMean <- if (quick) 0.12 else 0.07 # mean counts exp(m) per cell
tolGrid <- if (quick) 0.045 else 0.025 # the grid r posterior distribution

# ---- shipped constants (docs/design/negative-binomial.md sections 1, 3) ----

k <- 2
numTrees <- 1L
nodeScale <- 3 # nbinom's log-mean anchor
tau <- nodeScale / (k * sqrt(numTrees)) # leaf-prior sd
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
# the estimated arm's exposure: alternating 1 and 2 within each cell
offsetEst <- rep(c(0, log(2)), length.out = length(y))
offsetFixed <- rep(0, length(y))

shiftOf <- function(offset) log(max(sum(y), 0.5) / sum(exp(offset)))

# ---- exact posterior by structure enumeration + nested quadrature ----

nGrid <- length(grid)

# The exact posterior for one arm: for a cell with counts cnt, offsets off and
# dispersion r, integrate over the cell log mean m ~ N(c, tau^2): g = marginal
# likelihood, mc = its mean-count numerator E[exp(m)].
exactArm <- function(offset) {
  c0 <- shiftOf(offset)
  mGrid <- c0 + seq(-10, 10, by = 0.01)
  wM <- dnorm(mGrid, c0, tau) * (mGrid[2L] - mGrid[1L])
  cellIntegral <- function(rows, r) {
    cnt <- y[rows]
    off <- offset[rows]
    nobs <- length(cnt)
    loglik <- sum(lgamma(cnt + r)) -
      nobs * lgamma(r) -
      sum(lgamma(cnt + 1)) +
      sum(cnt) * mGrid +
      sum(cnt * off) +
      nobs * r * log(r)
    for (i in seq_len(nobs)) {
      loglik <- loglik - (cnt[i] + r) * log(r + exp(mGrid + off[i]))
    }
    wl <- wM * exp(loglik)
    list(g = sum(wl), mc = sum(wl * exp(mGrid)))
  }
  rowsA <- which(cell == 0L)
  rowsB <- which(cell == 1L)
  out <- list(
    rootG = numeric(nGrid),
    rootMC = numeric(nGrid),
    gA = numeric(nGrid),
    gB = numeric(nGrid),
    mcA = numeric(nGrid),
    mcB = numeric(nGrid)
  )
  for (ki in seq_len(nGrid)) {
    r <- grid[ki]
    cr <- cellIntegral(c(rowsA, rowsB), r)
    ca <- cellIntegral(rowsA, r)
    cb <- cellIntegral(rowsB, r)
    out$rootG[ki] <- cr$g
    out$rootMC[ki] <- cr$mc
    out$gA[ki] <- ca$g
    out$gB[ki] <- cb$g
    out$mcA[ki] <- ca$mc
    out$mcB[ki] <- cb$mc
  }
  out
}

# tree prior: a single binary predictor exhausts its one cut, so a split's
# children are forced leaves - root has prior 1 - base, split has prior base
priorRoot <- 1 - base
priorSplit <- base

# ---- estimated arm: sum over the r grid ----

ex <- exactArm(offsetEst)
rootMass <- priorRoot * priorW * ex$rootG # per grid point
splitMass <- priorSplit * priorW * ex$gA * ex$gB
den <- sum(rootMass) + sum(splitMass)
gridPost <- (rootMass + splitMass) / den
# root's mean count is shared by both cells; split's is per-cell (the other
# leaf's marginal integrates out to its g factor)
mcRootNum <- priorRoot * sum(priorW * ex$rootMC)
exactMeanA <- (mcRootNum + priorSplit * sum(priorW * ex$mcA * ex$gB)) / den
exactMeanB <- (mcRootNum + priorSplit * sum(priorW * ex$mcB * ex$gA)) / den

# ---- fixed arm: r pinned at rFixed, no grid sum ----

fx0 <- exactArm(offsetFixed)
kf <- which(grid == rFixed)
denF <- priorRoot * fx0$rootG[kf] + priorSplit * fx0$gA[kf] * fx0$gB[kf]
mcRootNumF <- priorRoot * fx0$rootMC[kf]
exactFixedA <- (mcRootNumF + priorSplit * fx0$mcA[kf] * fx0$gB[kf]) / denF
exactFixedB <- (mcRootNumF + priorSplit * fx0$mcB[kf] * fx0$gA[kf]) / denF

# ---- sampler fits: single tree, per-draw r and mean counts from the state ----

fitSeed <- function(seed, dispersion, offset) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = numTrees,
    updateState = FALSE,
    proposal.probs = c(
      birth_death = 0.5,
      swap = 0.1,
      change = 0.4,
      birth = 0.5
    )
  )
  sampler <- dbarts(
    data.frame(x1 = factor(cell)), # the cell predictor, categorical
    y,
    offset = offset,
    family = nbinom(dispersion = dispersion),
    control = control,
    tree.prior = cgm(power, base),
    leaf.prior = normal(k),
    verbose = FALSE
  )
  shiftGap <- abs(sampler$getLeafPrior()$response.shift - shiftOf(offset))
  if (!(shiftGap <= 1e-12)) {
    cat(sprintf(
      "engine log-mean shift off the formula by %g <- FAIL\n",
      shiftGap
    ))
    quit(status = 1L)
  }

  iA <- which(cell == 0L)[1L]
  iB <- which(cell == 1L)[1L]
  meanA <- 0
  meanB <- 0
  gridCounts <- numeric(nGrid)
  for (s in seq_len(ndpost)) {
    r <- sampler$run(if (s == 1L) nburn else 0L, 1L)
    rDraw <- sampler$getDispersion()
    gridCounts[match(rDraw, grid)] <- gridCounts[match(rDraw, grid)] + 1
    meanA <- meanA + exp(r$train[iA, 1L] - offset[iA])
    meanB <- meanB + exp(r$train[iB, 1L] - offset[iB])
  }
  c(meanA / ndpost, meanB / ndpost, gridCounts / ndpost)
}

runArm <- function(dispersion, offset) {
  rows <- do.call(
    rbind,
    lapply(seq_len(nSeeds), function(sd) fitSeed(sd, dispersion, offset))
  )
  colMeans(rows)
}

# estimated arm
est <- runArm(NULL, offsetEst)
fitMeanA <- est[1L]
fitMeanB <- est[2L]
fitGrid <- est[-(1:2)]

# fixed arm
fx <- runArm(rFixed, offsetFixed)
fitFixedA <- fx[1L]
fitFixedB <- fx[2L]

gapMeanEst <- max(abs(c(fitMeanA - exactMeanA, fitMeanB - exactMeanB)))
gapGrid <- max(abs(fitGrid - gridPost))
gapMeanFixed <- max(abs(c(fitFixedA - exactFixedA, fitFixedB - exactFixedB)))

cat("Negative-binomial exact-posterior gate (single tree, two cells):\n")
cat("--- estimated r (grid posterior), exposure offset ---\n")
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
