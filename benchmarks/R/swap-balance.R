#!/usr/bin/env Rscript

# Permanent exact-posterior gate for the SWAP move's detailed balance. The
# swap gives a node a child's split rule and the swapped child(ren) the node's;
# its acceptance is a pure Metropolis pi-ratio (symmetric proposal, no Hastings
# term) - moves.hpp swapMove scores exp(yLogPi + yLogL - xLogPi - xLogL). A
# defective swap that dropped the tree-prior ratio and accepted on the
# likelihood alone would target the wrong stationary distribution; this gate
# isolates that with a case where the likelihood ratio is EXACTLY one, so the
# whole signal lives in the prior term the defect would omit.
#
# Design (pure-interaction XOR on a 2-D ordinal grid):
#   x1 in 1:6 (K1-1 = 5 cuts), x2 in 1:4 (K2-1 = 3 cuts) - UNEQUAL cut counts.
#   y = M*sign + noise, sign = +1 iff (x1 > 3) == (x2 > 2). Every marginal is
#   flat; only a depth-2 crossed tree separates the +/-M cells. Two depth-2
#   arrangements give the IDENTICAL 4-cell partition:
#     A: root splits x1 @3.5, both children split x2 @2.5;
#     B: root splits x2 @2.5, both children split x1 @3.5.
#   A single root swap converts A<->B (the "children share a rule" branch:
#   both children carry the same rule, so both are swapped). Because A and B
#   induce the same partition their integrated likelihood is EXACTLY equal, so
#     pi(A)/pi(B) = prior(A)/prior(B),
#   a PRIOR-ONLY, calibration-independent ratio. Under the CGM rule prior each
#   split contributes -log(numAvailVars) - log(availCuts); every node here has
#   numAvailVars = 2, the growth/shape factors are identical between A and B,
#   and the only surviving difference is the chosen variable's cut count, so
#     prior(A)/prior(B) = (K1-1)/(K2-1) = 5/3,  P(A | A or B) = 5/8 = 0.625.
#   A prior-dropped / likelihood-only swap accepts A<->B with ratio 1 both ways
#   (equal likelihood) -> a symmetric two-state chain with fixed point 0.5. With
#   the swap dominating the A<->B flow (birth_death only 0.1, change 0) that 0.5
#   is the defect's realized target; the 0.625-vs-0.5 gap is the gate's power.
#
# The blocker this design had to crack: a single unregularized tree overfits the
# pure XOR blocks (redundant same-mean sub-splits -> ~23-node trees), so clean
# 7-node depth-2 trees never occur. The fix is a LARGER fixed residual variance:
# fixed(6) here makes the internal sigma large enough that within-block noise is
# not worth fitting (no redundant splits) yet small enough that the +/-M jump
# still forces the full depth-2 crossed structure. At this setting ~85% of draws
# are clean 7-node trees and every one lands on the informative cuts (x1 @3.5,
# x2 @2.5), so the A/B classes are exactly {A,B} and the 5/8 target is exact.
#
# Two arms:
#   analytic : P(A | A or B) = 5/8, derived from the CGM prior (asserted below
#              against an explicit two-tree prior computation, with lik(A) ==
#              lik(B) asserted as a guard).
#   engine   : a long single-tree, swap-dominant dbarts run; the frequency of A
#              among clean crossed 7-node draws at the informative cuts.
# The engine must match 5/8 (small |z|) and sit many sigma from the 0.5
# wrong-target (demonstrated power).
#
# Usage: Rscript swap-balance.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

# ---- fixed design ----

K1 <- 6L # x1 in 1:6 -> 5 cuts
K2 <- 4L # x2 in 1:4 -> 3 cuts
xd1 <- 3.5 # x1 XOR boundary cut (between values 3 and 4)
xd2 <- 2.5 # x2 XOR boundary cut (between values 2 and 3)
M <- 2.0
noiseSD <- 0.2
nPer <- 4L # replicates per grid cell
fixedVar <- 6.0 # fixed residual variance - the concentration lever (see header)
base <- 0.95
power <- 2.0
kLeaf <- 2.0
nodeScale <- 0.5

dataSeed <- 20260723L
engineSeed <- 20260723L
nKept <- if (quick) 60000L else 300000L
batchSize <- if (quick) 20000L else 50000L
nThin <- 10L
nBurn <- 4000L
proposalProbs <- c(birth_death = 0.1, swap = 0.9, change = 0, birth = 0.5)

zBound <- 4 # |z| vs exact must be within this
powerBound <- 20 # z vs wrong-target must exceed this (gate must have power)
concBound <- 0.99 # clean A/B must concentrate on the informative cuts

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# ---- analytic target from the CGM tree prior ----
#
# Log-prior of a specific depth-2 tree under cgm(power, base): each internal
# node contributes log(psplit_depth) + rule, each leaf log(1 - psplit_depth),
# with psplit_depth = base / (1 + depth)^power and the rule term
# -log(numAvailVars) - log(availCuts_of_chosen_var). A and B share the shape
# (root at depth 0, two internals at depth 1, four leaves at depth 2) so every
# growth/leaf factor cancels; only the chosen variables' cut counts differ.

psplit <- function(depth) base / (1 + depth)^power
ruleLP <- function(numAvail, availCuts) -log(numAvail) - log(availCuts)

# node availability for A and B (verified numAvail == 2 everywhere below):
#   A root: split x1 over full grid   -> numAvail 2, availCuts K1-1
#   A childL (x1<=3): split x2 full   -> numAvail 2, availCuts K2-1
#   A childR (x1>3):  split x2 full   -> numAvail 2, availCuts K2-1
#   B root: split x2 over full grid   -> numAvail 2, availCuts K2-1
#   B childL (x2<=2): split x1 full   -> numAvail 2, availCuts K1-1
#   B childR (x2>2):  split x1 full   -> numAvail 2, availCuts K1-1
availA <- list(
  c(numAvail = 2L, cuts = K1 - 1L), # root
  c(numAvail = 2L, cuts = K2 - 1L), # childL
  c(numAvail = 2L, cuts = K2 - 1L) # childR
)
availB <- list(
  c(numAvail = 2L, cuts = K2 - 1L),
  c(numAvail = 2L, cuts = K1 - 1L),
  c(numAvail = 2L, cuts = K1 - 1L)
)
stopifnot(all(vapply(c(availA, availB), function(a) a["numAvail"] == 2L, NA)))

logPriorTree <- function(avail) {
  # depths: root 0, two internals 1, four leaves 2
  internalDepths <- c(0L, 1L, 1L)
  lp <- 0
  for (i in seq_along(avail)) {
    lp <- lp +
      log(psplit(internalDepths[i])) +
      ruleLP(avail[[i]]["numAvail"], avail[[i]]["cuts"])
  }
  lp + 4 * log(1 - psplit(2L)) # four depth-2 leaves stop
}

logPriorA <- unname(logPriorTree(availA))
logPriorB <- unname(logPriorTree(availB))
priorRatioAB <- exp(logPriorA - logPriorB)
targetExact <- 1 / (1 + exp(logPriorB - logPriorA)) # P(A | A or B) = piA/(piA+piB)
targetWrong <- 0.5 # prior-dropped / likelihood-only swap fixed point

# cross-checks: the ratio must reduce to (K1-1)/(K2-1) and the target to 5/8
stopifnot(abs(priorRatioAB - (K1 - 1) / (K2 - 1)) < 1e-12)
stopifnot(abs(targetExact - (K1 - 1) / ((K1 - 1) + (K2 - 1))) < 1e-12)
stopifnot(abs(targetExact - 0.625) < 1e-12)

# ---- lik(A) == lik(B) guard ----
#
# A and B induce the identical 4-cell partition, so their integrated likelihood
# is exactly equal. Build the data, compute the constant-leaf conjugate marginal
# (verbatim from ConstantGaussianLeaf, matched calibration) for each tree's four
# leaves, and assert equality. The target above is calibration-independent, so
# this only guards the "equal partition" premise.

set.seed(dataSeed)
grid <- expand.grid(x1 = seq_len(K1), x2 = seq_len(K2))
grid <- grid[rep(seq_len(nrow(grid)), each = nPer), ]
x1 <- as.double(grid$x1)
x2 <- as.double(grid$x2)
sgn <- ifelse((x1 > 3) == (x2 > 2), 1, -1)
y <- M * sgn + rnorm(length(x1), sd = noiseSD)
n <- length(y)

yRange <- max(y) - min(y)
yScaled <- (y - min(y)) / yRange - 0.5
residVar <- fixedVar / yRange^2 # fixed(v) -> internal residual variance v/range^2
priorPrecision <- (kLeaf / nodeScale)^2

logIL <- function(idx) {
  z <- yScaled[idx]
  nLeaf <- length(z)
  if (nLeaf == 0L) {
    return(0)
  }
  posteriorPrecision <- nLeaf / residVar
  mean <- sum(z) / nLeaf
  centeredSumOfSquares <- sum(z * z) - sum(z) * mean
  0.5 *
    log(priorPrecision / (priorPrecision + posteriorPrecision)) -
    0.5 * centeredSumOfSquares / residVar -
    0.5 *
      ((priorPrecision * mean) * (posteriorPrecision * mean)) /
      (priorPrecision + posteriorPrecision)
}

# A leaves: {x1<=3,x2<=2},{x1<=3,x2>2},{x1>3,x2<=2},{x1>3,x2>2}
# B leaves: {x2<=2,x1<=3},{x2<=2,x1>3},{x2>2,x1<=3},{x2>2,x1>3}  (same 4 cells)
likA <-
  logIL(x1 <= 3 & x2 <= 2) +
  logIL(x1 <= 3 & x2 > 2) +
  logIL(x1 > 3 & x2 <= 2) +
  logIL(x1 > 3 & x2 > 2)
likB <-
  logIL(x2 <= 2 & x1 <= 3) +
  logIL(x2 <= 2 & x1 > 3) +
  logIL(x2 > 2 & x1 <= 3) +
  logIL(x2 > 2 & x1 > 3)
stopifnot(abs(likA - likB) < 1e-9)

# ---- tree classifiers ----
#
# getTrees rows are preorder; a balanced 7-node depth-2 tree is
# 1=root,2=Lchild,3=leaf,4=leaf,5=Rchild,6=leaf,7=leaf. classifyShape is the
# variable-only lump (root=x1/children=x2 -> A, mirror -> B). classifyExact
# additionally pins the informative cuts, so exact A/B are precisely {A,B} and
# the 5/8 target carries no approximation. Concentration of shape onto exact is
# reported as a health check.

classifyShape <- function(sub) {
  if (nrow(sub) != 7L) {
    return("other")
  }
  isLeaf <- sub$var == -1L
  if (!all(!isLeaf[c(1L, 2L, 5L)]) || !all(isLeaf[c(3L, 4L, 6L, 7L)])) {
    return("other")
  }
  rv <- sub$var[1L]
  lv <- sub$var[2L]
  rrv <- sub$var[5L]
  if (rv == 1L && lv == 2L && rrv == 2L) {
    return("A")
  }
  if (rv == 2L && lv == 1L && rrv == 1L) {
    return("B")
  }
  "other"
}

classifyExact <- function(sub, shape) {
  if (shape == "A") {
    if (sub$value[1L] == xd1 && sub$value[2L] == xd2 && sub$value[5L] == xd2) {
      return("A")
    }
  } else if (shape == "B") {
    if (sub$value[1L] == xd2 && sub$value[2L] == xd1 && sub$value[5L] == xd1) {
      return("B")
    }
  }
  "other"
}

# ---- engine arm: swap-dominant single-tree kernel ----

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  useQuantiles = TRUE,
  keepTrees = TRUE,
  n.samples = batchSize,
  n.burn = nBurn,
  n.thin = nThin,
  updateState = TRUE,
  rngSeed = engineSeed,
  n.cuts = 100L
)
sampler <- dbarts(
  cbind(x1 = x1, x2 = x2),
  y,
  control = ctl,
  tree.prior = cgm(power, base),
  node.prior = normal(kLeaf),
  resid.prior = fixed(fixedVar),
  proposal.probs = proposalProbs
)
stopifnot(is.null(sampler$data@offset))

nBatch <- as.integer(ceiling(nKept / batchSize))
shapeCls <- character(0L) # per-sample shape class (loose)
exactCls <- character(0L) # per-sample exact class (informative cuts)
first <- TRUE
for (b in seq_len(nBatch)) {
  if (first) {
    invisible(sampler$run(nBurn, batchSize))
    first <- FALSE
  } else {
    invisible(sampler$run(0L, batchSize))
  }
  tr <- sampler$getTrees()
  bySample <- split(seq_len(nrow(tr)), tr$sample)
  sh <- character(length(bySample))
  ex <- character(length(bySample))
  for (i in seq_along(bySample)) {
    sub <- tr[bySample[[i]], ]
    sh[i] <- classifyShape(sub)
    ex[i] <- if (sh[i] %in% c("A", "B")) classifyExact(sub, sh[i]) else "other"
  }
  shapeCls <- c(shapeCls, sh)
  exactCls <- c(exactCls, ex)
}
nDraws <- length(exactCls)

# ---- batch-means MC error over the A/B subsequence ----

batchMeanSE <- function(v, nBatches = 400L) {
  len <- (length(v) %/% nBatches) * nBatches
  bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
  sd(bm) / sqrt(nBatches)
}

isAB <- exactCls %in% c("A", "B")
indA <- as.integer(exactCls[isAB] == "A") # ordered A/B subsequence
nAB <- length(indA)
pEngine <- mean(indA)
se <- batchMeanSE(indA)
zExact <- (pEngine - targetExact) / se
zWrong <- (pEngine - targetWrong) / se

# concentration of loose-clean shape onto the informative cuts
shapeAB <- shapeCls %in% c("A", "B")
concFrac <- sum(isAB) / sum(shapeAB)
cleanFrac <- mean(shapeAB)

# mixing health: A<->B transitions along the A/B subsequence
seqAB <- exactCls[isAB]
transitions <- sum(seqAB[-1L] != seqAB[-length(seqAB)])

# ---- report ----

cat(sprintf(
  "swap-balance gate  (n=%d obs, x1 %d cuts, x2 %d cuts, draws=%d, thin=%d)\n",
  n,
  K1 - 1L,
  K2 - 1L,
  nDraws,
  nThin
))
cat(sprintf(
  "  prior(A)/prior(B) = %.6f (= (K1-1)/(K2-1) = %d/%d);  lik(A)-lik(B) = %.2e\n",
  priorRatioAB,
  K1 - 1L,
  K2 - 1L,
  likA - likB
))
cat(sprintf(
  "  clean 7-node crossed fraction = %.3f;  on informative cuts = %.4f (of clean)\n",
  cleanFrac,
  concFrac
))
cat(sprintf(
  "  exact A/B draws = %d;  A<->B transitions = %d (%.1f%% of A/B steps)\n",
  nAB,
  transitions,
  100 * transitions / max(1L, nAB - 1L)
))
cat("\n")
cat(sprintf("%-16s %10s %10s %10s\n", "quantity", "value", "MCse", "z"))
cat(sprintf(
  "%-16s %10.4f %10.5f %10s\n",
  "P(A|A+B) engine",
  pEngine,
  se,
  ""
))
cat(sprintf(
  "%-16s %10.4f %10s %+10.1f\n",
  "  vs exact 5/8",
  targetExact,
  "",
  zExact
))
cat(sprintf(
  "%-16s %10.4f %10s %+10.1f\n",
  "  vs wrong 1/2",
  targetWrong,
  "",
  zWrong
))

# ---- verdict ----

passTarget <- is.finite(zExact) && abs(zExact) < zBound
passPower <- is.finite(zWrong) && zWrong > powerBound
passConc <- is.finite(concFrac) && concFrac > concBound
passClean <- cleanFrac > 0.3
pass <- passTarget && passPower && passConc && passClean

cat("\n================ VERDICT ================\n")
cat(sprintf(
  "  matches exact 5/8   : %s (|z| = %.1f < %g)\n",
  if (passTarget) "PASS" else "FAIL",
  abs(zExact),
  zBound
))
cat(sprintf(
  "  power vs wrong 1/2  : %s (z = %.1f > %g)\n",
  if (passPower) "PASS" else "FAIL",
  zWrong,
  powerBound
))
cat(sprintf(
  "  cut concentration   : %s (%.4f > %g)\n",
  if (passConc) "PASS" else "FAIL",
  concFrac,
  concBound
))
cat(sprintf(
  "  clean fraction      : %s (%.3f > 0.30)\n",
  if (passClean) "PASS" else "FAIL",
  cleanFrac
))
cat(sprintf(
  "\nSWAP BALANCE GATE: %s\n",
  if (pass) "PASS" else "FAIL"
))

# Exit nonzero on failure so a runner/CI can catch it (matches bd-balance.R).
if (!pass) {
  quit(status = 1L, save = "no")
}
