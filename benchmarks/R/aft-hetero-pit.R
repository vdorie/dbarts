#!/usr/bin/env Rscript

# Validation of the PER-OBSERVATION scale the heteroscedastic AFT sampler
# draws its censored latents at (docs/design/aft-variance-forest.md).
#
# The gate has two legs. RECOVERY: under a two-level true s(x) and heavy
# right-censoring the posterior surface must separate the two levels and the
# mean surface must stay unbiased. PIT: on the same fit, after burn-in, each
# censored row's redraw is a lower-truncated normal at (mu_i, s_i, b_i), so its
# probability integral transform
#
#   v = (Phi((z - mu)/s) - Phi((b - mu)/s)) / (1 - Phi((b - mu)/s))
#
# is U(0, 1) whenever s is the scale the draw actually used. Conditional on the
# chain state every v is exactly uniform and the v's are conditionally
# independent across rows, and because that conditional law does not depend on
# the state the pooled collection is unconditionally i.i.d. uniform - so a
# per-cell ks.test against punif is an exact test, not an asymptotic one.
# Pooling by x-cell is what makes it sensitive to the SCALE rather than to the
# mean: a redraw run at a scale constant across rows piles the low-s cell's v
# at 0 and the high-s cell's at 1, which the per-cell KS rejects by many orders
# of magnitude. It is the only gate sensitive to a wrong per-observation scale
# installed by the chain, and it needs no oracle.
#
# The pairing is load-bearing and is itself part of what is tested: sweep t's
# redraw runs against the surface the variance forest left at the END of sweep
# t - 1, and against the mean fits recorded at sweep t (the variance sweep does
# not move them). A one-sweep misalignment is a non-uniform v.
#
# Usage: Rscript aft-hetero-pit.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

numBurnIn <- if (quick) 500L else 1500L
numKept <- if (quick) 60L else 200L
thin <- 5L
# The null is exact, so any nominal level is valid; this one is set far below
# 0.05 so the gate does not false-alarm on a host whose draws differ - it can
# afford to, the poison arm rejecting by hundreds of orders of magnitude rather
# than marginally.
alpha <- 1e-5

# ---- a two-level variance surface under heavy right-censoring ----

set.seed(2027L)
n <- 600L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
cell <- ifelse(x[, 1L] < 0.5, 1L, 2L)
sTrue <- c(0.3, 1.2)[cell]
fTrue <- 1.5 * x[, 2L]
logTtrue <- fTrue + sTrue * rnorm(n)
logC <- fTrue + 0.05 + 0.8 * rnorm(n)
status <- as.numeric(logTtrue <= logC)
obsLogT <- ifelse(status == 1, logTtrue, logC)
censored <- status == 0
cat(sprintf(
  "censoring rate: %.2f (cell 1 %.2f, cell 2 %.2f)\n",
  mean(censored),
  mean(censored[cell == 1L]),
  mean(censored[cell == 2L])
))

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  seed = 91L
)
sampler <- dbarts(
  x,
  obsLogT,
  control = control,
  variance = varianceForest(n.trees = 20L)
)
ctrl <- sampler$control
attr(ctrl, "bartcore.survival") <- status
sampler$control <- ctrl
bc <- dbarts:::bartcoreSampler(sampler, family = "aft")

invisible(bartcoreRun(bc, numBurnIn, 0L))

# ---- drive the chain one thinned block at a time ----

sSum <- numeric(n)
fitSum <- numeric(n)
pit <- vector("list", numKept)
atBound <- 0L
for (t in seq_len(numKept)) {
  r <- bartcoreRun(bc, 0L, thin)
  z <- bartcoreGetLatents(bc)
  mu <- r$train[, thin]
  # the surface the LAST sweep's redraw ran against: the variance forest sweeps
  # after the latent refresh, so sweep `thin`'s draw saw sweep `thin - 1`'s
  s <- sqrt(r$variance[, thin - 1L])
  sSum <- sSum + rowMeans(sqrt(r$variance))
  fitSum <- fitSum + rowMeans(r$train)

  zc <- (z[censored] - mu[censored]) / s[censored]
  bc.std <- (obsLogT[censored] - mu[censored]) / s[censored]
  # the upper-tail ratio, in logs: stable where the truncation is deep, which
  # the naive difference of two CDFs is not
  v <- 1 -
    exp(
      pnorm(zc, lower.tail = FALSE, log.p = TRUE) -
        pnorm(bc.std, lower.tail = FALSE, log.p = TRUE)
    )
  atBound <- atBound + sum(zc <= bc.std)
  pit[[t]] <- pmin(pmax(v, 0), 1)
}

sHat <- sSum / numKept
fitHat <- fitSum / numKept

# ---- recovery leg ----

sLow <- mean(sHat[cell == 1L])
sHigh <- mean(sHat[cell == 2L])
bias <- mean(fitHat - fTrue)
fitCor <- cor(fitHat, fTrue)
cat(sprintf(
  "surface: cell 1 %.3f (truth 0.30), cell 2 %.3f (truth 1.20), ratio %.2f\n",
  sLow,
  sHigh,
  sHigh / sLow
))
cat(sprintf("mean surface: bias %.3f, correlation %.3f\n", bias, fitCor))

recoveryFailed <- sHigh < 2 * sLow ||
  sLow < 0.15 ||
  sLow > 0.6 ||
  sHigh < 0.7 ||
  sHigh > 1.8 ||
  abs(bias) > 0.15 ||
  fitCor < 0.7

# ---- PIT leg, per cell ----

cellOfCensored <- cell[censored]
pitAll <- unlist(pit, use.names = FALSE)
cellAll <- rep(cellOfCensored, times = numKept)

pValue <- numeric(2L)
for (cc in 1:2) {
  u <- pitAll[cellAll == cc]
  pValue[cc] <- suppressWarnings(ks.test(u, "punif")$p.value)
  cat(sprintf(
    "cell %d: %d draws, mean PIT %.4f, KS p %.3g%s\n",
    cc,
    length(u),
    mean(u),
    pValue[cc],
    if (pValue[cc] < alpha) " <- FAIL" else ""
  ))
}
if (atBound > 0L) {
  cat(sprintf("latents at the truncation bound: %d\n", atBound))
}

if (recoveryFailed) {
  cat("recovery leg FAILED\n")
}
if (recoveryFailed || any(pValue < alpha)) {
  quit(status = 1L)
}
cat("\nOK: censored latents are redrawn at each row's own s(x)\n")
