# grouped-mixing.R
#
# Permanent mixing harness for the in-engine grouped random-effect sampler
# (rbart_vi with a built-in tau prior -> the GroupedResponse path). It
# institutionalizes the throwaway measurements of docs/plans/tau-slice-review.md
# (whose scripts lived in a since-deleted tmp dir) so any future grouped-mixing
# work has a fixed gate to score against. It reproduces two findings:
#
#   Part A (review 3a) -- the tau IACT/ESS grid over K (few..many groups) x
#     signal strength (weak/strong group effect). tau mixes fine with a strong
#     group signal and POORLY in the weak-signal / small-K corner.
#
#   Part B (review 3d) -- THE ATTRIBUTION. The same weak-signal fits WITH a
#     strong mean forest f vs WITHOUT it (f identically 0, intercept + ranef
#     only), to isolate the forest-ranef confounding contribution. Removing f
#     collapses the tau IACT; the ratio quantifies how much of the weak-signal
#     bottleneck is f-b confounding (the review measured ~25x at K=3, ~3x at
#     K=10). See docs/design/forest-ranef-interweaving.md for the analysis.
#
# METRIC IS LOAD-INDEPENDENT. IACT/ESS depend only on the chain's statistical
# autocorrelation, not on wall-clock, so -- unlike benchmarks/R/bench-sampler.R
# -- this needs NO quiet machine and may run alongside other load.
#
# The single-seed IACT of the weak-signal tau chain is very noisy (the tau
# posterior is heavy-tailed and barely identified at small K -- the review and
# the bcf-ridge landing both flag single-run IACT as unstable), so every cell
# is averaged over several independent-seed replicates; the spread is printed.
#
# Measures the INSTALLED package (tinytest/benchmark convention); rebuild with
# R CMD INSTALL before scoring a change. HEAD's default cauchy prior draws tau
# by the exact Makalic-Schmidt inverse-gamma conjugate step (103a9ef), not the
# slice sampler -- the review established these two mix identically, so the
# attribution is unchanged by that landing.
#
# Self-contained, seed-fixed. Run: Rscript benchmarks/R/grouped-mixing.R
# Tunables below (sizes, seed count) trade runtime for IACT stability.

suppressMessages(library(dbarts))
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("grouped-mixing.R needs the 'coda' package for effectiveSize")
}

## ---- tunables -------------------------------------------------------------
N_OBS <- 1000L # observations (balanced across groups)
N_TREES <- 50L # mean-forest size
N_SAMPLES <- 6000L # kept tau draws per chain (n.thin = 1)
N_BURN <- 2000L # warmup sweeps
SIGMA_EPS <- 1.0 # residual sd (the signal scale reference)
TAU_WEAK <- 0.2 # weak group-effect sd (<< sd(f), << sigma)
TAU_STRONG <- 2.0 # strong group-effect sd
S_WITHIN <- 0.8 # within-group spread of the group-aliased predictor (below)
SEEDS_GRID <- 1:6 # replicate seeds for Part A
SEEDS_ATTR <- 1:8 # replicate seeds for Part B (ratio needs more stability)
BASE_SEED <- 20260720L

## ---- Friedman (1991) mean surface, the canonical BART test function -------
## sd ~ 4.7 with U(0,1) inputs, i.e. a STRONG mean forest relative to a weak
## group signal of sd 0.2 -- the regime where f and b compete hardest.
friedman <- function(X) {
  10 *
    sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] +
    5 * X[, 5]
}

## ---- tau chain diagnostics ------------------------------------------------
## coda::effectiveSize (spectral / AR); IACT = N / ESS; lag-1 autocorrelation.
diag_tau <- function(tau) {
  tau <- as.numeric(tau)
  ess <- as.numeric(coda::effectiveSize(tau))
  list(
    postmean = mean(tau),
    ess = ess,
    iact = length(tau) / ess,
    lag1 = as.numeric(acf(tau, lag.max = 1L, plot = FALSE)$acf[2L])
  )
}

## ---- one grouped fit ------------------------------------------------------
## The confounding regime needs GENUINE identifiability aliasing: group
## membership must correlate with a predictor the forest splits on, so f's
## contribution to a group's mean competes with b_g to explain it. Predictor
## 4 is clustered by group (x4 = plogis(group_center + N(0, S_WITHIN)), group
## centers spread across the latent scale); f's 10*x4 term is then a near-per-
## group constant the forest fights b over. At small K there are few group
## levels so the aliasing swings tau (sd of few b_j) hard; at large K the
## leverage per group falls. (Groups INDEPENDENT of x exercise only the milder
## fit-variance-leakage mode and do NOT reproduce the review's small-K
## dominance -- see docs/design/forest-ranef-interweaving.md.)
##
## S_WITHIN is deliberately MODERATE: tighter clustering makes x4 predict the
## group so sharply that the no-f forest fits b_g directly through x4, dirtying
## the control and collapsing the contrast. withf = FALSE zeroes the mean
## surface (y = b_g + eps): with no y-signal BART shrinks to near-stumps and
## does not split on x4, so the ranef stays cleanly identified -- the review's
## "NO f" arm. Same seed drives both arms of an attribution pair.
fit_tau <- function(K, tau_true, withf, seed) {
  set.seed(BASE_SEED + seed)
  grp <- rep(seq_len(K), length.out = N_OBS)
  groupCenter <- qnorm((seq_len(K) - 0.5) / K)
  x4 <- plogis(groupCenter[grp] + rnorm(N_OBS, 0, S_WITHIN))
  X <- cbind(runif(N_OBS), runif(N_OBS), runif(N_OBS), x4, runif(N_OBS))
  b <- rnorm(K, 0, tau_true)
  fx <- if (withf) friedman(X) else rep(0.0, N_OBS)
  y <- fx + b[grp] + rnorm(N_OBS, 0, SIGMA_EPS)
  fit <- rbart_vi(
    y ~ .,
    data = data.frame(y = y, X),
    group.by = grp,
    prior = cauchy,
    n.trees = N_TREES,
    n.samples = N_SAMPLES,
    n.burn = N_BURN,
    n.thin = 1L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = seed
  )
  diag_tau(fit$tau)
}

## mean over replicate seeds of a per-seed statistic
cell <- function(K, tau_true, withf, seeds) {
  d <- lapply(seeds, function(s) fit_tau(K, tau_true, withf, s))
  iacts <- vapply(d, `[[`, numeric(1), "iact")
  ess <- vapply(d, `[[`, numeric(1), "ess")
  list(
    postmean = mean(vapply(d, `[[`, numeric(1), "postmean")),
    ess = mean(ess),
    iact = mean(iacts),
    iact_lo = min(iacts),
    iact_hi = max(iacts),
    lag1 = mean(vapply(d, `[[`, numeric(1), "lag1"))
  )
}

t_start <- proc.time()[3]
cat(sprintf(
  "grouped-mixing: n=%d, n.trees=%d, n.samples=%d, n.burn=%d, sigma=%.1f\n",
  N_OBS,
  N_TREES,
  N_SAMPLES,
  N_BURN,
  SIGMA_EPS
))
cat(sprintf(
  "prior=cauchy (in-engine exact-IG), 1 chain; IACT averaged over seeds\n\n"
))

## ---- Part A: the K x signal grid (review 3a) ------------------------------
cat("== Part A: tau IACT/ESS grid (K x signal), all WITH forest f ==\n")
cat(sprintf(
  "%4s %8s %10s %8s %8s %6s %6s\n",
  "K",
  "signal",
  "postmean",
  "ESS",
  "IACT",
  "lag1",
  "seeds"
))
gridK <- c(3L, 10L, 40L)
for (sig in c("weak", "strong")) {
  tau_true <- if (sig == "weak") TAU_WEAK else TAU_STRONG
  for (K in gridK) {
    r <- cell(K, tau_true, TRUE, SEEDS_GRID)
    cat(sprintf(
      "%4d %8s %10.3f %8.0f %8.1f %6.3f %6d\n",
      K,
      sprintf("%s(%.1f)", sig, tau_true),
      r$postmean,
      r$ess,
      r$iact,
      r$lag1,
      length(SEEDS_GRID)
    ))
  }
}

## ---- Part B: the attribution (review 3d) ----------------------------------
cat("\n== Part B: forest-ranef confounding attribution (weak signal) ==\n")
cat("   with f = strong Friedman mean; no f = intercept-only (f == 0)\n")
cat(sprintf(
  "%4s %14s %14s %8s   %s\n",
  "K",
  "IACT_with_f",
  "IACT_no_f",
  "ratio",
  "IACT range (with f | no f)"
))
for (K in gridK) {
  wf <- cell(K, TAU_WEAK, TRUE, SEEDS_ATTR)
  nf <- cell(K, TAU_WEAK, FALSE, SEEDS_ATTR)
  cat(sprintf(
    "%4d %14.1f %14.1f %8.1f   [%.0f-%.0f | %.0f-%.0f]\n",
    K,
    wf$iact,
    nf$iact,
    wf$iact / nf$iact,
    wf$iact_lo,
    wf$iact_hi,
    nf$iact_lo,
    nf$iact_hi
  ))
}

cat(sprintf("\ntotal runtime: %.0f s\n", proc.time()[3] - t_start))
