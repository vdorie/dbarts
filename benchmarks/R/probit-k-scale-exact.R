#!/usr/bin/env Rscript

# Exact-posterior gate for the probit rescaling step (dbartsControl's
# probitRescaleForest; docs/design/probit-k-scale-move.md). The step multiplies
# the probit latents and the leaves by one factor and divides k by it, the
# factor drawn from its exact conditional, so it must leave the posterior of
# (k, leaves) exactly where it was while moving along the direction k mixes
# slowest in.
#
# THE CONFIGURATION. One tree whose structure is frozen (every structural
# proposal probability zero), drawn once from the prior at a fixed seed to hold
# three or four leaves; probit; k ~ chi(1.5, 2); n = 150. Given the tree the
# leaves are independent given k, so the posterior of k is exact by one
# one-dimensional integral a leaf,
#
#     p(k | y) propto p(k) prod_l int prod_{i in l} Phi(s_i (o_i + mu))
#                                        N(mu; 0, (c / k)^2) dmu,
#
# s_i = 2 y_i - 1 over the active rows, c the leaf scale (3 for one tree), and
# E[mu_l | y] by the same integrals. A leaf holding one class has no finite
# likelihood mode and is integrated by adaptive quadrature on z = k mu / c; a
# leaf holding both is summed on a fine grid in mu. The law of log k is then
# integrated by the trapezoid rule: counting whole grid cells biased every
# arm's quantile probabilities by 0.001 to 0.003, which put correct samplers
# at z -3 to -5.
#
# Four arms, each exercising one term of the step's conditional:
#   pure    every leaf separated: the small-k mass the step exists to reach;
#   mixed   no leaf separated;
#   offset  mixed, o = 0.4 on every row: the offset's cross term;
#   mask    mixed, a fifth of the rows inactive: the active-row count and the
#           latents left unscaled.
#
# STATISTICS. P(k < q_j) at the exact deciles q_j, j = 0.1, 0.25, 0.5, 0.75,
# 0.9, and E[mu_l] on each unseparated leaf, each a z against a 50-batch
# batch-means error, failing above 4.5. And one mixing bound: the pure arm's
# largest batch-means error of P(k < q_j). With the step that error is 0.0039
# at 4e6 sweeps and without it 0.0194, so a step silently switched off fails
# here rather than passing more loosely. Full mode bounds it at 0.01; quick
# mode at QUICK_MIXING_BOUND below, set between the two sides' quick values
# (0.0064 with the step, 0.0190 without, at the landing).
#
# The `never` argument runs with the step off, which must fail: without the
# step the pure arm is the sampler the step was built to replace.
#
# Usage: Rscript probit-k-scale-exact.R [quick] [never] [pure] [mixed]
#   [offset] [mask]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
rescale <- !("never" %in% args)
armNames <- intersect(c("pure", "mixed", "offset", "mask"), args)
if (length(armNames) == 0L) {
  armNames <- c("pure", "mixed", "offset", "mask")
}

n <- 150L
thin <- 10L
numSweeps <- if (quick) 1e6 else 4e6
numBatches <- 50L
zBound <- 4.5
QUICK_MIXING_BOUND <- 0.012
mixingBound <- if (quick) QUICK_MIXING_BOUND else 0.01
deciles <- c(0.1, 0.25, 0.5, 0.75, 0.9)
leafScale <- 3
kPrior <- c(df = 1.5, scale = 2)

set.seed(401L)
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, paste0("x", 1:3)))
drawnK <- dbartsPriors$normal(dbartsPriors$chi(
  kPrior[["df"]],
  kPrior[["scale"]]
))
frozen <- c(birth_death = 0, swap = 0, change = 0, perturb = 0, rule_gibbs = 0)
control <- dbartsControl(
  n.trees = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = as.integer(numSweeps / thin),
  n.thin = thin,
  updateState = FALSE,
  verbose = FALSE,
  seed = 20261008L,
  probitRescaleForest = rescale,
  proposal.probs = c(frozen, birth = 0.5)
)
sampler <- dbarts(
  x,
  as.double(rep_len(0:1, n)),
  leaf.prior = drawnK,
  control = control,
  family = "probit"
)
repeat {
  sampler$sampleTreesFromPrior()
  sampler$setLeafPrior(dbartsPriors$normal(k = 1))
  sampler$sampleLeafParametersFromPrior()
  fits <- as.numeric(sampler$predict(x))
  numLeaves <- length(unique(fits))
  if (numLeaves >= 3L && numLeaves <= 4L) break
}
leaf <- match(fits, unique(fits))

# ---- the exact law ----

logKGrid <- seq(log(1e-3), log(50), length.out = 4001L)
muStep <- 0.004
muGrid <- seq(-80, 80, by = muStep)

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# log p(k) + log k, the chi(df, scale) density of k on the log-k scale
logKPrior <- function(k) {
  u <- k / kPrior[["scale"]]
  (kPrior[["df"]] - 1) * log(u) - 0.5 * u^2 + log(k)
}

exactLaw <- function(y, o, active) {
  perLeaf <- lapply(seq_len(numLeaves), function(l) {
    rows <- active > 0 & leaf == l
    s <- 2 * y[rows] - 1
    offsets <- o[rows]
    if (length(unique(s)) == 1L && all(offsets == 0)) {
      # separated: the likelihood tends to one, so integrate on z = k mu / c
      count <- length(s)
      logIntegral <- vapply(
        exp(logKGrid),
        function(k) {
          log(
            integrate(
              function(z) {
                exp(count * pnorm(leafScale * z / k, log.p = TRUE)) * dnorm(z)
              },
              -Inf,
              Inf,
              rel.tol = 1e-12,
              subdivisions = 2000L
            )$value
          )
        },
        0
      )
      return(list(
        logIntegral = logIntegral,
        mean = rep(NA_real_, length(logKGrid))
      ))
    }
    logLik <- numeric(length(muGrid))
    for (key in unique(paste(s, offsets))) {
      group <- paste(s, offsets) == key
      j <- which(group)[1L]
      logLik <- logLik +
        sum(group) * pnorm(s[j] * (offsets[j] + muGrid), log.p = TRUE)
    }
    byK <- vapply(
      exp(logKGrid),
      function(k) {
        logWeight <- logLik + dnorm(muGrid, 0, leafScale / k, log = TRUE)
        weight <- exp(logWeight - max(logWeight))
        c(
          logSumExp(logWeight) + log(muStep),
          sum(weight * muGrid) / sum(weight)
        )
      },
      c(0, 0)
    )
    list(logIntegral = byK[1L, ], mean = byK[2L, ])
  })
  logPosterior <- logKPrior(exp(logKGrid)) +
    Reduce(`+`, lapply(perLeaf, `[[`, "logIntegral"))
  density <- exp(logPosterior - max(logPosterior))
  cumulative <- c(0, cumsum((density[-1L] + density[-length(density)]) / 2))
  cumulative <- cumulative / cumulative[length(cumulative)]
  weight <- density / sum(density)
  # the tails underflow to flat runs of 0 and 1, which the inverse skips
  rising <- !duplicated(cumulative)
  list(
    quantiles = exp(
      approx(cumulative[rising], logKGrid[rising], xout = deciles)$y
    ),
    leafMean = vapply(perLeaf, function(p) sum(weight * p$mean), 0)
  )
}

batchMeansError <- function(v) {
  batches <- colMeans(matrix(
    v[seq_len(numBatches * (length(v) %/% numBatches))],
    ncol = numBatches
  ))
  sd(batches) / sqrt(numBatches)
}

# ---- the arms ----

anyFailure <- FALSE
for (arm in armNames) {
  set.seed(match(arm, c("pure", "mixed", "offset", "mask")) + 500L)
  o <- if (arm == "offset") rep(0.4, n) else rep(0, n)
  active <- if (arm == "mask") as.double(seq_len(n) %% 5L != 0L) else rep(1, n)
  probability <- if (arm == "pure") {
    c(1, 0, 1, 0)[leaf]
  } else {
    pnorm(o + c(0.8, -0.5, 0.3, -1.2)[leaf])
  }
  y <- as.double(rbinom(n, 1L, probability))
  exact <- exactLaw(y, o, active)

  sampler$setOffset(if (arm == "offset") o else NULL)
  sampler$setLeafPrior(dbartsPriors$normal(k = 2))
  sampler$sampleLeafParametersFromPrior()
  sampler$setLeafPrior(drawnK)
  sampler$setResponse(y)
  sampler$setActiveRows(if (arm == "mask") active else NULL)
  started <- proc.time()[[3L]]
  run <- sampler$run(2000L, as.integer(numSweeps / thin))
  elapsed <- proc.time()[[3L]] - started

  k <- as.numeric(run$k)
  leafDraws <- vapply(
    seq_len(numLeaves),
    function(l) colMeans(run$train[leaf == l, , drop = FALSE]),
    numeric(ncol(run$train))
  ) -
    o[1L]
  indicator <- lapply(exact$quantiles, function(q) as.numeric(k < q))
  errorP <- vapply(indicator, batchMeansError, 0)
  zP <- (vapply(indicator, mean, 0) - deciles) / errorP
  unseparated <- which(!is.na(exact$leafMean))
  zMu <- vapply(
    unseparated,
    function(l) {
      (mean(leafDraws[, l]) - exact$leafMean[l]) /
        batchMeansError(leafDraws[, l])
    },
    0
  )
  failed <- anyNA(c(zP, zMu)) || any(abs(c(zP, zMu)) > zBound)
  mixingFailed <- arm == "pure" && !(max(errorP) < mixingBound)
  cat(sprintf(
    "%-6s %d leaves, %.0e sweeps in %.0f s | z(P(k < q)) %s | z(E mu) %s | max se(P) %.4f%s%s\n",
    arm,
    numLeaves,
    numSweeps,
    elapsed,
    paste(sprintf("%+.1f", zP), collapse = " "),
    if (length(zMu) > 0L) paste(sprintf("%+.1f", zMu), collapse = " ") else "-",
    max(errorP),
    if (failed) " <- FAIL" else "",
    if (mixingFailed) {
      sprintf(" <- FAIL mixing (bound %.4f)", mixingBound)
    } else {
      ""
    }
  ))
  anyFailure <- anyFailure || failed || mixingFailed
}

if (anyFailure) {
  cat(
    "\nFAIL: the probit k chain deviates from the exact posterior or mixes too slowly\n"
  )
  quit(status = 1L)
}
cat("\nOK: the probit k chain matches the exact posterior\n")
