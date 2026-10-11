#!/usr/bin/env Rscript

# Exact-posterior gate for the moves that rescale a binary forest against its
# drawn k: the probit rescaling step (docs/design/probit-k-scale-move.md),
# which multiplies the probit latents and the leaves by one factor and divides
# k by it, the factor drawn from its exact conditional, and the collapsed
# rescaling move, which draws the factor with the latents integrated out and
# then draws the latents again. Each must leave the posterior of (k, leaves)
# exactly where it was while moving along the direction k mixes slowest in.
#
# THE CONFIGURATION, per link. One tree whose structure is frozen (every
# structural proposal probability zero), drawn once from the prior at a fixed
# seed to hold three or four leaves; probit or logistic; k ~ chi(1.5, 2);
# n = 150. Given the tree the leaves are independent given k, so the posterior
# of k is exact by one one-dimensional integral a leaf,
#
#     p(k | y) propto p(k) prod_l int prod_{i in l} F(s_i (o_i + mu))^w_i
#                                        N(mu; 0, (c / k)^2) dmu,
#
# F the link's distribution function (Phi or the logistic), s_i = 2 y_i - 1
# over the active rows, w_i the row's count (1 but in the counts arm), c the
# leaf scale, read from the sampler (3 for one probit tree, pi sqrt(3) for one
# logistic tree), and
# E[mu_l | y] by the same integrals. A leaf holding one class has no finite
# likelihood mode and is integrated by adaptive quadrature on z = k mu / c; a
# leaf holding both is summed on a fine grid in mu. The posterior of log k is then
# integrated by the trapezoid rule: counting whole grid cells biased every
# arm's quantile probabilities by 0.001 to 0.003, which put correct samplers
# at z -3 to -5.
#
# The arms, each exercising one term of a move's conditional:
#   pure    every leaf separated: the small-k mass the moves exist to reach;
#   mixed   no leaf separated;
#   offset  mixed, o = 0.4 on every row: the offset's term;
#   mask    mixed, a fifth of the rows inactive: the rows in the likelihood
#           and the latents left alone;
#   counts  logistic only: mixed, counts 1 to 3 as weights.
#
# STATISTICS. P(k < q_j) at the exact deciles q_j, j = 0.1, 0.25, 0.5, 0.75,
# 0.9, and E[mu_l] on each unseparated leaf, each a z against a 50-batch
# batch-means error, failing above 4.5. And one mixing bound, under probit:
# the pure arm's
# largest batch-means error of P(k < q_j). With the step that error is 0.0039
# at 4e6 sweeps and without it 0.0194, so a step that silently declines fails
# here rather than passing more loosely. Full mode bounds it at 0.01; quick
# mode at QUICK_MIXING_BOUND below, set between the two sides' quick values
# at the landing: with the step 0.0064 on arm64 and 0.0099 on x86-64, whose
# chain takes another path; without it 0.0190 (measured while the step had a
# switch; it has none now). The logistic pure arm reports the same error and
# is held to no bound.
#
# The gate cannot see k left unscaled after the latents and leaves are scaled
# correctly: that mutant passes quick mode at a worst |z| of 4.3. The engine's
# mapping test (tests/cpp, runForestRescaleTests) is what catches it.
#
# Usage: Rscript probit-k-scale-exact.R [quick] [probit] [logistic]
#          [pure] [mixed] [offset] [mask] [counts]
# With no link named both run; with no arm named, every arm of each link.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
allArms <- c("pure", "mixed", "offset", "mask", "counts")
armNames <- intersect(allArms, args)
if (length(armNames) == 0L) {
  armNames <- allArms
}
links <- intersect(c("probit", "logistic"), args)
if (length(links) == 0L) {
  links <- c("probit", "logistic")
}

n <- 150L
thin <- 10L
numSweeps <- if (quick) 1e6 else 4e6
numBatches <- 50L
zBound <- 4.5
QUICK_MIXING_BOUND <- 0.014
mixingBound <- if (quick) QUICK_MIXING_BOUND else 0.01
deciles <- c(0.1, 0.25, 0.5, 0.75, 0.9)
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
  proposal.probs = c(frozen, birth = 0.5)
)

# ---- the exact posterior ----

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

# tree is the frozen tree as list(leaf, numLeaves, leafScale); logCdf the
# link's log distribution function; w the counts
exactPosterior <- function(tree, logCdf, y, o, active, w) {
  leafScale <- tree$leafScale
  perLeaf <- lapply(seq_len(tree$numLeaves), function(l) {
    rows <- active > 0 & tree$leaf == l
    s <- 2 * y[rows] - 1
    offsets <- o[rows]
    counts <- w[rows]
    if (length(unique(s)) == 1L && all(offsets == 0)) {
      # separated: the likelihood tends to one, so integrate on z = k mu / c
      count <- sum(counts)
      logIntegral <- vapply(
        exp(logKGrid),
        function(k) {
          log(
            integrate(
              function(z) {
                exp(count * logCdf(leafScale * z / k)) * dnorm(z)
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
        sum(counts[group]) * logCdf(s[j] * (offsets[j] + muGrid))
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

# ---- the arms, per link ----

anyFailure <- FALSE
for (link in links) {
  logCdf <- if (link == "probit") {
    function(q) pnorm(q, log.p = TRUE)
  } else {
    function(q) plogis(q, log.p = TRUE)
  }
  cdf <- if (link == "probit") pnorm else plogis
  sampler <- dbarts(
    x,
    as.double(rep_len(0:1, n)),
    leaf.prior = drawnK,
    control = control,
    family = link
  )
  repeat {
    sampler$sampleTreesFromPrior()
    sampler$setLeafPrior(dbartsPriors$normal(k = 1))
    sampler$sampleLeafParametersFromPrior()
    fits <- as.numeric(sampler$predict(x))
    numLeaves <- length(unique(fits))
    if (numLeaves >= 3L && numLeaves <= 4L) break
  }
  tree <- list(
    leaf = match(fits, unique(fits)),
    numLeaves = numLeaves,
    # the family's own scale for the one tree, not a constant written here
    leafScale = sampler$getLeafPrior()$k.scale
  )
  leaf <- tree$leaf

  for (arm in setdiff(armNames, if (link == "probit") "counts")) {
    set.seed(match(arm, allArms) + 500L)
    o <- if (arm == "offset") rep(0.4, n) else rep(0, n)
    active <- if (arm == "mask") {
      as.double(seq_len(n) %% 5L != 0L)
    } else {
      rep(1, n)
    }
    w <- if (arm == "counts") as.double(rep_len(1:3, n)) else rep(1, n)
    probability <- if (arm == "pure") {
      c(1, 0, 1, 0)[leaf]
    } else {
      cdf(o + c(0.8, -0.5, 0.3, -1.2)[leaf])
    }
    y <- as.double(rbinom(n, 1L, probability))
    exact <- exactPosterior(tree, logCdf, y, o, active, w)

    sampler$setOffset(if (arm == "offset") o else NULL)
    sampler$setLeafPrior(dbartsPriors$normal(k = 2))
    sampler$sampleLeafParametersFromPrior()
    sampler$setLeafPrior(drawnK)
    sampler$setResponse(y)
    if (link == "logistic") {
      sampler$setWeights(w)
    }
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
    mixingFailed <- link == "probit" &&
      arm == "pure" &&
      !(max(errorP) < mixingBound)
    cat(sprintf(
      "%-8s %-6s %d leaves, %.0e sweeps in %.0f s | z(P(k < q)) %s | z(E mu) %s | max se(P) %.4f%s%s\n",
      link,
      arm,
      numLeaves,
      numSweeps,
      elapsed,
      paste(sprintf("%+.1f", zP), collapse = " "),
      if (length(zMu) > 0L) {
        paste(sprintf("%+.1f", zMu), collapse = " ")
      } else {
        "-"
      },
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
}

if (anyFailure) {
  cat(
    "\nFAIL: a k chain deviates from the exact posterior or mixes too slowly\n"
  )
  quit(status = 1L)
}
cat("\nOK: the k chains match the exact posterior\n")
