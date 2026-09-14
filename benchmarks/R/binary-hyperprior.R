#!/usr/bin/env Rscript

# Re-evaluates the binary (probit) end-node hyperprior default, chi(1.5, 2).
#
# The July 2026 study (docs/plans/archive/chi-default-research.md) held the
# degrees of freedom at 1.5, varied only the scale, and ran four simulated
# data-generating processes. This harness varies BOTH parameters and widens
# the case set: six simulated processes crossed with sample size, predictor
# count and base rate, plus twenty-two real datasets scored by repeated
# splits: six from R and its recommended packages, and sixteen from the UCI
# Machine Learning Repository loaded by benchmarks/R/uci-binary.R.
#
# REAL DATA. The real-data column is the one that could move the package
# default, and the six R datasets are too few and too alike to settle it, so
# the sixteen UCI datasets restore the breadth the original k-sensitivity work
# had: 208 to 48,842 rows, 3 to 60 predictors, positive rates from 0.085 to
# 0.65, and seven datasets with factor predictors. A split is 80/20 of the
# whole dataset, as before, except that a dataset of more than 5,000 rows
# draws 4,000 training rows and 1,000 held-out rows per split, so that the
# cost of a fit stays bounded and the large datasets differ from the small
# ones in their predictors rather than in what a forest can afford to see.
# The recorded nTrain column says which rule a row got. The UCI files are
# downloaded on first use into $DBARTS_BENCH_DATA, or into
# tools::R_user_dir("dbarts", "cache") when that is unset, and verified
# against a recorded sha256; nothing is fetched until a UCI block runs.
#
# ARMS. Twenty-four hyperprior arms, chi(df, scale) for df in
# {1, 1.25, 1.5, 2, 3} crossed with scale in {1, 2, 5, Inf} plus df in
# {1.5, 3} crossed with scale in {0.5, 0.25}, and four fixed-k arms, k in
# {1, 1.5, 2, 3}. k = 2 is the fixed value BayesTree used and the value the
# chi(1.5, 2) prior is centered near; the scales below 1 are there because
# coverage rises monotonically in interval width across the rest of the grid,
# so an interior optimum, if there is one, is below scale 1. Arms are paired:
# within a case and a repetition every arm sees the same data and starts from
# the same MCMC seed, so a difference between arms is a difference in the
# prior, not in the draw. BINARY_HYPERPRIOR_ARMS narrows the set to a
# semicolon-separated list of arm names or to a named subset ("mixing").
#
# SCORES, all on held-out rows. Log score and Brier score against the held-out
# outcome; for the simulated cases, where the true probability is known, the
# coverage and mean width of the 90 percent posterior interval for that
# probability, and the root mean squared error of the posterior mean
# probability. Every arm also reports the posterior median, 90th percentile
# and maximum of the sampled k (a fixed arm reports its own k) and the fit's
# elapsed seconds.
#
# MIXING. Every fit reports split-Rhat and effective sample size for two
# scalar series: the sampled k, which is what the prior is about and which a
# fixed arm leaves missing, and the held-out mean probability, which the
# forest moves whether or not k is sampled. They say whether a score is a
# property of the posterior or of the chain length, which is the whole point
# of running the same cells at two lengths and calling `compare`.
#
# BLOCKS. One block is one invocation and writes one rds, so a full run splits
# across a session. A simulated block is a data-generating process at one
# sample size ("sim:friedman:500"), which is nine cells - three predictor
# counts by three base rates - at the repetition count below; appending a
# predictor count ("sim:friedman:2000:50") narrows it to that count's three
# cells and its own rds, which is how the large-n blocks are kept short. A
# real block is one dataset ("real:biopsy", "real:spambase"). `blocks` lists
# them, `all` runs every one in sequence, and `summarize` reads a directory of
# rds files and prints the tables the plan doc reports.
#
# Usage:
#   Rscript benchmarks/R/binary-hyperprior.R blocks
#   Rscript benchmarks/R/binary-hyperprior.R sim:friedman:500 [outdir] [quick]
#   Rscript benchmarks/R/binary-hyperprior.R real:biopsy [outdir]
#   Rscript benchmarks/R/binary-hyperprior.R all [outdir]
#   Rscript benchmarks/R/binary-hyperprior.R summarize [outdir]
#   Rscript benchmarks/R/binary-hyperprior.R compare [shortdir] [longdir]
#
# A simulated block also narrows to one base rate ("sim:separable:100:50:0.2"),
# which is how a single cell is refit at a length the grid cannot afford.
# `compare` pairs two directories on (cell, repetition, arm) and reports the
# second minus the first, so with a short directory first and a long one
# second it reports what the chain length does.
#
# outdir defaults to benchmarks/results/binary-hyperprior. `quick` cuts the
# repetitions and the draws for a smoke run; a quick rds is marked as such and
# summarize refuses to mix the two. BINARY_HYPERPRIOR_CORES sets the worker
# count (default 4), BINARY_HYPERPRIOR_REPS the simulated repetitions,
# BINARY_HYPERPRIOR_SPLITS the real-data splits, and BINARY_HYPERPRIOR_BURN,
# _DRAWS and _CHAINS the MCMC length a fit gets - the last three are how the
# convergence check reported in the plan doc was run, and they ride the saved
# settings so summarize says which length produced a directory. No baseline
# and no pass/fail exit status: this is a measurement whose verdict a person
# writes, in docs/plans/binary-hyperprior.md.
#
# Findings live in docs/plans/binary-hyperprior.md. Nothing here changes a
# package default.

suppressPackageStartupMessages(library(dbarts))
suppressPackageStartupMessages(library(parallel))

## ---------------------------------------------------------------- settings

nTrees <- 75L
nBurn <- 500L
nSamples <- 500L
nChains <- 1L
nTestSim <- 1000L
simReps <- 8L
realSplits <- 60L
probClamp <- 1e-6

quickSettings <- function() {
  nBurn <<- 100L
  nSamples <<- 100L
  nTestSim <<- 400L
  simReps <<- 2L
  realSplits <<- 3L
}

# Repetition counts are overridable so an expensive stratum can be run at
# fewer repetitions than the rest; summarize reports the count per sample
# size, so a run that mixes them says so.
envCount <- function(name, fallback) {
  value <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(value) || value < 1L) fallback else value
}

applyEnvironmentOverrides <- function() {
  simReps <<- envCount("BINARY_HYPERPRIOR_REPS", simReps)
  realSplits <<- envCount("BINARY_HYPERPRIOR_SPLITS", realSplits)
  nBurn <<- envCount("BINARY_HYPERPRIOR_BURN", nBurn)
  nSamples <<- envCount("BINARY_HYPERPRIOR_DRAWS", nSamples)
  nChains <<- envCount("BINARY_HYPERPRIOR_CHAINS", nChains)
}

defaultCores <- function() {
  envCount("BINARY_HYPERPRIOR_CORES", 4L)
}

## -------------------------------------------------------------------- arms

chiArm <- function(df, scale) {
  list(
    name = sprintf("chi(%s, %s)", format(df), format(scale)),
    kind = "chi",
    df = df,
    scale = scale,
    fixed.k = NA_real_
  )
}

fixedArm <- function(k) {
  list(
    name = sprintf("k = %s", format(k)),
    kind = "fixed",
    df = NA_real_,
    scale = NA_real_,
    fixed.k = k
  )
}

makeArms <- function() {
  arms <- list()
  for (scale in c(1, 2, 5, Inf)) {
    for (df in c(1, 1.25, 1.5, 2, 3)) {
      arms[[length(arms) + 1L]] <- chiArm(df, scale)
    }
  }
  # Scales below 1 at two degrees of freedom. chi(3, 0.5) and chi(3, 0.25)
  # sit near the prior medians of chi(1.5, 1) and chi(1.5, 0.5), so the pairs
  # test whether df and scale stay redundant through the prior median once
  # the median drops below one.
  for (scale in c(0.5, 0.25)) {
    for (df in c(1.5, 3)) {
      arms[[length(arms) + 1L]] <- chiArm(df, scale)
    }
  }
  for (k in c(1, 1.5, 2, 3)) {
    arms[[length(arms) + 1L]] <- fixedArm(k)
  }
  arms
}

incumbent <- "chi(1.5, 2)"

# A named arm subset, for legs too expensive to run at the full width. The
# incumbent is always kept, since summarize pairs every arm against it.
armSubsets <- list(
  mixing = c(
    "chi(1.5, 2)",
    "chi(1.5, 0.5)",
    "chi(1.5, 0.25)",
    "chi(1.25, 1)",
    "chi(1.5, Inf)",
    "k = 2"
  )
)

selectArms <- function(all) {
  request <- Sys.getenv("BINARY_HYPERPRIOR_ARMS")
  if (!nzchar(request) || request == "all") {
    return(all)
  }
  wanted <- if (!is.null(armSubsets[[request]])) {
    armSubsets[[request]]
  } else {
    trimws(strsplit(request, ";", fixed = TRUE)[[1L]])
  }
  wanted <- union(wanted, incumbent)
  names <- vapply(all, function(a) a$name, character(1L))
  unknown <- setdiff(wanted, names)
  if (length(unknown) > 0L) {
    stop("unknown arm(s): ", paste(unknown, collapse = ", "))
  }
  all[names %in% wanted]
}

arms <- selectArms(makeArms())

nodePriorFor <- function(arm) {
  if (arm$kind == "chi") {
    dbartsPriors$normal(dbartsPriors$chi(arm$df, arm$scale))
  } else {
    dbartsPriors$normal(arm$fixed.k)
  }
}

## ------------------------------------------------------- simulated designs

# Latent surfaces. Each takes a matrix of uniform predictors and returns an
# unstandardized latent value; only the first few columns are active, so the
# rest are noise at p = 20 and p = 50. The four July processes are friedman,
# linear, strong and weak; separable and interaction are new. strong and weak
# share one smooth additive surface and differ only in the signal standard
# deviation below, which is the axis their names name.
latentSurfaces <- list(
  friedman = function(x) {
    10 *
      sin(pi * x[, 1L] * x[, 2L]) +
      20 * (x[, 3L] - 0.5)^2 +
      10 * x[, 4L] +
      5 * x[, 5L]
  },
  linear = function(x) 2 * x[, 1L] + x[, 2L] - x[, 3L],
  strong = function(x) {
    sin(2 * pi * x[, 1L]) + 2 * (x[, 2L] - 0.5)^2 + 0.5 * x[, 3L]
  },
  weak = function(x) {
    sin(2 * pi * x[, 1L]) + 2 * (x[, 2L] - 0.5)^2 + 0.5 * x[, 3L]
  },
  separable = function(x) x[, 1L] + x[, 2L],
  interaction = function(x) {
    3 *
      (x[, 1L] > 0.5) *
      (x[, 2L] > 0.5) -
      3 * (x[, 3L] > 0.5) * (x[, 4L] > 0.5) +
      2 * (2 * x[, 1L] - 1) * (2 * x[, 5L] - 1)
  }
)

# Latent standard deviation each surface is scaled to. The signal-to-noise
# ratio is what separates strong, weak and separable from each other and from
# the two shape processes.
signalSd <- c(
  friedman = 1.5,
  linear = 1.5,
  strong = 3,
  weak = 0.35,
  separable = 6,
  interaction = 1.5
)

dgpNames <- names(latentSurfaces)
sampleSizes <- c(100L, 500L, 2000L)
predictorCounts <- c(5L, 20L, 50L)
baseRates <- c(0.05, 0.2, 0.5)

# Reference moments and the intercept that puts a process at a base rate, both
# computed once per (dgp, p, rate) from a large fixed pool so the process is a
# fixed truth rather than something that moves with the training draw.
calibrationCache <- new.env(parent = emptyenv())

calibrationFor <- function(dgp, p, rate) {
  key <- sprintf("%s|%d|%s", dgp, p, format(rate))
  if (!is.null(calibrationCache[[key]])) {
    return(calibrationCache[[key]])
  }
  set.seed(90000L + p)
  pool <- matrix(runif(20000L * p), 20000L, p)
  raw <- latentSurfaces[[dgp]](pool)
  center <- mean(raw)
  spread <- sd(raw)
  mu <- signalSd[[dgp]] * (raw - center) / spread
  intercept <- uniroot(
    function(a) mean(pnorm(mu + a)) - rate,
    interval = c(-30, 30),
    tol = 1e-8
  )$root
  value <- list(center = center, spread = spread, intercept = intercept)
  calibrationCache[[key]] <- value
  value
}

# A probit fit needs both outcome classes present, and at n = 100 with a base
# rate of 0.05 an all-zero training draw is not rare. Such a draw is discarded
# and redrawn from the next seed in a stated sequence; the attempt count rides
# the result so the affected cells are visible rather than silently dropped.
simulatedCase <- function(dgp, n, p, rate, seed) {
  cal <- calibrationFor(dgp, p, rate)
  probabilityOf <- function(x) {
    raw <- latentSurfaces[[dgp]](x)
    pnorm(signalSd[[dgp]] * (raw - cal$center) / cal$spread + cal$intercept)
  }
  for (attempt in seq_len(100L)) {
    set.seed(seed + 1000000L * (attempt - 1L))
    xTrain <- matrix(runif(n * p), n, p)
    xTest <- matrix(runif(nTestSim * p), nTestSim, p)
    pTrain <- probabilityOf(xTrain)
    pTest <- probabilityOf(xTest)
    yTrain <- rbinom(n, 1L, pTrain)
    if (length(unique(yTrain)) == 2L) break
  }
  if (length(unique(yTrain)) < 2L) {
    stop("no two-class training draw after 100 attempts")
  }
  colnames(xTrain) <- colnames(xTest) <- paste0("x", seq_len(p))
  list(
    train = data.frame(xTrain, y = yTrain),
    test = data.frame(xTest),
    yTest = rbinom(nTestSim, 1L, pTest),
    pTest = pTest,
    attempts = attempt
  )
}

## ------------------------------------------------------------- real data

# Six binary-outcome datasets from R and its recommended packages: small to
# moderate n, base rates from 0.21 to 0.40, numeric and factor predictors,
# one dichotomized survival outcome. The sixteen UCI datasets are appended
# below, after these, so the seed of an existing dataset does not move.
realDatasets <- list(
  pima = function() {
    data("Pima.tr", package = "MASS", envir = environment())
    data("Pima.te", package = "MASS", envir = environment())
    d <- rbind(Pima.tr, Pima.te)
    list(x = d[, setdiff(names(d), "type")], y = as.integer(d$type == "Yes"))
  },
  biopsy = function() {
    data("biopsy", package = "MASS", envir = environment())
    d <- biopsy[complete.cases(biopsy), ]
    list(
      x = d[, paste0("V", 1:9)],
      y = as.integer(d$class == "malignant")
    )
  },
  infert = function() {
    d <- datasets::infert
    list(
      x = d[, c(
        "education",
        "age",
        "parity",
        "induced",
        "spontaneous"
      )],
      y = as.integer(d$case)
    )
  },
  kyphosis = function() {
    d <- rpart::kyphosis
    list(
      x = d[, c("Age", "Number", "Start")],
      y = as.integer(d$Kyphosis == "present")
    )
  },
  pbc = function() {
    d <- survival::pbc
    keep <- c(
      "age",
      "sex",
      "ascites",
      "hepato",
      "spiders",
      "edema",
      "bili",
      "chol",
      "albumin",
      "copper",
      "alk.phos",
      "ast",
      "trig",
      "platelet",
      "protime",
      "stage"
    )
    d <- d[!is.na(d$trt), c("status", keep)]
    d <- d[complete.cases(d), ]
    list(x = d[, keep], y = as.integer(d$status == 2L))
  },
  birthwt = function() {
    data("birthwt", package = "MASS", envir = environment())
    d <- birthwt
    d$race <- factor(d$race)
    d$ftv <- factor(pmin(d$ftv, 2L))
    list(
      x = d[, c("age", "lwt", "race", "smoke", "ptl", "ht", "ui", "ftv")],
      y = as.integer(d$low)
    )
  }
)

# The sixteen UCI datasets, in descending size. benchmarks/R/uci-binary.R
# holds the urls, the sha256 of every file and the per-dataset cleaning; it is
# sourced the first time a UCI block asks for data, so a simulated block or a
# block on one of the six above needs no cache and no network. The names are
# listed here rather than read from that file because `blocks` and the seed of
# a real block are both functions of this order, and neither should depend on
# whether anything has been downloaded yet.
uciDatasetNames <- c(
  "adult",
  "bank",
  "magic",
  "mushroom",
  "spambase",
  "banknote",
  "german",
  "tictactoe",
  "transfusion",
  "creditapproval",
  "wdbc",
  "climate",
  "ionosphere",
  "haberman",
  "cleveland",
  "sonar"
)

uciLoaded <- FALSE

harnessDir <- function() {
  flag <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(flag) == 0L) {
    "benchmarks/R"
  } else {
    dirname(sub("^--file=", "", flag[1L]))
  }
}

ensureUci <- function() {
  if (!uciLoaded) {
    source(file.path(harnessDir(), "uci-binary.R"))
    uciLoaded <<- TRUE
  }
}

uciLoader <- function(name) {
  force(name)
  function() {
    ensureUci()
    d <- get("uciBinary")[[name]]()
    list(x = d[, setdiff(names(d), "y"), drop = FALSE], y = d$y)
  }
}

for (uciName in uciDatasetNames) {
  realDatasets[[uciName]] <- uciLoader(uciName)
}
rm(uciName)

# Above this many rows a split takes a fixed 4,000 training rows and 1,000
# held-out rows instead of 80/20 of everything, which is what keeps adult and
# bank inside the same per-fit budget as the rest.
subsampleAbove <- 5000L
subsampleTrain <- 4000L
subsampleTest <- 1000L

## ----------------------------------------------------- mixing diagnostics

# Split-Rhat and effective sample size, on the usual definitions (Vehtari et
# al. 2021), written out here so the harness needs no package beyond dbarts:
# split every chain in half, compare the between-half and within-half
# variances, and sum the averaged autocorrelations with Geyer's initial
# monotone rule. `draws` is iterations by chains.
splitSequences <- function(draws) {
  draws <- as.matrix(draws)
  n <- nrow(draws)
  half <- n %/% 2L
  if (half < 8L) {
    return(NULL)
  }
  cbind(
    draws[seq_len(half), , drop = FALSE],
    draws[seq.int(n - half + 1L, n), , drop = FALSE]
  )
}

# Autocovariance at every lag by FFT, divisor n, mean removed.
autocovariance <- function(x) {
  n <- length(x)
  padded <- 2^ceiling(log2(2 * n))
  spectrum <- fft(c(x - mean(x), rep(0, padded - n)))
  Re(fft(spectrum * Conj(spectrum), inverse = TRUE))[seq_len(n)] / (padded * n)
}

mixingDiagnostics <- function(draws) {
  none <- c(rhat = NA_real_, ess = NA_real_)
  seqs <- splitSequences(draws)
  if (is.null(seqs) || any(!is.finite(seqs))) {
    return(none)
  }
  n <- nrow(seqs)
  m <- ncol(seqs)
  within <- apply(seqs, 2L, var)
  if (any(!is.finite(within)) || min(within) <= 0) {
    return(none)
  }
  W <- mean(within)
  B <- if (m > 1L) n * var(colMeans(seqs)) else 0
  varPlus <- ((n - 1) * W + B) / n
  rhat <- sqrt(varPlus / W)

  acov <- vapply(seq_len(m), function(j) autocovariance(seqs[, j]), numeric(n))
  rho <- 1 - (W - rowMeans(acov)) / varPlus
  rho[1L] <- 1
  # Successive pairs, cut at the first non-positive one and then forced
  # decreasing; this is what keeps the autocorrelation sum from running into
  # the noise at long lags.
  nPairs <- (n - 1L) %/% 2L
  pairs <- rho[2L * seq_len(nPairs) - 1L] + rho[2L * seq_len(nPairs)]
  nonPositive <- which(pairs <= 0)
  if (length(nonPositive) > 0L && nonPositive[1L] > 1L) {
    pairs <- pairs[seq_len(nonPositive[1L] - 1L)]
  } else if (length(nonPositive) > 0L) {
    pairs <- pairs[1L]
  }
  pairs <- cummin(pairs)
  tau <- -1 + 2 * sum(pairs)
  total <- n * m
  ess <- if (tau > 0) min(total / tau, total * log10(total)) else NA_real_
  c(rhat = rhat, ess = ess)
}

## ------------------------------------------------------------ fit and score

# Pools the chain margin of a draws array into a draws-by-column matrix, with
# the draw index running within a chain, so a reshape by column recovers the
# iterations-by-chains layout the diagnostics want.
poolChains <- function(x) {
  if (length(dim(x)) < 3L) {
    return(x)
  }
  matrix(aperm(x, c(2L, 1L, 3L)), dim(x)[2L] * dim(x)[1L], dim(x)[3L])
}

fitAndScore <- function(arm, case, mcmcSeed) {
  elapsed <- system.time(
    fit <- bart(
      y ~ .,
      data = case$train,
      test = case$test,
      family = "probit",
      n.trees = nTrees,
      n.samples = nSamples,
      n.burn = nBurn,
      n.chains = nChains,
      n.threads = 1L,
      keepTrees = FALSE,
      verbose = FALSE,
      seed = mcmcSeed,
      combineChains = FALSE,
      node.prior = nodePriorFor(arm)
    )
  )[["elapsed"]]

  probabilities <- pnorm(poolChains(fit$yhat.test))
  pHat <- colMeans(probabilities)
  clamped <- pmin(pmax(pHat, probClamp), 1 - probClamp)
  y <- case$yTest
  logScore <- -mean(y * log(clamped) + (1 - y) * log(1 - clamped))
  brier <- mean((pHat - y)^2)

  if (is.null(case$pTest)) {
    coverage <- NA_real_
    width <- NA_real_
    probRmse <- NA_real_
  } else {
    bounds <- apply(probabilities, 2L, quantile, probs = c(0.05, 0.95))
    coverage <- mean(case$pTest >= bounds[1L, ] & case$pTest <= bounds[2L, ])
    width <- mean(bounds[2L, ] - bounds[1L, ])
    probRmse <- sqrt(mean((pHat - case$pTest)^2))
  }

  # Two scalar series per fit carry the mixing report: the sampled k, which is
  # the quantity the prior is about, and the held-out mean probability, which
  # is the forest's own summary and moves even when k is fixed.
  kMatrix <- if (is.null(fit$k)) {
    NULL
  } else {
    matrix(as.vector(t(as.matrix(fit$k))), nSamples, nChains)
  }
  kDraws <- if (is.null(kMatrix)) {
    rep(arm$fixed.k, nSamples * nChains)
  } else {
    as.vector(kMatrix)
  }
  kMixing <- if (is.null(kMatrix)) {
    c(rhat = NA_real_, ess = NA_real_)
  } else {
    mixingDiagnostics(kMatrix)
  }
  pMixing <- mixingDiagnostics(matrix(
    rowMeans(probabilities),
    nSamples,
    nChains
  ))

  data.frame(
    arm = arm$name,
    df = arm$df,
    scale = arm$scale,
    fixed.k = arm$fixed.k,
    logScore = logScore,
    brier = brier,
    coverage = coverage,
    width = width,
    probRmse = probRmse,
    kMedian = median(kDraws),
    kQ90 = unname(quantile(kDraws, 0.9)),
    kMax = max(kDraws),
    kRhat = unname(kMixing[["rhat"]]),
    kEss = unname(kMixing[["ess"]]),
    pRhat = unname(pMixing[["rhat"]]),
    pEss = unname(pMixing[["ess"]]),
    seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

runCell <- function(cellLabel, caseFor, reps, seedBase, cores, extra) {
  rows <- list()
  for (rep in seq_len(reps)) {
    case <- caseFor(seedBase + rep)
    mcmcSeed <- seedBase + 500000L + rep
    scored <- mclapply(
      arms,
      function(arm) fitAndScore(arm, case, mcmcSeed),
      mc.cores = cores,
      mc.preschedule = TRUE
    )
    failed <- vapply(scored, function(s) !is.data.frame(s), logical(1L))
    if (any(failed)) {
      stop(sprintf(
        "%s rep %d: %d arm(s) failed; first: %s",
        cellLabel,
        rep,
        sum(failed),
        as.character(scored[[which(failed)[1L]]])
      ))
    }
    scored <- do.call(rbind, scored)
    rows[[length(rows) + 1L]] <- cbind(
      cell = cellLabel,
      rep = rep,
      nTrain = nrow(case$train),
      nPositive = sum(case$train$y),
      attempts = if (is.null(case$attempts)) 1L else case$attempts,
      extra,
      scored,
      row.names = NULL
    )
  }
  do.call(rbind, rows)
}

## ------------------------------------------------------------------ blocks

simBlockNames <- function() {
  as.vector(t(outer(
    dgpNames,
    sampleSizes,
    function(d, n) sprintf("sim:%s:%d", d, n)
  )))
}

realBlockNames <- function() paste0("real:", names(realDatasets))

allBlockNames <- function() c(simBlockNames(), realBlockNames())

# A stable integer per case, so every seed in the run is a stated function of
# the case and the repetition and nothing depends on evaluation order.
simCaseIndex <- function(dgp, n, p, rate) {
  1000L *
    match(dgp, dgpNames) +
    100L * match(n, sampleSizes) +
    10L * match(p, predictorCounts) +
    match(rate, baseRates)
}

runSimBlock <- function(
  dgp,
  n,
  cores,
  counts = predictorCounts,
  rates = baseRates
) {
  rows <- list()
  for (p in counts) {
    for (rate in rates) {
      label <- sprintf("%s|%d|%d|%s", dgp, n, p, format(rate))
      seedBase <- 10L * simCaseIndex(dgp, n, p, rate)
      cat(sprintf("  %-28s ", label))
      elapsed <- system.time(
        block <- runCell(
          label,
          function(seed) simulatedCase(dgp, n, p, rate, seed),
          simReps,
          seedBase,
          cores,
          data.frame(
            kind = "sim",
            dgp = dgp,
            n = n,
            p = p,
            rate = rate,
            stringsAsFactors = FALSE
          )
        )
      )[["elapsed"]]
      cat(sprintf("%6.1f s\n", elapsed))
      rows[[length(rows) + 1L]] <- block
    }
  }
  do.call(rbind, rows)
}

runRealBlock <- function(name, cores) {
  raw <- realDatasets[[name]]()
  n <- nrow(raw$x)
  subsampled <- n > subsampleAbove
  nTrain <- if (subsampled) subsampleTrain else floor(0.8 * n)
  label <- sprintf(
    "%s|%d|%d|%s",
    name,
    n,
    ncol(raw$x),
    format(round(mean(raw$y), 3))
  )
  cat(sprintf("  %-28s ", label))
  # A large dataset draws train and test together out of one permutation, so
  # the held-out rows are a fixed 1,000 rather than the other 80 percent; a
  # small one keeps the 80/20 rule the six R datasets were run under, down to
  # the order of the random draws, so their earlier output still pairs.
  caseFor <- function(seed) {
    for (attempt in seq_len(100L)) {
      set.seed(seed + 1000000L * (attempt - 1L))
      if (subsampled) {
        drawn <- sample.int(n, nTrain + subsampleTest)
        trainRows <- drawn[seq_len(nTrain)]
        testRows <- drawn[nTrain + seq_len(subsampleTest)]
      } else {
        trainRows <- sample.int(n, nTrain)
        testRows <- setdiff(seq_len(n), trainRows)
      }
      if (length(unique(raw$y[trainRows])) == 2L) break
    }
    list(
      train = data.frame(
        raw$x[trainRows, , drop = FALSE],
        y = raw$y[trainRows]
      ),
      test = raw$x[testRows, , drop = FALSE],
      yTest = raw$y[testRows],
      pTest = NULL,
      attempts = attempt
    )
  }
  elapsed <- system.time(
    block <- runCell(
      label,
      caseFor,
      realSplits,
      10L * (8000L + match(name, names(realDatasets))),
      cores,
      data.frame(
        kind = "real",
        dgp = name,
        n = n,
        p = ncol(raw$x),
        rate = mean(raw$y),
        stringsAsFactors = FALSE
      )
    )
  )[["elapsed"]]
  cat(sprintf("%6.1f s\n", elapsed))
  block
}

runBlock <- function(block, outDir, cores, quick) {
  parts <- strsplit(block, ":", fixed = TRUE)[[1L]]
  cat(sprintf("%s (%d arms, %d cores)\n", block, length(arms), cores))
  result <- if (parts[1L] == "sim") {
    counts <- if (length(parts) >= 4L) {
      as.integer(parts[4L])
    } else {
      predictorCounts
    }
    rates <- if (length(parts) >= 5L) {
      matched <- baseRates[
        which.min(abs(baseRates - as.numeric(parts[5L])))
      ]
      if (!isTRUE(all.equal(matched, as.numeric(parts[5L])))) {
        stop("no such base rate: ", parts[5L])
      }
      matched
    } else {
      baseRates
    }
    runSimBlock(parts[2L], as.integer(parts[3L]), cores, counts, rates)
  } else if (parts[1L] == "real") {
    runRealBlock(parts[2L], cores)
  } else {
    stop("unknown block: ", block)
  }
  attr(result, "quick") <- quick
  attr(result, "settings") <- list(
    nArms = length(arms),
    nTrees = nTrees,
    nBurn = nBurn,
    nSamples = nSamples,
    nChains = nChains,
    nTestSim = nTestSim,
    simReps = simReps,
    realSplits = realSplits
  )
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  file <- file.path(outDir, paste0(gsub(":", "-", block, fixed = TRUE), ".rds"))
  saveRDS(result, file)
  cat(sprintf("  wrote %s (%d rows)\n", file, nrow(result)))
  invisible(result)
}

## --------------------------------------------------------------- summarize

standardError <- function(x) sd(x) / sqrt(length(x))

# Mean and standard error of one score over a set of rows, plus the paired
# difference from the incumbent arm: within a (cell, rep) both arms saw the
# same data and the same seed, so the paired spread is the one that answers
# "is this arm better", and it is much tighter than the unpaired one.
armSummary <- function(rows, score) {
  rows <- rows[!is.na(rows[[score]]), ]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  base <- rows[rows$arm == incumbent, c("cell", "rep", score)]
  names(base)[3L] <- "incumbent"
  merged <- merge(rows, base, by = c("cell", "rep"))
  parts <- split(merged, merged$arm)
  out <- do.call(
    rbind,
    lapply(names(parts), function(a) {
      d <- parts[[a]]
      data.frame(
        arm = a,
        mean = mean(d[[score]]),
        se = standardError(d[[score]]),
        diff = mean(d[[score]] - d$incumbent),
        diffSe = standardError(d[[score]] - d$incumbent),
        stringsAsFactors = FALSE
      )
    })
  )
  out[order(out$mean), ]
}

# Worst-cell regret: per cell average each arm over its repetitions, subtract
# the best arm's cell average, and report the largest shortfall any cell
# imposes on an arm. This is the robustness measure the July study used.
worstCell <- function(rows, score) {
  rows <- rows[!is.na(rows[[score]]), ]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  cellMeans <- aggregate(
    rows[[score]],
    by = list(cell = rows$cell, arm = rows$arm),
    FUN = mean
  )
  names(cellMeans)[3L] <- "value"
  best <- aggregate(
    cellMeans$value,
    by = list(cell = cellMeans$cell),
    FUN = min
  )
  names(best)[2L] <- "best"
  cellMeans <- merge(cellMeans, best, by = "cell")
  cellMeans$regret <- cellMeans$value - cellMeans$best
  parts <- split(cellMeans, cellMeans$arm)
  out <- do.call(
    rbind,
    lapply(names(parts), function(a) {
      d <- parts[[a]]
      worst <- which.max(d$regret)
      data.frame(
        arm = a,
        meanRegret = mean(d$regret),
        worstRegret = d$regret[worst],
        worstCell = d$cell[worst],
        stringsAsFactors = FALSE
      )
    })
  )
  out[order(out$worstRegret), ]
}

printTable <- function(x, digits = 4L) {
  if (is.null(x)) {
    cat("  (no rows)\n")
    return(invisible(NULL))
  }
  print(format(x, digits = digits), row.names = FALSE)
  cat("\n")
}

# Blocks written before the mixing columns existed are missing them, so the
# union of the column names is filled in rather than demanded.
bindBlocks <- function(pieces) {
  columns <- unique(unlist(lapply(pieces, names)))
  do.call(
    rbind,
    lapply(pieces, function(p) {
      attributes(p)[c("quick", "settings")] <- NULL
      for (missing in setdiff(columns, names(p))) {
        p[[missing]] <- NA
      }
      p[, columns, drop = FALSE]
    })
  )
}

readRun <- function(outDir) {
  files <- list.files(outDir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("no rds files in ", outDir)
  }
  pieces <- lapply(files, readRDS)
  quickFlags <- vapply(
    pieces,
    function(p) isTRUE(attr(p, "quick")),
    logical(1L)
  )
  if (length(unique(quickFlags)) > 1L) {
    stop("directory mixes quick and full blocks; separate them")
  }
  # A directory is one leg. Blocks may differ in their repetition counts, an
  # expensive stratum having been run at fewer, but a directory that mixed
  # chain lengths would have summarize report the first block's length for
  # all of them, so that is refused rather than labelled.
  lengths <- unique(vapply(
    pieces,
    function(p) {
      settings <- attr(p, "settings")
      sprintf(
        "%s chain(s) of %s after %s, %s trees",
        settings$nChains,
        settings$nSamples,
        settings$nBurn,
        settings$nTrees
      )
    },
    character(1L)
  ))
  if (length(lengths) > 1L) {
    stop(
      "directory mixes chain lengths; separate them: ",
      paste(lengths, collapse = " / ")
    )
  }
  rows <- bindBlocks(pieces)
  attr(rows, "settings") <- attr(pieces[[1L]], "settings")
  attr(rows, "quick") <- any(quickFlags)
  rows
}

summarizeRun <- function(outDir) {
  rows <- readRun(outDir)
  settings <- attr(rows, "settings")
  quickFlags <- attr(rows, "quick")

  cat(sprintf(
    "%d rows, %d arms, %d cells (%d simulated, %d real), %d trees, %d chain(s) of %d draws after %d burn%s\n\n",
    nrow(rows),
    length(unique(rows$arm)),
    length(unique(rows$cell)),
    length(unique(rows$cell[rows$kind == "sim"])),
    length(unique(rows$cell[rows$kind == "real"])),
    settings$nTrees,
    if (is.null(settings$nChains)) 1L else settings$nChains,
    settings$nSamples,
    settings$nBurn,
    if (isTRUE(quickFlags)) "  [QUICK]" else ""
  ))

  repCounts <- aggregate(rows$rep, by = list(n = rows$n), FUN = max)
  cat("repetitions per cell by sample size: ")
  cat(paste(sprintf("n=%d: %d", repCounts$n, repCounts$x), collapse = ", "))
  cat("\n\n")

  rows$covGap <- abs(rows$coverage - 0.9)
  sim <- rows[rows$kind == "sim", ]
  real <- rows[rows$kind == "real", ]

  cat("=== simulated cases: held-out log score ===\n")
  printTable(armSummary(sim, "logScore"))
  cat("=== simulated cases: held-out Brier score ===\n")
  printTable(armSummary(sim, "brier"))
  cat("=== simulated cases: 90 percent interval coverage of true p(x) ===\n")
  printTable(armSummary(sim, "coverage"))
  cat("=== simulated cases: distance of coverage from the nominal 0.9 ===\n")
  printTable(armSummary(sim, "covGap"))
  cat("=== simulated cases: 90 percent interval width ===\n")
  printTable(armSummary(sim, "width"))
  cat("=== simulated cases: RMSE of posterior mean probability ===\n")
  printTable(armSummary(sim, "probRmse"))
  cat("=== real datasets: held-out log score ===\n")
  printTable(armSummary(real, "logScore"))
  cat("=== real datasets: held-out Brier score ===\n")
  printTable(armSummary(real, "brier"))

  cat("=== worst-cell regret: log score, all cases ===\n")
  printTable(worstCell(rows, "logScore"))
  cat("=== worst-cell regret: Brier score, all cases ===\n")
  printTable(worstCell(rows, "brier"))
  cat(
    "=== worst-cell regret: distance of coverage from 0.9, simulated cases ===\n"
  )
  printTable(worstCell(sim, "covGap"))
  cat("=== worst-cell regret: probability RMSE, simulated cases ===\n")
  printTable(worstCell(sim, "probRmse"))

  cat("=== sampled k ===\n")
  kRows <- aggregate(
    rows[, c("kMedian", "kQ90", "kMax", "seconds")],
    by = list(arm = rows$arm),
    FUN = mean
  )
  kMax <- aggregate(rows$kMax, by = list(arm = rows$arm), FUN = max)
  names(kMax)[2L] <- "kMaxOverall"
  printTable(merge(kRows, kMax, by = "arm"))

  if ("kRhat" %in% names(rows) && any(!is.na(rows$kRhat))) {
    cat("=== mixing: sampled k and held-out mean probability ===\n")
    mix <- rows[!is.na(rows$pRhat), ]
    diagnostics <- do.call(
      rbind,
      lapply(split(mix, mix$arm), function(d) {
        data.frame(
          arm = d$arm[1L],
          kRhatMean = mean(d$kRhat),
          kRhatMax = suppressWarnings(max(d$kRhat)),
          kRhatOver1.05 = mean(d$kRhat > 1.05),
          kEssMedian = suppressWarnings(median(d$kEss)),
          pRhatMean = mean(d$pRhat),
          pRhatMax = max(d$pRhat),
          pEssMedian = median(d$pEss),
          stringsAsFactors = FALSE
        )
      })
    )
    printTable(diagnostics[order(diagnostics$pRhatMean), ])

    cat("=== mixing by sample size: incumbent and the fixed arms ===\n")
    focus <- mix[mix$arm %in% c(incumbent, "k = 2"), ]
    printTable(aggregate(
      focus[, c("kRhat", "kEss", "pRhat", "pEss", "coverage")],
      by = list(arm = focus$arm, n = focus$n),
      FUN = function(x) mean(x, na.rm = TRUE)
    ))

    cat("=== coverage against how well the fit mixed, simulated cases ===\n")
    simMix <- mix[mix$kind == "sim" & !is.na(mix$coverage), ]
    if (nrow(simMix) > 0L) {
      simMix$band <- cut(
        simMix$pRhat,
        breaks = c(0, 1.01, 1.05, 1.2, Inf),
        labels = c("<=1.01", "1.01-1.05", "1.05-1.2", ">1.2")
      )
      printTable(aggregate(
        simMix[, c("coverage", "covGap", "width", "kRhat")],
        by = list(band = simMix$band),
        FUN = function(x) mean(x, na.rm = TRUE)
      ))
      kBand <- simMix[!is.na(simMix$kRhat), ]
      if (nrow(kBand) > 0L) {
        kBand$band <- cut(
          kBand$kRhat,
          breaks = c(0, 1.01, 1.05, 1.2, Inf),
          labels = c("<=1.01", "1.01-1.05", "1.05-1.2", ">1.2")
        )
        cat("  banded on the sampled k instead:\n")
        printTable(aggregate(
          kBand[, c("coverage", "covGap", "width", "kEss")],
          by = list(band = kBand$band),
          FUN = function(x) mean(x, na.rm = TRUE)
        ))
      }
    }
  }

  cat("=== coverage by base rate, selected arms ===\n")
  byRate <- aggregate(
    sim[, c("coverage", "width", "logScore")],
    by = list(arm = sim$arm, rate = sim$rate),
    FUN = mean
  )
  printTable(byRate[order(byRate$arm, byRate$rate), ])

  cat("=== log score by sample size, all arms ===\n")
  bySize <- aggregate(
    sim[, c("logScore", "coverage", "kMedian")],
    by = list(arm = sim$arm, n = sim$n),
    FUN = mean
  )
  printTable(bySize[order(bySize$arm, bySize$n), ])

  cat("=== degrees of freedom at fixed scale: simulated averages ===\n")
  chiRows <- sim[!is.na(sim$df), ]
  byDf <- aggregate(
    chiRows[, c(
      "logScore",
      "brier",
      "coverage",
      "width",
      "probRmse",
      "kMedian"
    )],
    by = list(scale = chiRows$scale, df = chiRows$df),
    FUN = mean
  )
  printTable(byDf[order(byDf$scale, byDf$df), ])

  invisible(rows)
}

## ----------------------------------------------------------------- compare

settingsLine <- function(settings) {
  sprintf(
    "%d chain(s) of %d draws after %d burn",
    if (is.null(settings$nChains)) 1L else settings$nChains,
    settings$nSamples,
    settings$nBurn
  )
}

# Pairs two directories on (cell, repetition, arm). The data, the prior and
# the seed are the same on both sides and only the chain length differs, so
# the difference is what the length does and not what the prior does. Second
# directory minus first, so a positive coverage difference is the second
# configuration covering more.
compareRuns <- function(dirA, dirB) {
  a <- readRun(dirA)
  b <- readRun(dirB)
  scores <- c(
    "logScore",
    "brier",
    "coverage",
    "width",
    "probRmse",
    "kMedian",
    "kRhat",
    "kEss",
    "pRhat",
    "pEss"
  )
  scores <- intersect(scores, intersect(names(a), names(b)))
  keep <- c("cell", "rep", "arm", "kind", "n", scores)
  merged <- merge(
    a[, keep],
    b[, keep],
    by = c("cell", "rep", "arm"),
    suffixes = c(".a", ".b")
  )
  if (nrow(merged) == 0L) {
    stop("the two directories share no (cell, repetition, arm)")
  }
  # The pairing is an intersection, so a leg run over fewer cells, arms or
  # repetitions than the other contributes nothing for the rows it lacks.
  # The unpaired counts are reported so a reader sees how much was dropped.
  cat(sprintf(
    "A %s: %s
B %s: %s
%d paired fits over %d cells and %d arms
%d of A's %d rows and %d of B's %d rows had no partner

",
    dirA,
    settingsLine(attr(a, "settings")),
    dirB,
    settingsLine(attr(b, "settings")),
    nrow(merged),
    length(unique(merged$cell)),
    length(unique(merged$arm)),
    nrow(a) - nrow(merged),
    nrow(a),
    nrow(b) - nrow(merged),
    nrow(b)
  ))

  armTable <- function(rows, columns) {
    parts <- split(rows, rows$arm)
    out <- do.call(
      rbind,
      lapply(names(parts), function(armName) {
        d <- parts[[armName]]
        piece <- data.frame(arm = armName, stringsAsFactors = FALSE)
        for (column in columns) {
          x <- d[[paste0(column, ".a")]]
          y <- d[[paste0(column, ".b")]]
          ok <- !is.na(x) & !is.na(y)
          piece[[paste0(column, ".a")]] <- mean(x[ok])
          piece[[paste0(column, ".b")]] <- mean(y[ok])
          piece[[paste0(column, ".diff")]] <- mean(y[ok] - x[ok])
          piece[[paste0(column, ".diffSe")]] <- standardError(y[ok] - x[ok])
        }
        piece
      })
    )
    out
  }

  sim <- merged[merged$kind.a == "sim", ]
  real <- merged[merged$kind.a == "real", ]

  if (nrow(sim) > 0L) {
    cat("=== what the chain length does to coverage, simulated cases ===\n")
    printTable(armTable(sim, "coverage"))
    cat("=== and to interval width ===\n")
    printTable(armTable(sim, "width"))
    cat("=== and to log score, simulated cases ===\n")
    printTable(armTable(sim, "logScore"))
    cat("=== coverage by sample size ===\n")
    bySize <- do.call(
      rbind,
      lapply(split(sim, sim$n.a), function(d) {
        cbind(n = d$n.a[1L], armTable(d, "coverage"))
      })
    )
    printTable(bySize[order(bySize$arm, bySize$n), ])
    if ("kMedian" %in% scores) {
      cat("=== and to the sampled k ===\n")
      printTable(armTable(sim, "kMedian"))
    }
  }
  if (nrow(real) > 0L) {
    cat("=== what the chain length does to log score, real datasets ===\n")
    printTable(armTable(real, "logScore"))
  }
  invisible(merged)
}

## -------------------------------------------------------------------- main

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    cat(
      "usage: binary-hyperprior.R <block|blocks|all|summarize|compare> [outdir] [outdir2] [quick]\n"
    )
    quit(status = 1L)
  }
  command <- args[1L]
  quick <- "quick" %in% args
  if (quick) {
    quickSettings()
  }
  applyEnvironmentOverrides()
  rest <- setdiff(args[-1L], "quick")
  outDir <- if (length(rest) >= 1L) {
    rest[1L]
  } else {
    "benchmarks/results/binary-hyperprior"
  }
  cores <- defaultCores()

  if (command == "blocks") {
    cat(paste(allBlockNames(), collapse = "\n"), "\n", sep = "")
  } else if (command == "summarize") {
    summarizeRun(outDir)
  } else if (command == "compare") {
    if (length(rest) < 2L) {
      stop("compare needs two directories")
    }
    compareRuns(rest[1L], rest[2L])
  } else if (command == "all") {
    for (block in allBlockNames()) {
      runBlock(block, outDir, cores, quick)
    }
  } else {
    runBlock(command, outDir, cores, quick)
  }
}

main()
