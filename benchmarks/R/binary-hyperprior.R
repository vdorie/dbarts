#!/usr/bin/env Rscript

# Re-evaluates the binary (probit) end-node hyperprior default, chi(1.5, 2).
#
# The July 2026 study (docs/plans/archive/chi-default-research.md) held the
# degrees of freedom at 1.5, varied only the scale, and ran four simulated
# data-generating processes. This harness varies BOTH parameters and widens
# the case set: six simulated processes crossed with sample size, predictor
# count and base rate, plus six real datasets from R and its recommended
# packages scored by repeated 80/20 splits.
#
# ARMS. Twenty hyperprior arms, chi(df, scale) for df in {1, 1.25, 1.5, 2, 3}
# crossed with scale in {1, 2, 5, Inf}, and three fixed-k arms, k in {1, 2, 3}.
# k = 2 is the fixed value BayesTree used and the value the chi(1.5, 2) prior
# is centered near. Arms are paired: within a case and a repetition every arm
# sees the same data and starts from the same MCMC seed, so a difference
# between arms is a difference in the prior, not in the draw.
#
# SCORES, all on held-out rows. Log score and Brier score against the held-out
# outcome; for the simulated cases, where the true probability is known, the
# coverage and mean width of the 90 percent posterior interval for that
# probability, and the root mean squared error of the posterior mean
# probability. Every arm also reports the posterior median, 90th percentile
# and maximum of the sampled k (a fixed arm reports its own k) and the fit's
# elapsed seconds.
#
# BLOCKS. One block is one invocation and writes one rds, so a full run splits
# across a session. A simulated block is a data-generating process at one
# sample size ("sim:friedman:500"), which is nine cells - three predictor
# counts by three base rates - at the repetition count below; appending a
# predictor count ("sim:friedman:2000:50") narrows it to that count's three
# cells and its own rds, which is how the large-n blocks are kept short. A
# real block is one dataset ("real:biopsy"). `blocks` lists them, `all` runs
# every one in sequence, and `summarize` reads a directory of rds files and
# prints the tables the plan doc reports.
#
# Usage:
#   Rscript benchmarks/R/binary-hyperprior.R blocks
#   Rscript benchmarks/R/binary-hyperprior.R sim:friedman:500 [outdir] [quick]
#   Rscript benchmarks/R/binary-hyperprior.R real:biopsy [outdir]
#   Rscript benchmarks/R/binary-hyperprior.R all [outdir]
#   Rscript benchmarks/R/binary-hyperprior.R summarize [outdir]
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

makeArms <- function() {
  dfs <- c(1, 1.25, 1.5, 2, 3)
  scales <- c(1, 2, 5, Inf)
  arms <- list()
  for (scale in scales) {
    for (df in dfs) {
      arms[[length(arms) + 1L]] <- list(
        name = sprintf("chi(%s, %s)", format(df), format(scale)),
        kind = "chi",
        df = df,
        scale = scale,
        fixed.k = NA_real_
      )
    }
  }
  for (k in c(1, 2, 3)) {
    arms[[length(arms) + 1L]] <- list(
      name = sprintf("k = %s", format(k)),
      kind = "fixed",
      df = NA_real_,
      scale = NA_real_,
      fixed.k = k
    )
  }
  arms
}

arms <- makeArms()
incumbent <- "chi(1.5, 2)"

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
# one dichotomized survival outcome.
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

## ------------------------------------------------------------ fit and score

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
      node.prior = nodePriorFor(arm)
    )
  )[["elapsed"]]

  probabilities <- pnorm(fit$yhat.test)
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

  kDraws <- if (is.null(fit$k)) {
    rep(arm$fixed.k, nSamples * nChains)
  } else {
    as.vector(fit$k)
  }
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

runSimBlock <- function(dgp, n, cores, counts = predictorCounts) {
  rows <- list()
  for (p in counts) {
    for (rate in baseRates) {
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
  nTrain <- floor(0.8 * n)
  label <- sprintf(
    "%s|%d|%d|%s",
    name,
    n,
    ncol(raw$x),
    format(round(mean(raw$y), 3))
  )
  cat(sprintf("  %-28s ", label))
  caseFor <- function(seed) {
    for (attempt in seq_len(100L)) {
      set.seed(seed + 1000000L * (attempt - 1L))
      trainRows <- sample.int(n, nTrain)
      if (length(unique(raw$y[trainRows])) == 2L) break
    }
    testRows <- setdiff(seq_len(n), trainRows)
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
    runSimBlock(parts[2L], as.integer(parts[3L]), cores, counts)
  } else if (parts[1L] == "real") {
    runRealBlock(parts[2L], cores)
  } else {
    stop("unknown block: ", block)
  }
  attr(result, "quick") <- quick
  attr(result, "settings") <- list(
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

summarizeRun <- function(outDir) {
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
  rows <- do.call(
    rbind,
    lapply(pieces, function(p) {
      attributes(p)[c("quick", "settings")] <- NULL
      p
    })
  )
  settings <- attr(pieces[[1L]], "settings")

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
    if (any(quickFlags)) "  [QUICK]" else ""
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

## -------------------------------------------------------------------- main

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    cat(
      "usage: binary-hyperprior.R <block|blocks|all|summarize> [outdir] [quick]\n"
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
  } else if (command == "all") {
    for (block in allBlockNames()) {
      runBlock(block, outDir, cores, quick)
    }
  } else {
    runBlock(command, outDir, cores, quick)
  }
}

main()
