#!/usr/bin/env Rscript

# Stage 2 of change-move-fix (docs/plans/change-move-fix.md): mixing cost of the
# two exact change kernels ("prior", "hybrid") against the shipped defect
# ("current"). The balance gate (benchmarks/R/change-balance.R) must pass for
# both repairs before this runs; here we measure effective sample size per
# second and per sweep on the stage-1 config grid, three kernels per config,
# same seeds, strictly serial single chain/thread.
#
# Kernel is selected via DBARTS_CHANGE_KERNEL, read once when dbarts() builds
# the C++ sampler. Change-move acceptance and no-op counts come from the
# env-gated engine counters (DBARTS_CHANGE_STATS -> a per-run CSV, written at
# the end of each $run); they are zero-overhead and bit-identical when unset.
#
# ESS estimator: coda::effectiveSize throughout (one consistent single-chain
# spectral estimator). ESS per sweep = ESS / nSamples (n.thin = 1, so one kept
# draw per sweep); ESS per second uses the wall time of the sampling run only.
#
# Usage: Rscript change-fix-stage2.R [quick]

suppressPackageStartupMessages({
  library(dbarts)
  library(coda)
})

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nBurn <- if (quick) 300L else 1000L
nSamples <- if (quick) 800L else 3000L
seed <- 20260708L
kernels <- c("current", "prior", "hybrid")

# ---------------------------------------------------------------------------
# data generators (verbatim from change-fix-instrumentation.R)
# ---------------------------------------------------------------------------

makeFriedman <- function(n, p = 10L, s = 1L) {
  set.seed(s)
  x <- matrix(runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5] + rnorm(n)
  list(x = x, y = y)
}

makeMixed <- function(n, s = 2L) {
  set.seed(s)
  cont <- as.data.frame(matrix(runif(n * 5L), n, 5L))
  names(cont) <- paste0("c", 1:5)
  c4 <- factor(sample.int(4L, n, replace = TRUE))
  c8 <- factor(sample.int(8L, n, replace = TRUE))
  c16 <- factor(sample.int(16L, n, replace = TRUE))
  # continuous signal plus categorical effects that reward same-variable stacking
  eff8 <- (as.integer(c8) - 4.5) * 0.6
  eff16 <- (as.integer(c16) %% 4L) * 1.0
  y <- 10 * sin(pi * cont$c1 * cont$c2) + 10 * cont$c3 +
    (as.integer(c4) == 1L) * 4 + eff8 + eff16 + rnorm(n)
  x <- cbind(cont, c4 = c4, c8 = c8, c16 = c16)
  list(x = x, y = y)
}

# ---------------------------------------------------------------------------
# ESS helpers
# ---------------------------------------------------------------------------

# effectiveSize of one trace, guarding constant / degenerate series
essOne <- function(v) {
  if (length(v) < 8L || !is.finite(stats::var(v)) || stats::var(v) <= 0)
    return(NA_real_)
  as.numeric(coda::effectiveSize(v))
}

# mean ESS over the rows of a (quantity x nSamples) matrix, dropping NAs
essRowsMean <- function(m) {
  vals <- apply(m, 1L, essOne)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) NA_real_ else mean(vals)
}

essRowsMin <- function(m) {
  vals <- apply(m, 1L, essOne)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) NA_real_ else min(vals)
}

# ---------------------------------------------------------------------------
# one (config, kernel) fit
# ---------------------------------------------------------------------------

baseControl <- function(ntree) {
  dbartsControl(
    n.chains = 1L, n.threads = 1L, n.trees = ntree,
    n.samples = nSamples, n.burn = nBurn, n.thin = 1L,
    keepTrainingFits = TRUE, keepTrees = FALSE, updateState = FALSE,
    rngSeed = seed
  )
}

runFit <- function(name, x, y, ntree, kernel, treePriorExpr = NULL) {
  statsPath <- tempfile("changefix2-stats-", fileext = ".csv")
  Sys.setenv(DBARTS_CHANGE_KERNEL = kernel)
  Sys.setenv(DBARTS_CHANGE_STATS = statsPath)
  on.exit({
    Sys.unsetenv("DBARTS_CHANGE_KERNEL")
    Sys.unsetenv("DBARTS_CHANGE_STATS")
    unlink(statsPath)
  }, add = TRUE)

  ctl <- baseControl(ntree)
  s <- if (is.null(treePriorExpr)) {
    dbarts(x, y, control = ctl)
  } else {
    eval(bquote(dbarts(x, y, control = ctl, tree.prior = .(treePriorExpr))))
  }

  invisible(s$run(nBurn, 0L))            # burn-in, not timed; stats overwritten
  t0 <- Sys.time()
  samp <- s$run(0L, nSamples)            # sampling run; its stats are dumped
  seconds <- as.numeric(Sys.time() - t0, units = "secs")

  stats <- utils::read.csv(statsPath)
  proposals <- stats$proposals[1]
  acceptRate <- if (proposals > 0) stats$accepted[1] / proposals else NA_real_
  noopRate <- if (proposals > 0) stats$noop[1] / proposals else NA_real_

  # a fixed handful of train-fit coordinates, spread across the rows
  n <- nrow(samp$train)
  coords <- unique(round(seq(1, n, length.out = 6L)))
  trainMat <- samp$train[coords, , drop = FALSE]

  essSigma <- essOne(samp$sigma)
  essK <- if (!is.null(samp$k)) essOne(samp$k) else NA_real_
  essTrainMean <- essRowsMean(trainMat)
  essTrainMin <- essRowsMin(trainMat)
  essVarMean <- essRowsMean(samp$varcount)
  essVarMin <- essRowsMin(samp$varcount)

  perSec <- function(e) if (is.finite(e)) e / seconds else NA_real_
  perSweep <- function(e) if (is.finite(e)) e / nSamples else NA_real_

  cat(sprintf(
    "  %-8s %6.2fs  accept=%.3f noop=%.4f | ESS/s sigma=%.1f train=%.1f var=%.1f\n",
    kernel, seconds, acceptRate, noopRate,
    perSec(essSigma), perSec(essTrainMean), perSec(essVarMean)))

  data.frame(
    config = name, kernel = kernel, ntree = ntree,
    nBurn = nBurn, nSamples = nSamples, seconds = seconds,
    changeProposals = proposals, changeAcceptRate = acceptRate,
    noopRate = noopRate,
    essSigma = essSigma, essSigmaPerSec = perSec(essSigma),
    essSigmaPerSweep = perSweep(essSigma),
    essK = essK, essKPerSec = perSec(essK), essKPerSweep = perSweep(essK),
    essTrainMean = essTrainMean, essTrainMeanPerSec = perSec(essTrainMean),
    essTrainMeanPerSweep = perSweep(essTrainMean), essTrainMin = essTrainMin,
    essVarMean = essVarMean, essVarMeanPerSec = perSec(essVarMean),
    essVarMeanPerSweep = perSweep(essVarMean), essVarMin = essVarMin,
    stringsAsFactors = FALSE
  )
}

runConfig <- function(name, x, y, ntree, treePriorExpr = NULL) {
  cat(sprintf("\n=== %s (ntree=%d) ===\n", name, ntree))
  do.call(rbind, lapply(kernels, function(kn)
    runFit(name, x, y, ntree, kn, treePriorExpr)))
}

# ---------------------------------------------------------------------------
# config grid (stage-1 grid; generators reused verbatim)
# ---------------------------------------------------------------------------

rows <- list()
add <- function(df) rows[[length(rows) + 1L]] <<- df

# (1) Friedman-1, continuous-only, defaults
for (n in c(1000L, 10000L)) {
  fr <- makeFriedman(n)
  for (m in c(75L, 200L))
    add(runConfig(sprintf("Friedman n=%d m=%d", n, m), fr$x, fr$y, m))
}

# (2) mixed 5 continuous + 3 categorical (4,8,16 levels), n=1e4, m=75
mx <- makeMixed(10000L)
add(runConfig("Mixed n=1e4 m=75", mx$x, mx$y, 75L))

# (3) deep single-tree stress: ntree=1, power=1, n=1e3
fr3 <- makeFriedman(1000L)
add(runConfig("Deep single-tree n=1e3 power=1", fr3$x, fr3$y, 1L,
  treePriorExpr = quote(cgm(power = 1, base = 0.95))))

# (4) config (2) with DART enabled
add(runConfig("Mixed+DART n=1e4 m=75", mx$x, mx$y, 75L,
  treePriorExpr = quote(dart(power = 2, base = 0.95))))

# ---------------------------------------------------------------------------
# combined table + CSV
# ---------------------------------------------------------------------------

tab <- do.call(rbind, rows)
cat("\n================ STAGE 2 ESS SUMMARY (per second) ================\n")
compact <- tab[, c("config", "kernel", "seconds", "changeAcceptRate",
  "noopRate", "essSigmaPerSec", "essTrainMeanPerSec", "essVarMeanPerSec")]
print(compact, digits = 4, row.names = FALSE)

cat("\n================ per sweep ================\n")
compactSweep <- tab[, c("config", "kernel", "essSigmaPerSweep",
  "essTrainMeanPerSweep", "essVarMeanPerSweep")]
print(compactSweep, digits = 4, row.names = FALSE)

outCsv <- file.path("benchmarks", "results", "change-fix-stage2.csv")
dir.create(dirname(outCsv), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, outCsv, row.names = FALSE)
cat(sprintf("\nwrote %s\n", outCsv))
