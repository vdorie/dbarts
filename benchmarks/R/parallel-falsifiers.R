#!/usr/bin/env Rscript

# Drivers for the first three falsifiers of
# docs/design/parallel-bart-frontier.md section 5. These exercise an
# INSTRUMENTED build of the engine (src/bartcore counters that never land):
# the worktree engine reads BC_FALSIFIER_MODE and writes CSV rows keyed by
# BC_FALSIFIER_TAG to <BC_FALSIFIER_OUT>-e{1,2,3}.csv. On the shipped engine
# these env vars are inert (the instrumentation is reverted before commit), so
# this script only reproduces numbers against a checkout that still carries it.
#
#   E2  field-fraction profile: per-sweep DRAM-byte traffic split between
#       residual-field maintenance (removed by 3.4's atom map) and per-move
#       scan work (kept). Kills 3.4 if field is not the dominant (~90%) share.
#   E1  atom census: distinct per-observation leaf-assignment tuples over
#       consecutive tree blocks b in {4,8,16} on real fitted forests. Kills
#       3.4 if atoms approach n at small b.
#   E3  stale-residual agreement: per proposed move, the accept decision under
#       the frozen start-of-sweep residual vs the true rolled residual (same
#       proposal randomness). Kills 3.2a if they disagree often OR the frozen
#       (stage-1) survivor rate is near 1.
#
# Usage:  Rscript parallel-falsifiers.R {e1|e2|e3|all} [outdir]

suppressMessages(library(dbarts))

args    <- commandArgs(trailingOnly = TRUE)
which   <- if (length(args) >= 1L) args[1L] else "all"
outdir  <- if (length(args) >= 2L) args[2L] else tempfile("falsifier-")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outpref <- file.path(outdir, "fal")
Sys.setenv(BC_FALSIFIER_OUT = outpref)

## Friedman-1 with 5 noise columns; fixed seed per n.
makeFriedman <- function(n, p = 10L, seed = 11L) {
  set.seed(seed)
  x <- matrix(runif(n * p), n, p)
  Ey <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
  list(x = x, y = Ey + rnorm(n))
}

## config grid: default m = 75 and one deep-tree config (lower power => deeper)
configs <- list(
  list(tag = "default", ntree = 75L, k = 2.0, base = 0.95, power = 2.0),
  list(tag = "deep",    ntree = 75L, k = 2.0, base = 0.95, power = 1.0)
)
ns   <- c(1e4, 1e5)
seed <- 99L

runOne <- function(mode, cfg, n, ndpost, nskip, keepevery = 1L,
                   keeptrees = FALSE) {
  d <- makeFriedman(as.integer(n))
  Sys.setenv(BC_FALSIFIER_MODE = mode)
  Sys.setenv(BC_FALSIFIER_TAG  = sprintf("%s_n%g", cfg$tag, n))
  invisible(bart(d$x, d$y, ntree = cfg$ntree, k = cfg$k, base = cfg$base,
                 power = cfg$power, ndpost = ndpost, nskip = nskip,
                 keepevery = keepevery, keeptrees = keeptrees,
                 keeptrainfits = FALSE, verbose = FALSE, nthread = 1L,
                 seed = seed))
  Sys.unsetenv("BC_FALSIFIER_MODE")
}

run_e2 <- function() {
  cat("== E2 field-fraction profile ==\n")
  for (cfg in configs) for (n in ns)
    runOne("e2", cfg, n, ndpost = 200L, nskip = 200L)
  f <- paste0(outpref, "-e2.csv")
  if (!file.exists(f)) { cat("(no e2 output)\n"); return(invisible()) }
  h <- c("tag","n","m","iters","burn","sweeps","residualRoll","nodeAverages",
         "fitWrite","totalRebuild","response","bookkeep","scanPartRd",
         "scanPartWr","scanLeafStats","field","scanHi")
  e2 <- read.csv(f, header = FALSE, col.names = h)
  e2$scanLo <- e2$scanPartRd + e2$scanLeafStats
  e2$fieldPctHi <- 100 * e2$field / (e2$field + e2$scanHi)  # scan upper bnd
  e2$fieldPctLo <- 100 * e2$field / (e2$field + e2$scanLo)  # scan lower bnd
  for (i in seq_len(nrow(e2))) {
    r <- e2[i, ]
    cat(sprintf(
      "%-12s field%%=[%.1f,%.1f]  field/sweep=%.2e scan/sweep=[%.2e,%.2e]\n",
      r$tag, r$fieldPctHi, r$fieldPctLo, r$field / r$sweeps,
      r$scanLo / r$sweeps, r$scanHi / r$sweeps))
    cat(sprintf("             roll=%.2e nodeAvg=%.2e fit=%.2e rebuild=%.2e resp=%.2e book=%.2e | partRd=%.2e partWr=%.2e leafStat=%.2e\n",
      r$residualRoll/r$sweeps, r$nodeAverages/r$sweeps, r$fitWrite/r$sweeps,
      r$totalRebuild/r$sweeps, r$response/r$sweeps, r$bookkeep/r$sweeps,
      r$scanPartRd/r$sweeps, r$scanPartWr/r$sweeps, r$scanLeafStats/r$sweeps))
  }
  invisible(e2)
}

run_e1 <- function() {
  cat("== E1 atom census ==\n")
  for (cfg in configs) for (n in ns)
    runOne("e1", cfg, n, ndpost = 10L, nskip = 200L, keepevery = 20L,
           keeptrees = TRUE)
  f <- paste0(outpref, "-e1.csv")
  if (!file.exists(f)) { cat("(no e1 output)\n"); return(invisible()) }
  e1 <- read.csv(f, header = FALSE,
                 col.names = c("tag","n","m","sample","start","b","atoms"))
  agg <- aggregate(atoms ~ tag + n + b, data = e1, FUN = mean)
  agg <- agg[order(agg$tag, agg$n, agg$b), ]
  for (i in seq_len(nrow(agg))) {
    r <- agg[i, ]
    cat(sprintf("%-12s n=%-7g b=%-2d  atoms=%.0f  atoms/n=%.3f  occ=n/atoms=%.1f\n",
      r$tag, r$n, r$b, r$atoms, r$atoms / r$n, r$n / r$atoms))
  }
  invisible(e1)
}

run_e3 <- function() {
  cat("== E3 stale-residual agreement ==\n")
  for (cfg in configs) for (n in ns)
    runOne("e3", cfg, n, ndpost = 100L, nskip = 100L)
  f <- paste0(outpref, "-e3.csv")
  if (!file.exists(f)) { cat("(no e3 output)\n"); return(invisible()) }
  e3 <- read.csv(f, header = FALSE,
                 col.names = c("tag","n","m","move","proposed","degenerate",
                               "agree","survivor","trueAccept","survTrue"))
  ## collapse move types per config for an overall line, and per-type
  for (tg in unique(e3$tag)) {
    s <- e3[e3$tag == tg, ]
    tot <- sum(s$proposed)
    if (tot == 0) next
    cat(sprintf("-- %s (proposed=%g, degenerate=%g)\n", tg, tot,
                sum(s$degenerate)))
    cat(sprintf("   OVERALL agree=%.3f survivor=%.3f trueAccept=%.3f\n",
        sum(s$agree)/tot, sum(s$survivor)/tot, sum(s$trueAccept)/tot))
    for (i in seq_len(nrow(s))) {
      r <- s[i, ]
      if (r$proposed == 0) next
      cat(sprintf("   %-7s prop=%-8g agree=%.3f survivor=%.3f trueAcc=%.3f\n",
          r$move, r$proposed, r$agree/r$proposed, r$survivor/r$proposed,
          r$trueAccept/r$proposed))
    }
  }
  invisible(e3)
}

cat("output prefix:", outpref, "\n")
if (which %in% c("e2","all")) run_e2()
if (which %in% c("e1","all")) run_e1()
if (which %in% c("e3","all")) run_e3()
