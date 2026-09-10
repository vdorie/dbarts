#!/usr/bin/env Rscript

# Probes the xint_t caps (src/bartcore/data.hpp): maxNumCutsRepresentable
# (65533 cuts on an ordinal column), maxCategories (65535 categorical levels)
# and maxLevelsForKind, which stops an ordered factor one level lower because
# its K - 1 midpoint grid spends a code per cut. All three fall out of the
# 16-bit predictor code, whose widening the plan puts out of scope.
#
# This is a PROBE, not a timing sweep: it asks what each user-reachable
# channel does at, and one past, its cap - a silent clamp, a named refusal or
# a wrong answer - and then asks whether a design anyone writes gets near one.
# The last question is the one that decides whether the caps bind: the
# per-column cut count is min(n.cuts, distinct values - 1) and n.cuts
# defaults to 100, so a column reaches 65533 cuts only when a caller asks for
# it on a design with more rows than that. The two fit-time rows at the end
# say what asking costs.
#
# Usage: Rscript benchmarks/R/constant-xint-caps.R [quick]
# No baseline, no pass/fail exit status: an informational probe, run by hand
# when the caps are revisited.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

maxNumCutsRepresentable <- 65533L
maxCategories <- 65535L

report <- function(channel, value, outcome) {
  cat(sprintf("%-34s %10s  %s\n", channel, format(value), outcome))
}

outcomeOf <- function(expr) {
  tryCatch(
    {
      force(expr)
      "accepted"
    },
    error = function(e) paste0("REFUSED: ", conditionMessage(e))
  )
}

n <- 70000L
set.seed(43)
x <- matrix(runif(n * 2L), n, 2L)
colnames(x) <- c("x1", "x2")
y <- x[, 1L] + rnorm(n)
baseControl <- function(nCuts) {
  dbartsControl(
    n.cuts = nCuts,
    useQuantiles = TRUE,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    n.burn = 0L,
    updateState = FALSE,
    verbose = FALSE
  )
}

cat("channel                              request  outcome\n")
for (nCuts in c(65532L, 65533L, 65534L, 100000L)) {
  report(
    "control n.cuts",
    nCuts,
    outcomeOf(dbarts(x, y, control = baseControl(nCuts)))
  )
}

sampler <- dbarts(x, y, control = baseControl(100L))
for (nCuts in c(65533L, 65534L)) {
  cuts <- seq.int(0, 1, length.out = nCuts + 2L)[-c(1L, nCuts + 2L)]
  report(
    "sampler$setCutPoints length",
    nCuts,
    outcomeOf(sampler$setCutPoints(list(cuts), 1L))
  )
}
rm(sampler)

for (K in c(65534L, 65535L, 65536L)) {
  levelValues <- as.character(seq_len(K))
  g <- factor(
    levelValues[c(seq_len(K), rep.int(1L, 10L))],
    levels = levelValues
  )
  yg <- rnorm(length(g))
  report(
    "categorical predictor levels",
    K,
    outcomeOf(dbartsData(yg ~ g, data.frame(g = g, yg = yg)))
  )
  rm(g, yg, levelValues)
}
for (K in c(65533L, 65534L, 65535L)) {
  levelValues <- as.character(seq_len(K))
  g <- factor(
    levelValues[c(seq_len(K), rep.int(1L, 10L))],
    levels = levelValues,
    ordered = TRUE
  )
  yg <- rnorm(length(g))
  report(
    "ordered factor predictor levels",
    K,
    outcomeOf(dbartsData(yg ~ g, data.frame(g = g, yg = yg)))
  )
  rm(g, yg, levelValues)
}

if (!quick) {
  cat("\nwhat asking for the widest grid costs (n = 70000, 5 trees):\n")
  cat(sprintf("%12s %14s\n", "n.cuts", "msec/iter"))
  for (nCuts in c(100L, maxNumCutsRepresentable)) {
    fitSampler <- dbarts(x, y, control = baseControl(nCuts))
    invisible(fitSampler$run(5L, 0L))
    start <- Sys.time()
    invisible(fitSampler$run(0L, 20L))
    elapsed <- as.numeric(Sys.time() - start, units = "secs")
    cat(sprintf("%12d %14.3f\n", nCuts, elapsed / 20 * 1000))
    rm(fitSampler)
  }
}

cat(
  "\nnote: a cut request past its cap is accepted and silently clamped where",
  "\nthe store sizes its grid (ColumnStore, maxNumCutsRepresentable), while a",
  "\nlevel count past its cap is a named refusal at the bridge",
  "(refuseLevelCountPastCeiling).\n"
)
cat(sprintf(
  "\ncaps: %d cuts, %d categorical levels, %d ordered-factor levels\n",
  maxNumCutsRepresentable,
  maxCategories,
  maxNumCutsRepresentable + 1L
))
