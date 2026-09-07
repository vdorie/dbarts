#!/usr/bin/env Rscript

# Stage 0 of the tree-mixing falsifier: the move census
# (docs/design/tree-mixing-proposals.md section 6.1). Pure measurement - no
# kill criterion, no engine behaviour change - answering the three questions
# the veto-scaffold counts left open: the distribution of the log-likelihood
# difference among REJECTED structural proposals, change-move proposals and
# acceptances by target node depth, and the cut displacement that puts a
# same-variable cut move at a workable acceptance rate.
#
# The numbers come from the engine's own acceptance terms, appended one line
# per structural proposal by scaffolding in the move kernels that is compiled
# only under -DBARTCORE_MOVE_CENSUS. THIS SCRIPT DOES NOTHING against an
# ordinary build: the library must be built with the flag.
#
# Building the instrumented library. src/Makevars sets PKG_CPPFLAGS with '=',
# so the flag rides CPPFLAGS instead, which R composes into ALL_CPPFLAGS and
# no package makefile touches. R CMD INSTALL reads the user Makevars LAST, so
# an append there survives:
#
#   printf 'CPPFLAGS += -DBARTCORE_MOVE_CENSUS\n' > /tmp/census-makevars
#   mkdir -p /tmp/census-lib
#   R_MAKEVARS_USER=/tmp/census-makevars \
#     R CMD INSTALL --preclean --library=/tmp/census-lib .
#   R_LIBS=/tmp/census-lib Rscript benchmarks/R/move-census.R
#
# A private --library keeps the instrumented build off the ordinary one; the
# flag also reaches misc/external/rc, which ignore it. The run mode stops with
# a diagnostic if the first cell writes no census file, which is what an
# uninstrumented library looks like from R.
#
# Four cells, 200 burn-in plus 500 sampled sweeps, one chain, one thread,
# fixed seeds:
#
#   default   dbarts defaults: Friedman, n = 5000, p = 10, m = 75, sigma = 1
#   lownoise  the same at sigma^2 = 0.1, where structure freezing is
#             established (section 3.2)
#   wide      Friedman p = 50 with 45 noise columns, the regime that killed
#             the warm-start default
#   bcf       the causal-forest strong-scale cell: two forests, prognostic
#             amplitude 8 sigma, which is docs/design/bcf.md's strong-|a|
#             regime (|a|/sigma large, where the mu forest's structure mixes
#             slowly at high SNR)
#
# Usage:
#   Rscript move-census.R                          run every cell, summarize
#   Rscript move-census.R run [dir] [cell ...]     run cells, write census
#   Rscript move-census.R summarize [dir]          summarize existing files
# Append 'quick' for a smoke test (fewer sweeps, smaller n; not comparable).
#
# One cell runs per R process, spawned by the run mode: the engine opens the
# census file once, on the first record, so a second cell in the same process
# would append to the first cell's file.
#
# Summaries are taken over the SAMPLED sweeps only; burn-in acceptance is
# reported beside them for context.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
modes <- c("run", "summarize", "runcell")
mode <- if (length(args) >= 1L && args[[1L]] %in% modes) args[[1L]] else "both"
args <- setdiff(args, modes)
outputDir <- if (length(args) >= 1L) args[[1L]] else "benchmarks/census"
cellArgs <- if (length(args) >= 2L) args[-1L] else character()

nBurn <- if (quick) 20L else 200L
nSamples <- if (quick) 50L else 500L
nObservations <- if (quick) 500L else 5000L
nTrees <- 75L
dataSeed <- 20260907L
samplerSeed <- 7L

# ---------------------------------------------------------------- generators

genFriedman <- function(n, p, sigma) {
  x <- matrix(runif(n * p), n, p)
  f <- 10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
  list(x = x, y = f + rnorm(n, sd = sigma))
}

# Hahn, Murray and Carvalho's shape with the prognostic surface standardized
# and rescaled to strength * sigma, so the cell's |a|/sigma is set explicitly;
# the propensity is a function of mu, so the confounding rides the amplitude.
genCausal <- function(n, p, sigma, strength) {
  x <- matrix(runif(n * p), n, p)
  mu <- 2 * sin(pi * x[, 1L] * x[, 2L]) + 2 * (x[, 3L] - 0.5)^2 + x[, 4L]
  mu <- strength * sigma * (mu - mean(mu)) / sd(mu)
  tau <- 1 + 2 * x[, 3L]
  propensity <- 0.05 + 0.8 * pnorm(mu / (strength * sigma) - 0.5 * x[, 1L])
  z <- rbinom(n, 1L, propensity)
  list(x = x, z = z, y = mu + z * tau + rnorm(n, sd = sigma))
}

censusControl <- function() {
  dbartsControl(
    verbose = FALSE,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = nSamples,
    updateState = FALSE,
    seed = samplerSeed
  )
}

# ---------------------------------------------------------------- the cells

cells <- list(
  default = function() {
    data <- genFriedman(nObservations, 10L, 1)
    dbarts(data$x, data$y, control = censusControl())
  },
  lownoise = function() {
    data <- genFriedman(nObservations, 10L, sqrt(0.1))
    dbarts(data$x, data$y, control = censusControl())
  },
  wide = function() {
    data <- genFriedman(nObservations, 50L, 1)
    dbarts(data$x, data$y, control = censusControl())
  },
  bcf = function() {
    data <- genCausal(nObservations, 10L, 1, strength = 8)
    z <- data$z
    dbarts(
      data$x,
      data$y,
      forests = list(forest(), forest(basis = ~ factor(z))),
      control = censusControl()
    )
  }
)

censusFile <- function(dir, cell) {
  file.path(dir, sprintf("move-census-%s.csv", cell))
}

runCell <- function(dir, cell) {
  file <- censusFile(dir, cell)
  unlink(file)
  Sys.setenv(BARTCORE_MOVE_CENSUS_FILE = file)
  on.exit(Sys.unsetenv("BARTCORE_MOVE_CENSUS_FILE"), add = TRUE)
  set.seed(dataSeed)
  sampler <- cells[[cell]]()
  elapsed <- system.time(invisible(sampler$run(nBurn, nSamples)))[["elapsed"]]
  rm(sampler)
  invisible(gc(FALSE))
  # the engine's stream is still buffered here, so existence is the test: an
  # uninstrumented library never opens the file at all
  if (!file.exists(file)) {
    stop(
      "no census records for cell '",
      cell,
      "': the loaded dbarts was not built with -DBARTCORE_MOVE_CENSUS ",
      "(see this file's header)"
    )
  }
  cat(sprintf(
    "%-9s %6.1f s  %s (%.1f MB)\n",
    cell,
    elapsed,
    file,
    file.size(file) / 1024^2
  ))
}

# ------------------------------------------------------------- the summaries

proposalNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "move",
  "noop",
  "accepted",
  "nodeDepth",
  "treeDepth",
  "interior",
  "logLik",
  "logPrior",
  "logCorr"
)
cutNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "nodeDepth",
  "displacement",
  "logRatio"
)

readCensus <- function(file) {
  lines <- readLines(file)
  kind <- substr(lines, 1L, 1L)
  list(
    proposals = read.csv(
      text = lines[kind == "p"],
      header = FALSE,
      col.names = proposalNames
    ),
    cuts = read.csv(
      text = lines[kind == "d"],
      header = FALSE,
      col.names = cutNames
    )
  )
}

# proposals, no-op rate, and the two acceptance denominators the veto scaffold
# separated: per proposal made and per proposal that reached a score
moveTable <- function(p) {
  by <- split(p, p$move)
  rows <- lapply(names(by), function(move) {
    m <- by[[move]]
    scored <- sum(m$noop == 0L)
    data.frame(
      move = move,
      proposals = nrow(m),
      noop.pct = 100 * mean(m$noop == 1L),
      accept.pct = 100 * mean(m$accepted == 1L),
      accept.scored.pct = if (scored > 0L) {
        100 * sum(m$accepted) / scored
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  scored <- sum(p$noop == 0L)
  rbind(
    out,
    data.frame(
      move = "all",
      proposals = nrow(p),
      noop.pct = 100 * mean(p$noop == 1L),
      accept.pct = 100 * mean(p$accepted == 1L),
      accept.scored.pct = if (scored > 0L) {
        100 * sum(p$accepted) / scored
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  )
}

# The fork the whole program turns on: are rejections close calls (scale, so
# the temperature family is live) or wrong proposals (only better-aimed
# proposals help)?
rejectionTable <- function(p) {
  r <- p[p$noop == 0L & p$accepted == 0L, ]
  rows <- lapply(c(sort(unique(r$move)), "all"), function(move) {
    d <- if (move == "all") r$logLik else r$logLik[r$move == move]
    finite <- d[is.finite(d)]
    q <- quantile(finite, c(0.01, 0.05, 0.25, 0.5, 0.75), names = FALSE)
    data.frame(
      move = move,
      rejected = length(d),
      vetoed.pct = 100 * mean(!is.finite(d)),
      q01 = q[1L],
      q05 = q[2L],
      q25 = q[3L],
      q50 = q[4L],
      q75 = q[5L],
      within1.pct = 100 * mean(abs(finite) <= 1),
      within2.pct = 100 * mean(abs(finite) <= 2),
      within5.pct = 100 * mean(abs(finite) <= 5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Section 3.3's never-measured claim: does change acceptance collapse with the
# depth of the node whose rule is redrawn?
changeDepthTable <- function(p) {
  c0 <- p[p$move == "change" & p$noop == 0L, ]
  if (nrow(c0) == 0L) {
    return(NULL)
  }
  by <- split(c0, c0$nodeDepth)
  do.call(
    rbind,
    lapply(names(by), function(depth) {
      d <- by[[depth]]
      data.frame(
        depth = as.integer(depth),
        proposals = nrow(d),
        accepted = sum(d$accepted),
        accept.pct = 100 * mean(d$accepted == 1L),
        stringsAsFactors = FALSE
      )
    })
  )
}

# The window-width input: the acceptance a same-variable cut move would have
# had at each displacement, min(1, exp(logRatio)) averaged over the interior
# nodes a change proposal visited. Signed displacements are pooled by
# magnitude after reporting both signs' counts.
cutTable <- function(d) {
  if (nrow(d) == 0L) {
    return(NULL)
  }
  d$magnitude <- abs(d$displacement)
  by <- split(d, d$magnitude)
  do.call(
    rbind,
    lapply(names(by), function(magnitude) {
      m <- by[[magnitude]]
      alpha <- pmin(1, exp(m$logRatio))
      data.frame(
        displacement = as.integer(magnitude),
        probes = nrow(m),
        left = sum(m$displacement < 0L),
        accept.pct = 100 * mean(alpha),
        vetoed.pct = 100 * mean(!is.finite(m$logRatio)),
        median.logratio = median(m$logRatio[is.finite(m$logRatio)]),
        stringsAsFactors = FALSE
      )
    })
  )
}

roundFrame <- function(x, digits = 3L) {
  numeric <- vapply(x, is.numeric, logical(1L))
  x[numeric] <- lapply(x[numeric], round, digits)
  x
}

summarizeCell <- function(dir, cell) {
  file <- censusFile(dir, cell)
  if (!file.exists(file)) {
    cat("\n== ", cell, ": no census file at ", file, "\n", sep = "")
    return(invisible(NULL))
  }
  census <- readCensus(file)
  p <- census$proposals
  d <- census$cuts
  sampled <- p$sweep >= nBurn
  cat(
    "\n== ",
    cell,
    " (",
    format(nrow(p), big.mark = ","),
    " proposals, ",
    length(unique(p$forest)),
    " forest(s), burn-in ",
    "acceptance ",
    round(100 * mean(p$accepted[!sampled] == 1L), 2),
    "%)\n",
    sep = ""
  )
  p <- p[sampled, ]
  d <- d[d$sweep >= nBurn, ]

  cat("\nper move, sampled sweeps:\n")
  print(roundFrame(moveTable(p), 2L), row.names = FALSE)
  cat("\nlog-likelihood difference among rejected proposals:\n")
  print(roundFrame(rejectionTable(p), 2L), row.names = FALSE)
  cat("\nchange proposals by target node depth:\n")
  print(roundFrame(changeDepthTable(p), 2L), row.names = FALSE)
  cat("\nsame-variable cut move, acceptance against displacement:\n")
  print(roundFrame(cutTable(d), 2L), row.names = FALSE)
  invisible(NULL)
}

# ------------------------------------------------------------------- driver

selected <- if (length(cellArgs) > 0L) {
  intersect(names(cells), cellArgs)
} else {
  names(cells)
}
if (length(selected) == 0L) {
  stop(
    "no known cell selected; cells are: ",
    paste(names(cells), collapse = ", ")
  )
}

# one cell per process; see the header
spawnCell <- function(dir, cell) {
  file <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  if (length(file) != 1L) {
    stop("run mode needs the script's path; invoke it with Rscript")
  }
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(file, "runcell", dir, cell, if (quick) "quick")
  )
  if (status != 0L) {
    stop("cell '", cell, "' failed")
  }
}

if (mode %in% c("run", "both")) {
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  cat(
    "writing census files to ",
    outputDir,
    " (",
    nBurn,
    " burn-in + ",
    nSamples,
    " sampled sweeps, n = ",
    nObservations,
    ")\n",
    sep = ""
  )
  for (cell in selected) {
    spawnCell(outputDir, cell)
  }
}

if (mode == "runcell") {
  if (length(selected) != 1L) {
    stop("runcell takes exactly one cell")
  }
  runCell(outputDir, selected)
}

if (mode %in% c("summarize", "both")) {
  for (cell in selected) {
    summarizeCell(outputDir, cell)
  }
}
