#!/usr/bin/env Rscript

# P1's structural acceptance readout: per-move acceptance for one seed of
# each (rung, arm) of P1-friedman.R, beside the acceptance rates Pratola
# (2016) section 2.2 published for the same cell under birth/death only,
# around 18% at sigma^2 = 1 and around 4% at sigma^2 = 0.1.
#
# Acceptance is not a quantity an ordinary build reports, so this script needs
# the instrumented library that benchmarks/R/move-census.R documents, built
# with -DBARTCORE_MOVE_CENSUS appended to CPPFLAGS through a user Makevars:
#
#   printf 'CPPFLAGS += -DBARTCORE_MOVE_CENSUS\n' > /tmp/census-makevars
#   mkdir -p /tmp/census-lib
#   R_MAKEVARS_USER=/tmp/census-makevars \
#     R CMD INSTALL --preclean --library=/tmp/census-lib .
#   R_LIBS=/tmp/census-lib Rscript benchmarks/R/surfaces/P1-friedman-census.R
#
# It is kept apart from P1-friedman.R so that the cell itself runs against an
# ordinary build. The chains are shorter than the cell's: acceptance is a
# per-proposal rate read over sampled sweeps, and the instrumented build
# prices a schedule of cut displacements on every change proposal, which
# costs far more than the sweep it instruments.
#
# One cell runs per R process, as the move census does: the engine opens the
# census file once, on the first record, so a second cell in the same process
# would append to the first cell's file.
#
# Usage:
#   Rscript P1-friedman-census.R [outputDir] [cell ...]
#   Rscript P1-friedman-census.R run [outputDir] [cell ...]
#   Rscript P1-friedman-census.R summarize [outputDir] [cell ...]
# Append 'quick' for a smoke test at reduced sweeps and n.

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
modes <- c("run", "summarize", "runcell")
mode <- if (length(args) >= 1L && args[[1L]] %in% modes) args[[1L]] else "both"
args <- setdiff(args, modes)

nTrees <- 200L
p <- 10L
replicate <- 1L
nSamples <- if (quick) 50L else 200L

arms <- list(
  default = NULL,
  birthdeath = c(birth_death = 1, swap = 0, change = 0, birth = 0.5),
  swap = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
)

# rung, sigma and frequency are P1-friedman.R's, so a census cell fits the
# same data one of its replicates fits; n and the sweep counts are not.
designs <- list(
  house = list(
    n = if (quick) 500L else 2000L,
    nBurn = if (quick) 100L else 1000L,
    sigma = 0.25,
    frequency = 1,
    seedKey = "housesigma025"
  ),
  variance1 = list(
    n = if (quick) 500L else 5000L,
    nBurn = if (quick) 100L else 5000L,
    sigma = 1,
    frequency = 2,
    seedKey = "pratolavariance1"
  ),
  variance01 = list(
    n = if (quick) 500L else 5000L,
    nBurn = if (quick) 100L else 5000L,
    sigma = sqrt(0.1),
    frequency = 2,
    seedKey = "pratolavariance01"
  )
)

cells <- character(0)
for (designName in names(designs)) {
  cells <- c(cells, paste(designName, names(arms), sep = "-"))
}

published <- c(variance1 = 18, variance01 = 4)

censusFile <- function(dir, cell) {
  file.path(dir, sprintf("P1-move-census-%s.csv", cell))
}

runCell <- function(dir, cell) {
  parts <- strsplit(cell, "-", fixed = TRUE)[[1L]]
  design <- designs[[parts[1L]]]
  arm <- arms[[parts[2L]]]
  file <- censusFile(dir, cell)
  unlink(file)
  Sys.setenv(BARTCORE_MOVE_CENSUS_FILE = file)
  on.exit(Sys.unsetenv("BARTCORE_MOVE_CENSUS_FILE"), add = TRUE)

  set.seed(surfacesDataSeed("P1", design$seedKey, replicate))
  data <- surfacesFriedman(
    design$n,
    nTest = 1000L,
    p,
    sigma = design$sigma,
    frequency = design$frequency
  )
  # The sampler's own run, not bart2's: the engine stamps each record with the
  # sweep index of the call that produced it, so burn-in and sampling have to
  # share one call for `sweep >= nBurn` to separate them.
  call <- list(
    data$x,
    data$y,
    control = dbartsControl(
      verbose = FALSE,
      n.trees = nTrees,
      n.chains = 1L,
      n.threads = 1L,
      n.samples = nSamples,
      updateState = FALSE,
      seed = surfacesSamplerSeed(replicate)
    )
  )
  if (!is.null(arm)) {
    call$proposal.probs <- arm
  }
  sampler <- do.call(dbarts, call)
  elapsed <- system.time(
    invisible(sampler$run(design$nBurn, nSamples))
  )[["elapsed"]]
  rm(sampler)
  invisible(gc(FALSE))
  if (!file.exists(file)) {
    stop(
      "no census records for cell '",
      cell,
      "': the loaded dbarts was not built with -DBARTCORE_MOVE_CENSUS ",
      "(see this file's header)"
    )
  }
  cat(sprintf(
    "%-22s %6.1f s  %s (%.1f MB)\n",
    cell,
    elapsed,
    file,
    file.size(file) / 1024^2
  ))
}

# Only the proposal records; the cut-displacement records the instrumented
# build also writes belong to the move census, not to this readout.
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

readProposals <- function(file) {
  lines <- readLines(file)
  read.csv(
    text = lines[substr(lines, 1L, 1L) == "p"],
    header = FALSE,
    col.names = proposalNames
  )
}

# Acceptance per proposal MADE, which is the denominator Pratola's figures
# use, and per proposal that reached a score, which separates a move refused
# on its merits from one that never had a candidate.
moveTable <- function(p) {
  rows <- lapply(c(sort(unique(p$move)), "all"), function(move) {
    m <- if (move == "all") p else p[p$move == move, ]
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
  numeric <- vapply(out, is.numeric, logical(1L))
  out[numeric] <- lapply(out[numeric], round, 2L)
  out
}

summarizeCell <- function(dir, cell) {
  file <- censusFile(dir, cell)
  if (!file.exists(file)) {
    cat("\n== ", cell, ": no census file at ", file, "\n", sep = "")
    return(invisible(NULL))
  }
  proposals <- readProposals(file)
  designName <- strsplit(cell, "-", fixed = TRUE)[[1L]][1L]
  sampled <- proposals$sweep >= designs[[designName]]$nBurn
  reference <- if (designName %in% names(published)) {
    sprintf(", published birth/death acceptance %g%%", published[[designName]])
  } else {
    ""
  }
  cat(sprintf(
    "\n== %s (%s proposals, %d burn-in + %d sampled sweeps%s)\n",
    cell,
    format(nrow(proposals), big.mark = ","),
    designs[[designName]]$nBurn,
    nSamples,
    reference
  ))
  print(moveTable(proposals[sampled, ]), row.names = FALSE)
  invisible(NULL)
}

selected <- intersect(cells, args)
if (length(selected) == 0L) {
  selected <- cells
}
outputDir <- surfacesOutputDir(args, flags = cells)

spawnCell <- function(dir, cell) {
  file <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(file, "runcell", dir, cell, if (quick) "quick")
  )
  if (status != 0L) {
    stop("cell '", cell, "' failed")
  }
}

if (mode %in% c("run", "both")) {
  cat(sprintf(
    "writing census files to %s (%d sampled sweeps after each rung's burn-in)\n",
    outputDir,
    nSamples
  ))
  surfacesUptime("uptime before")
  for (cell in selected) {
    spawnCell(outputDir, cell)
  }
  surfacesUptime("uptime after")
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
