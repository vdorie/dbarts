#!/usr/bin/env Rscript

# The checkpoint after the corrected monotone move
# (docs/plans/monotone-exact-birth-death.md, Staging): monotone-order-size.R's
# fits on the corrected engine under prior = "leaf", each run as a child
# Rscript under caps, reading per-move count time, memory and the Barker
# hybrid's quantities from the move census.
#
# The library must be built with the census:
#
#   printf 'CPPFLAGS += -DBARTCORE_MOVE_CENSUS\n' > /tmp/census-makevars
#   R_MAKEVARS_USER=/tmp/census-makevars \
#     R CMD INSTALL --preclean --library=/tmp/census-lib .
#   R_LIBS=/tmp/census-lib Rscript benchmarks/R/monotone-checkpoint.R \
#     run out=/tmp/checkpoint
#
# Caps: macOS has no timeout and refuses ulimit -v, so each fit is a child
# Rscript started with system2(wait = FALSE) that writes its pid first; the
# driver polls `ps -o rss=` every second and kills the child with
# tools::pskill past 60 minutes of wall time or 8 GB resident. A killed fit is
# recorded as capped, with the cap it hit.
#
# Usage:
#   run out=<dir> [trees=1,5,10,20] [designs=1:1,2:1,1:2,1:3] [seeds=1]
#       [n=5000] [burn=1000] [kept=200] [minutes=60] [gb=8] [count.all]
#     designs are constrained:free axis counts (p the sum); count.all sets
#     BARTCORE_MOVE_CENSUS_COUNT_ALL, so moves the free bound decides are
#     counted too (the decision is unchanged; their time is not a cost the
#     sampler pays, and the report keeps it apart)
#   report out=<dir>       summarize the census files already in <dir>
#   child <file.rds>       one fit (the driver's own entry)
#
# Per fit: the largest count time and memory among counts the sampler needed,
# the share of sweeps with such a count over 1 s, accepted deaths leaving a
# merged component past 2^20 down-sets, and, over counted moves, the share
# the hybrid would switch at B = 2^22 and a_B / a_MH on those. A move the
# bound decided with m <= 17 cannot pass B (W <= 2^m m), so it counts as
# unswitched; the report gives how many moves stay unknown.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else "run"
option <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", name, "="), "", hit[1L]) else default
}
script <- local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  normalizePath(file[1L])
})

# ---- child: one fit ------------------------------------------------------------

runChild <- function(specFile) {
  spec <- readRDS(specFile)
  writeLines(as.character(Sys.getpid()), spec$pidFile)
  suppressPackageStartupMessages(library(dbarts))
  p <- spec$nc + spec$nf
  set.seed(spec$seed)
  x <- matrix(
    runif(spec$n * p),
    spec$n,
    p,
    dimnames = list(NULL, paste0("x", seq_len(p)))
  )
  f <- 2 * rowSums(x[, seq_len(spec$nc), drop = FALSE])
  if (spec$nf > 0L) {
    f <- f + 2 * sin(2 * pi * x[, spec$nc + 1L]) * (1 + x[, 1L])
  }
  y <- f + rnorm(spec$n, 0, 0.3)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = spec$trees,
    n.samples = spec$kept,
    updateState = FALSE,
    seed = spec$seed + 1L
  )
  directions <- setNames(rep(1L, spec$nc), colnames(x)[seq_len(spec$nc)])
  sampler <- dbarts(
    x,
    y,
    control = control,
    monotone = dbartsForests$monotone(directions, prior = "leaf")
  )
  burn <- system.time(invisible(sampler$run(spec$burn, 0L)))[["elapsed"]]
  kept <- system.time(invisible(sampler$run(0L, spec$kept)))[["elapsed"]]
  saveRDS(
    list(burnSeconds = burn, keptSeconds = kept),
    spec$doneFile
  )
}

# ---- driver: fits under caps ---------------------------------------------------

rssBytes <- function(pid) {
  out <- suppressWarnings(system2(
    "ps",
    c("-o", "rss=", "-p", pid),
    stdout = TRUE,
    stderr = FALSE
  ))
  out <- trimws(out)
  if (!length(out) || !nzchar(out[1L])) NA_real_ else 1024 * as.numeric(out[1L])
}

runFitCapped <- function(spec, capSeconds, capBytes, countAll) {
  unlink(c(spec$pidFile, spec$doneFile, spec$censusFile))
  specFile <- file.path(spec$dir, "spec.rds")
  saveRDS(spec, specFile)
  env <- paste0("BARTCORE_MOVE_CENSUS_FILE=", shQuote(spec$censusFile))
  if (countAll) {
    env <- c(env, "BARTCORE_MOVE_CENSUS_COUNT_ALL=1")
  }
  started <- Sys.time()
  system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), "child", shQuote(specFile)),
    env = env,
    wait = FALSE,
    stdout = file.path(spec$dir, "child.log"),
    stderr = file.path(spec$dir, "child.log")
  )
  pid <- NA_integer_
  capped <- NA_character_
  peak <- 0
  repeat {
    Sys.sleep(1)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (is.na(pid) && file.exists(spec$pidFile)) {
      pid <- as.integer(readLines(spec$pidFile, warn = FALSE)[1L])
    }
    if (is.na(pid)) {
      if (elapsed > 120) {
        capped <- "no pid"
        break
      }
      next
    }
    rss <- rssBytes(pid)
    if (is.na(rss)) {
      break # the child has exited
    }
    peak <- max(peak, rss)
    if (elapsed > capSeconds) {
      capped <- "time"
    } else if (rss > capBytes) {
      capped <- "memory"
    }
    if (!is.na(capped)) {
      tools::pskill(pid, tools::SIGKILL)
      break
    }
  }
  done <- if (file.exists(spec$doneFile)) readRDS(spec$doneFile) else NULL
  if (is.null(done) && is.na(capped)) {
    capped <- "failed"
  }
  list(
    capped = capped,
    wallSeconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    peakRssBytes = peak,
    secondsPerKeptSweep = if (is.null(done)) {
      NA_real_
    } else {
      done$keptSeconds / spec$kept
    }
  )
}

# ---- report: the census, per fit -----------------------------------------------

zColumns <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "move",
  "accepted",
  "m",
  "pairComponents",
  "needed",
  "counted",
  "logNormalizer",
  "seconds",
  "downSets",
  "work",
  "peakBytes",
  "switched",
  "barkerShare",
  "mergedDownSets"
)

readCensus <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }
  # a killed child can leave a partial last line
  all <- readLines(file, warn = FALSE)
  lines <- grep("^z,", all, value = TRUE)
  lines <- lines[lengths(strsplit(lines, ",")) == length(zColumns)]
  if (!length(lines)) {
    return(NULL)
  }
  z <- utils::read.csv(
    text = lines,
    header = FALSE,
    col.names = zColumns,
    na.strings = "NA"
  )
  sweeps <- unique(as.integer(sub(
    "^[a-z],([0-9-]+),.*$",
    "\\1",
    grep("^t,[0-9]+,", all, value = TRUE)
  )))
  attr(z, "numSweeps") <- length(sweeps)
  z
}

summarizeFit <- function(z) {
  if (is.null(z)) {
    return(data.frame(moves = 0L))
  }
  needed <- z[z$needed == 1L & z$counted == 1L, ]
  slowSweeps <- unique(needed$sweep[needed$seconds > 1])
  deaths <- z[z$move == "death" & z$accepted == 1L, ]
  counted <- z[z$counted == 1L, ]
  known <- z$counted == 1L | z$m <= 17
  switched <- counted[counted$switched == 1L, ]
  data.frame(
    moves = nrow(z),
    neededShare = mean(z$needed == 1L),
    maxCountSeconds = if (nrow(needed)) max(needed$seconds) else 0,
    maxCountBytes = if (nrow(needed)) max(needed$peakBytes) else 0,
    maxDownSets = if (nrow(needed)) max(needed$downSets) else 0,
    slowSweepShare = length(slowSweeps) / max(1L, attr(z, "numSweeps")),
    acceptedDeaths = nrow(deaths),
    mergedPast2e20 = sum(deaths$mergedDownSets > 2^20, na.rm = TRUE),
    switchedShare = sum(counted$switched == 1L) / max(1L, sum(known)),
    switchUnknown = sum(!known),
    barkerShareSwitched = if (nrow(switched)) {
      mean(switched$barkerShare, na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

# ---- main ----------------------------------------------------------------------

if (mode == "child") {
  runChild(args[2L])
  quit(status = 0L)
}

out <- option("out", NA_character_)
if (is.na(out)) {
  stop("give out=<dir>")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (mode == "report") {
  rows <- list()
  for (dir in list.dirs(out, recursive = FALSE)) {
    summary <- summarizeFit(readCensus(file.path(dir, "census.csv")))
    info <- if (file.exists(file.path(dir, "result.rds"))) {
      readRDS(file.path(dir, "result.rds"))
    } else {
      list()
    }
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(fit = basename(dir)),
      as.data.frame(info[setdiff(names(info), "spec")]),
      summary
    )
  }
  table <- do.call(rbind, rows)
  print(table, digits = 3L)
  utils::write.csv(table, file.path(out, "summary.csv"), row.names = FALSE)
  quit(status = 0L)
}

trees <- as.integer(strsplit(option("trees", "1,5,10,20"), ",")[[1L]])
designs <- strsplit(
  strsplit(option("designs", "1:1,2:1,1:2,1:3"), ",")[[1L]],
  ":"
)
seeds <- as.integer(strsplit(option("seeds", "1"), ",")[[1L]])
n <- as.integer(option("n", "5000"))
burn <- as.integer(option("burn", "1000"))
kept <- as.integer(option("kept", "200"))
capSeconds <- 60 * as.numeric(option("minutes", "60"))
capBytes <- 2^30 * as.numeric(option("gb", "8"))
countAll <- "count.all" %in% args

for (design in designs) {
  for (nTrees in trees) {
    for (seed in seeds) {
      nc <- as.integer(design[1L])
      nf <- as.integer(design[2L])
      label <- sprintf("t%02d-c%d-f%d-s%d", nTrees, nc, nf, seed)
      dir <- file.path(normalizePath(out), label)
      dir.create(dir, showWarnings = FALSE)
      spec <- list(
        trees = nTrees,
        nc = nc,
        nf = nf,
        seed = seed,
        n = n,
        burn = burn,
        kept = kept,
        dir = dir,
        pidFile = file.path(dir, "pid"),
        doneFile = file.path(dir, "done.rds"),
        censusFile = file.path(dir, "census.csv")
      )
      result <- runFitCapped(spec, capSeconds, capBytes, countAll)
      saveRDS(result, file.path(dir, "result.rds"))
      cat(sprintf(
        "%s: %s, %.0f s wall, peak RSS %.0f MB, %.4g s per kept sweep\n",
        label,
        if (is.na(result$capped)) "done" else paste("capped:", result$capped),
        result$wallSeconds,
        result$peakRssBytes / 2^20,
        result$secondsPerKeptSweep
      ))
    }
  }
}
