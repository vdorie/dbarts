#!/usr/bin/env Rscript

# The checkpoint after the corrected monotone move
# (docs/plans/monotone-exact-birth-death.md, Staging): monotone-order-size.R's
# fits on the corrected engine under prior = "leaf", each run as a child
# Rscript under caps, reading per-move count time, memory and the Barker
# hybrid's quantities from the move census.
#
# The library must be built with the census; the driver refuses to run when a
# probe fit writes no census file:
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
  # a fit shorter than the driver's poll is never sampled there: this is the
  # resident size at exit, a floor on the peak
  saveRDS(
    list(
      burnSeconds = burn,
      keptSeconds = kept,
      finalRssBytes = rssBytes(Sys.getpid())
    ),
    spec$doneFile
  )
}

# ---- census probe: refuse a library built without it ---------------------------

# The engine opens its census file at the first move, so a tiny monotone fit
# that leaves no file ran against a library built without the census, whose
# report would read as zero moves.
checkCensusBuild <- function(dir) {
  # per driver, so drivers sharing an out= directory do not race
  file <- tempfile("census-probe-", tmpdir = dir, fileext = ".csv")
  code <- paste(
    "suppressPackageStartupMessages(library(dbarts));",
    "set.seed(1); x <- matrix(runif(200), 100, 2,",
    "dimnames = list(NULL, c('x1', 'x2')));",
    "y <- x[, 1] + rnorm(100, 0, 0.1);",
    "s <- dbarts(x, y, control = dbartsControl(n.chains = 1L,",
    "n.threads = 1L, n.trees = 1L, updateState = FALSE),",
    "monotone = dbartsForests$monotone(c(x1 = 1L), prior = 'leaf'));",
    "invisible(s$run(5L, 5L))"
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(code)),
    env = paste0("BARTCORE_MOVE_CENSUS_FILE=", shQuote(file))
  )
  if (!identical(status, 0L) || !file.exists(file) || !file.size(file)) {
    stop(
      "the dbarts on the library path was not built with ",
      "-DBARTCORE_MOVE_CENSUS (a probe fit wrote no census; see this file's ",
      "header): the report would read as zero moves"
    )
  }
  unlink(file)
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
  polls <- 0L
  # a poll gap far past the 1 s interval is the machine sleeping: the wall
  # clock ran on while the child did not, so the time cap counts awake time
  # only, and a fit with any gap is marked (its per-move times may span it)
  lastPoll <- started
  slept <- 0
  repeat {
    Sys.sleep(1)
    now <- Sys.time()
    gap <- as.numeric(difftime(now, lastPoll, units = "secs"))
    if (gap > 30) {
      slept <- slept + gap - 1
    }
    lastPoll <- now
    elapsed <- as.numeric(difftime(now, started, units = "secs")) - slept
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
    polls <- polls + 1L
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
  if (!is.null(done$finalRssBytes) && !is.na(done$finalRssBytes)) {
    peak <- max(peak, done$finalRssBytes)
  }
  list(
    capped = capped,
    started = format(started, "%Y-%m-%d %H:%M:%S"),
    ended = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    # nonzero: the fit spanned a sleep and its census times are suspect; rerun
    sleptSeconds = slept,
    wallSeconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    peakRssBytes = if (polls > 0L || peak > 0) peak else NA_real_,
    # 0: the fit ended within one poll, and the peak is the exit sample
    rssPolls = polls,
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
  all <- all[grepl("^[a-z],[0-9]+,", all)]
  # each run() call numbers its sweeps from 0, so burn-in and kept sweeps
  # share indices: a drop in the index starts a new call
  index <- as.integer(sub("^[a-z],([0-9]+),.*$", "\\1", all))
  call <- cumsum(c(0L, diff(index) < 0L))
  sweep <- index + c(0L, cumsum(tapply(index, call, max) + 1L))[call + 1L]
  isZ <- startsWith(all, "z,") &
    lengths(strsplit(all, ",")) == length(zColumns)
  if (!any(isZ)) {
    return(NULL)
  }
  z <- utils::read.csv(
    text = all[isZ],
    header = FALSE,
    col.names = zColumns,
    na.strings = "NA"
  )
  z$sweep <- sweep[isZ]
  attr(z, "numSweeps") <- length(unique(sweep[startsWith(all, "t,")]))
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
      list(capped = "no result (running or killed with its driver)")
    }
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(fit = basename(dir)),
      as.data.frame(info[setdiff(names(info), "spec")]),
      summary
    )
  }
  columns <- unique(unlist(lapply(rows, names)))
  table <- do.call(
    rbind,
    lapply(rows, function(row) {
      row[setdiff(columns, names(row))] <- NA
      row[columns]
    })
  )
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
checkCensusBuild(normalizePath(out))

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
        "%s: %s%s, %.0f s wall, peak RSS %s, %.4g s per kept sweep\n",
        label,
        if (is.na(result$capped)) "done" else paste("capped:", result$capped),
        if (result$sleptSeconds > 0) {
          sprintf(" (SPANNED A %.0f s SLEEP: rerun)", result$sleptSeconds)
        } else {
          ""
        },
        result$wallSeconds,
        if (is.na(result$peakRssBytes)) {
          "unsampled (< poll)"
        } else if (result$rssPolls == 0L) {
          sprintf(
            ">= %.0f MB (< poll; exit sample)",
            result$peakRssBytes / 2^20
          )
        } else {
          sprintf("%.0f MB", result$peakRssBytes / 2^20)
        },
        result$secondsPerKeptSweep
      ))
    }
  }
}
