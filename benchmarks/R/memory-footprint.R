#!/usr/bin/env Rscript

# Peak-resident-set audit of a fit against the closed-form model in
# docs/design/memory-footprint.md. One subprocess per grid cell under
# /usr/bin/time (-l on macOS, -v on Linux); its maximum resident set size is
# the measurement, and a baseline subprocess that loads the package and fits
# nothing is subtracted, so what remains is what the fit allocated.
#
# Usage:
#   Rscript memory-footprint.R                  run the grid and report
#   Rscript memory-footprint.R record [out.csv] run and write a baseline CSV
#   Rscript memory-footprint.R fit base.csv     re-report a recorded CSV
# Append 'quick' for a three-cell smoke test of the plumbing; its cells are
# far too small to discriminate the model, so quick mode reports but does not
# score.
#
# Exits non-zero when the residual misses the note's tolerance: every cell
# within max(10 pct, 20 MB) and a median absolute relative residual under
# 5 pct over the grid.
#
# Three terms are measured rather than derived, each by the same subprocess
# method as the cells themselves and each reported in its own column beside
# the closed form, so the gate stays auditable:
#   warmup_mb     the fit path's one-off session growth - byte-compiling the
#                 fit closures and populating the S4 dispatch tables
#   churn_mb      the training-fit mean's collector churn: apply allocates
#                 two small R objects per observation and leaves them
#                 resident until the collector's next cycle
#   ingest_mb     the predictor copies dbartsData makes over and above the
#                 caller's matrix - one persistent subset copy and one
#                 transient complete-cases copy, of which only as much
#                 reaches the peak as the collector has not reclaimed
# All three are host-dependent; none is a byte count the source fixes.
#
# Every cell supplies 'sigest'. Left unset, a gaussian fit estimates the
# starting sigma with an lm() over the whole design, whose model frame, na
# filter, model matrix and QR dominate a short run's peak and are base R's
# allocations, not the sampler's. The sigma-estimate excursion below prices
# them on their own rows.
#
# Maintainer-run on an otherwise-idle machine: a peak RSS measured beside
# other load is not evidence. One subprocess runs at a time.

suppressPackageStartupMessages(library(dbarts))

MB <- 1e6

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
mode <- if (length(args) >= 1L) args[[1L]] else "print"

# Mean live node count of a tree: measured at 5.73 by measureDuplicates()
# below (n = 2e4, T = 75, C = 2, n.burn = 0) and carried here. It enters only
# the live-tree and saved-tree rows, both small at the grid's draw counts.
MEAN.NODES <- 5.73

# Mean number of populated arena depth levels per tree in the designated-leaf
# statistics cache: 12.9 measured bytes per n*T*C over 4 bytes per index_t.
# See the leaf statistics cache row of the design note for why a level count
# rather than a node count is the right multiplier.
CACHE.LEVELS <- 3.2

baseCell <- function() {
  list(
    n = 1e5,
    p = 20L,
    n.trees = 200L,
    n.chains = 1L,
    n.samples = 10L,
    n.test = 0,
    family = "gaussian",
    leaf.columns = NULL,
    keepTrees = FALSE,
    sigest = 1
  )
}

withCell <- function(...) {
  cell <- baseCell()
  changes <- list(...)
  for (field in names(changes)) {
    cell[[field]] <- changes[[field]]
  }
  cell
}

cellName <- function(cell) {
  name <- sprintf(
    "n%.0f-p%d-t%d-c%d-s%d",
    cell$n,
    cell$p,
    cell$n.trees,
    cell$n.chains,
    cell$n.samples
  )
  if (cell$n.test > 0) {
    name <- paste0(name, "-test")
  }
  if (cell$keepTrees) {
    name <- paste0(name, "-keeptrees")
  }
  if (!is.null(cell$leaf.columns)) {
    name <- paste0(name, "-linear")
  }
  if (!identical(cell$family, "gaussian")) {
    name <- paste0(name, "-", cell$family)
  }
  name
}

# A base cell plus one-axis excursions identifies every coefficient of an
# additive separable model; the crossing at the smallest n is the
# interaction check.
buildGrid <- function(quick) {
  if (quick) {
    return(list(
      withCell(n = 1e4),
      withCell(n = 1e4, n.trees = 75L),
      withCell(n = 1e4, n.chains = 2L)
    ))
  }
  cells <- list(baseCell())
  add <- function(...) cells[[length(cells) + 1L]] <<- withCell(...)
  for (n in c(1e4, 1e6)) {
    add(n = n)
  }
  for (p in c(10L, 50L)) {
    add(p = p)
  }
  add(n.trees = 75L)
  for (n.chains in c(2L, 4L)) {
    add(n.chains = n.chains)
  }
  add(keepTrees = TRUE)
  add(n.test = 2e4)
  add(leaf.columns = c("x1", "x2", "x3"))
  add(leaf.columns = c("x1", "x2", "x3"), n.trees = 75L)
  add(family = "probit")
  add(n.samples = 200L, keepTrees = TRUE)
  for (p in c(10L, 50L)) {
    for (n.trees in c(75L, 200L)) {
      for (n.chains in c(1L, 2L)) {
        for (keepTrees in c(FALSE, TRUE)) {
          add(
            n = 1e4,
            p = p,
            n.trees = n.trees,
            n.chains = n.chains,
            keepTrees = keepTrees
          )
        }
      }
    }
  }
  cells[!duplicated(vapply(cells, cellName, character(1L)))]
}

# ---------------------------------------------------------------------------
# The model. Every term of the closed form is a row of the design note's
# table, read off the element type and the allocation site; the three
# measured allowances are returned beside it, never folded into it.
# ---------------------------------------------------------------------------

predictBytes <- function(
  cell,
  warmup = 0,
  churn = function(rows) 0,
  ingest = function(n, p) 0
) {
  n <- cell$n
  p <- cell$p
  n.test <- cell$n.test
  n.trees <- cell$n.trees
  draws <- cell$n.samples * cell$n.chains
  m <- MEAN.NODES
  constant.leaf <- is.null(cell$leaf.columns)
  binary <- identical(cell$family, "probit")

  # once per sampler: packed codes, cut grid, column metadata, the bridge's
  # owned conditioning vectors, and a designated leaf's gathered raw columns
  sampler <- 2 * n * p + 2 * n.test * p + 800 * p + 70 * p + 24 * n + 8 * n.test
  if (!constant.leaf) {
    # the gathered raw columns and the standardized copy the leaf keeps, plus
    # the leaf's sufficient-statistic cache. Leaf memberships partition the
    # observations, so the member lists LIVE in one tree sum to at most 4*n -
    # but the cache is arena-indexed and never pruned, and each slot's member
    # vector keeps the capacity of the largest membership that slot ever held,
    # so the resident total is 4*n per populated arena depth level per tree
    # per chain. The 256 MiB budget bounds the tracked live bytes only, and
    # at every cell here those stay well under it.
    sampler <- sampler +
      16 * n * length(cell$leaf.columns) +
      4 * n * n.trees * cell$n.chains * CACHE.LEVELS
  }

  # per chain: the dominant n*T pair, the total-fit, residual and move
  # scratch slabs, the family's own channels, the live trees, the test slabs
  per.chain <- 4 *
    n *
    n.trees +
    (if (constant.leaf) 4 else 8) * n * n.trees +
    20 * n +
    (if (binary) 16 else 8) * n +
    56 * n.trees * m +
    296 * n.trees +
    (if (constant.leaf) 8 * n.trees * m + 24 * n.trees else 0) +
    16 * n.test

  saved.trees <- if (cell$keepTrees) draws * n.trees * (24 * m + 40) else 0

  # the R vectors the bridge allocates and the engine writes into
  channels <- 8 * draws * (n + n.test) + 4 * p * draws + 8 * draws

  # R-layer transient copies of the prediction arrays, live at packaging.
  # Gaussian reaches three (the engine's array, the reshape, apply's aperm);
  # a binary fit takes no mean, so it reaches three only when the
  # combineChains reshape is the two-step matrix()/t() of a 3-D channel.
  copies <- if (!binary || cell$n.chains > 1L) 3 else 2
  transients <- (copies - 1) * 8 * draws * (n + n.test)

  # the caller's own predictor matrix and response; the copies dbartsData
  # makes on top of them are the measured ingestion allowance
  predictors <- 8 * n * p + 24 * n.test * p + 8 * n

  # the training-fit mean's collector churn, absent on a binary fit
  mean.churn <- if (binary) 0 else churn(n) + churn(n.test)

  closed.form <- sampler +
    cell$n.chains * per.chain +
    saved.trees +
    channels +
    transients +
    predictors
  list(
    closed.form = closed.form,
    warmup = warmup,
    churn = mean.churn,
    ingest = ingest(n, p),
    total = closed.form + warmup + mean.churn + ingest(n, p)
  )
}

# ---------------------------------------------------------------------------
# Measurement.
# ---------------------------------------------------------------------------

WORKER.SOURCE <- '
suppressPackageStartupMessages(library(dbarts))
cell <- readRDS(commandArgs(trailingOnly = TRUE)[[1L]])
if (identical(cell$what, "churn")) {
  m <- matrix(rnorm(cell$rows * cell$draws), cell$draws, cell$rows)
  invisible(apply(m, 2L, mean))
} else if (identical(cell$what, "ingest")) {
  set.seed(99L)
  n <- as.integer(cell$n)
  x <- runif(n * cell$p)
  dim(x) <- c(n, cell$p)
  colnames(x) <- paste0("x", seq_len(cell$p))
  y <- runif(n)
  invisible(list(x, y, dbarts:::dbartsData(x, y)))
} else if (identical(cell$what, "cell")) {
  set.seed(99L)
  n <- as.integer(cell$n)
  # dim<- rather than matrix(), which would double the predictor block for
  # the length of the call and put the peak somewhere the model does not
  # describe
  x <- runif(n * cell$p)
  dim(x) <- c(n, cell$p)
  colnames(x) <- paste0("x", seq_len(cell$p))
  f <- 10 * sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
  y <- if (identical(cell$family, "probit")) {
    rbinom(n, 1L, pnorm(f - mean(f)))
  } else {
    f + rnorm(n)
  }
  rm(f)
  call.args <- list(
    formula = x,
    data = y,
    n.samples = cell$n.samples,
    n.burn = 0L,
    n.trees = cell$n.trees,
    n.chains = cell$n.chains,
    n.threads = 1L,
    keepTrees = cell$keepTrees,
    family = cell$family,
    verbose = FALSE,
    seed = 1L
  )
  if (!is.null(cell$sigest)) call.args$sigest <- cell$sigest
  if (cell$n.test > 0) {
    call.args$test <- x[seq_len(as.integer(cell$n.test)), , drop = FALSE]
  }
  if (!is.null(cell$leaf.columns)) {
    call.args$node.prior <- dbarts::dbartsPriors$linear(cell$leaf.columns)
  }
  invisible(do.call(dbarts::bart, call.args))
}
'

# Parses maximum resident set size out of /usr/bin/time's report: bytes on
# macOS (-l), kilobytes on Linux (-v).
parseMaxRss <- function(lines) {
  hit <- grep("maximum resident set size", lines, ignore.case = TRUE)
  if (length(hit) == 0L) {
    stop("no maximum resident set size in the /usr/bin/time report")
  }
  line <- lines[[hit[[1L]]]]
  value <- as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", trimws(line)))
  if (grepl("kbytes", line, fixed = TRUE)) value * 1024 else value
}

measure <- function(spec, worker.file) {
  spec.file <- tempfile(fileext = ".rds")
  err.file <- tempfile()
  on.exit(unlink(c(spec.file, err.file)), add = TRUE)
  saveRDS(spec, spec.file)
  status <- system2(
    "/usr/bin/time",
    c(
      if (Sys.info()[["sysname"]] == "Darwin") "-l" else "-v",
      shQuote(file.path(R.home("bin"), "Rscript")),
      shQuote(worker.file),
      shQuote(spec.file)
    ),
    stdout = FALSE,
    stderr = err.file
  )
  report <- readLines(err.file, warn = FALSE)
  if (!identical(status, 0L)) {
    stop(paste(c("worker failed:", report), collapse = "\n"))
  }
  parseMaxRss(report)
}

runGrid <- function(quick) {
  worker.file <- tempfile(fileext = ".R")
  on.exit(unlink(worker.file), add = TRUE)
  writeLines(WORKER.SOURCE, worker.file)

  baseline <- measure(list(what = "none"), worker.file)
  warmupCell <- withCell(n = 1e3, p = 5L, n.trees = 1L, n.samples = 1L)
  warmupCell$what <- "cell"
  warmup <- measure(warmupCell, worker.file) - baseline
  cat(sprintf(
    "baseline (package loaded, no fit) %.1f MB; fit-path warm-up %.1f MB\n",
    baseline / MB,
    warmup / MB
  ))

  # the allowances, each measured once per distinct shape and cached
  allowance.cache <- new.env(parent = emptyenv())
  cached <- function(key, compute) {
    if (is.null(allowance.cache[[key]])) {
      allowance.cache[[key]] <- compute()
    }
    allowance.cache[[key]]
  }
  churnFor <- function(draws) {
    function(rows) {
      if (rows == 0) {
        return(0)
      }
      cached(paste("churn", rows, draws), function() {
        peak <- measure(
          list(what = "churn", rows = rows, draws = draws),
          worker.file
        )
        # the probe's own matrix and apply's aperm copy are modelled rows
        max(0, peak - baseline - 16 * rows * draws)
      })
    }
  }
  ingestFor <- function(n, p) {
    cached(paste("ingest", n, p), function() {
      peak <- measure(list(what = "ingest", n = n, p = p), worker.file)
      # the probe's own caller-side matrix and response are modelled rows
      max(0, peak - baseline - 8 * n * p - 8 * n)
    })
  }

  rows <- data.frame()
  addRow <- function(scenario, metric, value) {
    rows <<- rbind(
      rows,
      data.frame(scenario = scenario, metric = metric, value = value)
    )
  }
  addRow("_baseline", "peak_rss_mb", baseline / MB)
  addRow("_warmup", "peak_rss_mb", warmup / MB)

  cat(sprintf(
    "%-34s %9s %9s %7s %7s %7s %9s\n",
    "cell",
    "measured",
    "predicted",
    "closed",
    "warmup",
    "churn",
    "ingest"
  ))
  for (cell in buildGrid(quick)) {
    spec <- cell
    spec$what <- "cell"
    measured <- (measure(spec, worker.file) - baseline) / MB
    draws <- cell$n.samples * cell$n.chains
    parts <- predictBytes(cell, warmup, churnFor(draws), ingestFor)
    predicted <- parts$total / MB
    name <- cellName(cell)
    cat(sprintf(
      "%-34s %9.1f %9.1f %7.1f %7.1f %7.1f %9.1f\n",
      name,
      measured,
      predicted,
      parts$closed.form / MB,
      parts$warmup / MB,
      parts$churn / MB,
      parts$ingest / MB
    ))
    addRow(name, "peak_rss_mb", measured)
    addRow(name, "predicted_mb", predicted)
    addRow(name, "residual_mb", measured - predicted)
    addRow(name, "closed_form_mb", parts$closed.form / MB)
    addRow(name, "churn_mb", parts$churn / MB)
    addRow(name, "ingest_mb", parts$ingest / MB)
  }

  # the starting-sigma estimate, priced against the same cell with sigest
  # supplied: base R's lm() over the whole design, not a sampler allocation
  for (p in if (quick) 20L else c(10L, 20L, 50L)) {
    paired <- withCell(p = p, n = if (quick) 1e4 else 1e5)
    paired$what <- "cell"
    without <- measure(paired, worker.file)
    paired$sigest <- NULL
    with.lm <- measure(paired, worker.file)
    addRow(
      sprintf("sigma-estimate-n%.0f-p%d", paired$n, p),
      "extra_peak_mb",
      (with.lm - without) / MB
    )
  }

  rows$value <- round(rows$value, 3L)
  rows$host <- paste(Sys.info()[c("sysname", "machine")], collapse = " ")
  rows$rev <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
  rows$date <- format(Sys.Date())
  rows$quick <- quick
  rows
}

reportResiduals <- function(rows, score = TRUE) {
  rows <- rows[
    rows$metric %in%
      c("peak_rss_mb", "predicted_mb", "residual_mb") &
      !startsWith(rows$scenario, "_"),
  ]
  wide <- reshape(
    rows[c("scenario", "metric", "value")],
    direction = "wide",
    idvar = "scenario",
    timevar = "metric"
  )
  names(wide) <- sub("^value[.]", "", names(wide))
  wide$tolerance_mb <- pmax(0.10 * wide$predicted_mb, 20)
  wide$rel <- abs(wide$residual_mb) / wide$predicted_mb
  wide$flag <- ifelse(abs(wide$residual_mb) > wide$tolerance_mb, "MISS", "")

  cat("\n")
  print(
    wide[c("scenario", "peak_rss_mb", "predicted_mb", "residual_mb", "flag")],
    row.names = FALSE
  )
  worst <- wide[which.max(abs(wide$residual_mb) / wide$tolerance_mb), ]
  cat(sprintf(
    "\nworst cell: %s, residual %.1f MB against a tolerance of %.1f MB\n",
    worst$scenario,
    worst$residual_mb,
    worst$tolerance_mb
  ))
  cat(sprintf(
    "median absolute relative residual %.1f pct over %d cell(s) (limit 5.0)\n",
    100 * median(wide$rel),
    nrow(wide)
  ))
  if (!score) {
    cat("\nquick mode: a plumbing smoke test, not scored\n")
    return(TRUE)
  }
  misses <- sum(wide$flag == "MISS")
  if (misses > 0L || median(wide$rel) >= 0.05) {
    cat(sprintf("\nFAIL: %d cell(s) outside tolerance\n", misses))
    return(FALSE)
  }
  cat("\nOK: every cell within tolerance\n")
  TRUE
}

# ---------------------------------------------------------------------------
# The R-layer duplicates: gc()'s max-used column around each reshape and the
# apply mean, object.size on what a fit retains, and the mean live node
# count the saved-tree rows need.
# ---------------------------------------------------------------------------

# The high-water of the R heap over expr, above what was already live when it
# started. gc(reset = TRUE) reports the current usage and sets the high-water
# mark to it, so the difference is what expr itself put on the heap.
maxUsedMb <- function(expr) {
  before <- gc(reset = TRUE, full = TRUE)
  force(expr)
  after <- gc(full = TRUE)
  (sum(after[, "max used"] * c(56, 8)) - sum(before[, "used"] * c(56, 8))) / MB
}

measureDuplicates <- function(quick) {
  n <- if (quick) 2000L else 20000L
  n.samples <- if (quick) 20L else 100L
  set.seed(99L)
  x <- runif(n * 10L)
  dim(x) <- c(n, 10L)
  colnames(x) <- paste0("x", seq_len(10L))
  y <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 5 * x[, 5L] + rnorm(n)

  rows <- data.frame()
  addRow <- function(scenario, metric, value) {
    rows <<- rbind(
      rows,
      data.frame(scenario = scenario, metric = metric, value = value)
    )
  }

  full.size <- 8 * n * n.samples * 2L / MB
  # the engine's own layout, n x n.samples x n.chains
  samples <- array(rnorm(n * n.samples * 2L), c(n, n.samples, 2L))
  sizeOf <- function(object) as.numeric(object.size(object)) / MB
  addRow("yhat.train-before-convert", "object_size_mb", sizeOf(samples))
  combined <- NULL
  used <- maxUsedMb(
    combined <- dbarts:::convertSamplesFromDbartsToBart(samples, 2L, TRUE)
  )
  addRow("reshape-combined", "peak_heap_ratio", used / full.size)
  addRow("yhat.train-after-convert", "object_size_mb", sizeOf(combined))
  used <- maxUsedMb(invisible(apply(combined, length(dim(combined)), mean)))
  addRow("apply-mean-combined", "peak_heap_ratio", used / full.size)
  used <- maxUsedMb(
    invisible(dbarts:::convertSamplesFromDbartsToBart(samples, 2L, FALSE))
  )
  addRow("reshape-uncombined", "peak_heap_ratio", used / full.size)
  addRow("prediction-array", "full_size_mb", full.size)
  rm(combined, samples)

  fit <- dbarts::bart(
    x,
    y,
    sigest = 1,
    n.samples = n.samples,
    n.burn = 0L,
    n.trees = 75L,
    n.chains = 2L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 1L
  )
  addRow("x", "object_size_mb", sizeOf(x))
  addRow("y", "object_size_mb", sizeOf(y))
  addRow("yhat.train", "object_size_mb", sizeOf(fit$yhat.train))
  # object.size sees only the external pointer, so the sampler the fit
  # retains is priced from the model instead
  addRow(
    "retained-sampler",
    "engine_total_mb",
    predictBytes(list(
      n = n,
      p = 10L,
      n.trees = 75L,
      n.chains = 2L,
      n.samples = n.samples,
      n.test = 0,
      family = "gaussian",
      leaf.columns = NULL,
      keepTrees = TRUE
    ))$closed.form /
      MB
  )

  sizes <- unlist(lapply(fit$fit$state, function(chain) {
    unlist(lapply(chain$forests, function(forest) forest$tree.sizes))
  }))
  addRow("mean-live-nodes", "nodes_per_tree", mean(sizes))
  # what storeState writes for one draw, against the same run's whole
  # prediction array
  addRow(
    "stored-state-ratio",
    "ratio",
    13 * sum(sizes) / (8 * n * n.samples * 2L)
  )
  rows$value <- round(rows$value, 4L)
  rows$host <- paste(Sys.info()[c("sysname", "machine")], collapse = " ")
  rows$rev <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
  rows$date <- format(Sys.Date())
  rows$quick <- quick
  rows
}

if (mode == "fit") {
  if (length(args) < 2L) {
    stop("usage: memory-footprint.R fit baseline.csv")
  }
  if (!reportResiduals(read.csv(args[[2L]]), !quick)) quit(status = 1L)
} else {
  results <- runGrid(quick)
  duplicates <- measureDuplicates(quick)
  cat("\nR-layer duplicates:\n")
  print(duplicates[c("scenario", "metric", "value")], row.names = FALSE)
  if (mode == "record") {
    out.file <- if (length(args) >= 2L) args[[2L]] else "memory-footprint.csv"
    write.csv(rbind(results, duplicates), out.file, row.names = FALSE)
    cat("\nwrote", nrow(results) + nrow(duplicates), "rows to", out.file, "\n")
  }
  if (!reportResiduals(results, !quick)) quit(status = 1L)
}
