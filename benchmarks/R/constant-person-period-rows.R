#!/usr/bin/env Rscript

# Measures the person-period row cap, the max.rows argument of hazard()
# (R/family.R), which the expander expandDiscreteTimeHazard (R/dbarts.R)
# enforces before it allocates: a discrete-time hazard fit whose grid would
# turn n subjects into more than 1e7 at-risk rows is refused, naming the
# coarsening levers. The cap is a guard against an over-fine time grid, not a
# statement about what the sampler can fit.
#
# Reported per expanded row count N': the EXPANSION time, and the peak
# resident set the expansion reaches, since the expanded design is a full
# n' x (p + 1) double matrix and the guard exists to stop exactly that
# allocation. Each row count runs in a FRESH R process, so the resident-set
# figure is the whole cost of holding one expansion rather than a delta
# against an already-grown heap. The cap BINDS if a design a user would
# plausibly write - the row counts here are 10 to 20 columns wide - reaches
# 1e7 rows at a memory cost the machine could actually absorb.
#
# Usage: Rscript benchmarks/R/constant-person-period-rows.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cap is revisited. Run it on a quiet machine.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
childArg <- grep("^child=", args, value = TRUE)

p <- 10L
targets <- if (quick) c(1e5, 1e6) else c(1e6, 1e7, 3e7)

residentMB <- function() {
  as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
}

measure <- function(target) {
  # a subject observed to period t_i contributes t_i rows; drawing t_i
  # uniformly on 1..K gives a mean of (K + 1) / 2 rows per subject
  periods <- 20L
  n <- as.integer(ceiling(target / ((periods + 1) / 2)))
  set.seed(37)
  x <- matrix(runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  time <- sample.int(periods, n, replace = TRUE)
  status <- rbinom(n, 1L, 0.4)

  invisible(gc(FALSE))
  before <- residentMB()
  start <- Sys.time()
  expanded <- dbarts:::expandDiscreteTimeHazard(
    x,
    time,
    status,
    max.rows = 1e9
  )
  elapsed <- as.numeric(Sys.time() - start, units = "secs")
  after <- residentMB()
  rows <- nrow(expanded$x)
  c(subjects = n, rows = rows, secs = elapsed, mb = after - before)
}

if (length(childArg) > 0L) {
  values <- measure(as.numeric(sub("^child=", "", childArg[[1L]])))
  cat("CHILD", paste0(values, collapse = " "), "\n")
} else {
  scriptPath <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
  )
  cat(sprintf("%d predictor columns, 20 periods, cap = 1e7 rows\n", p))
  cat(sprintf(
    "%12s %12s %12s %12s %12s\n",
    "target",
    "subjects",
    "rows",
    "seconds",
    "peak_MB"
  ))
  for (target in targets) {
    output <- system2(
      file.path(R.home("bin"), "Rscript"),
      c(scriptPath, paste0("child=", format(target, scientific = FALSE))),
      stdout = TRUE,
      stderr = FALSE
    )
    values <- as.numeric(strsplit(
      grep("^CHILD ", output, value = TRUE)[[1L]],
      " +"
    )[[1L]][-1L])
    cat(sprintf(
      "%12.1e %12.0f %12.0f %12.3f %12.1f\n",
      target,
      values[[1L]],
      values[[2L]],
      values[[3L]],
      values[[4L]]
    ))
  }
  # what the guard itself does at the boundary
  set.seed(41)
  small <- matrix(runif(100L * p), 100L, p)
  smallTime <- sample.int(20L, 100L, replace = TRUE)
  refusal <- tryCatch(
    {
      dbarts:::expandDiscreteTimeHazard(
        small,
        smallTime,
        rbinom(100L, 1L, 0.5),
        max.rows = 100
      )
      "no refusal"
    },
    error = function(e) conditionMessage(e)
  )
  cat(sprintf(
    "\nguard at max.rows = 100 on a %d-row expansion:\n  %s\n",
    sum(smallTime),
    refusal
  ))
}
