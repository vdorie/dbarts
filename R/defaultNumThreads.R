## the thread count the fitting functions default to: the cores on the
## machine, no more than the tasks to run, and no more than two where R CMD
## check asks for that. The variable is read as R's own check reads it
## (parallel:::.check_ncores): set, non-empty and not "false" in any case. NA, where the
## cores cannot be counted, passes through for fallBackToOneThread.
defaultNumThreads <- function(n.tasks = Inf) {
  ## a default argument is evaluated before its door has checked n.chains, so
  ## a value that is no task count is left to that door's own refusal
  if (
    !is.numeric(n.tasks) ||
      length(n.tasks) != 1L ||
      is.na(n.tasks) ||
      n.tasks < 1
  ) {
    n.tasks <- Inf
  }
  cores <- guessNumCores()
  if (is.na(cores)) {
    return(NA_integer_)
  }
  limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(limit) && limit != "false") {
    n.tasks <- min(n.tasks, 2)
  }
  as.integer(min(cores, n.tasks))
}
