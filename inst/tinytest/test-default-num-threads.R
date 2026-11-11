# defaultNumThreads caps the default thread count at two under R CMD check's
# variable, read as R's own check reads it: set, non-empty, not "false". The
# variable is set and restored in process.
original <- Sys.getenv("_R_CHECK_LIMIT_CORES_", NA_character_)
originalMc <- options(mc.cores = 8L)
restore <- function() {
  if (is.na(original)) {
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  } else {
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = original)
  }
  options(originalMc)
}
cores <- dbarts::guessNumCores()
tasks <- function(n) {
  if (is.na(cores)) NA_integer_ else as.integer(min(cores, n))
}
capped <- function(n) {
  if (is.na(cores)) NA_integer_ else as.integer(min(cores, n, 2))
}
tryCatch(
  {
    for (value in c(NA, "", "false", "FALSE")) {
      if (is.na(value)) {
        Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
      } else {
        Sys.setenv("_R_CHECK_LIMIT_CORES_" = value)
      }
      expect_identical(dbarts::defaultNumThreads(100L), tasks(100L))
      expect_identical(dbarts::defaultNumThreads(), tasks(Inf))
    }
    expect_identical(dbarts::defaultNumThreads(1L), tasks(1L))
    for (value in c("TRUE", "warn", "true", "WARN")) {
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = value)
      expect_identical(dbarts::defaultNumThreads(100L), capped(100L))
      expect_identical(dbarts::defaultNumThreads(), capped(Inf))
      expect_identical(dbarts::defaultNumThreads(1L), capped(1L))
    }
    # the defaults call it; a stated count is never capped, mc.cores is no
    # input to it
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
    expect_identical(
      dbarts::dbartsControl(n.chains = 4L)@n.threads,
      if (is.na(cores)) 1L else capped(4L)
    )
    expect_identical(
      dbarts::dbartsControl(n.chains = 4L, n.threads = 4L)@n.threads,
      4L
    )
    expect_identical(dbarts::defaultNumThreads(4L), capped(4L))
    # no task count: the door of the caller refuses it, not this default
    expect_identical(dbarts::defaultNumThreads(0), capped(Inf))
  },
  finally = restore()
)
expect_identical(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA_character_), original)
