# A real SIGINT stops a loop of fits wrapped in try(): the interrupt is not an
# error, so try() lets it through to the handler outside the loop, and the
# fits after the stopped one never start. A child process, at home only, on
# unix, with one thread, loading the dbarts under test from its own library.
# The signal comes from a shell the child starts, since a process started in
# the background has SIGINT ignored. What this shows that the in-process tests
# cannot is R's own handling of a signal the poll consumes; the unhandled path
# (restarts, the options hooks, the blank line) is pinned in test-monotone.R.
if (at_home() && .Platform$OS.type == "unix") {
  dir <- tempfile("interrupt-signal")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  script <- file.path(dir, "child.R")
  lib <- dirname(system.file(package = "dbarts"))
  writeLines(
    c(
      sprintf("suppressMessages(library(dbarts, lib.loc = %s))", deparse(lib)),
      "set.seed(1L)",
      "x <- matrix(rnorm(100L), 50L, 2L)",
      "y <- rnorm(50L)",
      "fit <- function() {",
      "  dbarts(y ~ x, control = dbartsControl(n.threads = 1L, n.chains = 1L,",
      "    n.trees = 10L, n.samples = 1L, updateState = FALSE))",
      "}",
      "started <- 0L",
      "outcome <- tryCatch(",
      "  {",
      "    for (i in 1:3) {",
      "      sampler <- fit()",
      "      started <- started + 1L",
      "      if (i == 1L) {",
      "        system2('sh', c('-c', shQuote(paste('sleep 1; kill -INT',",
      "          Sys.getpid()))), wait = FALSE)",
      "      }",
      "      try(sampler$run(1000000L, 1L), silent = TRUE)",
      "    }",
      "    'finished'",
      "  },",
      "  interrupt = function(cond) 'interrupted'",
      ")",
      "writeLines(c(outcome, started), file.path(commandArgs(TRUE)[1L], 'out'))"
    ),
    script
  )
  # timeout kills the child when it does not end by itself
  system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(dir)),
    stdout = FALSE,
    stderr = FALSE,
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = ":")),
    timeout = 60
  )
  expect_true(file.exists(file.path(dir, "out")))
  if (file.exists(file.path(dir, "out"))) {
    expect_equal(readLines(file.path(dir, "out")), c("interrupted", "1"))
  }
}
