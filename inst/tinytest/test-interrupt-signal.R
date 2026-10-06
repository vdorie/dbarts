# A real SIGINT stops a loop of fits wrapped in try(): the interrupt is not an
# error, so try() lets it through to the handler outside the loop, and the
# fits after the stopped one never start. A child process, at home only, on
# unix, with one thread. The signal comes from a shell the child starts, since
# a process started in the background has SIGINT ignored.
if (at_home() && .Platform$OS.type == "unix") {
  dir <- tempfile("interrupt-signal")
  dir.create(dir)
  script <- file.path(dir, "child.R")
  writeLines(
    c(
      "suppressMessages(library(dbarts))",
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
      "      try(sampler$run(5000000L, 1L), silent = TRUE)",
      "    }",
      "    'finished'",
      "  },",
      "  interrupt = function(cond) 'interrupted'",
      ")",
      "writeLines(c(outcome, started), file.path(commandArgs(TRUE)[1L], 'out'))"
    ),
    script
  )
  system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(dir)),
    stdout = FALSE,
    stderr = FALSE,
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = ":")),
    timeout = 120
  )
  expect_true(file.exists(file.path(dir, "out")))
  if (file.exists(file.path(dir, "out"))) {
    expect_equal(readLines(file.path(dir, "out")), c("interrupted", "1"))
  }
  unlink(dir, recursive = TRUE)
}
