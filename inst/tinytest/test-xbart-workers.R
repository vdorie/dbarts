source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

run <- function(...) {
  dbarts::xbart(
    x,
    y,
    n.samples = 6L,
    n.burn = c(5L, 3L),
    n.reps = 2L,
    n.trees = 5L,
    k = c(1, 2),
    seed = 12L,
    ...
  )
}

windows <- .Platform$OS.type == "windows"
ref <- run(n.threads = 1L)

if (!windows) {
  expect_equal(run(n.threads = 2L, parallel = "fork"), ref)
}
expect_equal(run(n.threads = 2L, parallel = "socket"), ref)

cl <- parallel::makeCluster(2L)
expect_equal(run(n.threads = 4L, cl = cl), ref)
expect_equal(run(n.threads = 2L, cl = cl), ref)
expect_equal(run(n.threads = 1L, cl = cl), ref)
expect_true(
  identical(unlist(parallel::clusterEvalQ(cl, 1L)), c(1L, 1L))
)
parallel::stopCluster(cl)

if (windows) {
  expect_error(run(n.threads = 2L, parallel = "fork"), "socket")
}

expect_error(run(parallel = "bogus"), "arg")
expect_error(run(cl = 1L), "'cl' must be")
expect_error(run(cl = structure(list(), class = "cluster")), "'cl' must be")

if (!windows) {
  boom <- function(y.test, y.test.hat, weights) stop("loss went boom")
  expect_error(
    run(n.threads = 2L, parallel = "fork", loss = boom),
    "loss went boom"
  )
}

old <- options(dbarts.parallel = "bogus")
expect_error(run(n.threads = 2L), "arg")
options(old)

verboseLine <- function(...) {
  out <- capture.output(run(n.threads = 2L, verbose = TRUE, ...))
  grep("running", out, value = TRUE)
}
old <- options(dbarts.parallel = "socket")
expect_true(grepl("socket worker", verboseLine()))
options(old)
if (!windows) {
  expect_true(grepl("fork worker", verboseLine(parallel = "fork")))
}
oldEnv <- Sys.getenv("RSTUDIO", NA)
Sys.setenv(RSTUDIO = "1")
expect_true(grepl("socket worker", verboseLine(parallel = "auto")))
if (is.na(oldEnv)) {
  Sys.unsetenv("RSTUDIO")
} else {
  Sys.setenv(RSTUDIO = oldEnv)
}

warner <- function(y.test, y.test.hat, weights) {
  warning("loss warned")
  0
}
oldWarn <- options(warn = 2)
expect_error(run(n.threads = 1L, loss = warner), "loss warned")
expect_error(
  run(n.threads = 2L, parallel = "socket", loss = warner),
  "loss warned"
)
if (!windows) {
  expect_error(
    run(n.threads = 2L, parallel = "fork", loss = warner),
    "loss warned"
  )
}
options(oldWarn)

if (!windows) {
  # exactly one of the two children kills itself without raising an error
  lock <- tempfile()
  dier <- function(y.test, y.test.hat, weights) {
    if (dir.create(lock, showWarnings = FALSE)) {
      tools::pskill(Sys.getpid(), tools::SIGKILL)
      Sys.sleep(10)
    }
    0
  }
  expect_error(
    suppressWarnings(run(n.threads = 2L, parallel = "fork", loss = dier)),
    "exited without a result"
  )
  unlink(lock, recursive = TRUE)
}
