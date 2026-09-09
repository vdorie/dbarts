source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that works with fixed seed; the fold splits are drawn via R's
# sample() even with a chain seed set, so the sampling kind is pinned
# against leakage from other test files
oldSampleKind <- RNGkind()[3L]
suppressWarnings(RNGkind(sample.kind = "Rejection"))

x <- testData$x
y <- testData$y

k <- c(4, 8)

runXval <- function(n.threads, ...) {
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 4L,
    n.samples = 20L,
    n.burn = c(10L, 5L),
    n.test = 5,
    k = k,
    n.threads = n.threads,
    seed = 0L,
    ...
  )
}

xval.1 <- runXval(1L)
xval.2 <- runXval(1L)

expect_true(all(!is.na(xval.1)))
expect_equal(dim(xval.1), c(4L, length(k)))
expect_equal(xval.1, xval.2)

# a seed reproduces at ANY thread count: work is distributed over
# (replication, fold) units and each unit draws its split and its fits from
# seeds derived from the call's seed and its own index, so which worker ran a
# unit - and how many there were - reaches no draw. Under R CMD check
# --as-cran more than two simultaneous worker processes are refused, so the
# four-thread arm runs everywhere else.
maxWorkers <- if (
  nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) ||
    parallel::detectCores() < 4L
) {
  2L
} else {
  4L
}

for (n.threads in c(2L, maxWorkers)) {
  xval.threaded <- runXval(n.threads)
  expect_true(all(!is.na(xval.threaded)))
  expect_equal(dim(xval.threaded), c(4L, length(k)))
  expect_identical(
    xval.1,
    xval.threaded,
    info = paste("n.threads =", n.threads)
  )
  # and a rerun at the same count is itself stable
  expect_identical(xval.threaded, runXval(n.threads))
}

# the fold count exceeding the worker count changes nothing either: 10 folds
# over one replication is 10 units, which no thread count divides evenly
runFolds <- function(n.threads) {
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 1L,
    n.samples = 15L,
    n.burn = c(10L, 5L),
    n.test = 10,
    k = k,
    n.trees = 10L,
    n.threads = n.threads,
    seed = 12L
  )
}
folds.1 <- runFolds(1L)
expect_true(all(!is.na(folds.1)))
for (n.threads in c(2L, 3L, maxWorkers)) {
  expect_identical(
    folds.1,
    runFolds(n.threads),
    info = paste("10 folds at n.threads =", n.threads)
  )
}

# an unseeded run is not reproducible; the pin above is the seed's doing and
# not an accident of the grid being deterministic
unseeded <- function() {
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 2L,
    n.samples = 15L,
    n.burn = c(10L, 5L),
    n.test = 5,
    k = k,
    n.threads = 1L
  )
}
expect_true(any(unseeded() != unseeded()))

rm(
  unseeded,
  folds.1,
  runFolds,
  maxWorkers,
  n.threads,
  xval.threaded,
  xval.2,
  xval.1,
  runXval,
  k,
  y,
  x
)

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind)

rm(testData)
