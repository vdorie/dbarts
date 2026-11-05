# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

# The pinned values are exact, so they are checked only where exactness is
# the contract: the reference build (--enable-reference-build), whose draw
# path stays scalar and fixed-order should the shipped build ever vectorize it
# (today both builds compile the same kernels). They were recorded on arm64
# macOS; across hosts draws are never bitwise, but on x86-64 Linux they land
# within 2e-13 of these values, well inside the tolerance. CI installs that
# build and runs these files on both; locally, install it the same way and run
# them with tinytest::run_test_file(). This exit guards test runs only;
# tools/regenerate-snapshots.R evaluates this file outside tinytest, where
# exit_file() returns its message and stops nothing, so that tool carries its
# own refusal.
if (!identical(dbarts:::buildInfo()$mode, "reference")) {
  exit_file("seeded-drift snapshots are pinned to the reference build")
}

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# basic probit example
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  k = 4.5,
  verbose = FALSE
)

referenceBase <- list(
  yhatTrain = c(
    0.0576368530059772,
    0.184665389804701,
    0.208272585469654,
    0.200604543271416,
    0.0420271289864169
  ),
  varcount = c(15L, 20L, 21L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceBase$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceBase$varcount)
expect_equal(extract(bartFit), pnorm(bartFit$yhat.train))
rm(bartFit, n.sims, n.burn)

# basic probit example with a non-zero binary offset
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  k = 4.5,
  binaryOffset = 0.1,
  verbose = FALSE
)

n.sims <- nrow(bartFit$yhat.train)

referenceOffset <- list(
  yhatTrain = c(
    -0.147850430756404,
    0.553946087227406,
    0.296439194431698,
    0.29244544138844,
    -0.344340170266188
  ),
  varcount = c(18L, 22L, 24L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceOffset$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceOffset$varcount)
rm(bartFit, n.sims, n.burn)

rm(testData)
