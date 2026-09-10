# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

# The pinned values only mean anything on the reference build: its draw path
# is scalar and fixed-order, where the shipped build's is free to reassociate
# and lands elsewhere from the same seed. This exit guards test runs only;
# tools/regenerate-snapshots.R evaluates this file outside tinytest, where
# exit_file() returns its message and stops nothing, so that tool carries its
# own refusal.
if (!identical(dbarts:::buildInfo()$mode, "reference")) {
  exit_file("seeded-drift snapshots are pinned to the reference build")
}

source(
  system.file("common", "multithreadData.R", package = "dbarts"),
  local = TRUE
)

n.sims <- 5L
n.burn <- 0L
n.tree <- 3L
weights <- c(
  rep(1, floor(.9 * nrow(testData$x))),
  rep(2, nrow(testData$x) - floor(.9 * nrow(testData$x)))
)

set.seed(99L)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  weights = weights,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = n.tree,
  verbose = FALSE,
  nthread = 2L
)

reference <- list(
  sigest = 1.08105693951868,
  sigma = c(
    1.2372451390987,
    1.14277060463474,
    1.14262103884284,
    1.12976218154258,
    1.1262386359894
  ),
  yhatTrain = c(
    0.167526877778454,
    0.167526877778454,
    0.167526877778455,
    0.167526877778454,
    0.167526877778454
  )
)

expect_equal(bartFit$sigest, reference$sigest)
expect_equal(bartFit$sigma, reference$sigma)
expect_equal(bartFit$yhat.train[n.sims, 1L:5L], reference$yhatTrain)

rm(testData, bartFit, weights, n.sims, n.burn, n.tree)
