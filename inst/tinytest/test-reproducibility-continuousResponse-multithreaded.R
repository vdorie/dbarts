# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

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
    1.12824224176549
  ),
  yhatTrain = c(
    0.165454366580266,
    0.165454366580266,
    0.165454366580266,
    0.165454366580266,
    0.165454366580266
  )
)

expect_equal(bartFit$sigest, reference$sigest)
expect_equal(bartFit$sigma, reference$sigma)
expect_equal(bartFit$yhat.train[n.sims, 1L:5L], reference$yhatTrain)

rm(testData, bartFit, weights, n.sims, n.burn, n.tree)
