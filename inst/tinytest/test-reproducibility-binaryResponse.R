# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

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
    -0.135254389899304,
    0.305296648320592,
    0.667004345507011,
    0.627483394905247,
    -0.152297510506958
  ),
  varcount = c(17L, 21L, 27L)
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
    0.283145614996885,
    0.00876684849456576,
    0.471386423757719,
    0.484111685661707,
    0.0461271758169887
  ),
  varcount = c(21L, 27L, 26L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceOffset$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceOffset$varcount)
rm(bartFit, n.sims, n.burn)

rm(testData)
