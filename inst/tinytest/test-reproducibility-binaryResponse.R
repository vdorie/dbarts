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
    0.103021425098473,
    0.467048094413091,
    0.514231624032179,
    0.329182361876733,
    -0.343432894333086
  ),
  varcount = c(23L, 23L, 25L)
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
    -0.0621219392664718,
    0.232716454999861,
    0.939205045377481,
    0.348071560452961,
    -0.033666929068414
  ),
  varcount = c(25L, 19L, 28L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceOffset$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceOffset$varcount)
rm(bartFit, n.sims, n.burn)

rm(testData)
