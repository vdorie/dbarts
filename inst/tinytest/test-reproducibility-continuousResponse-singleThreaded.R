# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

# The pinned values only mean anything on the reference build: its draw path
# is scalar and fixed-order, where the shipped build's is free to reassociate
# and lands elsewhere from the same seed. Regenerate on the reference build
# too - tools/regenerate-snapshots.R evaluates this file top to bottom, and
# the exit below is not a function it defines.
if (!identical(dbarts:::buildInfo()$mode, "reference")) {
  exit_file("seeded-drift snapshots are pinned to the reference build")
}

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# basic Friedman example
n.burn <- 100L
n.sims <- 3000L

set.seed(99L)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = 50L,
  verbose = FALSE
)

burnRange <- -4L:0L + n.burn
simRange <- -4L:0L + n.sims

referenceBasic <- list(
  sigest = 2.75657293556356,
  firstSigma = c(
    1.12293994498356,
    1.15396446900544,
    1.0294546509834,
    1.09067737916805,
    0.935362507635552
  ),
  sigma = c(
    0.759129374624587,
    0.838878898485429,
    0.701691862434377,
    0.721724940811735,
    0.779734901422799
  ),
  yhatTrain = c(
    7.06413329929853,
    17.4337419205153,
    16.7082841802654,
    4.8703219286024,
    20.0268128639406
  ),
  yhatTrainMean = c(
    6.9797170241275,
    17.1201525933919,
    16.248609827515,
    3.58302464724473,
    19.431144046293
  ),
  varcount = c(8L, 8L, 9L, 8L, 5L, 5L, 4L, 5L, 9L, 4L)
)

expect_equal(bartFit$sigest, referenceBasic$sigest)
expect_equal(bartFit$first.sigma[burnRange], referenceBasic$firstSigma)
expect_equal(bartFit$sigma[simRange], referenceBasic$sigma)
expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceBasic$yhatTrain)
expect_equal(bartFit$yhat.train.mean[1L:5L], referenceBasic$yhatTrainMean)
expect_null(bartFit$yhat.test)
expect_null(bartFit$yhat.test.mean)
expect_equal(bartFit$varcount[n.sims, ], referenceBasic$varcount)
expect_equal(bartFit$y, testData$y)

rm(bartFit, n.burn, n.sims, burnRange, simRange)

# weighted Friedman example: 10 observations duplicated with weight 2
n.burn <- 100L
n.sims <- 3000L
weights <- c(rep(1, 90), rep(2, 10))

set.seed(99L)
sampler <- dbarts::dbarts(
  y ~ x,
  testData,
  weights = weights,
  n.samples = n.sims,
  control = dbarts::dbartsControl(
    n.tree = 50L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
samples <- sampler$run(n.burn)

simRange <- -4L:0L + n.sims

referenceWeighted <- list(
  sigma = c(
    1.17035584829227,
    1.15883218115005,
    1.22170718852339,
    1.03614059922771,
    1.10927630240842
  ),
  train = c(
    6.56651165927518,
    18.3155589369094,
    17.0943087002347,
    1.77370147913461,
    19.2840333732315
  ),
  trainMean = c(
    7.05489423546351,
    17.2644434863267,
    16.5933476250534,
    3.6358581065713,
    19.466296865403
  ),
  varcount = c(6L, 10L, 8L, 12L, 11L, 4L, 6L, 10L, 1L, 5L)
)

expect_equal(samples$sigma[simRange], referenceWeighted$sigma)
expect_equal(samples$train[1L:5L, n.sims], referenceWeighted$train)
expect_equal(
  apply(samples$train, 1L, mean)[1L:5L],
  referenceWeighted$trainMean
)
expect_null(samples$test)
expect_equal(samples$varcount[, n.sims], referenceWeighted$varcount)

rm(sampler, samples, n.burn, n.sims, simRange, weights, testData)
