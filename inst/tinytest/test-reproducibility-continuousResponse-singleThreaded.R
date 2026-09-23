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
    1.22035209496754,
    1.24648756689355,
    1.10421701104289,
    1.37899283411952,
    1.28272279671776
  ),
  sigma = c(
    0.762183027792739,
    0.838079808724326,
    0.733699418433162,
    0.827060546096168,
    0.89624589303342
  ),
  yhatTrain = c(
    6.96486740992556,
    16.8139211255427,
    16.2952179267739,
    3.72657275123245,
    19.3773544687538
  ),
  yhatTrainMean = c(
    6.94984991421709,
    16.9631263634629,
    16.3837459684629,
    3.64535500201471,
    19.5771644621008
  ),
  varcount = c(7L, 11L, 6L, 14L, 5L, 7L, 6L, 8L, 6L, 2L)
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
    1.17035584829226,
    1.15883218115004,
    1.22170718852337,
    1.03614059922771,
    1.10927630240841
  ),
  train = c(
    6.5665116592752,
    18.3155589369094,
    17.0943087002348,
    1.77370147913458,
    19.2840333732315
  ),
  trainMean = c(
    7.05489423546351,
    17.2644434863267,
    16.5933476250534,
    3.63585810657127,
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
