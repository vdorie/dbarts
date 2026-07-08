# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

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
    0.76218302779274,
    0.838079808724327,
    0.73369941843316,
    0.827060546096167,
    0.89624589303342
  ),
  yhatTrain = c(
    6.96486740992556,
    16.8139211255427,
    16.2952179267739,
    3.72657275123249,
    19.3773544687538
  ),
  yhatTrainMean = c(
    6.94984991421711,
    16.9631263634629,
    16.3837459684629,
    3.64535500201474,
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
    0.841804052178214,
    0.939948827039949,
    0.90045532341112,
    0.842498335510937,
    0.778535847501784
  ),
  train = c(
    7.32492731937251,
    18.2449254335483,
    16.9077260101996,
    3.67917481598543,
    18.8310103164836
  ),
  trainMean = c(
    6.84567942087165,
    16.9728459619374,
    16.4521278192749,
    3.59229995799396,
    19.5219252778269
  ),
  varcount = c(12L, 11L, 10L, 9L, 8L, 7L, 3L, 12L, 4L, 5L)
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
