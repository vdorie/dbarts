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
    1.27799889304745,
    1.16830649427195,
    1.04690003561717,
    0.985013989929164,
    0.980069275377834
  ),
  sigma = c(
    0.897514659173368,
    0.873520798380078,
    1.01218368081941,
    1.0375340741489,
    1.12939178334252
  ),
  yhatTrain = c(
    6.49024934155385,
    17.2619800575445,
    16.8192337581922,
    4.48570325819818,
    18.491528865972
  ),
  yhatTrainMean = c(
    7.10159282464669,
    16.9865298978667,
    16.2756168982109,
    3.62427158501144,
    19.6012840157722
  ),
  varcount = c(8L, 13L, 10L, 10L, 7L, 4L, 6L, 7L, 2L, 15L)
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
    1.01951051483557,
    0.969417528916808,
    0.911990709420783,
    0.891616199703843,
    0.958402465287379
  ),
  train = c(
    6.68921014951903,
    18.102633322784,
    16.3311173703059,
    3.06424981681408,
    20.3337068253742
  ),
  trainMean = c(
    6.82900204406982,
    17.4291752920266,
    16.4014347835872,
    3.42275934393803,
    19.4338261577946
  ),
  varcount = c(12L, 10L, 8L, 12L, 4L, 4L, 4L, 8L, 4L, 5L)
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
