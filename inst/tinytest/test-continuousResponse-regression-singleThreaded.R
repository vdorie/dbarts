source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test_that basic Friedman example passes regression test
set.seed(99L)
n.burn <- 100L
n.sims <- 3000L
bartFit <- dbarts::bart(
  testData$x, testData$y,
  ndpost = n.sims, nskip = n.burn, ntree = 50L,
  verbose = FALSE
)

burnRange <- -4L:0L + n.burn
simRange <- -4L:0L + n.sims

## values used to be nabbed from BayesTree but since default compilation no longer suffices
## we just hope for the best
expect_equal(bartFit$sigest, 2.75657293556356)
expect_equal(
  bartFit$first.sigma[burnRange],
  c(
    0.978613806496775, 0.966751172964945, 0.971960034686007,
    1.13080177483346, 1.16526455713829
  )
)
expect_equal(
  bartFit$sigma[simRange],
  c(
    0.820050134071397, 0.6995317490121, 0.657370394544752,
    0.741948236947241, 0.790859875668413
  ),
)
expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    6.9004600366795, 19.428366367298, 16.7515843829575,
    4.03664814000156, 19.2115096206922
  )
)
expect_equal(
  bartFit$yhat.train.mean[1L:5L],
  c(
    6.94987391847224, 17.1189458912545, 16.4845933421212,
    3.59659979897874, 19.402121358429
  )
)
expect_null(bartFit$yhat.test)
expect_null(bartFit$yhat.test.mean)
expect_equal(
  bartFit$varcount[n.sims,],
  c(11L, 8L, 8L, 10L, 5L, 6L, 5L, 7L, 5L, 4L)
)
expect_equal(bartFit$y, testData$y)



# test that weighted Friedman example passes regression test
## We would run this in BayesTree to get the numbers, but it has
## some pecularities with end nodes that end up with less than 5 observations.
##
## x <- rbind(testData$x, testData$x[91:100,])
## y <- c(testData$y, testData$y[91:100])
## set.seed(99)
## bartFit <- bart(x, y, ndpost = 3000L, ntree = 50L, verbose = FALSE, sigest = 2.96994035586992)
n.burn <- 100L
n.sims <- 3000L

weights <- c(rep(1, 90), rep(2, 10))
set.seed(99L)
sampler <- dbarts::dbarts(
  y ~ x, testData, weights = weights, n.samples = n.sims,
  control = dbarts::dbartsControl(
    n.tree = 50L, n.chains = 1L, n.threads = 1L, updateState = FALSE
  )
)
samples <- sampler$run(n.burn)

simRange <- -4L:0L + n.sims

expect_equal(
  samples$sigma[simRange],
  c(
    0.799590220807799, 0.798519164113922, 0.795786934536967,
    0.819845412348088, 0.803202410155125
  )
)
expect_equal(
  samples$train[1L:5L, n.sims],
  c(
    6.41247200504924, 17.779541101313, 16.6384682149653,
    2.99182248714801, 20.0367716413078
  )
)
expect_equal(
  apply(samples$train, 1L, mean)[1L:5L],
  c(
    6.95393593883588, 16.8725634668488, 16.2790155289519,
    3.56431296797806, 19.505179739922
  )
)
expect_null(samples$test)
expect_equal(
  samples$varcount[, n.sims],
  c(13L, 10L, 10L, 10L, 6L, 6L, 5L, 2L, 3L, 7L)
)

rm(testData)

