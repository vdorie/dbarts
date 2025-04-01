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
    1.1012953868359, 1.11674494632181, 1.19220453402753,
    0.991117100284716, 1.1414230266799
  )
)
expect_equal(
  bartFit$sigma[simRange],
  c(
    0.791065582262117, 0.779960242055518, 0.765895229637772,
    0.680435570569646, 0.715885883433391
  ),
)
expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    5.90394273643642, 17.3027315068948, 17.1364922710491,
    4.88896148963325, 18.5629371958413
  )
)
expect_equal(
  bartFit$yhat.train.mean[1L:5L],
  c(
    7.05955386454589, 17.1465372981429, 16.2879432909552,
    3.61515808656625, 19.6911646008683
  )
)
expect_null(bartFit$yhat.test)
expect_null(bartFit$yhat.test.mean)
expect_equal(
  bartFit$varcount[n.sims,],
  c(15L, 16L, 3L, 9L, 4L, 8L, 6L, 5L, 4L, 5L)
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
    0.710459609625003, 0.766615204051126, 0.77320226463128,
    0.795150560866139, 0.91496203135795
  )
)
expect_equal(
  samples$train[1L:5L, n.sims],
  c(
    7.52723993609673, 15.9289672965162, 17.3166342084901,
    4.38749632007886, 18.9820324697806
  )
)
expect_equal(
  apply(samples$train, 1L, mean)[1L:5L],
  c(
    6.94519043557289, 16.9761581257151, 16.5971535791954,
    3.55388932832975, 19.4361777198272
  )
)
expect_null(samples$test)
expect_equal(
  samples$varcount[, n.sims],
  c(13L, 7L, 10L, 11L, 9L, 5L, 8L, 5L, 3L, 5L)
)

rm(testData)

