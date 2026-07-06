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
    1.27799889304745, 1.16830649427195, 1.04690003561717,
    0.985013989929163, 0.980069275377835
  )
)
expect_equal(
  bartFit$sigma[simRange],
  c(
    0.89751465917337, 0.873520798380078, 1.01218368081942,
    1.0375340741489, 1.12939178334252
  ),
)
expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    6.49024934155389, 17.2619800575445, 16.8192337581922,
    4.48570325819815, 18.491528865972
  )
)
expect_equal(
  bartFit$yhat.train.mean[1L:5L],
  c(
    7.10159282464671, 16.9865298978667, 16.2756168982109,
    3.62427158501144, 19.6012840157722
  )
)
expect_null(bartFit$yhat.test)
expect_null(bartFit$yhat.test.mean)
expect_equal(
  bartFit$varcount[n.sims,],
  c(8L, 13L, 10L, 10L, 7L, 4L, 6L, 7L, 2L, 15L)
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
    1.01951051483556, 0.969417528916804, 0.911990709420778,
    0.891616199703838, 0.958402465287374
  )
)
expect_equal(
  samples$train[1L:5L, n.sims],
  c(
    6.68921014951904, 18.1026333227841, 16.3311173703059,
    3.06424981681408, 20.3337068253742
  )
)
expect_equal(
  apply(samples$train, 1L, mean)[1L:5L],
  c(
    6.82900204406983, 17.4291752920266, 16.4014347835872,
    3.42275934393804, 19.4338261577946
  )
)
expect_null(samples$test)
expect_equal(
  samples$varcount[, n.sims],
  c(12L, 10L, 8L, 12L, 4L, 4L, 4L, 8L, 4L, 5L)
)

rm(testData)

