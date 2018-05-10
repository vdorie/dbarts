context("bart w/continuous response")

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("basic Friedman example passes regression test", {
  set.seed(99)
  n.burn <- 100L
  n.sims <- 3000L
  bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
  
  burnRange <- -4L:0L + n.burn
  simRange <- -4L:0L + n.sims
  
  ## values used to be nabbed from BayesTree but since default compilation no longer suffices
  ## we just hope for the best
  expect_equal(bartFit$sigest, 2.75657293556356)
  expect_equal(bartFit$first.sigma[burnRange], c(1.1012953868359, 1.11674494632181, 1.19220453402753, 0.991117100284716, 1.1414230266799))
  expect_equal(bartFit$sigma[simRange], c(0.791065582262117, 0.779960242055518, 0.765895229637772, 0.680435570569646, 0.715885883433391))
  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(5.90394273643642, 17.3027315068948, 17.1364922710491, 4.88896148963325, 18.5629371958413))
  expect_equal(bartFit$yhat.train.mean[1:5], c(7.05955386454589, 17.1465372981429, 16.2879432909552, 3.61515808656625, 19.6911646008683))
  expect_identical(bartFit$yhat.test, NULL)
  expect_identical(bartFit$yhat.test.mean, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(15, 16, 3, 9, 4, 8, 6, 5, 4, 5))
  expect_equal(bartFit$y, testData$y)
})

test_that("weighted Friedman example passes regression test", {
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
  set.seed(99)
  sampler <- dbarts(y ~ x, testData, weights = weights, n.samples = n.sims,
                    control = dbartsControl(n.tree = 50L, n.chains = 1L, n.threads = 1L, updateState = FALSE))
  samples <- sampler$run(n.burn)

  simRange <- -4L:0L + n.sims
  
  expect_equal(samples$sigma[simRange], c(0.710459609625003, 0.766615204051126, 0.77320226463128, 0.795150560866139, 0.91496203135795))
  expect_equal(samples$train[1:5, n.sims], c(7.52723993609673, 15.9289672965162, 17.3166342084901, 4.38749632007886, 18.9820324697806))
  expect_equal(apply(samples$train, 1, mean)[1:5], c(6.94519043557289, 16.9761581257151, 16.5971535791954, 3.55388932832975, 19.4361777198272))
  expect_identical(samples$test, NULL)
  expect_equal(samples$varcount[, n.sims], c(13, 7, 10, 11, 9, 5, 8, 5, 3, 5))
})

test_that("Friedman example with test data passes regression test", {
  n.test <- 25
  set.seed(99)
  testData$x.test <- matrix(runif(n.test * 10), n.test, 10)

  n.burn <- 100L
  n.sims <- 3000L
  set.seed(99)
  bartFit <- bart(testData$x, testData$y, testData$x.test, ndpost = n.sims, ntree = 50L, verbose = FALSE)

  burnRange <- -4L:0L + n.burn
  simRange <- -4L:0L + n.sims
  
  expect_equal(bartFit$sigest, 2.75657293556356)
  expect_equal(bartFit$first.sigma[burnRange], c(1.1012953868359, 1.11674494632181, 1.19220453402753, 0.991117100284716, 1.1414230266799))
  expect_equal(bartFit$sigma[simRange], c(0.791065582262117, 0.779960242055518, 0.765895229637772, 0.680435570569646, 0.715885883433391))
  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(5.90394273643642, 17.3027315068948, 17.1364922710491, 4.88896148963325, 18.5629371958413))
  expect_equal(bartFit$yhat.train.mean[1:5], c(7.05955386454589, 17.1465372981429, 16.2879432909552, 3.61515808656625, 19.6911646008683))
  expect_equal(bartFit$yhat.test[n.sims, 1:5], c(8.90623977366087, 8.54395198404081, 13.4788133890168, 10.5518618865746, 13.2281193102996))
  expect_equal(bartFit$yhat.test.mean[1:5], c(8.47539157876457, 8.75335221762624, 14.5533301222872, 12.9291317904526, 14.4509807756208))
  expect_equal(bartFit$varcount[n.sims,], c(15, 16, 3, 9, 4, 8, 6, 5, 4, 5))
  expect_equal(bartFit$y, testData$y)
})
