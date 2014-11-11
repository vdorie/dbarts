context("multithreaded bart")

source(system.file("common", "multithreadData.R", package = "dbarts"))

test_that("multithreaded matches single threaded", {
  ## something weak so that it runs quickly w/500k observations
  n.sims <- 15L
  n.burn <- 0L
  n.tree <- 3L

  set.seed(99)
  singleThreadedFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = n.tree, verbose = FALSE,
                            nthread = 1L)

  set.seed(99)
  multiThreadedFit <-  bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = n.tree, verbose = FALSE,
                            nthread = 2L)

  expect_equal(singleThreadedFit$sigma, multiThreadedFit$sigma)
  expect_equal(singleThreadedFit$yhat.train[n.sims, 1:5], multiThreadedFit$yhat.train[n.sims, 1:5])
  expect_equal(singleThreadedFit$yhat.train.mean[1:5], multiThreadedFit$yhat.train.mean[1:5])
  expect_identical(multiThreadedFit$yhat.test, NULL)
  expect_identical(multiThreadedFit$yhat.test.mean, NULL)
  expect_equal(singleThreadedFit$varcount[n.sims,], multiThreadedFit$varcount[n.sims,])
})
