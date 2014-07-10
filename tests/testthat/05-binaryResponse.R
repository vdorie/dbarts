context("bart w/binary response")

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("basic probit example matches previous runs", {
  n.sims <- 1500L

  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, k = 4.5)

  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(0.013067286960227, 0.328688722224162, 0.779692569544106, 0.0328631379576858, -0.142278282618683))
  expect_identical(bartFit$yhat.test, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(25, 11, 19))
}
  
  
