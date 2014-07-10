source(system.file("common", "friedmanData.R", package = "dbarts"))

context("bart w/continuous response")

test_that("basic Friedman example matches BayesTree", {
  set.seed(99)
  bartFit <- bart(testData$x, testData$y, ndpost = 3000L, ntree = 50L, verbose = FALSE)

  ## values nabbed from BayesTree
  expect_equal(bartFit$sigest, 2.75657293556356)
  expect_equal(bartFit$first.sigma[96:100], c(1.13037519988882, 1.23646726844162, 1.23697074878426, 1.43717531142355, 1.2966006609167))
  expect_equal(bartFit$sigma[2996:3000], c(1.13250597992587, 1.04402216223967, 1.13484065640547, 1.13361395130178, 1.17303068865435))
  expect_equal(bartFit$yhat.train[3000, 1:5], c(6.87144867255257, 17.7282195714198, 15.3745425088033, 2.35166037611414, 17.513954638809))
  expect_equal(bartFit$yhat.train.mean[1:5], c(7.10888053702327, 17.1287345028716, 16.4373304388333, 3.43030743438885, 19.4871611291663))
  expect_identical(bartFit$yhat.test, NULL)
  expect_identical(bartFit$yhat.test.mean, NULL)
  expect_equal(bartFit$varcount[3000,], c(4, 11, 13, 7, 4, 4, 3, 9, 6, 0))
  expect_equal(bartFit$y, testData$y)
})

test_that("weighted Friedman example matches BayesTree (almost)", {
  ## We would run this in BayesTree to get the numbers, but it has
  ## some pecularities with end nodes that end up with less than 5 observations.
  ##
  ## x <- rbind(testData$x, testData$x[91:100,])
  ## y <- c(testData$y, testData$y[91:100])
  ## set.seed(99)
  ## bartFit <- bart(x, y, ndpost = 3000L, ntree = 50L, verbose = FALSE, sigest = 2.96994035586992)

  weights <- c(rep(1, 90), rep(2, 10))
  sampler <- dbarts(y ~ x, testData, weights = weights, n.samples = 3000L, control = dbartsControl(n.tree = 50L))
  set.seed(99)
  samples <- sampler$run()
  
  
  expect_equal(samples$sigma[2996:3000], c(1.54667729319934, 1.68262681627311, 1.81265277638574, 1.76616264026719, 1.66945397458847))
  expect_equal(samples$train[1:5, 3000], c(6.06186554505718, 15.5789410273734, 14.7430324722169, 4.11101178490336, 20.7592870835837))
  expect_equal(apply(samples$train, 1, mean)[1:5], c(7.85410930377738, 16.4362260369858, 15.1945202440969, 4.33967778959696, 20.1430332315336))
  expect_identical(samples$test, NULL)
  expect_equal(samples$varcount[,3000], c(5, 12, 9, 5, 5, 5, 2, 5, 4, 5))
})
