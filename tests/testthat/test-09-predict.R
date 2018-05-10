context("predict")

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("predict fails if sampler not saved", {
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE)
  expect_error(predict(bartFit, testData$x))
})

test_that("predict gives same result as x_train", {
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x)
  expect_equal(predictions, bartFit$yhat.train)
  
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, nchain = 4L, nthread = 1L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x)
  expect_equal(predictions, bartFit$yhat.train)
})

test_that("fixed sample mode when run sequentially gives same predictions as sequential updates mode", {
  set.seed(0)
  pred.bart <- bart2(testData$x, testData$y, testData$x, n.samples = 5, n.burn = 0L, n.trees = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE)$yhat.test
  
  set.seed(0)
  sampler <- dbarts(testData$x, testData$y,
                    control = dbartsControl(n.samples = 5, n.burn = 0L, n.trees = 4L, n.chains = 1L, n.threads = 1L, keepTrees = TRUE, updateState = FALSE))
  sampler$sampleTreesFromPrior()
  for (i in seq_len(5L))
    invisible(sampler$run(0L, 1L))
  pred.dbarts <- sampler$predict(testData$x)
  
  expect_equal(pred.bart, t(pred.dbarts))
})

test_that("sequentially running samples don't overflow with fixed trees", {
  sampler <- dbarts(testData$x, testData$y,
                    control = dbartsControl(n.samples = 5, n.burn = 0L, n.trees = 4L, n.chains = 1L, n.threads = 1L, keepTrees = TRUE, updateState = FALSE))
  sampler$sampleTreesFromPrior()
  for (i in seq_len(6L))
    invisible(sampler$run(0L, 1L))
  
  expect_is(sampler, "dbartsSampler")
})
