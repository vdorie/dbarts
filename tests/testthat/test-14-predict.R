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

