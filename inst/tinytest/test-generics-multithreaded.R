source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that predict gives same result when single or multi-threaded
bartFit <- dbarts::bart(
  testData$x, testData$y,
  ndpost = 20L, nskip = 5L, ntree = 5L,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 2L)
expect_equal(predictions, bartFit$yhat.train)

bartFit <- dbarts::bart(
  testData$x, testData$y,
  ndpost = 20L, nskip = 5L, ntree = 5L, nchain = 4L, nthread = 1L,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 2L)
expect_equal(predictions, bartFit$yhat.train)

rm(bartFit)

rm(testData)

