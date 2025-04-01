source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that predict fails if sampler not saved
bartFit <- dbarts::bart(
  testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE
)
expect_error(predict(bartFit, testData$x, n.threads = 1L))

rm(bartFit)

rm(testData)

