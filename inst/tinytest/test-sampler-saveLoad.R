source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test thanks to Jeremy Coyle
# test that sampler saves/loads correctly
set.seed(99L)
bartFit <- dbarts::bart(
  testData$x, testData$y,
  ntree = 3L, ndpost = 7L, nskip = 0L,
  keeptrees = TRUE, verbose = FALSE
)

preds.old <- predict(bartFit, testData$x)

invisible(bartFit$fit$state)

tempFile <- tempfile()
saveRDS(bartFit, file = tempFile)
rm(bartFit)
bartFit <- readRDS(tempFile)
unlink(tempFile)

preds.new <- predict(bartFit, testData$x)
expect_equal(preds.old, preds.new)

