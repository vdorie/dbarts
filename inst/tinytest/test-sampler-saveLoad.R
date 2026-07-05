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

# saving without forcing $state cannot restore the sampler, but forcing the
# serialized promise against the dead pointer must yield the stored-state
# error, not a C-level null-pointer one
set.seed(99L)
bartFit.untouched <- dbarts::bart(
  testData$x, testData$y,
  ntree = 3L, ndpost = 7L, nskip = 0L,
  keeptrees = TRUE, verbose = FALSE
)
tempFile <- tempfile()
saveRDS(bartFit.untouched, file = tempFile)
bartFit.untouched <- readRDS(tempFile)
unlink(tempFile)
expect_error(predict(bartFit.untouched, testData$x),
             pattern = "stored state")

