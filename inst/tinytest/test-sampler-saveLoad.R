source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test thanks to Jeremy Coyle
# test that a keepTrees fit saves/loads correctly with no manual state ritual:
# keepTrees now stores the sampler state automatically at fit time
set.seed(99L)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ntree = 3L,
  ndpost = 7L,
  nskip = 0L,
  keeptrees = TRUE,
  verbose = FALSE
)

preds.old <- predict(bartFit, testData$x)

tempFile <- tempfile()
saveRDS(bartFit, file = tempFile)
rm(bartFit)
bartFit <- readRDS(tempFile)
unlink(tempFile)

preds.new <- predict(bartFit, testData$x)
expect_equal(preds.old, preds.new)

# a low-level sampler with updateState = FALSE stores no state, so serializing
# it and forcing the promise against the dead pointer must yield the stored-
# state error naming storeState(), not a C-level null-pointer one
control <- dbarts::dbartsControl(
  updateState = FALSE,
  n.burn = 0L,
  n.chains = 1L,
  n.trees = 3L,
  verbose = FALSE
)
sampler <- dbarts::dbarts(testData$x, testData$y, control = control)
invisible(sampler$run(0L, 7L))

tempFile <- tempfile()
saveRDS(sampler, file = tempFile)
sampler <- readRDS(tempFile)
unlink(tempFile)
expect_error(sampler$predict(testData$x), pattern = "storeState")
