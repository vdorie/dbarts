source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that predict fails if sampler not saved
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20,
  nskip = 5,
  ntree = 5L,
  verbose = FALSE
)
# the pattern matters: bartFit is built with a namespace-qualified call, whose
# call[[1L]] is a `::` call rather than a bare symbol. An unanchored
# expect_error passed even while the guard died with "the condition has length
# > 1" instead of naming 'keeptrees'.
expect_error(
  predict(bartFit, testData$x, n.threads = 1L),
  pattern = "keeptrees"
)
expect_error(extract(bartFit, type = "trees"), pattern = "keeptrees")

bart2Fit <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 20L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = FALSE,
  keepTrainingFits = FALSE
)
expect_error(
  predict(bart2Fit, testData$x, n.threads = 1L),
  pattern = "keepTrees"
)
expect_error(extract(bart2Fit, type = "trees"), pattern = "keepTrees")
expect_error(extract(bart2Fit, sample = "train"), pattern = "keepTrainingFits")

rm(bartFit, bart2Fit)

rm(testData)
