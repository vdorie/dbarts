source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that formula specification creates valid objects
trainData_df <- as.data.frame(testData)
set.seed(0)
trainData_df$weights <- runif(nrow(trainData_df))
trainData_df$offset  <- rnorm(nrow(trainData_df))

modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

expect_inherits(
  dbarts::dbartsData(modelFormula, trainData_df),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(modelFormula, trainData_df, weights = weights),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(modelFormula, trainData_df, offset = offset),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(modelFormula, trainData_df, weights = weights, offset = offset),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(
    modelFormula,
    trainData_df,
    subset = 1:10,
    weights = weights,
    offset = offset
  ),
  "dbartsData"
)

testData_df <- trainData_df[1:20,]
expect_inherits(
  dbarts::dbartsData(
    modelFormula,
    trainData_df,
    test = testData_df,
    weights = weights
  ),
  "dbartsData"
)

rm(testData_df, trainData_df)


# test that test argument creates valid objects
## test when is embedded in passed data
testData$test <- testData$x[11:20,]
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, test),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, testData$test),
  "dbartsData"
)

## test when is in environment of formula
test <- testData$test
testData$test <- NULL
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, test),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, testData$x[11:20,]),
  "dbartsData"
)
rm(test)


# test that test weights are created correctly
trainData <- as.data.frame(testData)
trainData$weights <- runif(nrow(trainData))

modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

testData_df <- trainData[1:20,]
data <- dbarts::dbartsData(
  modelFormula,
  trainData,
  test = testData_df,
  weights = weights
)
expect_inherits(data, "dbartsData")
expect_equal(data@weights, trainData$weights)
expect_equal(data@weights.test, testData_df$weights)
rm(data, testData_df, modelFormula, trainData)

rm(testData)

