source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that formula specification raises errors
expect_error(
  dbarts::dbartsData("not-a-formula", testData),
  "unrecognized 'formula' type"
)
expect_error(
  dbarts::dbartsData(y ~ x),
  "object 'y' not found"
)
expect_error(
  dbarts::dbartsData(NULL, testData),
  "unrecognized 'formula' type"
)
expect_error(
  dbarts::dbartsData(y ~ 0, testData),
  "predictors must be specified"
)


# test that extra arguments for formula specification raises errors
set.seed(0)
testData_df <- as.data.frame(testData)
testData_df$weights <- runif(nrow(testData_df))
testData_df$offset  <- rnorm(nrow(testData_df))

modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

expect_error(
  dbarts::dbartsData(modelFormula, testData_df, subset = "not-a-number"),
  "empty 'subset' specified"
)
expect_error(
  dbarts::dbartsData(modelFormula, testData_df, weights = "not-a-number"),
  "variable lengths differ \\(found for '\\(weights\\)'\\)"
)
expect_error(
  dbarts::dbartsData(
    modelFormula,
    testData_df,
    weights = rep("not-a-number", nrow(testData_df))
  ),
  "'weights' must be of type numeric"
)
expect_error(
  dbarts::dbartsData(modelFormula, testData_df, weights = offset),
  " 'weights' must all be non-negative"
)

rm(modelFormula, testData_df)


# test that compatibility specification raises errors
expect_error(
  dbarts::dbartsData(testData$x, "not-a-number"),
  "'data' must be numeric as well"
)
expect_error(
  dbarts::dbartsData(testData$x, testData$y[1]),
  "'x' must have the same number of observations as 'y'"
)


# test that test argument raises errors
expect_error(
  dbarts::dbartsData(y ~ x, testData, testData$x[11:20, 1:9]),
  "number of columns in 'test' must be equal to that of 'x'"
)
expect_error(
  dbarts::dbartsData(y ~ x, testData, "not-a-matrix"),
  "test matrix must be numeric"
)
expect_error(
  dbarts::dbartsData(y ~ x, testData, outOfScope),
  "object 'outOfScope' not found"
)

test <- testData$x[11:20,]
colnames(test) <- paste0("x", c(1:9, 11))
expect_warning(
  dbarts::dbartsData(y ~ x, testData, test),
  "column names of 'test' does not equal that of 'x'"
)
rm(test)

rm(testData)
