source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)
testX <- testData$X
colnames(testX) <- paste0("X.", seq_len(ncol(testX)))

# test that test offset fills in control logicals depending on specification
data <- dbarts::dbartsData(Z ~ X, testData, testX)

expect_null(data@offset)
expect_null(data@offset.test)
expect_identical(data@testUsesRegularOffset, NA)


data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = 0.2)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data)


otherOffset <- 0.2 + 0.1
data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = otherOffset)

expect_equal(data@offset[1:5], rep(0.3, 5))
expect_equal(data@offset.test[1:5], rep(0.3, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data, otherOffset)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = NULL
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_null(data@offset.test)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = 0.1
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.1, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = offset + 0.1
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.3, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = offset
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data)


set.seed(0)
otherOffset <- runif(nrow(testData$X))
data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = otherOffset)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, TRUE)


expect_error(
  dbarts::dbartsData(Z ~ X, testData, testX[-1, ], offset = otherOffset),
  "'offset' was given as 'otherOffset', which cannot be evaluated on the rows of 'test'"
)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = NULL
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_null(data@offset.test)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = 0.2
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = offset + 0.1
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], otherOffset[1:5] + 0.1)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = otherOffset
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = NULL,
  offset.test = otherOffset
)

expect_null(data@offset)
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data, otherOffset)

rm(testData)
