source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that compatibility specification creates valid objects
testData$weights <- runif(length(testData$y))
testData$offset  <- rnorm(length(testData$y))

attach(testData)
expect_inherits(dbarts::dbartsData(x, y), "dbartsData")
expect_inherits(
  dbarts::dbartsData(x, y, weights = weights),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(x, y, offset = offset),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(x, y, weights = weights, offset = offset),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(x, y, subset = 1:10, weights = weights, offset = offset),
  "dbartsData"
)
detach(testData)

testData$weights <- NULL
testData$offset <- NULL


# test that compatibility specification works with dimnames
x <- testData$x
y <- testData$y

colnames(x) <- paste0("x.", seq_len(ncol(x)))
expect_inherits(dbarts::dbartsData(x, y), "dbartsData")

x <- x[,1L,drop=FALSE]
expect_inherits(dbarts::dbartsData(x, y), "dbartsData")
rm(x, y)


# test that compatibility specification works with duplicated dimnames
x <- testData$x
y <- testData$y

colnames(x) <- c(paste0("x.", seq_len(ncol(x) - 2L)), "", "")
x.test <- x
x.test[,ncol(x.test)] <- x.test[,ncol(x.test)] + 1
 
expect_equal(dbarts::dbartsData(x, y, x.test)@x.test, x.test)
rm(x.test, y, x)

rm(testData)
