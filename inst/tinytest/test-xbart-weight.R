source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that runs correctly with weighted input
x <- testData$x
y <- testData$y
weights <- rep(1, length(y))

xval <- dbarts::xbart(
  x, y, weights = weights,
  n.samples = 6L, n.burn = c(5L, 3L, 1L),
  n.reps = 3, n.test = 5,
  k = 2,
  n.threads = 2L
)

expect_equal(length(xval), 3L)

rm(weights, y, x, testData)

