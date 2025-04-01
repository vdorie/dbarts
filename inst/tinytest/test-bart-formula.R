source(system.file("common", "almostLinearBinaryData.R", package = "dbarts"), local = TRUE)

# test that bart creates viable sampler with formula, data specification"
data <- data.frame(y = testData$y, x = testData$x)
modelFormula <- y ~ x.1 + x.2 + x.3

expect_inherits(
  dbarts::bart(modelFormula, data, nskip = 0L, ndpost = 1L, verbose = FALSE),
  "bart"
)
expect_inherits(
  dbarts::bart(
    modelFormula,
    data[1L:100L,],
    data[101L:200L,],
    nskip = 0L,
    ndpost = 1L,
    verbose = FALSE
  ),
  "bart"
)
rm(modelFormula, data)

rm(testData)

