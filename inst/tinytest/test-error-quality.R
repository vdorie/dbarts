source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# (a) a non-finite (Inf/-Inf) response is named, alongside the pre-existing
# missing-value message that NaN still receives
yInf <- testData$y
yInf[1L] <- Inf
expect_error(
  dbarts::dbartsData(testData$x, yInf),
  "response contains non-finite values"
)
yNegInf <- testData$y
yNegInf[1L] <- -Inf
expect_error(
  dbarts::dbartsData(testData$x, yNegInf),
  "response contains non-finite values"
)
yNaN <- testData$y
yNaN[1L] <- NaN
expect_error(
  dbarts::dbartsData(testData$x, yNaN),
  "response contains missing values"
)

# (b) zero-row training data is named rather than faulting deeper
expect_error(
  dbarts::dbartsData(testData$x[0L, , drop = FALSE], testData$y[0L]),
  "data has zero rows"
)

# (c) a zero (or thinned-to-zero) sample count in bart2 is named
expect_error(
  dbarts::bart2(
    testData$x,
    testData$y,
    n.samples = 0L,
    n.burn = 0L,
    verbose = FALSE
  ),
  "'n.samples' must be a positive integer"
)

# (d) a newdata variable missing from a formula fit is named, even when it
# collides with a base object (here 'c') that the model-frame eval would
# otherwise silently resolve
df <- data.frame(y = testData$y, a = testData$x[, 1L], c = testData$x[, 2L])
newdf <- data.frame(a = testData$x[1:5, 1L])
expect_error(
  dbarts::dbartsData(y ~ a + c, df, test = newdf),
  "missing variable required by the model: 'c'"
)

# (e) an unnamed-matrix test set against a named-predictor fit still matches
# by position (a warning, not an error, so legitimately-unnamed use is not
# broken), but the warning now spells out the assumed column mapping
xNamed <- testData$x[, 1:3]
colnames(xNamed) <- c("aa", "bb", "cc")
testUnnamed <- unname(testData$x[1:5, 1:3])
expect_warning(
  dbarts::dbartsData(xNamed, testData$y, test = testUnnamed),
  "column 1 = 'aa'"
)

rm(testData)
