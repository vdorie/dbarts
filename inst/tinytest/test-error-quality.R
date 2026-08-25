source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
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

# (c) a zero (or thinned-to-zero) sample count is named, by one shared
# message across bart2/xbart/rbart_vi (dbarts()/dbartsControl() accept it -
# a host-loop-driven sampler - so the message states the split too)
zeroSampleMessage <- paste0(
  "'n.samples' must leave at least one draw after thinning \\(n.samples ",
  "%/% n.thin = 0\\)"
)
expect_error(
  dbarts::bart2(
    testData$x,
    testData$y,
    n.samples = 0L,
    n.burn = 0L,
    verbose = FALSE
  ),
  zeroSampleMessage
)
expect_error(
  dbarts::xbart(testData$x, testData$y, n.samples = 0L, n.burn = c(0L, 0L)),
  zeroSampleMessage
)
groupBy <- factor(rep(seq_len(4L), length.out = length(testData$y)))
expect_error(
  dbarts::rbart_vi(
    testData$x,
    testData$y,
    group.by = groupBy,
    n.samples = 0L,
    n.burn = 0L,
    verbose = FALSE
  ),
  zeroSampleMessage
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
warnings.unnamedTest <- captureWarnings(
  dbarts::dbartsData(xNamed, testData$y, test = testUnnamed)
)
expect_equal(length(warnings.unnamedTest), 1L)
expect_match(conditionMessage(warnings.unnamedTest[[1L]]), "column 1 = 'aa'")

# (f) an offset naming the kind actually supplied: a matrix offset is only
# ever a multinomial category shift (paired with 'counts'), and a flat
# offset is refused on a fit that carries 'counts'
matOff <- matrix(0.0, nrow(testData$x), 3L)
expect_error(
  dbarts::dbartsData(testData$x, offset = matOff),
  "which only family = \"multinomial\" accepts"
)
oneHotY <- matrix(0L, nrow(testData$x), 3L)
oneHotY[cbind(seq_len(nrow(testData$x)), rep_len(1:3, nrow(testData$x)))] <- 1L
expect_error(
  dbarts::dbartsData(testData$x, counts = oneHotY, offset = testData$y),
  "softmax's own null direction"
)

# (g) xbart names a non-scalar n.threads by the argument, not R's raw
# coercion error, and resolves 'family' before the response is ingested
expect_error(
  dbarts::xbart(testData$x, testData$y, n.threads = c(1L, 2L)),
  "'n.threads' must be of length 1"
)
badSurvY <- cbind(
  time = abs(testData$y),
  status = rep(c(0L, 1L), length.out = length(testData$y))
)
expect_error(
  dbarts::xbart(testData$x, badSurvY, family = "zzz"),
  "'arg' should be one of"
)

rm(testData)
