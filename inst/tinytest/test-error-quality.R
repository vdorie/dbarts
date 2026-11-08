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
# NaN is a missing response, so the default na.action drops its row; under
# na.pass the completeness check still names it
expect_equal(
  length(dbarts::dbartsData(testData$x, yNaN)@y),
  length(yNaN) - 1L
)
expect_error(
  dbarts::dbartsData(testData$x, yNaN, na.action = na.pass),
  "response contains missing values"
)

# (b) zero-row training data is named rather than faulting deeper
expect_error(
  dbarts::dbartsData(testData$x[0L, , drop = FALSE], testData$y[0L]),
  "data has zero rows"
)

# (c) a zero (or thinned-to-zero) sample count is named, by one shared
# message across bart/xbart (dbarts()/dbartsControl() accept it -
# a host-loop-driven sampler - so the message states the split too)
zeroSampleMessage <- paste0(
  "'n.samples' must leave at least one draw after thinning \\(n.samples ",
  "%/% n.thin = 0\\)"
)
expect_error(
  dbarts::bart(
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
  "'family' should be one of"
)

# a quoted prior for k (0.9-34's bart2 took quote(chi)) is refused with the
# forms that do work named, at every door that reads k; other wrong values
# keep their own message
quotedChi <- quote(chi)
quotedChiCall <- quote(chi(1.25, Inf))
quotedFits <- list(
  bart = function(k) {
    dbarts::bart(
      testData$x,
      testData$y,
      k = k,
      n.samples = 3L,
      n.burn = 3L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE
    )
  },
  bartBT = function(k) {
    dbarts::bartBT(
      testData$x,
      testData$y,
      k = k,
      nskip = 3L,
      ndpost = 3L,
      nchain = 1L,
      nthread = 1L,
      verbose = FALSE
    )
  },
  xbart = function(k) {
    dbarts::xbart(testData$x, testData$y, k = k, n.reps = 2L, n.threads = 1L)
  }
)
for (door in names(quotedFits)) {
  for (quoted in list(quotedChi, quotedChiCall)) {
    quotedMsg <- tryCatch(quotedFits[[door]](quoted), error = conditionMessage)
    expect_true(grepl("quoted expression", quotedMsg, fixed = TRUE))
    expect_true(grepl("chi(1.25, Inf) directly", quotedMsg, fixed = TRUE))
    expect_true(grepl("dbartsPriors$chi(", quotedMsg, fixed = TRUE))
    expect_true(grepl("the string \"chi(", quotedMsg, fixed = TRUE))
  }
  # each form the message names fits
  for (works in list(
    dbarts::dbartsPriors$chi(1.25, Inf),
    "chi(1.25, Inf)"
  )) {
    expect_silent(quotedFits[[door]](works))
  }
}
expect_error(
  quotedFits$bart(-1),
  "'k' must be a positive scalar or a hyperprior specification",
  fixed = TRUE
)
rm(quotedChi, quotedChiCall, quotedFits, door, quoted, quotedMsg, works)

# (g) argument mistakes that reached base R's own errors, or named an argument
# the caller did not write, are refused in the package's words
quickArgs <- list(
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
quickBart <- function(...) {
  do.call(dbarts::bart, c(list(testData$x, testData$y, ...), quickArgs))
}
for (flag in c("keepCall", "keepSampler", "samplerOnly")) {
  values <- if (flag == "samplerOnly") {
    list(NA, "yes")
  } else {
    list(NA, "yes", NULL)
  }
  for (value in values) {
    expect_error(
      do.call(quickBart, structure(list(value), names = flag)),
      paste0("'", flag, "' must be TRUE or FALSE"),
      fixed = TRUE,
      info = flag
    )
  }
}
# a number is read as as.logical reads it, as 0.9-34's bart2 read one: 1 is
# TRUE and 0 is FALSE
expect_null(quickBart(keepCall = 0)$call)
expect_false(is.null(quickBart(keepCall = 1)$call))
expect_inherits(quickBart(samplerOnly = 1), "dbartsSampler")
expect_inherits(quickBart(samplerOnly = 0), "bart")
expect_inherits(quickBart(keepSampler = 1)$fit, "dbartsSampler")
expect_null(quickBart(keepSampler = 0)$fit)
sigestWords <- paste0(
  "'sigest' must be a positive finite number, or \"auto\", \"dense\" ",
  "or \"sparse\""
)
expect_error(quickBart(sigest = Inf), sigestWords, fixed = TRUE)
expect_error(
  dbarts::dbarts(testData$x, testData$y, sigest = Inf),
  sigestWords,
  fixed = TRUE
)
quickBT <- function(...) {
  dbarts::bartBT(
    testData$x,
    testData$y,
    ...,
    ntree = 3L,
    nskip = 2L,
    verbose = FALSE
  )
}
expect_error(
  quickBT(binaryOffset = "a", ndpost = 2L),
  "'binaryOffset' must be numeric",
  fixed = TRUE
)
expect_error(
  quickBT(binaryOffset = c(1, 2), ndpost = 2L),
  "'binaryOffset' must be one value or one per row of 'y.train'",
  fixed = TRUE
)
expect_error(
  quickBT(x.test = testData$x[1:5, 1:2], ndpost = 2L),
  "'x.test' must have as many columns as 'x.train'",
  fixed = TRUE
)
expect_error(
  dbarts::bartBT(testData$x, testData$y, verbose = NA),
  "'verbose' must be TRUE or FALSE",
  fixed = TRUE
)
expect_error(
  quickBT(ndpost = 2L, keepevery = 5L),
  "'ndpost' must leave at least one draw after thinning by 'keepevery'",
  fixed = TRUE
)
expect_error(
  quickBT(ndpost = 2L, seed = c(1, 2)),
  "'seed' must be a single whole number",
  fixed = TRUE
)
expect_error(
  dbarts::xbart(testData$x, testData$y, seed = c(1, 2)),
  "'seed' must be a single whole number",
  fixed = TRUE
)
expect_error(
  dbarts::dbartsFamilies$hazard(breaks = numeric(0)),
  "hazard 'breaks' must be NULL or hold at least one time",
  fixed = TRUE
)
cutsData <- data.frame(y = testData$y, testData$x)
expect_error(
  do.call(
    dbarts::bart,
    c(list(y ~ ., cutsData, n.cuts = c(X1 = 10L, 20L)), quickArgs)
  ),
  "'n.cuts' names some values and not others; name every value, or none",
  fixed = TRUE
)
expect_error(
  dbarts::pdbart(testData$x, testData$y, pl = NA),
  "'pl' must be TRUE or FALSE",
  fixed = TRUE
)
expect_error(
  dbarts::pd2bart(testData$x, testData$y, pl = NA),
  "'pl' must be TRUE or FALSE",
  fixed = TRUE
)
expect_error(
  dbarts::dbartsControl(proposal.probs = list(swap = numeric(0))),
  "'proposal.probs' entries must be single numbers; 'swap' has 0",
  fixed = TRUE
)
expect_identical(
  dbarts::dbartsControl(
    proposal.probs = list(birth_death = 0.5, swap = 0.1, change = 0.4)
  )@proposal.probs,
  dbarts::dbartsControl(
    proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4)
  )@proposal.probs
)
rm(quickArgs, quickBart, flag, values, value, sigestWords, quickBT, cutsData)

rm(testData)
