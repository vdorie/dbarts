source(
  system.file("common", "almostLinearBinaryData.R", package = "dbarts"),
  local = TRUE
)

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
    data[1L:100L, ],
    data[101L:200L, ],
    nskip = 0L,
    ndpost = 1L,
    verbose = FALSE
  ),
  "bart"
)
rm(modelFormula)

# a one-sided formula (~ x.1 + x.2) names no response: dbartsData() defaults
# it to zeros for a composed sampler, but a fitting entry point used to
# return that all-zero fit silently instead of refusing it by name
expect_error(
  dbarts::bart(~ x.1 + x.2, data, nskip = 0L, ndpost = 1L, verbose = FALSE),
  "bart() requires a two-sided formula",
  fixed = TRUE
)
expect_error(
  dbarts::bart2(
    ~ x.1 + x.2,
    data,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "bart2() requires a two-sided formula",
  fixed = TRUE
)
expect_error(
  dbarts::xbart(~ x.1 + x.2, data, n.samples = 3L, verbose = FALSE),
  "xbart() requires a two-sided formula",
  fixed = TRUE
)
expect_error(
  dbarts::rbart_vi(
    ~ x.1 + x.2,
    data,
    group.by = rep(1:2, length.out = nrow(data)),
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "rbart_vi() requires a two-sided formula",
  fixed = TRUE
)

# a ':'/'*' term expands to a term.labels entry the model frame never
# carries, dying inside makeModelMatrix with R's own "undefined columns
# selected"; refused by name instead, naming what IS supported
interactionReason <- "':' and '*' terms are not supported in 'formula'"
expect_error(
  dbarts::bart2(
    y ~ x.1 * x.2,
    data,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  interactionReason,
  fixed = TRUE
)
expect_error(
  dbarts::bart2(
    y ~ x.1:x.2 + x.3,
    data,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  interactionReason,
  fixed = TRUE
)
rm(interactionReason)

# poly()/log()/offset() terms are unaffected - each still fits
expect_inherits(
  dbarts::bart2(
    y ~ poly(x.1, 2) + log(abs(x.2) + 1) + offset(x.3),
    data,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "bart"
)

rm(data)

rm(testData)
