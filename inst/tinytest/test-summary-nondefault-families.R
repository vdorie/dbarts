# summary() methods for the three bart2() families with no earlier method
# (bartOrdinal, bartNegbin, bartHurdle): before these existed, summary()
# fell through to summary.default and printed a raw Length/Class/Mode table.
# Pinned here: each returns a non-default summary object whose print output
# names the family-specific rows (thresholds, dispersion, or both hurdle
# components), reusing summary.bart's own row-table machinery rather than
# duplicating it.

set.seed(61)
n <- 50L
x <- matrix(runif(n * 2L), n, 2L)

seededControl <- list(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 8L,
  n.samples = 15L,
  n.burn = 5L,
  verbose = FALSE
)

## --- ordinal: thresholds, not summary.default -------------------------------
yOrdinal <- ordered(
  sample(c("lo", "mid", "hi"), n, replace = TRUE),
  levels = c("lo", "mid", "hi")
)
fitOrdinal <- do.call(
  bart2,
  c(list(x, yOrdinal), seededControl)
)
expect_false(identical(class(summary(fitOrdinal)), "summary.default"))
summaryOrdinal <- summary(fitOrdinal)
expect_true(inherits(summaryOrdinal, "summary.bart"))
ordinalOutput <- capture.output(print(summaryOrdinal))
expect_true(any(grepl("threshold\\[1\\]", ordinalOutput)))
expect_true(any(grepl("threshold\\[2\\]", ordinalOutput)))

## --- nbinom: dispersion, not summary.default -------------------------------
yCount <- rpois(n, lambda = 4)
fitNegbin <- do.call(
  bart2,
  c(list(x, yCount, family = "nbinom"), seededControl)
)
expect_false(identical(class(summary(fitNegbin)), "summary.default"))
summaryNegbin <- summary(fitNegbin)
expect_true(inherits(summaryNegbin, "summary.bart"))
negbinOutput <- capture.output(print(summaryNegbin))
expect_true(any(grepl("^1 dispersion", negbinOutput)))

## --- hurdle: both components, not summary.default --------------------------
yHurdle <- ifelse(runif(n) < 0.3, 0, rexp(n, rate = 0.5))
fitHurdle <- do.call(
  bart2,
  c(list(x, yHurdle, family = "hurdle.lognormal"), seededControl)
)
expect_false(identical(class(summary(fitHurdle)), "summary.default"))
summaryHurdle <- summary(fitHurdle)
expect_true(inherits(summaryHurdle, "summary.bartHurdle"))
hurdleOutput <- capture.output(print(summaryHurdle))
expect_true(any(grepl("Occupancy component", hurdleOutput)))
expect_true(any(grepl("Positive-part component", hurdleOutput)))
# each component's own summary.bart object is reused, not recomputed
expect_identical(summaryHurdle$occupancy, summary(fitHurdle$occupancy))
expect_identical(summaryHurdle$positive, summary(fitHurdle$positive))

## --- print() already has a method for all three; no fallback to print.default
expect_false(identical(class(fitOrdinal)[1L], "list"))
expect_true(exists("print.bartOrdinal", where = asNamespace("dbarts")))
expect_true(exists("print.bartNegbin", where = asNamespace("dbarts")))
expect_true(exists("print.bartHurdle", where = asNamespace("dbarts")))
