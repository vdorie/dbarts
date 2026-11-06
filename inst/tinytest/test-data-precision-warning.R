# dbartsData warns when the response is precision-degenerate - a huge offset
# with a tiny range quantizes to (near-)identical doubles before the engine
# ever sees it - detected by range(y)/max(abs(y)) < 1e-10. There is
# deliberately no distinct-value-count check: a genuine rounding collapse
# always trips the ratio (a few surviving doubles means a range of a few
# ulps, ~1e-15 relative), while low-cardinality data the ratio does not
# flag is legitimately discrete and fits cleanly after standardization.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(0)
n <- 200L
x <- matrix(runif(n * 2), n)

# the review-6 case: y in [1e15, 1e15 + 1e-3] collapses to a single
# representable double
y.huge <- 1e15 + runif(n) * 1e-3
warnings.huge <- captureWarnings(dbarts::dbartsData(x, y.huge))
expect_equal(length(warnings.huge), 1L)
expect_match(conditionMessage(warnings.huge[[1L]]), "indistinguishable")
expect_true(length(unique(y.huge)) < n)

# the threshold itself, from below. A response whose range is 1e-14 of its
# scale keeps plenty of distinct doubles, so nothing about it is degenerate to
# a cardinality check - it is the RATIO that is small, and the 1e-10 constant
# is what has to catch it
y.near <- 1e15 + runif(n) * 10
expect_true(length(unique(y.near)) > 50L)
warnings.near <- captureWarnings(dbarts::dbartsData(x, y.near))
expect_equal(length(warnings.near), 1L)
expect_true(any(grepl(
  "indistinguishable",
  vapply(warnings.near, conditionMessage, "")
)))

# a merely large, but not precision-degenerate, response does not warn:
# range is still a healthy fraction of scale
y.large.ok <- 1e6 + rnorm(n, 0, 100)
expect_silent(dbarts::dbartsData(x, y.large.ok))

# few distinct values atop a moderate offset stay silent: the 3 values are
# separated by ~1e10 ulps - perfectly distinguishable doubles that the
# pre-fit standardization maps cleanly - so this is discrete data, not a
# precision artifact
y.fewDistinct <- 1e6 + sample(c(0, 1, 2), n, replace = TRUE)
expect_equal(diff(range(y.fewDistinct)), 2)
expect_silent(dbarts::dbartsData(x, y.fewDistinct))

# ordinary continuous data never warns
y.ok <- rnorm(n)
expect_silent(dbarts::dbartsData(x, y.ok))

# a constant, but not all-zero, response is caught by the range/scale ratio
y.const <- rep_len(5, n)
warnings.constResp <- captureWarnings(dbarts::dbartsData(x, y.const))
expect_equal(length(warnings.constResp), 1L)
expect_match(conditionMessage(warnings.constResp[[1L]]), "indistinguishable")

# an all-zero response is guarded against the 0/0 division and is not
# itself a precision artifact (it is exactly constant, not approximately
# so), so it does not trip this warning
y.allZero <- rep_len(0, n)
expect_silent(dbarts::dbartsData(x, y.allZero))

# regression guards against any future cardinality-based check: binary
# responses (the one family dbarts formally treats as discrete) and
# small-n categorical responses stay silent
y.binary <- rep_len(c(0, 1), n)
expect_silent(dbarts::dbartsData(x, y.binary))

x.small <- matrix(runif(10), 5L)
y.binary.small <- c(0, 1, 0, 1, 0)
expect_silent(dbarts::dbartsData(x.small, y.binary.small))

x.smallN <- matrix(runif(80), 40L)
y.categorical.smallN <- sample(1:3, 40L, replace = TRUE)
expect_silent(dbarts::dbartsData(x.smallN, y.categorical.smallN))

rm(
  n,
  x,
  y.huge,
  y.near,
  y.large.ok,
  y.fewDistinct,
  y.ok,
  y.const,
  y.allZero,
  y.binary,
  x.small,
  y.binary.small,
  x.smallN,
  y.categorical.smallN
)

# A count or a category has nothing to rescale, so no fit of one passes the
# warning on: an ordinal or multinomial fit to one observed category gets its
# own one warning at every level, a constant count is silent, and a fit its
# family then refuses (a hurdle response with no zero, a probit response not
# coded 0/1) gets the refusal alone. A gaussian fit to a constant still warns.
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)
set.seed(9)
n <- 40L
x <- matrix(runif(n * 2), n)
singleCategory <- paste0(
  "is fitted to a response with a single observed category; the prior alone ",
  "drives the fitted probabilities"
)
fitQuietly <- function(...) {
  bart(
    x,
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}
for (level in c("a", "b", "c")) {
  y.ordinal <- factor(rep(level, n), levels = c("a", "b", "c"), ordered = TRUE)
  for (family in c("ordinal", "auto")) {
    warned <- captureWarnings(fitQuietly(y.ordinal, family = family))
    expect_equal(length(warned), 1L, info = paste(family, level))
    expect_identical(
      conditionMessage(warned[[1L]]),
      paste0("family \"ordinal\" ", singleCategory),
      info = paste(family, level)
    )
  }
  warned <- captureWarnings(dbarts::dbarts(
    x,
    y.ordinal,
    family = "ordinal",
    n.trees = 5L
  ))
  expect_equal(length(warned), 1L, info = level)
  y.multinomial <- factor(rep(level, n), levels = c("a", "b", "c"))
  warned <- captureWarnings(fitQuietly(y.multinomial, family = "multinomial"))
  expect_equal(length(warned), 1L, info = level)
  expect_identical(
    conditionMessage(warned[[1L]]),
    paste0("family \"multinomial\" ", singleCategory),
    info = level
  )
}
# a count matrix with one category ever counted is the same response
counts <- cbind(a = 0L, b = rpois(n, 2) + 1L, c = 0L)
warned <- captureWarnings(fitQuietly(counts, family = "multinomial"))
expect_equal(length(warned), 1L)
expect_match(conditionMessage(warned[[1L]]), singleCategory, fixed = TRUE)
# two observed categories of three are an ordinary fit
expect_silent(fitQuietly(
  factor(rep(c("a", "b"), length.out = n), levels = c("a", "b", "c")),
  family = "multinomial"
))
expect_equal(
  length(captureWarnings(fitQuietly(rep(3L, n), family = "nbinom"))),
  0L
)
warned <- captureWarnings(expect_error(
  fitQuietly(rep(2.5, n), family = "hurdle.lognormal"),
  "at least one exact zero"
))
expect_equal(length(warned), 0L)
warned <- captureWarnings(expect_error(
  fitQuietly(rep(2, n), family = "probit"),
  "requires a response coded 0/1"
))
expect_equal(length(warned), 0L)
warned <- captureWarnings(fitQuietly(rep(5, n), family = "gaussian"))
expect_equal(length(warned), 1L)
expect_match(conditionMessage(warned[[1L]]), "indistinguishable")
rm(
  n,
  x,
  singleCategory,
  fitQuietly,
  level,
  y.ordinal,
  family,
  warned,
  y.multinomial,
  counts
)
