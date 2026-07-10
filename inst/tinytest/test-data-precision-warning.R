# dbartsData warns when the response is precision-degenerate - a huge offset
# with a tiny range quantizes to (near-)identical doubles before the engine
# ever sees it - detected by range(y)/max(abs(y)) < 1e-10. There is
# deliberately no distinct-value-count check: a genuine rounding collapse
# always trips the ratio (a few surviving doubles means a range of a few
# ulps, ~1e-15 relative), while low-cardinality data the ratio does not
# flag is legitimately discrete and fits cleanly after standardization.

set.seed(0)
n <- 200L
x <- matrix(runif(n * 2), n)

# the review-6 case: y in [1e15, 1e15 + 1e-3] collapses to a single
# representable double
y.huge <- 1e15 + runif(n) * 1e-3
expect_warning(
  dbarts::dbartsData(x, y.huge),
  pattern = "indistinguishable"
)
expect_true(length(unique(y.huge)) < n)

# a merely large, but not precision-degenerate, response does not warn:
# range is still a healthy fraction of scale
y.large.ok <- 1e6 + rnorm(n, 0, 100)
expect_silent(dbarts::dbartsData(x, y.large.ok))

# few distinct values atop a moderate offset stay silent: the 3 values are
# separated by ~1e10 ulps - perfectly distinguishable doubles that the
# pre-fit standardization maps cleanly - so this is discrete data, not a
# precision artifact
y.fewDistinct <- 1e6 + sample(c(0, 1, 2), n, replace = TRUE)
expect_true(diff(range(y.fewDistinct)) / max(abs(y.fewDistinct)) >= 1e-10)
expect_silent(dbarts::dbartsData(x, y.fewDistinct))

# ordinary continuous data never warns
y.ok <- rnorm(n)
expect_silent(dbarts::dbartsData(x, y.ok))

# a constant, but not all-zero, response is caught by the range/scale ratio
y.const <- rep_len(5, n)
expect_warning(
  dbarts::dbartsData(x, y.const),
  pattern = "indistinguishable"
)

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
