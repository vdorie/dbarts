# missing = "incorporate" via the raw matrix (x.train, y.train) interface:
# previously the matrix branch's complete-cases filter ran unconditionally,
# silently dropping NA rows and defeating "incorporate" (only the formula
# interface, which uses na.action = na.pass, could actually reach it).

set.seed(7)
n <- 300L
x1 <- rnorm(n)
x2 <- rnorm(n)
isMissing <- runif(n) < 0.2
x1[isMissing] <- NA_real_
y <- x2 + 2 * isMissing + rnorm(n, 0, 0.3)

xMat <- cbind(x1 = x1, x2 = x2)
df <- data.frame(x1 = x1, x2 = x2, y = y)

# the raw matrix interface now keeps NA-bearing rows under "incorporate"
# (the default), with no spurious "row(s) dropped" warning
expect_silent(data.matrix.inc <- dbarts::dbartsData(xMat, y))
expect_equal(nrow(data.matrix.inc@x), n)
expect_equal(length(data.matrix.inc@y), n)
expect_true(anyNA(data.matrix.inc@x[, "x1"]))
expect_equal(which(is.na(data.matrix.inc@x[, "x1"])), which(isMissing))
expect_equal(data.matrix.inc@missing, "incorporate")

# missing = "error" is unaffected: it still rejects NA-bearing predictors
expect_error(
  dbarts::dbartsData(xMat, y, missing = "error"),
  pattern = "predictors contain missing"
)

# a missing response is never incorporated, matrix interface included - it
# is rejected outright rather than silently dropped, matching the formula
# interface's response handling
y.badY <- y
y.badY[1L] <- NA_real_
expect_error(
  dbarts::dbartsData(xMat, y.badY),
  pattern = "response contains missing"
)

# the formula interface on the identical data, for comparison
data.formula.inc <- dbarts::dbartsData(y ~ x1 + x2, df)
expect_equal(nrow(data.formula.inc@x), n)
expect_true(anyNA(data.formula.inc@x[, "x1"]))

# the two interfaces now agree on what is retained and modeled
expect_equal(
  data.matrix.inc@x[, c("x1", "x2")],
  data.formula.inc@x[, c("x1", "x2")]
)
expect_equal(unname(data.matrix.inc@y), unname(data.formula.inc@y))
expect_equal(data.matrix.inc@varTypes, data.formula.inc@varTypes)

# end to end: a matrix-interface incorporate fit reproduces a formula-
# interface fit on the same data bit for bit, given the same seed - the
# downstream sampler machinery does not distinguish how the data arrived
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  seed = 5L
)
sampler.matrix <- dbarts::dbarts(xMat, y, control = control)
sampler.formula <- dbarts::dbarts(y ~ x1 + x2, df, control = control)
samples.matrix <- sampler.matrix$run(100L, 50L)
samples.formula <- sampler.formula$run(100L, 50L)
expect_identical(samples.matrix$train, samples.formula$train)
expect_identical(samples.matrix$sigma, samples.formula$sigma)
expect_true(any(is.na(sampler.matrix$data@x[, "x1"])))

rm(
  n,
  x1,
  x2,
  isMissing,
  y,
  xMat,
  df,
  data.matrix.inc,
  data.formula.inc,
  y.badY,
  control,
  sampler.matrix,
  sampler.formula,
  samples.matrix,
  samples.formula
)
