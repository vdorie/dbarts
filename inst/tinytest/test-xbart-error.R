source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that fails with invalid inputs
x <- testData$x
y <- testData$y
expect_error(
  dbarts::xbart(y ~ x, method = "not-a-method"),
  "method must be in 'k-fold', 'random subsample'"
)
expect_error(
  dbarts::xbart(y ~ x, method = NULL),
  "method must be in 'k-fold', 'random subsample'"
)
expect_error(
  dbarts::xbart(y ~ x, method = NA_character_),
  "method must be in 'k-fold', 'random subsample'"
)

expect_error(
  dbarts::xbart(y ~ x, n.samples = 0L, n.threads = 1L),
  "'n.samples' must leave at least one draw after thinning"
)
expect_error(
  dbarts::xbart(y ~ x, n.samples = "not-a-integer"),
  "'n.samples' argument to xbart must be coercible to integer type"
)
expect_error(
  dbarts::xbart(y ~ x, n.samples = NULL),
  "'n.samples' must be of length 1"
)
expect_error(
  dbarts::xbart(y ~ x, n.samples = NA_integer_),
  "'n.samples' argument to xbart must be a non-negative integer"
)

expect_error(
  dbarts::xbart(y ~ x, method = "k-fold", n.test = 1),
  "for k-fold crossvalidation, 'n.test' must be an integer in"
)
expect_error(
  dbarts::xbart(y ~ x, method = "k-fold", n.test = length(testData$y) + 1),
  "for k-fold crossvalidation, 'n.test' must be an integer in"
)
expect_error(
  dbarts::xbart(y ~ x, method = "random subsample", n.test = 0),
  "for random subsample crossvalidation, 'n.test' must be in"
)
expect_error(
  dbarts::xbart(
    y ~ x,
    method = "random subsample",
    n.test = length(testData$y) + 1
  ),
  "for random subsample crossvalidation, 'n.test' must be in"
)
expect_error(
  dbarts::xbart(y ~ x, n.test = "not-a-numeric"),
  "'n.test' must be coercible to type: integer"
)
expect_error(dbarts::xbart(y ~ x, n.test = NULL), "'n.test' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, n.test = NA_real_),
  "missing value where TRUE/FALSE needed"
)

expect_error(
  dbarts::xbart(y ~ x, n.reps = 0L),
  "'n.reps' must be a positive integer"
)
expect_error(
  dbarts::xbart(y ~ x, n.reps = "not-a-integer"),
  "'n.reps' must be coercible to type: integer"
)
expect_error(dbarts::xbart(y ~ x, n.reps = NULL), "'n.reps' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, n.reps = NA_integer_),
  "'n.reps' must be a positive integer"
)

expect_error(
  dbarts::xbart(y ~ x, n.burn = c(200L, -1L, 50L)),
  "'n.burn' must contain non-negative integers"
)
expect_error(
  dbarts::xbart(y ~ x, n.burn = "not-a-integer"),
  "'n.burn' must be coercible to type: integer"
)
expect_error(dbarts::xbart(y ~ x, n.burn = NULL), "'n.burn' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, n.burn = NA_integer_),
  "'n.burn' must contain non-negative integers"
)

expect_equal(eval(formals(dbarts::xbart)$n.burn), c(200L, 150L))
# a length-3 n.burn is still accepted; the third element is silently
# dropped, mirroring n.test's own truncate-to-length posture
expect_equal(
  dbarts::xbart(
    y ~ x,
    method = "k-fold",
    n.reps = 1L,
    n.samples = 4L,
    n.burn = c(2L, 1L),
    n.test = 5,
    n.threads = 1L,
    seed = 0L
  ),
  dbarts::xbart(
    y ~ x,
    method = "k-fold",
    n.reps = 1L,
    n.samples = 4L,
    n.burn = c(2L, 1L, 99L),
    n.test = 5,
    n.threads = 1L,
    seed = 0L
  )
)

expect_error(
  dbarts::xbart(y ~ x, loss = "unknown-loss"),
  "loss must be in 'rmse', 'log', 'mcr', or a function"
)
expect_error(
  dbarts::xbart(y ~ x, loss = 2),
  "loss must be in 'rmse', 'log', 'mcr', or a function"
)
expect_error(
  dbarts::xbart(y ~ x, loss = function(x) x),
  "supplied loss function must take exactly three arguments"
)
expect_error(
  dbarts::xbart(y ~ x, loss = list(function(x, y) NULL, "not-an-environment")),
  "supplied loss function must take exactly three arguments"
)

expect_error(
  dbarts::xbart(y ~ x, n.threads = 0L),
  "'n.threads' must be a positive integer"
)
expect_error(
  dbarts::xbart(y ~ x, n.threads = "not-a-integer"),
  "'n.threads' must be coercible to type: integer"
)
expect_error(
  dbarts::xbart(y ~ x, n.threads = NULL),
  "'n.threads' cannot be NULL"
)

expect_error(
  dbarts::xbart(y ~ x, n.trees = 0L),
  "'n.trees' must contain only positive integers"
)
expect_error(
  dbarts::xbart(y ~ x, n.trees = "not-a-integer"),
  "'n.trees' must be coercible to type: integer"
)
expect_error(dbarts::xbart(y ~ x, n.trees = NULL), "'n.trees' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, n.trees = NA_integer_),
  "'n.trees' must contain only positive integers"
)

expect_error(
  dbarts::xbart(y ~ x, k = c(-0.5, 1)),
  "'k' must contain only positive values"
)
expect_error(
  dbarts::xbart(y ~ x, k = "not-a-numeric"),
  "'k' must be coercible to type: numeric"
)
expect_error(
  dbarts::xbart(y ~ x, k = NA_real_),
  "'k' must contain only positive values"
)

expect_error(
  dbarts::xbart(y ~ x, power = c(0, 0.5)),
  "'power' must contain only positive values"
)
expect_error(
  dbarts::xbart(y ~ x, power = "not-a-numeric"),
  "'power' must be coercible to type: numeric"
)
expect_error(dbarts::xbart(y ~ x, power = NULL), "'power' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, power = NA_real_),
  "'power' must contain only positive values"
)

expect_error(
  dbarts::xbart(y ~ x, base = c(0.5, 1)),
  "'base' must contain only values in"
)
expect_error(
  dbarts::xbart(y ~ x, base = "not-a-numeric"),
  "'base' must be coercible to type: numeric"
)
expect_error(dbarts::xbart(y ~ x, base = NULL), "'base' cannot be NULL")
expect_error(
  dbarts::xbart(y ~ x, base = NA_real_),
  "'base' must contain only values in"
)

expect_error(
  dbarts::xbart(y ~ x, verbose = "not-a-logical"),
  "'verbose' argument to xbart must be TRUE/FALSE"
)
expect_error(
  dbarts::xbart(y ~ x, verbose = NULL),
  "'verbose' argument to xbart must be TRUE/FALSE"
)
expect_error(
  dbarts::xbart(y ~ x, verbose = NA),
  "'verbose' argument to xbart must be TRUE/FALSE"
)

expect_error(dbarts::xbart(y ~ x, resid.prior = NULL), "dbartsResidPrior")

expect_error(
  dbarts::xbart(y ~ x, sigest = -1),
  "'sigest' argument to xbart must be positive"
)
expect_error(
  dbarts::xbart(y ~ x, sigest = "not-a-numeric"),
  "'sigest' argument to xbart must be coercible to numeric type"
)


rm(testData)
