source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that fails with invalid inputs
x <- testData$x
y <- testData$y
expect_error(dbarts::xbart(y ~ x, method = "not-a-method"))
expect_error(dbarts::xbart(y ~ x, method = NULL))
expect_error(dbarts::xbart(y ~ x, method = NA_character_))

expect_error(dbarts::xbart(y ~ x, n.samples = 0L))
expect_error(dbarts::xbart(y ~ x, n.samples = "not-a-integer"))
expect_error(dbarts::xbart(y ~ x, n.samples = NULL))
expect_error(dbarts::xbart(y ~ x, n.samples = NA_integer_))

expect_error(dbarts::xbart(y ~ x, method = "k-fold", n.test = 1))
expect_error(
  dbarts::xbart(y ~ x, method = "k-fold", n.test = length(testData$y) + 1)
)
expect_error(
  dbarts::xbart(y ~ x, method = "random subsample", n.test = 0)
)
expect_error(
  dbarts::xbart(y ~ x, method = "random subsample", n.test = length(testData$y) + 1)
)
expect_error(dbarts::xbart(y ~ x, n.test = "not-a-numeric"))
expect_error(dbarts::xbart(y ~ x, n.test = NULL))
expect_error(dbarts::xbart(y ~ x, n.test = NA_real_))

expect_error(dbarts::xbart(y ~ x, n.reps = 0L))
expect_error(dbarts::xbart(y ~ x, n.reps = "not-a-integer"))
expect_error(dbarts::xbart(y ~ x, n.reps = NULL))
expect_error(dbarts::xbart(y ~ x, n.reps = NA_integer_))

expect_error(dbarts::xbart(y ~ x, n.burn = c(200L, -1L, 50L)))
expect_error(dbarts::xbart(y ~ x, n.burn = "not-a-integer"))
expect_error(dbarts::xbart(y ~ x, n.burn = NULL))
expect_error(dbarts::xbart(y ~ x, n.burn = NA_integer_))

expect_error(dbarts::xbart(y ~ x, loss = "unknown-loss"))
expect_error(dbarts::xbart(y ~ x, loss = 2))
expect_error(dbarts::xbart(y ~ x, loss = function(x) x))
expect_error(
  dbarts::xbart(y ~ x, loss = list(function(x, y) NULL, "not-an-environment"))
)

expect_error(dbarts::xbart(y ~ x, n.threads = 0L))
expect_error(dbarts::xbart(y ~ x, n.threads = "not-a-integer"))
expect_error(dbarts::xbart(y ~ x, n.threads = NULL))

expect_error(dbarts::xbart(y ~ x, n.trees = 0L))
expect_error(dbarts::xbart(y ~ x, n.trees = "not-a-integer"))
expect_error(dbarts::xbart(y ~ x, n.trees = NULL))
expect_error(dbarts::xbart(y ~ x, n.trees = NA_integer_))

expect_error(dbarts::xbart(y ~ x, k = c(-0.5, 1)))
expect_error(dbarts::xbart(y ~ x, k = "not-a-numeric"))
expect_error(dbarts::xbart(y ~ x, k = NA_real_))

expect_error(dbarts::xbart(y ~ x, power = c(0, 0.5)))
expect_error(dbarts::xbart(y ~ x, power = "not-a-numeric"))
expect_error(dbarts::xbart(y ~ x, power = NULL))
expect_error(dbarts::xbart(y ~ x, power = NA_real_))

expect_error(dbarts::xbart(y ~ x, base = c(0.5, 1)))
expect_error(dbarts::xbart(y ~ x, base = "not-a-numeric"))
expect_error(dbarts::xbart(y ~ x, base = NULL))
expect_error(dbarts::xbart(y ~ x, base = NA_real_))

expect_error(dbarts::xbart(y ~ x, verbose = "not-a-logical"))
expect_error(dbarts::xbart(y ~ x, verbose = NULL))
expect_error(dbarts::xbart(y ~ x, verbose = NA))

expect_error(dbarts::xbart(y ~ x, resid.prior = NULL))

expect_error(dbarts::xbart(y ~ x, sigma = -1))
expect_error(dbarts::xbart(y ~ x, sigma = "not-a-numeric"))


rm(testData)

