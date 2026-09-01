source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that non-prior, model arguments raise errors
expect_error(
  dbarts::dbarts(y ~ x, testData, verbose = NA),
  "'verbose' argument to dbarts must be TRUE/FALSE"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, verbose = "not-a-logical"),
  "'verbose' argument to dbarts must be TRUE/FALSE"
)

expect_error(
  dbarts::dbarts(y ~ x, testData, n.samples = -1L),
  "'n.samples' argument to dbarts must be a non-negative integer"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, n.samples = "not-an-integer"),
  "'n.samples' must be coercible to type: integer"
)

expect_error(
  dbarts::dbarts(y ~ x, testData, sigma = -1.0),
  "'sigma' argument to dbarts must be positive"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, sigma = "not-an-integer"),
  "'sigma' argument to dbarts must be coercible to numeric type"
)

# test that prior model arguments raise errors
expect_error(
  dbarts::dbarts(y ~ x, testData, tree.prior = normal),
  "'tree.prior' must be a tree prior specification"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, tree.prior = cgm(0, 0)),
  "'power' must be positive"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, tree.prior = cgm(1, 0)),
  "'base' must be in \\(0, 1\\)"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, tree.prior = cgm(1, 0, "extra")),
  "'base' must be in \\(0, 1\\)"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, tree.prior = cgm(1, 1)),
  "'base' must be in \\(0, 1\\)"
)

expect_error(
  dbarts::dbarts(y ~ x, testData, node.prior = cgm),
  "'node.prior' must be a node prior specification"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, node.prior = normal(0)),
  "'k' must be a positive scalar"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, node.prior = normal(normal)),
  "'k' must be a positive scalar"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, node.prior = normal(chi(scale = -1))),
  "'scale' must be positive"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, node.prior = normal(chi(not_n_arg = 2.0))),
  "unused argument \\(not_n_arg = 2\\)"
)

expect_error(
  dbarts::dbarts(y ~ x, testData, resid.prior = binomial),
  "'resid.prior' must be a residual prior specification"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, resid.prior = chisq(0, 0)),
  "'df' must be positive"
)
expect_error(
  dbarts::dbarts(y ~ x, testData, resid.prior = chisq(1, 0)),
  "'quantile' must be positive"
)

rm(testData)
