context("dbarts model arguments")

source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

test_that("non-prior, model arguments raise errors", {
  expect_error(dbarts(y ~ x, testData, verbose = NA))
  expect_error(dbarts(y ~ x, testData, verbose = "not-a-logical"))

  expect_error(dbarts(y ~ x, testData, n.samples = -1L))
  expect_error(dbarts(y ~ x, testData, n.samples = "not-an-integer"))

  expect_error(dbarts(y ~ x, testData, sigma = -1.0))
  expect_error(dbarts(y ~ x, testData, sigma = "not-an-integer"))
})

test_that("prior model arguments raise errors", {
  expect_error(dbarts(y ~ x, testData, tree.prior = normal))
  expect_error(dbarts(y ~ x, testData, tree.prior = cgm(0, 0)))
  expect_error(dbarts(y ~ x, testData, tree.prior = cgm(1, 0)))
  expect_error(dbarts(y ~ x, testData, tree.prior = cgm(1, 0, "extra")))
  expect_error(dbarts(y ~ x, testData, tree.prior = cgm(1, 1)))
  
  expect_error(dbarts(y ~ x, testData, node.prior = cgm))
  expect_error(dbarts(y ~ x, testData, node.prior = normal(0)))
  expect_error(dbarts(y ~ x, testData, node.prior = normal(normal)))
  expect_error(dbarts(y ~ x, testData, node.prior = normal(chi(scale = -1))))
  expect_error(dbarts(y ~ x, testData, node.prior = normal(chi(not_n_arg = 2.0))))

  expect_error(dbarts(y ~ x, testData, resid.prior = binomial))
  expect_error(dbarts(y ~ x, testData, resid.prior = chisq(0, 0)))
  expect_error(dbarts(y ~ x, testData, resid.prior = chisq(1, 0)))
})

test_that("prior model arguments create valid objects", {
  expect_is(dbarts(y ~ x, testData, verbose = FALSE, n.samples = 500,
                   tree.prior = cgm(0.75, 0.5), node.prior = normal(3.5),
                   resid.prior = chisq(5, 0.9), sigma = 1.0,
                   control = dbartsControl(n.threads = 1L, n.chains = 1L)),
            "dbartsSampler")
  expect_is(dbarts(y ~ x, testData, verbose = FALSE, n.samples = 500,
                   tree.prior = cgm(0.75, 0.5), node.prior = normal(chi(1.0, Inf)),
                   resid.prior = chisq(5, 0.9), sigma = 1.0,
                   control = dbartsControl(n.threads = 1L, n.chains = 1L)),
            "dbartsSampler")
  k <- "chi(1.0, 2.0)"
  expect_is(dbarts(y ~ x, testData, verbose = FALSE, n.samples = 500,
                   tree.prior = cgm(0.75, 0.5), node.prior = normal(k),
                   resid.prior = chisq(5, 0.9), sigma = 1.0,
                   control = dbartsControl(n.threads = 1L, n.chains = 1L)),
            "dbartsSampler")

})

