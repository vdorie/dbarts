source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that prior model arguments create valid objects
expect_inherits(
  dbarts::dbarts(
    y ~ x,
    testData,
    verbose = FALSE,
    n.samples = 500,
    tree.prior = cgm(0.75, 0.5),
    node.prior = normal(3.5),
    family = gaussian(sigma = chisq(5, 0.9)),
    sigest = 1.0,
    control = dbartsControl(n.threads = 1L, n.chains = 1L),
  ),
  "dbartsSampler"
)
expect_inherits(
  dbarts::dbarts(
    y ~ x,
    testData,
    verbose = FALSE,
    n.samples = 500,
    tree.prior = cgm(0.75, 0.5),
    node.prior = normal(chi(1.0, Inf)),
    family = gaussian(sigma = chisq(5, 0.9)),
    sigest = 1.0,
    control = dbartsControl(n.threads = 1L, n.chains = 1L)
  ),
  "dbartsSampler"
)
k <- "chi(1.0, 2.0)"
expect_inherits(
  dbarts(
    y ~ x,
    testData,
    verbose = FALSE,
    n.samples = 500,
    tree.prior = cgm(0.75, 0.5),
    node.prior = normal(k),
    family = gaussian(sigma = chisq(5, 0.9)),
    sigest = 1.0,
    control = dbartsControl(n.threads = 1L, n.chains = 1L)
  ),
  "dbartsSampler"
)
rm(k)

rm(testData)
