# A data object that holds a starting sigma, a sampler's data, hands it on to
# dbarts() and bart() as it does to dbartsSpec(): the number stays unless
# sigest gives another; a data object that holds none is regressed as ever.

set.seed(20261111L)
n <- 60L
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- d$x + rnorm(n)
ctl <- dbartsControl(n.chains = 1L, n.threads = 1L, verbose = FALSE)
bartArgs <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)

s <- dbarts(y ~ x, d, sigest = 0.37, n.trees = 5L, control = ctl)
expect_identical(s$data@sigma, 0.37)
regressed <- dbarts(y ~ x, d, n.trees = 5L, control = ctl)$data@sigma
expect_true(regressed != 0.37)

# both doors keep it
expect_identical(dbarts(s$data, n.trees = 5L, control = ctl)$data@sigma, 0.37)
expect_identical(
  do.call(bart, c(list(s$data), bartArgs))$fit$data@sigma,
  0.37
)
expect_identical(dbartsSpec(s$data)$data@sigma, 0.37)

# a stated sigest overrides it
expect_identical(
  dbarts(s$data, sigest = 0.5, n.trees = 5L, control = ctl)$data@sigma,
  0.5
)
expect_identical(
  do.call(bart, c(list(s$data, sigest = 0.5), bartArgs))$fit$data@sigma,
  0.5
)
# and a rule to estimate by asks for the estimate
expect_identical(
  dbarts(s$data, sigest = "dense", n.trees = 5L, control = ctl)$data@sigma,
  regressed
)

# a data object that holds no estimate is regressed
plain <- dbartsData(y ~ x, d)
expect_true(is.na(plain@sigma))
expect_identical(
  dbarts(plain, n.trees = 5L, control = ctl)$data@sigma,
  regressed
)
expect_identical(
  do.call(bart, c(list(plain), bartArgs))$fit$data@sigma,
  regressed
)
