source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart fails with invalid group.by
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = NA, n.threads = 1L),
  "'group.by' must be coercible to factor type"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = not_a_symbol, n.threads = 1L),
  "'group.by' not found"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = testData$g[-1L], n.threads = 1L),
  "group.by' not of length equal to that of data"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = "not a factor", n.threads = 1L),
  "'group.by' not of length equal to that of data"
)

rm(testData)

# a seed argument fixes .Random.seed only for the duration of the call: a
# custom tau prior routes rbart_vi through the R callback fit path, so
# making it throw errors out from inside that fixed-seed window and checks
# the caller's stream is still restored rather than left on the fixed seed
seedProbeData <- data.frame(
  x = rnorm(20L),
  y = rnorm(20L),
  g = factor(rep(1:2, 10L))
)
failingPrior <- function(x, rel.scale) stop("engineered failure")

set.seed(0L)
preCallState <- .GlobalEnv$.Random.seed

expect_error(
  dbarts::rbart_vi(
    y ~ x,
    seedProbeData,
    group.by = g,
    prior = failingPrior,
    n.trees = 3L,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 0L,
    n.samples = 1L,
    n.thin = 1L,
    keepTrees = FALSE,
    verbose = FALSE,
    seed = 99L
  ),
  "engineered failure"
)
expect_equal(.GlobalEnv$.Random.seed, preCallState)

rm(seedProbeData, failingPrior, preCallState)
