source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that basic probit example with flat hyperprior superior to default
n.sims <- 200L
n.burn <- 100L

set.seed(99L)
bartFit <- dbarts::bartBT(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  verbose = FALSE
)

set.seed(99L)
bartFit.flat <- dbarts::bartBT(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  k = chi(1, Inf),
  verbose = FALSE
)

expect_true(
  cor(qnorm(testData$p), colMeans(bartFit$yhat.train)) <
    cor(qnorm(testData$p), colMeans(bartFit.flat$yhat.train))
)
rm(bartFit.flat, bartFit, n.burn, n.sims)

# test_that binary model with k hyperprior is reproducible when multithreaded
fit1 <- dbarts::bart(
  testData$X[1L:100L, ],
  testData$Z[1L:100L],
  n.trees = 5L,
  n.samples = 100L,
  n.burn = 0L,
  n.threads = 2L,
  n.chains = 2L,
  seed = 99L,
  verbose = FALSE
)
fit2 <- dbarts::bart(
  testData$X[1L:100L, ],
  testData$Z[1L:100L],
  n.trees = 5L,
  n.samples = 100L,
  n.burn = 0L,
  n.threads = 2L,
  n.chains = 2L,
  seed = 99L,
  verbose = FALSE
)
expect_equal(fit1$yhat.train, fit2$yhat.train)
rm(fit2, fit1)

source(
  system.file("common", "almostLinearBinaryData.R", package = "dbarts"),
  local = TRUE
)

fitSubset <- 1L:100L
testSubset <- 101L:200L

fitData <- list(y = testData$y[fitSubset], x = testData$x[fitSubset, ])
mu <- testData$mu[testSubset]

glmFit <- stats::glm(y ~ x, fitData, family = binomial(link = "probit"))

predictData <- list(x = testData$x[testSubset, ])
mu.hat.glm <- predict(glmFit, newdata = predictData)


set.seed(99L)
bartFit <- dbarts::bartBT(
  testData$x[fitSubset, ],
  testData$y[fitSubset],
  testData$x[testSubset, ],
  binaryOffset = testData$offset,
  verbose = FALSE
)
mu.hat.bart <- colMeans(bartFit$yhat.test)

# test that binary example using close to linear function provides sensible results
expect_true(cor(mu, mu.hat.glm) < cor(mu, mu.hat.bart))
expect_true((range(mu.hat.bart) * 1.2)[1L] >= range(mu)[1L])
expect_true((range(mu.hat.bart) * 1.2)[2L] <= range(mu)[2L])

# test_that binary example using a flat prior is similar to default in tuned model
set.seed(99L)
bartFit.flat <- dbarts::bartBT(
  testData$x[fitSubset, ],
  testData$y[fitSubset],
  testData$x[testSubset, ],
  binaryOffset = testData$offset,
  verbose = FALSE,
  k = chi(1, Inf)
)
mu.hat.bart.flat <- colMeans(bartFit.flat$yhat.test)


expect_true(cor(mu.hat.bart, mu.hat.bart.flat) > 0.95)
# the flat-prior k posterior is heavy-tailed and its median swings widely
# across seeds; assert it mixes and concentrates at modest
# values rather than pinning a tight seed-locked bound
expect_true(length(unique(bartFit.flat$k)) > 100L)
expect_true(median(bartFit.flat$k) < 15)

# An improper scale (chi's scale = Inf) leaves k's posterior improper: with one
# tree and a large df the chain carries k to infinity, where every leaf is
# pinned at zero. The run finishes with the offset as its fit and no NaN, and
# the trees keep moving, now drawn from their prior.
x.runaway <- testData$x[fitSubset, ]
y.runaway <- testData$y[fitSubset]
offset.runaway <- rep(testData$offset, length(fitSubset))
fitRunaway <- function(family, n.burn = 100L, n.samples = 100L, ...) {
  dbarts::bart(
    x.runaway,
    y.runaway,
    offset = offset.runaway,
    family = family,
    k = chi(1000, Inf),
    n.trees = 1L,
    n.burn = n.burn,
    n.samples = n.samples,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    ...
  )
}
for (family in c("probit", "logistic")) {
  bartFit.runaway <- fitRunaway(family, seed = 99L)
  k.runaway <- bartFit.runaway$k
  atInfinity <- is.infinite(k.runaway)
  expect_false(anyNA(k.runaway))
  expect_true(sum(atInfinity) > 10L && atInfinity[length(k.runaway)])
  expect_false(anyNA(bartFit.runaway$yhat.train))
  expect_true(all(
    abs(bartFit.runaway$yhat.train[atInfinity, ] - testData$offset) < 1e-12
  ))
  expect_true(
    nrow(unique(bartFit.runaway$varcount[atInfinity, , drop = FALSE])) > 1L
  )
}

# an infinite k survives storeState and setState, saveRDS and readRDS, and a
# warm start
bartFit.runaway$fit$storeState()
state.runaway <- bartFit.runaway$fit$state
expect_true(is.infinite(state.runaway[[1L]]$forests[[1L]]$k))
bartFit.restored <- fitRunaway("logistic", n.burn = 0L, seed = 1L)
bartFit.restored$fit$setState(state.runaway)
bartFit.restored$fit$storeState()
expect_true(is.infinite(bartFit.restored$fit$state[[1L]]$forests[[1L]]$k))
samples.restored <- bartFit.restored$fit$run(0L, 5L)
expect_true(all(is.infinite(samples.restored$k)))
expect_true(all(abs(samples.restored$train - testData$offset) < 1e-12))

serialized <- tempfile(fileext = ".rds")
saveRDS(bartFit.runaway, serialized)
bartFit.loaded <- readRDS(serialized)
unlink(serialized)
expect_equal(
  predict(bartFit.loaded, x.runaway, offset = offset.runaway, type = "link"),
  predict(bartFit.runaway, x.runaway, offset = offset.runaway, type = "link")
)
bartFit.loaded$fit$storeState()
expect_true(is.infinite(bartFit.loaded$fit$state[[1L]]$forests[[1L]]$k))

bartFit.warm <- fitRunaway(
  "logistic",
  n.burn = 0L,
  n.samples = 5L,
  seed = 1L,
  warm.start = bartFit.runaway
)
expect_true(all(is.infinite(bartFit.warm$k)))
expect_true(all(abs(bartFit.warm$yhat.train - testData$offset) < 1e-12))
rm(
  bartFit.warm,
  bartFit.loaded,
  serialized,
  samples.restored,
  bartFit.restored,
  state.runaway,
  atInfinity,
  k.runaway,
  bartFit.runaway,
  family,
  fitRunaway,
  offset.runaway,
  y.runaway,
  x.runaway
)

rm(
  mu.hat.bart.flat,
  bartFit.flat,
  mu.hat.bart,
  bartFit,
  mu.hat.glm,
  predictData,
  glmFit,
  mu,
  fitData,
  testSubset,
  fitSubset
)

rm(testData)
