source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart works with keepTrainingFits = FALSE
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  testData,
  group.by = g,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 2L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 1L,
  keepTrainingFits = FALSE,
  verbose = FALSE
)
expect_inherits(rbartFit, "rbart")
expect_true(is.null(rbartFit$yhat.train))
expect_true(is.null(rbartFit$yhat.train.mean))

rm(rbartFit)

# test that rbart works with k as a variable
# test thanks to Bruno Tancredi

k <- 1.8
expect_inherits(
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    n.samples = 2L,
    n.burn = 0L,
    n.thin = 2L,
    n.chains = 2L,
    n.trees = 3L,
    n.threads = 1L,
    k = k,
    keepTrainingFits = FALSE,
    verbose = FALSE
  ),
  "rbart"
)
rm(k)


# A user prior is taken by name, through a wrapper's argument, or as a
# built-in's name; an unknown name is refused, listing the built-ins.
priorFrame <- data.frame(y = rnorm(60L), x = runif(60L))
priorGroups <- rep(1:6, 10L)
myPrior <- function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE)
fitPrior <- function(prior) {
  suppressWarnings(rbart_vi(
    y ~ x,
    priorFrame,
    group.by = priorGroups,
    prior = prior,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  ))
}
fit <- suppressWarnings(rbart_vi(
  y ~ x,
  priorFrame,
  group.by = priorGroups,
  prior = myPrior,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
))
expect_true(inherits(fit, "rbart"))
expect_true(inherits(fitPrior("gamma"), "rbart"))
expect_true(inherits(fitPrior(myPrior), "rbart"))
expect_error(fitPrior("nope"), "'cauchy', 'gamma'")
rm(priorFrame, priorGroups, myPrior, fitPrior, fit)

rm(testData)
