# Internal bartcore engine surface (R/bartcore.R, src/bartcore/); the
# statistical-equivalence gates live outside the package in benchmarks/.

set.seed(99)
n <- 200L
p <- 5L
x <- matrix(runif(n * p), n, p)
f <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 5 * x[, 4L]
y <- f + rnorm(n)
x.test <- matrix(runif(10L * p), 10L, p)

control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 50L,
                         updateState = FALSE)
sampler <- dbarts(x, y, test = x.test, control = control)
bcSampler <- dbarts:::bartcoreSampler(sampler)

result <- dbarts:::bartcoreRun(bcSampler, 100L, 200L)

expect_equal(dim(result$yhat.train), c(n, 200L))
expect_equal(dim(result$yhat.test), c(10L, 200L))
expect_equal(dim(result$varcount), c(p, 200L))
expect_true(all(result$sigma > 0))

# the fit should explain most of the signal
fitMean <- rowMeans(result$yhat.train)
expect_true(mean((fitMean - f)^2) < 0.25 * mean((mean(y) - f)^2))

# embedded-Gibbs pattern: mutate offset between single draws
dbarts:::bartcoreSetOffset(bcSampler, rep(0.5, n))
result.offset <- dbarts:::bartcoreRun(bcSampler, 0L, 1L)
expect_equal(dim(result.offset$yhat.train), c(n, 1L))

dbarts:::bartcoreSetResponse(bcSampler, y + 1)
result.response <- dbarts:::bartcoreRun(bcSampler, 0L, 1L)
expect_true(mean(result.response$yhat.train) > mean(result.offset$yhat.train))

# no latents for a continuous response
expect_null(dbarts:::bartcoreGetLatents(bcSampler))

# binary response: probit latents match the response's signs; k pinned so
# this exercises the fixed-k path (the default binary chi hyperprior is
# tested below)
y.binary <- rbinom(n, 1L, pnorm(scale(f)))
sampler.binary <- dbarts(x, y.binary, control = control,
                         node.prior = normal(2))
bcSampler.binary <- dbarts:::bartcoreSampler(sampler.binary)
invisible(dbarts:::bartcoreRun(bcSampler.binary, 50L, 1L))

latents <- dbarts:::bartcoreGetLatents(bcSampler.binary)
expect_equal(length(latents), n)
expect_true(all(latents[y.binary == 1L] > 0) && all(latents[y.binary == 0L] < 0))

# fixed-k runs return no k samples
expect_null(result$k)

# the default binary spec uses the chi hyperprior on k
sampler.chik <- dbarts(x, y.binary, control = control)
bcSampler.chik <- dbarts:::bartcoreSampler(sampler.chik)
result.chik <- dbarts:::bartcoreRun(bcSampler.chik, 50L, 30L)
expect_equal(length(result.chik$k), 30L)
expect_true(all(result.chik$k > 0) && sd(result.chik$k) > 0)

# predictor mutation; on dedicated copies of x, since column and
# per-observation updates write into the borrowed matrix in place
x.mut <- x + 0
sampler.mut <- dbarts(x.mut, y, control = control)
bcSampler.mut <- dbarts:::bartcoreSampler(sampler.mut)
invisible(dbarts:::bartcoreRun(bcSampler.mut, 100L, 1L))

# an identity swap is always accepted; degenerate predictors would empty a
# leaf in some tree and roll back
expect_true(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.mut + 0))
x.degenerate <- matrix(0.5, n, p)
expect_false(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate))
expect_false(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate,
                                           updateCutPoints = TRUE))
result.mut <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 5L)
expect_true(all(is.finite(result.mut$yhat.train)))

# column-subset update: tiny jitter is accepted, degenerate column rejected
x.jitter <- x.mut[, 2L] + rnorm(n, 0, 1e-4)
expect_true(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, x.jitter, 2L))
expect_false(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, rep(0.5, n), 2L))
expect_error(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, rep(0.5, n),
                                              p + 1L))

# per-observation update: extreme values install except where an observation
# is the last occupant of a leaf
installed <- dbarts:::bartcoreUpdatePredictorPerObservation(bcSampler.mut,
                                                            rep(10, n), 1L)
expect_equal(length(installed), n)
expect_true(any(installed) && any(!installed))
result.perobs <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 5L)
expect_true(all(is.finite(result.perobs$yhat.train)))

# forced degenerate update collapses emptied splits instead of rolling back
expect_true(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate,
                                          forceUpdate = TRUE,
                                          updateCutPoints = TRUE))
result.forced <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 2L)
expect_true(all(is.finite(result.forced$yhat.train)))

# joint per-observation update: one mask, all-or-none across samplers that
# share an index-aligned column
x.jointA <- x + 0
x.jointB <- x + 0
x.jointB[, 2L:p] <- matrix(runif(n * (p - 1L)), n, p - 1L)
y.jointB <- -2 * x.jointB[, 1L] + rnorm(n)
bcSampler.jointA <- dbarts:::bartcoreSampler(dbarts(x.jointA, y,
                                                    control = control))
bcSampler.jointB <- dbarts:::bartcoreSampler(dbarts(x.jointB, y.jointB,
                                                    control = control))
invisible(dbarts:::bartcoreRun(bcSampler.jointA, 100L, 1L))
invisible(dbarts:::bartcoreRun(bcSampler.jointB, 100L, 1L))

installed.joint <- dbarts:::bartcoreUpdatePredictorPerObservationJointly(
  list(bcSampler.jointA, bcSampler.jointB), rep(10, n), c(1L, 1L))
expect_equal(length(installed.joint), n)
expect_true(any(installed.joint) && any(!installed.joint))
expect_true(all(is.finite(
  dbarts:::bartcoreRun(bcSampler.jointA, 0L, 1L)$yhat.train)))
expect_true(all(is.finite(
  dbarts:::bartcoreRun(bcSampler.jointB, 0L, 1L)$yhat.train)))
