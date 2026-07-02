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

# binary response: probit latents match the response's signs; fixed k, as
# bartcore does not yet support the chi hyperprior dbarts defaults to for
# binary responses
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
