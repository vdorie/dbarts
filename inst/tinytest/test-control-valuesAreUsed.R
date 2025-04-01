source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that control argument works
n.samples <- 500L
n.trees <- 50L
n.cuts <- 50L
n.chains <- 2L
n.threads <- 2L
control <- dbarts::dbartsControl(
  n.samples = n.samples, verbose = FALSE, n.trees = n.trees, n.cuts = n.cuts,
  n.chains = n.chains, n.threads = n.threads
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)

expect_equal(sampler$control@n.trees, n.trees)
expect_equal(sampler$control@verbose, FALSE)
expect_equal(sampler$control@n.samples, n.samples)
expect_true(all(sampler$data@n.cuts == n.cuts))
expect_equal(sampler$control@n.chains, n.chains)
expect_equal(sampler$control@n.threads, n.threads)

rm(
  sampler, control,
  n.threads, n.chains, n.cuts, n.trees, n.samples
)


# test that keepevery behaves as it did in BayesTree"
n.burn <- 0L
n.sims <- 50L
keepevery <- 5L

bartFit <- bart(
  testData$x, testData$y,
  ndpost = n.sims, nskip = n.burn, ntree = 5L,
  verbose = FALSE, keepevery = keepevery
)

expect_equal(nrow(bartFit$yhat.train), n.sims %/% keepevery)
rm(bartFit, keepevery, n.sims, n.burn)

# test_that call slot of control is propagated
n.samples <- 5L
n.burn    <- 1L
n.trees   <- 5L

control <- dbarts::dbartsControl(
  n.samples = n.samples,
  verbose = FALSE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
expect_equal(sampler$control@call[[1L]], quote(dbarts::dbarts))

control@call <- call("NULL")
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
expect_equal(sampler$control@call, control@call)

bartFit <- dbarts::bart(y ~ x, testData, verbose = FALSE, ndpost = n.samples, nskip = n.burn, ntree = n.trees)
expect_equal(bartFit$call[[1L]], quote(dbarts::bart))

bartFit <- dbarts::bart(y ~ x, testData, verbose = FALSE, ndpost = n.samples, nskip = n.burn, ntree = n.trees, keepcall = FALSE)
expect_equal(bartFit$call, call("NULL"))

rm(bartFit, sampler, control, n.trees, n.burn, n.samples)

rm(testData)

