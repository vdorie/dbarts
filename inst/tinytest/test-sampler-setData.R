source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler setData yields valid model when redefining cut points
n <- 105L

set.seed(0L)
x <- runif(n, -10, 10)
x.test <- seq.int(min(x), max(x), length.out = 101L)

beta.1 <- -0.75
beta.2 <-  0.5

n.cuts <- 10L
cutoffs <- min(x) + seq_len(n.cuts) * (max(x) - min(x)) / (n.cuts + 1)

y <- ifelse(x <= cutoffs[1L] | x > cutoffs[n.cuts], beta.1, beta.2) +
     rnorm(n, 0L, 0.15)

control <- dbarts::dbartsControl(
  n.trees = 1L, n.cuts = n.cuts,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE
)
sampler <- dbarts::dbarts(y ~ x, control = control)

samples1 <- sampler$run(500L, 1000L)

x.new <- x + diff(cutoffs[1L:2L])

sampler$setData(dbarts::dbartsData(y ~ x.new))
samples2 <- sampler$run(500L, 1000L)

expect_equal(
  sd(apply(samples1$train, 1, mean) - y),
  sd(apply(samples2$train, 1, mean) - y),
  tol = 1.0e-2
)

rm(
  samples2, x.new, samples1, sampler, control, y, cutoffs, n.cuts, beta.2,
  beta.1, x.test, x, n, testData
)
