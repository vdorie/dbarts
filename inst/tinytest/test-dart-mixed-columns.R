# DART meeting the pooled categorical mask and NA (incorporate) columns: the
# varprobs indexing across those column mappings is otherwise unguarded. Assert
# the fit runs and the per-sample variable probabilities are a proper
# distribution over the model-matrix columns, whose count matches the design.

set.seed(3)
n <- 250L
x1 <- runif(n)
x2 <- runif(n, -1, 1)
g <- factor(sample(letters[1:4], n, replace = TRUE))
x2[sample(n, 20L)] <- NA # NA predictor values, incorporated
mu <- 2 * x1 + ifelse(g == "a", 1, -0.5)
y <- mu + rnorm(n, 0, 0.3)
df <- data.frame(x1, x2, g, y)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE
)

# factors = "categorical": one pooled column per predictor, so the factor g
# occupies a single model-matrix column alongside x1 and the NA-carrying x2
sampler.cat <- dbarts(
  y ~ x1 + x2 + g,
  df,
  tree.prior = dart(),
  factors = "categorical",
  missing = "incorporate",
  control = control
)
expect_equal(ncol(sampler.cat$data@x), 3L)
expect_true(anyNA(sampler.cat$data@x)) # NA is modeled, not dropped

samples.cat <- sampler.cat$run(100L, 100L)
expect_true(all(is.finite(samples.cat$train)))
expect_equal(dim(samples.cat$varprobs), c(3L, 100L))
expect_true(all(is.finite(samples.cat$varprobs)))
expect_true(all(samples.cat$varprobs >= 0))
expect_true(all(abs(colSums(samples.cat$varprobs) - 1) < 1e-8))

# factors = "indicators": the factor expands to one column per level, and the
# varprobs still form a distribution over the widened design
sampler.ind <- dbarts(
  y ~ x1 + x2 + g,
  df,
  tree.prior = dart(),
  factors = "indicators",
  missing = "incorporate",
  control = control
)
expect_equal(ncol(sampler.ind$data@x), 6L) # x1, x2, g.a..g.d

samples.ind <- sampler.ind$run(100L, 100L)
expect_true(all(is.finite(samples.ind$train)))
expect_equal(dim(samples.ind$varprobs), c(6L, 100L))
expect_true(all(is.finite(samples.ind$varprobs)))
expect_true(all(samples.ind$varprobs >= 0))
expect_true(all(abs(colSums(samples.ind$varprobs) - 1) < 1e-8))

rm(
  sampler.cat,
  sampler.ind,
  samples.cat,
  samples.ind,
  control,
  df,
  x1,
  x2,
  g,
  mu,
  y,
  n
)
