# boundary inputs the suite otherwise never exercises: a single predictor,
# n = 2, a constant response, and a constant designated column under linear/gp
# leaves (which the docs claim degrade to an intercept at sd = 0). Every case
# either produces a finite fit or refuses at the R level with a clear error;
# none crash or return NaN.

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE
)

# p = 1: a single-predictor gaussian fit runs and stays finite
set.seed(1)
n <- 100L
x.p1 <- matrix(runif(n), n, 1L)
y.p1 <- 2 * x.p1[, 1L] + rnorm(n, 0, 0.3)
sampler.p1 <- dbarts(y.p1 ~ x.p1, control = control)
expect_equal(ncol(sampler.p1$data@x), 1L)
samples.p1 <- sampler.p1$run(50L, 50L)
expect_true(all(is.finite(samples.p1$train)))
expect_true(all(is.finite(samples.p1$sigma)))

# n = 2: the linear-fit sigma estimate is degenerate (residual df 0), so the
# sampler refuses at creation with a clear R-level error rather than crashing
x.n2 <- matrix(c(0.2, 0.8, 0.3, 0.7), 2L, 2L)
y.n2 <- c(1.0, -1.0)
expect_error(dbarts(y.n2 ~ x.n2, control = control), pattern = "sigma estimate")

# a constant response: the fit runs, sigma collapses to ~0, and the train
# fits equal the constant. Finite throughout.
set.seed(1)
x.const <- matrix(runif(n * 2L), n, 2L)
y.const <- rep(3.0, n)
sampler.const <- dbarts(y.const ~ x.const, control = control)
samples.const <- sampler.const$run(30L, 30L)
expect_true(all(is.finite(samples.const$train)))
expect_true(all(is.finite(samples.const$sigma)))
expect_true(all(abs(samples.const$train - 3.0) < 1e-6))

# a constant designated column under linear and gp leaves: an sd-0 column is
# kept at sd 1 and the leaf degrades to an intercept. Finite fits, no NaN.
set.seed(1)
x1.c <- rep(0.5, n) # constant designated column
x2.c <- runif(n)
y.c <- 2 * x2.c + rnorm(n, 0, 0.3)
df.c <- data.frame(x1 = x1.c, x2 = x2.c, y = y.c)

sampler.lin <- dbarts(
  y ~ x1 + x2,
  df.c,
  node.prior = linear("x1"),
  control = control
)
samples.lin <- sampler.lin$run(40L, 40L)
expect_true(all(is.finite(samples.lin$train)))
expect_true(all(is.finite(samples.lin$sigma)))

sampler.gp <- dbarts(
  y ~ x1 + x2,
  df.c,
  node.prior = gp("x1"),
  control = control
)
samples.gp <- sampler.gp$run(40L, 40L)
expect_true(all(is.finite(samples.gp$train)))
expect_true(all(is.finite(samples.gp$sigma)))

rm(
  control,
  n,
  x.p1,
  y.p1,
  sampler.p1,
  samples.p1,
  x.n2,
  y.n2,
  x.const,
  y.const,
  sampler.const,
  samples.const,
  x1.c,
  x2.c,
  y.c,
  df.c,
  sampler.lin,
  samples.lin,
  sampler.gp,
  samples.gp
)
