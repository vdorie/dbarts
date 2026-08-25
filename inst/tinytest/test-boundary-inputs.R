# boundary inputs the suite otherwise never exercises: a single predictor,
# n = 2, a constant response, and a constant designated column under linear/gp
# leaves (which the docs claim degrade to an intercept at sd = 0). Every case
# either produces a finite fit or refuses at the R level with a clear error;
# none crash or return NaN.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

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
# a non-degenerate fit is unaffected by estimateSigmaFromLinearModel's floor
expect_equal(sampler.p1$data@sigma, summary(lm(y.p1 ~ x.p1))$sigma)
expect_equal(ncol(sampler.p1$data@x), 1L)
samples.p1 <- sampler.p1$run(50L, 50L)
expect_true(all(is.finite(samples.p1$train)))
expect_true(all(is.finite(samples.p1$sigma)))

# n = 2: the linear-fit sigma estimate is degenerate (residual df 0, so lm's
# sigma is exactly 0/0 = NaN, which is not an estimate at all).
# estimateSigmaFromLinearModel falls back to the marginal residual sd with a
# classed warning, exactly as the sparse-predictor branch above does, then
# applies the same relative-epsilon floor; sd(y.n2) is well above that floor.
x.n2 <- matrix(c(0.2, 0.8, 0.3, 0.7), 2L, 2L)
y.n2 <- c(1.0, -1.0)
warnings.n2 <- captureWarnings(
  sampler.n2 <- dbarts(y.n2 ~ x.n2, control = control)
)
expect_equal(length(warnings.n2), 1L)
expect_true(inherits(warnings.n2[[1L]], "dbartsSigmaFallbackWarning"))
expect_equal(sampler.n2$data@sigma, sd(y.n2), tolerance = 0)
samples.n2 <- sampler.n2$run(10L, 10L)
expect_true(all(is.finite(samples.n2$train)))
expect_true(all(is.finite(samples.n2$sigma)))

# a constant response: the fit runs, sigma is floored at a relative epsilon,
# and the train fits equal the constant. Finite throughout. A constant
# response is exactly the precision-degenerate case dbartsData now warns
# about (see test-data-precision-warning.R), so the warning is expected here
# too.
set.seed(1)
x.const <- matrix(runif(n * 2L), n, 2L)
y.const <- rep(3.0, n)
warnings.const <- captureWarnings(
  sampler.const <- dbarts(y.const ~ x.const, control = control)
)
# The starting-sigma linear fit's own "essentially perfect fit" warning is
# muffled at its source (see estimateSigmaFromLinearModel), so only the
# degenerate-response warning should reach here.
expect_equal(length(warnings.const), 1L)
expect_true(grepl("indistinguishable", conditionMessage(warnings.const[[1L]])))
# the constant response's own residual, not its weight-scaled counterpart,
# sets the floor: sqrt(eps) * max(1, |y|)
expect_equal(
  sampler.const$data@sigma,
  sqrt(.Machine$double.eps) * 3.0,
  tolerance = 0
)
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
