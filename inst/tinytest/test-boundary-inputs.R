# boundary inputs the suite otherwise never exercises: a single predictor,
# n = 2, a constant response, and a constant designated column under linear/gp
# leaves (refused, but under a gp leaf whose lengthscale is stated). Every
# case either produces a finite fit or refuses with a clear error; none crash
# or return NaN.

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
# estimateSigmaFromLinearModel takes the marginal residual sd with no
# warning, saying so under verbose, then applies the same relative-epsilon
# floor; sd(y.n2) is well above that floor.
x.n2 <- matrix(c(0.2, 0.8, 0.3, 0.7), 2L, 2L)
y.n2 <- c(1.0, -1.0)
warnings.n2 <- captureWarnings(
  sampler.n2 <- dbarts(y.n2 ~ x.n2, control = control)
)
expect_equal(length(warnings.n2), 0L)
expect_equal(sampler.n2$data@sigma, sd(y.n2), tolerance = 0)
verboseControl <- control
verboseControl@verbose <- TRUE
expect_message(
  invisible(capture.output(dbarts(y.n2 ~ x.n2, control = verboseControl))),
  "starting sigma is the sd of the response (1.414): the linear model on 2",
  fixed = TRUE
)
samples.n2 <- sampler.n2$run(10L, 10L)
expect_true(all(is.finite(samples.n2$train)))
expect_true(all(is.finite(samples.n2$sigma)))

# n = 1: the residual has no degrees of freedom for lm either, and its
# marginal-fallback sd() is itself NA (undefined, not merely small); the
# finiteness re-check after the fallback, not just before it, is what
# floors this instead of erroring on the NA comparison.
x.n1 <- matrix(c(0.2, 0.3), 1L, 2L)
y.n1 <- c(1.0)
warnings.n1 <- captureWarnings(
  sampler.n1 <- dbarts(y.n1 ~ x.n1, control = control)
)
expect_equal(length(warnings.n1), 1L)
expect_true(grepl("indistinguishable", conditionMessage(warnings.n1[[1L]])))
expect_equal(sampler.n1$data@sigma, sqrt(.Machine$double.eps), tolerance = 0)
samples.n1 <- sampler.n1$run(10L, 10L)
expect_true(all(is.finite(samples.n1$train)))
expect_true(all(is.finite(samples.n1$sigma)))

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
# The starting-sigma least squares raises no warning of its own - it takes
# summary.lm's expression for sigma without going through summary.lm, so that
# function's "essentially perfect fit" is not on the path - and only the
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
# the transform's window is centered on the constant, c - 0.5 to c + 0.5, so the
# leaf prior is centered on it too; it was c to c + 1, half a unit above
prior.const <- sampler.const$getLeafPrior()
expect_identical(prior.const$prior.mean, 3.0)
expect_identical(prior.const$response.shift, 3.0)
samples.const <- sampler.const$run(30L, 30L)
expect_true(all(is.finite(samples.const$train)))
expect_true(all(is.finite(samples.const$sigma)))
expect_true(all(abs(samples.const$train - 3.0) < 1e-6))

# a constant designated column under linear and gp leaves: one value leaves
# no center and scale to take, so the creation is refused under linear and
# under gp; a gp leaf whose lengthscale is stated takes the placeholder, the
# column at its value divided by 1, and runs finite.
set.seed(1)
x1.c <- rep(0.5, n) # constant designated column
x2.c <- runif(n)
y.c <- 2 * x2.c + rnorm(n, 0, 0.3)
df.c <- data.frame(x1 = x1.c, x2 = x2.c, y = y.c)

expect_error(
  dbarts(y ~ x1 + x2, df.c, leaf.prior = linear("x1"), control = control),
  paste0(
    'the leaf covariate "x1" holds one value over the rows in the ',
    "likelihood \\(100 of 100 observations are in it\\), so its center ",
    "and scale cannot be taken from them; drop it from linear\\(\\)"
  )
)
expect_error(
  dbarts(y ~ x1 + x2, df.c, leaf.prior = gp("x1"), control = control),
  "state its lengthscale in gp\\(\\), drop it from gp\\(\\)"
)

sampler.gp <- dbarts(
  y ~ x1 + x2,
  df.c,
  leaf.prior = gp("x1", lengthscale = 1),
  control = control
)
expect_identical(sampler.gp$data@standardization[[1L]]$center[["x1"]], 0.5)
expect_true(is.na(sampler.gp$data@standardization[[1L]]$scale[["x1"]]))
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
  sampler.gp,
  samples.gp
)
