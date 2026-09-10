# The starting sigma estimate is summary(lm(y ~ x))$sigma to the last bit, on
# every shape the estimate is taken over. dbarts reaches it through the QR
# routine lm itself calls rather than through lm, so that no model frame and
# no second copy of the design exist at the largest n; nothing about the
# NUMBER may move with that, and expect_identical (not expect_equal) is the
# whole point of this file.

set.seed(11)
n <- 200L
f.ss <- factor(sample(c("a", "b", "c"), n, TRUE))
df.ss <- data.frame(x1 = rnorm(n), x2 = runif(n), f = f.ss)
y.ss <- 2 * df.ss$x1 - df.ss$x2 + as.numeric(f.ss) + rnorm(n)
w.ss <- runif(n, 0.5, 2)
o.ss <- rnorm(n)

# the design the estimate actually sees: dbarts expands the factor itself, so
# read the expanded matrix back rather than re-deriving lm's own coding
control.ss <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
sampler.ss <- dbarts(y.ss ~ ., df.ss, control = control.ss)
x.ss <- dbarts::extract(sampler.ss, "predictors")

expect_identical(sampler.ss$data@sigma, summary(lm(y.ss ~ x.ss))$sigma)

sampler.ssw <- dbarts(y.ss ~ ., df.ss, weights = w.ss, control = control.ss)
expect_identical(
  sampler.ssw$data@sigma,
  summary(lm(y.ss ~ x.ss, weights = w.ss))$sigma
)

sampler.sso <- dbarts(y.ss ~ ., df.ss, offset = o.ss, control = control.ss)
expect_identical(
  sampler.sso$data@sigma,
  summary(lm(y.ss ~ x.ss, offset = o.ss))$sigma
)

# the estimator itself, over the argument crossing the samplers above do not
# reach in one call
expect_identical(
  dbarts:::residualStandardError(y.ss, x.ss, w.ss, o.ss),
  summary(lm(y.ss ~ x.ss, weights = w.ss, offset = o.ss))$sigma
)

# a rank-deficient design: the QR pivots, and the pivot order follows the
# intercept-first column order model.matrix would have produced
x.dup <- cbind(x.ss, x.ss[, 1L])
expect_identical(
  dbarts:::residualStandardError(y.ss, x.dup, w.ss, NULL),
  summary(lm(y.ss ~ x.dup, weights = w.ss))$sigma
)

# a missing response: lm's default na.action drops the row before fitting,
# and so does this
y.na <- y.ss
y.na[c(3L, 17L, 198L)] <- NA_real_
expect_identical(
  dbarts:::residualStandardError(y.na, x.ss, NULL, NULL),
  summary(lm(y.na ~ x.ss))$sigma
)
expect_identical(
  dbarts:::residualStandardError(y.na, x.ss, w.ss, o.ss),
  summary(lm(y.na ~ x.ss, weights = w.ss, offset = o.ss))$sigma
)

# a response the design reproduces exactly, where the residual sum of squares
# is rounding noise and the summation order is all there is
expect_identical(
  dbarts:::residualStandardError(x.ss[, 1L], x.ss, NULL, NULL),
  suppressWarnings(summary(lm(x.ss[, 1L] ~ x.ss))$sigma)
)
