# The public bart(family = "hurdle.lognormal") surface: family-token
# routing (both spellings), the y >= 0 / require-a-zero / require-a-positive
# validation errors, and the family-vector refusals shared with every
# composed family (dbarts() cannot express two samplers; xbart omits the
# token). The analytic combine/retransform oracle, predict-on-
# newdata, save/load, and the recovery smoke live in test-hurdle.R.

set.seed(4401L)
n <- 60L
p <- 2L
x <- matrix(runif(n * p), n, p)
pi.true <- pnorm(0.8 * x[, 1L] - 0.4)
mu.true <- 0.5 + 0.5 * x[, 2L]
occupied <- rbinom(n, 1L, pi.true) == 1L
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), mu.true[occupied], 0.4))

fitArgs <- list(
  n.samples = 8L,
  n.burn = 6L,
  n.trees = 6L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 91L
)

# --- family routing: both spellings reach a bartHurdle, printing/reading the
# canonical token ---
fit <- do.call(bart, c(list(x, y, family = "hurdle.lognormal"), fitArgs))
expect_inherits(fit, "bartHurdle")
expect_false(inherits(fit, "bart"))
expect_equal(fit$family, "hurdle.lognormal")

# "twopart" is a retired spelling of the same model: refused by name at both
# doors rather than folded, so the package carries one token per model
expect_error(
  do.call(bart, c(list(x, y, family = "twopart"), fitArgs)),
  "hurdle.lognormal"
)

printed <- capture.output(print(fit))
expect_true(any(grepl("hurdle.lognormal", printed, fixed = TRUE)))

# --- y >= 0 validation and the require-a-zero / require-a-positive edges
# (splitHurdleResponse) ---
expect_error(
  bart(x, y - 10, family = "hurdle.lognormal"),
  "non-negative"
)
expect_error(
  bart(x, c(NA_real_, y[-1L]), family = "hurdle.lognormal"),
  "non-negative"
)
expect_error(
  bart(x, rep(0, n), family = "hurdle.lognormal"),
  "at least one positive"
)
expect_error(
  bart(x, abs(rnorm(n)) + 0.1, family = "hurdle.lognormal"),
  "at least one exact zero"
)

# --- dbarts() cannot express the two-sampler composition; directs to the
# front door ---
expect_error(dbarts(x, y, family = "hurdle.lognormal"), "bart\\(x.train")
expect_error(dbarts(x, y, family = "twopart"), "hurdle.lognormal")

# --- xbart does not fit it (its family vector is the refusal, the
# nbinom/hazard precedent) ---
expect_error(
  xbart(x, y, family = "hurdle.lognormal", n.samples = 10L, n.reps = 1L),
  "should be one of"
)
