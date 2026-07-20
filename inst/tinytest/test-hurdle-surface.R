# The public bart2(family = "hurdle.lognormal") surface (docs/design/
# hurdle.md sections 0-1, 13; docs/plans/hurdle.md C3): family-token routing
# (both spellings), the y >= 0 / require-a-zero / require-a-positive
# validation errors, and the family-vector refusals shared with every
# composed family (dbarts() cannot express two samplers; xbart and rbart_vi
# omit the token). The analytic combine/retransform oracle, predict-on-
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
# canonical token (docs/design/hurdle.md section 13, NAMING) ---
fit <- do.call(bart2, c(list(x, y, family = "hurdle.lognormal"), fitArgs))
expect_inherits(fit, "bartHurdle")
expect_false(inherits(fit, "bart"))
expect_equal(fit$family, "hurdle.lognormal")

fitAlias <- do.call(bart2, c(list(x, y, family = "twopart"), fitArgs))
expect_inherits(fitAlias, "bartHurdle")
expect_equal(fitAlias$family, "hurdle.lognormal")
# the alias resolves to the canonical token before anything RNG-affecting
# runs, so a "twopart" fit at the same seed reproduces the
# "hurdle.lognormal" fit's two component fits bit for bit
expect_identical(fitAlias$occupancy$yhat.train, fit$occupancy$yhat.train)
expect_identical(fitAlias$positive$yhat.train, fit$positive$yhat.train)

printed <- capture.output(print(fit))
expect_true(any(grepl("hurdle.lognormal", printed, fixed = TRUE)))
printedAlias <- capture.output(print(fitAlias))
expect_true(any(grepl("hurdle.lognormal", printedAlias, fixed = TRUE)))

# --- y >= 0 validation and the require-a-zero / require-a-positive edges
# (docs/design/hurdle.md sections 0-1, splitHurdleResponse) ---
expect_error(
  bart2(x, y - 10, family = "hurdle.lognormal"),
  "non-negative"
)
expect_error(
  bart2(x, c(NA_real_, y[-1L]), family = "hurdle.lognormal"),
  "non-negative"
)
expect_error(
  bart2(x, rep(0, n), family = "hurdle.lognormal"),
  "at least one positive"
)
expect_error(
  bart2(x, abs(rnorm(n)) + 0.1, family = "hurdle.lognormal"),
  "at least one exact zero"
)

# --- dbarts() cannot express the two-sampler composition; directs to
# bart2() (both spellings hit the same redirect) ---
expect_error(dbarts(x, y, family = "hurdle.lognormal"), "bart2")
expect_error(dbarts(x, y, family = "twopart"), "bart2")

# --- xbart and rbart_vi do not fit it (their family vectors are the
# refusal, the nbinom/hazard precedent) ---
expect_error(
  xbart(x, y, family = "hurdle.lognormal", n.samples = 10L, n.reps = 1L),
  "should be one of"
)
expect_error(
  rbart_vi(
    y ~ x,
    family = "hurdle.lognormal",
    group.by = rep(1:4, length.out = n)
  ),
  "should be one of"
)
