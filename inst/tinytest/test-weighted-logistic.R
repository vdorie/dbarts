# Weighted logistic responses: weights are observation counts (positive
# integers), so a weighted fit reproduces a fit on the physically replicated
# data up to Monte Carlo error, and the count semantics are enforced.
# Probit weights are refused (intractable); gaussian weights are unrestricted.

# pin the sampler kind: an earlier suite file leaks sample.kind = "Rounding",
# which would shift the drawn data and the Monte Carlo comparison below
set.seed(2, sample.kind = "Rejection")
n <- 160L
x <- matrix(runif(n * 2L), n)
f <- 2.5 * x[, 1L] - 1.2
p <- plogis(f)
y <- rbinom(n, 1L, p)
w <- sample(1:3, n, replace = TRUE)

# an integer-weighted logistic model fits and recovers the signal
fit.w <- bart(
  y ~ x,
  weights = w,
  family = "logistic",
  n.samples = 200L,
  n.burn = 100L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(fit.w$family, "logistic")
phat.w <- fitted(fit.w)
expect_true(all(is.finite(phat.w) & phat.w >= 0 & phat.w <= 1))
expect_true(cor(phat.w, p) > 0.5)

# the count semantics: replicate each row w times, fit unweighted, and the
# per-point posterior mean probabilities agree within Monte Carlo error
idx <- rep(seq_len(n), times = w)
yr <- y[idx]
xr <- x[idx, , drop = FALSE]
fit.rep <- bart(
  yr ~ xr,
  family = "logistic",
  n.samples = 200L,
  n.burn = 100L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
phat.rep <- fitted(fit.rep)[!duplicated(idx)]
expect_true(cor(phat.w, phat.rep) > 0.9)
expect_true(mean(abs(phat.w - phat.rep)) < 0.06)

# non-integer weights are not counts and are refused
expect_error(
  bart(
    y ~ x,
    weights = runif(n, 0.5, 2),
    family = "logistic",
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "must be positive integers"
)

# a zero count is a dropped row, also refused
expect_error(
  suppressWarnings(bart(
    y ~ x,
    weights = c(0, w[-1L]),
    family = "logistic",
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  pattern = "must be positive integers"
)

# probit (the default binary family) refuses weights entirely
expect_error(
  bart(
    y ~ x,
    weights = w,
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "probit models do not support weights"
)

# predict's posterior predictive draws at new rows hold their weights to the
# counts the fit takes, refusing a fractional one by name rather than drawing
# NA; integer counts draw that many trials
newRows <- data.frame(x = I(x[1:3, , drop = FALSE]))
expect_error(
  predict(fit.w, newRows, type = "ppd", weights = c(1.5, 2, 3)),
  "the posterior predictive 'weights' of a logistic fit are trial counts and must be positive integers",
  fixed = TRUE
)
expect_error(
  predict(fit.w, newRows, type = "ppd", weights = c(0, 2, 3)),
  "the posterior predictive 'weights' of a logistic fit are trial counts and must be positive integers",
  fixed = TRUE
)
ppdCounts <- predict(fit.w, newRows, type = "ppd", weights = c(1, 2, 3))
expect_true(all(
  ppdCounts >= 0 & ppdCounts <= rep(c(1, 2, 3), each = nrow(ppdCounts))
))
# a probit fit takes no weights of its own, but its draws at new rows take
# trial counts under the same policy
fit.p <- bart(
  y ~ x,
  family = "probit",
  n.samples = 10L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  predict(fit.p, newRows, type = "ppd", weights = c(1.5, 2, 3)),
  "the posterior predictive 'weights' of a probit fit are trial counts and must be positive integers",
  fixed = TRUE
)
ppdProbit <- predict(fit.p, newRows, type = "ppd", weights = c(1, 2, 3))
expect_true(all(
  ppdProbit >= 0 & ppdProbit <= rep(c(1, 2, 3), each = nrow(ppdProbit))
))
rm(newRows, ppdCounts, fit.p, ppdProbit)

# the dbarts sampler surface accepts integer-count logistic weights directly
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE
)
sampler <- dbarts(y ~ x, weights = w, family = "logistic", control = control)
samples <- sampler$run(50L, 50L)
expect_equal(length(unique(samples$sigma)), 1L) # logistic sd is fixed, not drawn
expect_true(all(is.finite(samples$train)))

rm(
  fit.w,
  fit.rep,
  phat.w,
  phat.rep,
  sampler,
  samples,
  x,
  f,
  p,
  y,
  w,
  idx,
  yr,
  xr,
  n,
  control
)
