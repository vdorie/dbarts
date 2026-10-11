# Weighted logistic responses: weights are observation counts (positive
# integers), so a weighted fit reproduces a fit on the physically replicated
# data up to Monte Carlo error, and the count semantics are enforced.
# Probit weights are refused (intractable); gaussian weights are unrestricted.

# pin the sampler kind: an earlier suite file leaks sample.kind = "Rounding",
# which would shift the drawn data and the Monte Carlo comparison below
set.seed(2, sample.kind = "Rejection")
# Binomial counts: at each predictor point 1 to 3 independent trials, entered
# as a row of successes and a row of failures, a count of 0 left out. A count
# states that many independent trials sharing the row's outcome, so a row's
# count must not repeat one trial.
n.points <- 160L
x.points <- matrix(runif(n.points * 2L), n.points)
p.points <- plogis(2.5 * x.points[, 1L] - 1.2)
trials <- sample(1:3, n.points, replace = TRUE)
successes <- rbinom(n.points, trials, p.points)
counts <- c(successes, trials - successes)
rows <- which(counts > 0L)
point <- rep(seq_len(n.points), 2L)[rows]
x <- x.points[point, , drop = FALSE]
p <- p.points[point]
y <- rep(c(1L, 0L), each = n.points)[rows]
w <- counts[rows]
n <- length(rows)
# both outcomes at one point occur, and every count from 1 to 3
expect_true(anyDuplicated(point) > 0L)
expect_equal(sort(unique(w)), 1:3)
expect_equal(sum(w), sum(trials))

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
expect_true(cor(phat.w, p) > 0.8)

# the count semantics: replicate each row w times, fit unweighted, and the
# per-row posterior mean probabilities agree within Monte Carlo error. The
# bounds must stay below what a fit that ignores the counts gives against the
# replicated one: a mean absolute difference near 0.04 at a correlation near
# 0.96, against 0.01 and 0.996 for the weighted fit.
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
expect_true(cor(phat.w, phat.rep) > 0.99)
expect_true(mean(abs(phat.w - phat.rep)) < 0.02)

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
expect_error(
  predict(fit.w, newRows, type = "ppd", weights = c(Inf, 2, 3)),
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
  p,
  y,
  w,
  n.points,
  x.points,
  p.points,
  trials,
  successes,
  counts,
  rows,
  point,
  idx,
  yr,
  xr,
  n,
  control
)
