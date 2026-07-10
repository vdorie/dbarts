# extract(type = "loglik") returns the per-draw, per-observation
# log-likelihood of the training response: for draw s and observation i,
# gaussian fits evaluate dnorm(y_i, ev_si, sigma_s / sqrt(w_i), log = TRUE)
# and binary fits w_i * dbinom(y_i, 1, p_si, log = TRUE). Checks below are
# bit-agreement with those forms on the stored draws, not statistical. See
# docs/plans/pointwise-loglik.md.

# 1. gaussian, unweighted, single chain: entries match dnorm on the stored
# draws; the S x n shape matches "ev" and feeds loo/WAIC directly
set.seed(2, sample.kind = "Rejection")
n <- 100L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

fit <- bart2(
  y ~ x,
  n.samples = 80L,
  n.burn = 40L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
ll <- extract(fit, type = "loglik")
ev <- extract(fit, type = "ev")

expect_equal(dim(ll), dim(ev))
for (i in c(1L, 37L, 100L)) {
  expect_identical(ll[, i], dnorm(y[i], ev[, i], fit$sigma, log = TRUE))
}

# 7. the type is extract-only: predict and fitted reject it with their
# standard vocabulary error, and there is no test response to evaluate
expect_error(predict(fit, x, type = "loglik"), pattern = "type must be in")
expect_error(fitted(fit, type = "loglik"), pattern = "type must be in")
expect_error(
  extract(fit, type = "loglik", sample = "test"),
  pattern = "no test response"
)

rm(ll, ev, fit, i)

# 2. gaussian, weighted: the residual sd for observation i is
# sigma_s / sqrt(w_i)
w <- rep_len(c(1, 4), n)
fit.w <- bart2(
  y ~ x,
  weights = w,
  n.samples = 60L,
  n.burn = 40L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ll.w <- extract(fit.w, type = "loglik")
ev.w <- extract(fit.w, type = "ev")

j4 <- which(w == 4)[1L]
expect_identical(
  ll.w[, j4],
  dnorm(y[j4], ev.w[, j4], fit.w$sigma / sqrt(w[j4]), log = TRUE)
)
j1 <- which(w == 1)[1L]
expect_identical(
  ll.w[, j1],
  dnorm(y[j1], ev.w[, j1], fit.w$sigma, log = TRUE)
)

rm(ll.w, ev.w, fit.w, j4, j1, w)

# 3. probit: matches the bernoulli log-likelihood of the fitted probability
set.seed(3, sample.kind = "Rejection")
y.b <- rbinom(n, 1L, pnorm(0.8 * x[, 1L] - 0.4))
fit.p <- bart2(
  y.b ~ x,
  n.samples = 60L,
  n.burn = 40L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ll.p <- extract(fit.p, type = "loglik")
ev.p <- extract(fit.p, type = "ev")

expect_equal(dim(ll.p), dim(ev.p))
for (i in c(2L, 50L)) {
  expect_identical(ll.p[, i], dbinom(y.b[i], 1L, ev.p[, i], log = TRUE))
}

rm(ll.p, ev.p, fit.p, i)

# 4. weighted logistic: integer weights are observation counts, so the
# log-likelihood is w_i times the unweighted bernoulli form on its ev
w.l <- rep_len(c(1L, 3L, 5L), n)
fit.l <- bart2(
  y.b ~ x,
  weights = w.l,
  family = "logistic",
  n.samples = 60L,
  n.burn = 40L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ll.l <- extract(fit.l, type = "loglik")
ev.l <- extract(fit.l, type = "ev")

j5 <- which(w.l == 5L)[1L]
expect_identical(
  ll.l[, j5],
  w.l[j5] * dbinom(y.b[j5], 1L, ev.l[, j5], log = TRUE)
)
j1 <- which(w.l == 1L)[1L]
expect_identical(ll.l[, j1], dbinom(y.b[j1], 1L, ev.l[, j1], log = TRUE))

rm(ll.l, ev.l, fit.l, j5, j1, w.l, y.b)

# 5. rbart_vi: the location conditions on the drawn group intercepts
set.seed(5, sample.kind = "Rejection")
g <- factor(rep_len(1:5, n))
y.r <- y + rnorm(5L, 0, 1.5)[as.integer(g)]
fit.r <- rbart_vi(
  y.r ~ x,
  group.by = g,
  n.samples = 40L,
  n.burn = 20L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)
ll.r <- extract(fit.r, type = "loglik")
ev.r <- extract(fit.r, type = "ev")
bart.r <- extract(fit.r, type = "bart")
ranef.r <- extract(fit.r, type = "ranef")

expect_equal(dim(ll.r), dim(ev.r))
i <- 3L
expect_identical(
  ll.r[, i],
  dnorm(
    y.r[i],
    bart.r[, i] + ranef.r[, as.character(g[i])],
    fit.r$sigma,
    log = TRUE
  )
)
expect_error(
  extract(fit.r, type = "loglik", sample = "test"),
  pattern = "no test response"
)

rm(ll.r, ev.r, bart.r, ranef.r, fit.r, y.r, g, i)

# 6. multi-chain shape conventions follow "ev" in both chain layouts, each
# chain pairs with its own sigma draws, and the stored layout does not
# matter (combineChains at fit time only changes packaging)
set.seed(6, sample.kind = "Rejection")
fit.mc <- bart2(
  y ~ x,
  n.samples = 30L,
  n.burn = 20L,
  n.trees = 20L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE,
  seed = 77L
)
ll.unc <- extract(fit.mc, type = "loglik", combineChains = FALSE)
ev.unc <- extract(fit.mc, type = "ev", combineChains = FALSE)
expect_equal(dim(ll.unc), dim(ev.unc))
for (chain in 1:2) {
  expect_identical(
    ll.unc[chain, , 5L],
    dnorm(y[5L], ev.unc[chain, , 5L], fit.mc$sigma[chain, ], log = TRUE)
  )
}

ll.comb <- extract(fit.mc, type = "loglik")
ev.comb <- extract(fit.mc, type = "ev")
expect_equal(dim(ll.comb), dim(ev.comb))
# combined rows are chain blocks: row s is chain 1's sample s, row
# n.samples + s is chain 2's
expect_identical(ll.comb[7L, 5L], ll.unc[1L, 7L, 5L])
expect_identical(ll.comb[37L, 5L], ll.unc[2L, 7L, 5L])

fit.mc2 <- bart2(
  y ~ x,
  n.samples = 30L,
  n.burn = 20L,
  n.trees = 20L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = TRUE,
  seed = 77L
)
expect_identical(extract(fit.mc2, type = "loglik"), ll.comb)

rm(ll.unc, ev.unc, ll.comb, ev.comb, fit.mc, fit.mc2, chain)
rm(x, y, n)
