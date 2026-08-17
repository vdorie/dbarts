# bart() carries $weights/$weights.test on its returned fit, matching
# bart2's tail: a weighted fit's loglik/ppd otherwise silently computes
# unweighted quantities. Oracle follows test-pointwise-loglik.R's own
# gaussian-weighted form: residual sd for observation i is sigma_s / sqrt(w_i).
set.seed(6, sample.kind = "Rejection")
n <- 80L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)
w <- rep_len(c(1, 4), n)

fit <- bart(
  x,
  y,
  weights = w,
  ndpost = 60L,
  nskip = 40L,
  ntree = 20L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)
expect_identical(fit$weights, w)
expect_null(fit$weights.test)

ll <- extract(fit, type = "loglik")
ev <- extract(fit, type = "ev")
j4 <- which(w == 4)[1L]
expect_identical(
  ll[, j4],
  dnorm(y[j4], ev[, j4], fit$sigma / sqrt(4), log = TRUE)
)
j1 <- which(w == 1)[1L]
expect_identical(ll[, j1], dnorm(y[j1], ev[, j1], fit$sigma, log = TRUE))

# reachable cross-check: bart2's tail has always attached $weights this way
fit2 <- bart2(
  y ~ x,
  weights = w,
  n.samples = 60L,
  n.burn = 40L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(fit$weights, fit2$weights)

rm(x, y, w, fit, fit2, ll, ev, j4, j1, n)

# a draw count that thins to zero refuses by name instead of faulting deep
# in the result-array reshape (dim(X) has no positive length)
expect_error(
  bart(
    matrix(rnorm(20L), 10L, 2L),
    rnorm(10L),
    ndpost = 5L,
    keepevery = 10L,
    nskip = 2L,
    ntree = 5L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE
  ),
  pattern = "ndpost"
)
