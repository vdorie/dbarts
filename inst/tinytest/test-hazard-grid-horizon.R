# survivalProbabilities on a discrete-time hazard fit at horizons on and off
# the period grid, against the product of (1 - hazard) draws predict() gives
# for each period: a horizon on a grid point includes its own period, one
# between two grid points the last period it passed, one before the first
# period none

set.seed(5101L)
n <- 60L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
time <- as.double(sample.int(3L, n, replace = TRUE))
status <- as.double(runif(n) < 0.7)
fit <- bart(
  x,
  cbind(time, status),
  family = "hazard",
  n.trees = 10L,
  n.burn = 20L,
  n.samples = 25L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 3L
)
expect_equal(fit$periods, c(1, 2, 3))

xn <- x[1:3, ]
expected <- array(1, c(25L, 4L, 3L))
survival <- 1
for (k in 1:3) {
  survival <- survival * (1 - predict(fit, cbind(xn, period = k), type = "ev"))
  expected[, k + 1L, ] <- survival
}

# on the grid
expect_equal(
  survivalProbabilities(fit, c(1, 2, 3), newdata = xn),
  expected[, 2:4, ]
)
# before the first period, and between grid points
expect_equal(
  survivalProbabilities(fit, c(0.5, 1.5, 2.5, 4), newdata = xn),
  expected[, c(1L, 2L, 3L, 4L), ]
)

rm(fit, x, xn, time, status, n, expected, survival, k)
