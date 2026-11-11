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

# the design the estimate sees: the factor as its indicator columns, which
# the regression takes however the trees are given the factor
control.ss <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
sampler.ss <- dbarts(y.ss ~ ., df.ss, control = control.ss, n.trees = 5L)
x.ss <- dbarts::extract(sampler.ss, "predictors")
x.ind <- makeModelMatrixFromDataFrame(df.ss)

expect_identical(sampler.ss$data@sigma, summary(lm(y.ss ~ x.ind))$sigma)

sampler.ssw <- dbarts(
  y.ss ~ .,
  df.ss,
  weights = w.ss,
  control = control.ss,
  n.trees = 5L
)
expect_identical(
  sampler.ssw$data@sigma,
  summary(lm(y.ss ~ x.ind, weights = w.ss))$sigma
)

sampler.sso <- dbarts(
  y.ss ~ .,
  df.ss,
  offset = o.ss,
  control = control.ss,
  n.trees = 5L
)
expect_identical(
  sampler.sso$data@sigma,
  summary(lm(y.ss ~ x.ind, offset = o.ss))$sigma
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

# --- the forests' bases enter the regression ---

# Each forest's basis is signal the model holds, so the regression takes it
# beside the predictors: a forest of numbers as its column and a forest of
# levels as the indicator of its second level, the one column a two-level
# factor predictor enters as. A BCF fit's estimate is then lm's on x and z.
set.seed(17)
n.b <- 100L
d.b <- data.frame(x1 = runif(n.b), x2 = runif(n.b))
d.b$z <- rbinom(n.b, 1L, 0.5)
d.b$w <- rnorm(n.b)
d.b$y <- sin(3 * d.b$x1) +
  4 * d.b$z +
  0.5 * d.b$w * d.b$x2 +
  rnorm(n.b, sd = 0.5)
w.b <- runif(n.b, 0.5, 2)
o.b <- rnorm(n.b, sd = 0.2)
forest <- dbartsForests$forest
f.levels <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = factor(z))
f.numbers <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = w)
sigmaOf <- function(...) {
  do.call(dbarts, list(..., control = control.ss))$data@sigma
}
lm.levels <- function(...) {
  summary(lm(y ~ x1 + x2 + factor(z), d.b, ...))$sigma
}
lm.numbers <- function(...) summary(lm(y ~ x1 + x2 + w, d.b, ...))$sigma

# the dense route is lm's own routine on lm's own columns: to the last bit
expect_identical(sigmaOf(f.levels, d.b), lm.levels())
expect_identical(
  sigmaOf(f.levels, d.b, weights = w.b),
  lm.levels(weights = w.b)
)
expect_identical(sigmaOf(f.levels, d.b, offset = o.b), lm.levels(offset = o.b))
expect_identical(sigmaOf(f.numbers, d.b), lm.numbers())
expect_identical(
  sigmaOf(f.numbers, d.b, weights = w.b),
  lm.numbers(weights = w.b)
)
expect_identical(
  sigmaOf(f.numbers, d.b, offset = o.b),
  lm.numbers(offset = o.b)
)
# three forests, one of each kind beside the forest with no basis
expect_identical(
  sigmaOf(
    y ~ forest(x1 + x2) + forest(x1, basis = factor(z)) + forest(x2, basis = w),
    d.b
  ),
  summary(lm(y ~ x1 + x2 + factor(z) + w, d.b))$sigma
)
# a data object, and a bart fit
expect_identical(
  sigmaOf(
    dbartsData(y ~ x1 + x2, d.b),
    forests = list(forest(), forest(basis = factor(z)))
  ),
  lm.levels()
)
# the estimator's own default design reads the bases as the fit's does
expect_identical(
  dbarts:::estimateSigmaFromLinearModel(
    dbarts(
      dbartsData(y ~ x1 + x2, d.b, keep.data = FALSE),
      forests = list(forest(), forest(basis = factor(d.b$z))),
      control = control.ss
    )$data
  ),
  lm.levels()
)
fit.b <- bart(
  f.levels,
  d.b,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(fit.b$sigest, lm.levels())

# a basis already among the predictors adds nothing to their span: the
# estimate is the fit's with no basis, bitwise (bartCause's BCF puts z
# among x and factor(z) among the bases). The rank the routine returns, not
# the columns it was given, sets the residual degrees of freedom
expect_identical(
  sigmaOf(y ~ forest(x1 + x2 + z) + forest(x1 + x2, basis = factor(z)), d.b),
  sigmaOf(y ~ x1 + x2 + z, d.b)
)
expect_identical(
  sigmaOf(y ~ forest(x1 + x2 + w) + forest(x1 + x2, basis = w), d.b),
  sigmaOf(y ~ x1 + x2 + w, d.b)
)

# the sparse routine, reached by sigest = "sparse" (a fit with a basis
# refuses a sparse x), agrees with lm to rounding, with weights and with an
# offset, the aliased basis included
expect_equal(
  sigmaOf(f.levels, d.b, sigest = "sparse"),
  lm.levels(),
  tolerance = 1e-12
)
expect_equal(
  sigmaOf(f.levels, d.b, sigest = "sparse", weights = w.b),
  lm.levels(weights = w.b),
  tolerance = 1e-12
)
expect_equal(
  sigmaOf(f.numbers, d.b, sigest = "sparse", offset = o.b),
  lm.numbers(offset = o.b),
  tolerance = 1e-12
)
expect_equal(
  sigmaOf(
    y ~ forest(x1 + x2 + z) + forest(x1 + x2, basis = factor(z)),
    d.b,
    sigest = "sparse"
  ),
  summary(lm(y ~ x1 + x2 + z, d.b))$sigma,
  tolerance = 1e-12
)

# dbartsSpec installs a declaration's bases before the estimate: a forests
# list over a data object, its basis a value or a column of the data
if (requireNamespace("Matrix", quietly = TRUE)) {
  x.sparse <- Matrix::Matrix(
    cbind(x1 = d.b$x1 * (d.b$x1 > 0.3), x2 = d.b$x2),
    sparse = TRUE
  )
  expect_equal(
    dbartsSpec(
      dbartsData(x.sparse, d.b$y),
      forests = list(forest(), forest(basis = factor(d.b$z))),
      control = control.ss
    )$data@sigma,
    summary(lm(d.b$y ~ as.matrix(x.sparse) + factor(d.b$z)))$sigma,
    tolerance = 1e-12
  )
}
expect_identical(
  dbartsSpec(
    dbartsData(y ~ x1 + x2, d.b),
    forests = list(forest(), forest(basis = factor(d.b$z))),
    control = control.ss
  )$data@sigma,
  lm.levels()
)
expect_identical(
  dbartsSpec(
    dbartsData(y ~ x1 + x2, d.b),
    forests = list(forest(), forest(basis = w)),
    control = control.ss
  )$data@sigma,
  lm.numbers()
)

# a stated sigest is taken as given
expect_identical(sigmaOf(f.levels, d.b, sigest = 1.5), 1.5)

# verbose counts the basis among the regression's columns, and the fallback
# line names it
announced <- function(...) {
  lines <- character()
  withCallingHandlers(
    do.call(dbarts, list(..., control = control.ss, verbose = TRUE)),
    dbartsStartingSigmaMessage = function(m) {
      lines <<- c(lines, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  lines
}
invisible(capture.output(lines.b <- announced(f.levels, d.b)))
expect_identical(
  lines.b,
  paste0(
    "estimating the starting sigma by a dense linear regression on 100 ",
    "rows and 3 columns; see 'sigest'\n"
  )
)
d.small <- d.b[1:5, ]
d.small$x3 <- runif(5L)
d.small$z <- c(0, 1, 0, 1, 1)
invisible(capture.output(
  lines.small <- announced(
    y ~ forest(x1 + x2 + x3) + forest(x1, basis = factor(z)),
    d.small
  )
))
expect_identical(length(lines.small), 2L)
expect_true(grepl(
  "on 4 columns (factors as indicators, each forest's basis among them) leaves no residual degrees of freedom in 5 rows",
  lines.small[2L],
  fixed = TRUE
))

# The floor that keeps an exact fit's rounding noise off the estimate is
# 2^16 eps relative to the estimate's own scale: the largest |y| over the rows
# in, each times the square root of its weight under precision weights; 2^16
# eps itself only for an all-zero response. A sound estimate in any units is
# lm's.
set.seed(17)
x.fl <- matrix(rnorm(60L), 30L)
y.fl <- drop(x.fl %*% c(1, -2)) + rnorm(30L)
exact.fl <- 1 + 2 * x.fl[, 1L]
w.fl <- rep(c(1, 4, 9), 10L)
counts.fl <- rep(c(1, 2, 3), 10L)
sigma.fl <- function(y, ...) {
  dbarts(x.fl, y, ..., control = control.ss, n.trees = 5L)$data@sigma
}
eps.fl <- 2^16 * .Machine$double.eps
for (scale in c(1e-9, 1e-12)) {
  tiny <- scale * y.fl
  expect_identical(sigma.fl(tiny), summary(lm(tiny ~ x.fl))$sigma)
  expect_identical(
    sigma.fl(tiny, weights = w.fl),
    summary(lm(tiny ~ x.fl, weights = w.fl))$sigma
  )
  expect_identical(
    sigma.fl(scale * exact.fl),
    eps.fl * max(abs(scale * exact.fl))
  )
}
# the weighted scale differs from the unweighted one here
expect_true(max(sqrt(w.fl) * abs(exact.fl)) > 1.5 * max(abs(exact.fl)))
expect_identical(
  sigma.fl(exact.fl, weights = w.fl),
  eps.fl * max(sqrt(w.fl) * abs(exact.fl))
)
# counts write rows out and leave each row's scale as it is
expect_identical(
  sigma.fl(exact.fl, weights = counts.fl, weights.type = "frequency"),
  eps.fl * max(abs(exact.fl))
)
expect_identical(sigma.fl(rep(0, 30L)), eps.fl)
# a spread of about 1 beside a mean of 1.7e9 is a sound estimate, 47 times
# the floor and far above an exact fit's rounding noise: lm's sigma
y.large <- 1.7e9 + y.fl
expect_true(summary(lm(y.large ~ x.fl))$sigma < 2)
expect_identical(sigma.fl(y.large), summary(lm(y.large ~ x.fl))$sigma)

# A fold's floor reads its training rows on the same scale: every fold's
# regression is exact, and under precision weights each floor is the largest
# sqrt(w) |y| over the fold's training rows
foldFloors <- function(rule) {
  original <- dbarts:::floorSigmaEstimate
  calls <- list()
  recording <- function(sigma, residual, weights = NULL, frequency = FALSE) {
    result <- original(sigma, residual, weights, frequency)
    calls[[length(calls) + 1L]] <<- list(
      residual = residual,
      weights = weights,
      result = result
    )
    result
  }
  assignInNamespace("floorSigmaEstimate", recording, "dbarts")
  on.exit(assignInNamespace("floorSigmaEstimate", original, "dbarts"))
  invisible(xbart(
    x.fl,
    exact.fl,
    weights = w.fl,
    sigest = rule,
    n.reps = 1L,
    n.test = 3L,
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 5L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  ))
  Filter(function(call) length(call$residual) < 30L, calls)
}
for (rule in c("dense", "sparse")) {
  folds <- foldFloors(rule)
  expect_equal(length(folds), 3L, info = rule)
  for (fold in folds) {
    expect_identical(
      fold$result,
      eps.fl * max(sqrt(fold$weights) * abs(fold$residual)),
      info = rule
    )
  }
}
