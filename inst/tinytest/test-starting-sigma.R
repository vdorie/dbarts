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
# the data object's own bases, and a bart fit
expect_identical(
  sigmaOf(dbartsData(
    y ~ x1 + x2,
    d.b,
    bases = list(NULL, factor(d.b$z))
  )),
  lm.levels()
)
# the estimator's own default design reads the bases as the fit's does
expect_identical(
  dbarts:::estimateSigmaFromLinearModel(
    dbartsData(
      y ~ x1 + x2,
      d.b,
      bases = list(NULL, factor(d.b$z)),
      keep.data = FALSE
    )
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
# list over a data object with none, and one replacing the object's own
if (requireNamespace("Matrix", quietly = TRUE)) {
  x.sparse <- Matrix::Matrix(
    cbind(x1 = d.b$x1 * (d.b$x1 > 0.3), x2 = d.b$x2),
    sparse = TRUE
  )
  expect_equal(
    dbartsSpec(
      dbartsData(x.sparse, d.b$y, bases = list(NULL, factor(d.b$z))),
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
    dbartsData(
      y ~ x1 + x2,
      d.b,
      bases = list(NULL, factor(d.b$z))
    ),
    forests = list(forest(), forest(basis = d.b$w)),
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
