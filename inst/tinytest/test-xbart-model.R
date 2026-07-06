source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that works with non-standard models
x <- testData$x
y <- testData$y

k <- c(4, 8)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    k = k, n.threads = 1L, resid.prior = chisq(2.5, 0.9)
  )
)

expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    k = k, n.threads = 1L, resid.prior = fixed(2)
  )
)

n.trees <- c(5L, 10L)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    n.trees = n.trees, n.threads = 1L
  )
)

# split.probs and dart reach the tree prior, which used to be hardcoded to a
# uniform cgm; both fix or sample the split-variable probabilities while power
# and base stay the swept grid
p <- ncol(x)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 2L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    n.threads = 1L, split.probs = c(5, rep(1, p - 1L))
  )
)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 2L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    n.threads = 1L, dart = TRUE
  )
)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 2L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    n.threads = 1L, dart = dbarts::dbartsPriors$dart(a = 1)
  )
)
expect_error(
  dbarts::xbart(x, y, n.reps = 1L, n.threads = 1L,
                dart = TRUE, split.probs = c(5, rep(1, p - 1L))),
  pattern = "cannot be combined with 'dart'")
expect_error(
  dbarts::xbart(x, y, n.reps = 1L, n.threads = 1L, dart = "yes"),
  pattern = "must be TRUE, FALSE")
expect_error(
  dbarts::xbart(x, y, n.reps = 1L, n.threads = 1L, split.probs = c(1, 2)),
  pattern = "does not equal number of columns")

rm(p, n.trees, k, y, x)

rm(testData)

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that runs with binary data and k hyperprior
x <- testData$X
z <- testData$Z

n.reps  <- 3L
power   <- c(1.5, 2)

xval <- dbarts::xbart(
  x, z, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps,
  power = power,
  n.threads = 2L
)

expect_inherits(xval, "matrix")
expect_equal(dim(xval), c(n.reps, length(power)))
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    power   = as.character(power)
  )
)

rm(xval)

# family routes through to the folds: logistic and forced-gaussian fit 0/1
# responses, binary families reject continuous ones
xval.logistic <- dbarts::xbart(
  x, z, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps, n.threads = 1L, family = "logistic"
)
expect_equal(dim(xval.logistic), NULL)
expect_true(!anyNA(xval.logistic))

xval.gaussian <- dbarts::xbart(
  x, z, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps, n.threads = 1L, family = "gaussian"
)
expect_true(!anyNA(xval.gaussian))

expect_error(
  dbarts::xbart(x, rnorm(nrow(x)), n.reps = 1L, n.threads = 1L,
                family = "probit"),
  pattern = "requires a response coded 0/1"
)

rm(xval.logistic, xval.gaussian, power, n.reps, z, x)

rm(testData)

