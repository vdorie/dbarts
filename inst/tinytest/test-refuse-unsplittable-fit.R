# bart(), xbart() and rbart_vi() refuse a fit in which every predictor holds
# one value on the rows kept and no forest has a basis, naming the columns;
# dbarts() and new() take such a sampler, a fit with one splittable column or
# with a basis forest is fitted, and samplerOnly returns the sampler.

set.seed(23L)
n <- 40L
y <- rnorm(n)
cst <- rep(1, n)
x1 <- rnorm(n)
args <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
fitWith <- function(...) do.call(dbarts::bart, c(list(...), args))

# a matrix of constant columns, one or several, names reported
expect_error(fitWith(matrix(cst, ncol = 1L), y), "every predictor holds one")
cm <- cbind(a = cst, b = 2)
expect_error(fitWith(cm, y), "'a', 'b'")
expect_error(fitWith(y ~ cst, data.frame(y = y, cst = cst)), "'cst'")

# sparse storage: an all-zero column is one value
sp <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  x = numeric(0),
  dims = c(n, 2L)
)
expect_error(fitWith(sp, y), "every predictor holds one")

# one splittable column among constants is fitted
expect_inherits(fitWith(cbind(cst, x1), y), "bart")

# the other doors
expect_error(
  dbarts::xbart(
    y ~ cst,
    data.frame(y = y, cst = cst),
    n.trees = 5L,
    n.reps = 2L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "every predictor holds one"
)
expect_error(
  dbarts::rbart_vi(
    cm,
    y,
    group.by = factor(rep(1:4, length.out = n)),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "every predictor holds one"
)

# the sampler door keeps taking it, as does samplerOnly
sampler <- dbarts::dbarts(
  matrix(cst, ncol = 1L),
  y,
  control = dbarts::dbartsControl(
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_inherits(sampler, "dbartsSampler")
expect_inherits(
  fitWith(matrix(cst, ncol = 1L), y, samplerOnly = TRUE),
  "dbartsSampler"
)
