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
  n.trees = 5L,
  control = dbarts::dbartsControl(n.threads = 1L, verbose = FALSE)
)
expect_inherits(sampler, "dbartsSampler")
expect_inherits(
  fitWith(matrix(cst, ncol = 1L), y, samplerOnly = TRUE),
  "dbartsSampler"
)

# bartBT, the BayesTree-style door, refuses as bart does; sampleronly returns
legacy <- function(...) {
  dbarts::bartBT(
    matrix(cst, ncol = 1L),
    y,
    ntree = 5L,
    ndpost = 5L,
    nskip = 5L,
    verbose = FALSE,
    ...
  )
}
expect_error(legacy(), "every predictor holds one")
expect_inherits(legacy(sampleronly = TRUE), "dbartsSampler")

# the other response families refuse too
cats <- factor(rep(c("a", "b", "c"), length.out = n))
counts <- rpois(n, 3)
ordinal <- rep(1:3, length.out = n)
expect_error(
  fitWith(
    cats ~ cst,
    data.frame(cats = cats, cst = cst),
    family = "multinomial"
  ),
  "every predictor holds one"
)
expect_error(
  fitWith(
    ordinal ~ cst,
    data.frame(ordinal = ordinal, cst = cst),
    family = "ordinal"
  ),
  "every predictor holds one"
)
expect_error(
  fitWith(
    counts ~ cst,
    data.frame(counts = counts, cst = cst),
    family = "nbinom"
  ),
  "every predictor holds one"
)

# a sparse column holds its zero where it stores fewer entries than rows, so
# one stored value on some rows is two values; on every row it is one
some <- Matrix::sparseMatrix(i = 1:10, j = rep(1L, 10L), x = 1, dims = c(n, 1L))
expect_inherits(fitWith(some, y), "bart")
every <- Matrix::sparseMatrix(i = 1:n, j = rep(1L, n), x = 1, dims = c(n, 1L))
expect_error(fitWith(every, y), "every predictor holds one")

# a missing value is a value of its own: one value and a missing one can split
withNA <- cst
withNA[1:5] <- NA
expect_inherits(fitWith(matrix(withNA, ncol = 1L), y), "bart")

# a fit with a basis forest is not refused, whatever its predictors hold
d <- data.frame(y = y, cst = cst, z = rnorm(n))
expect_inherits(fitWith(y ~ forest(cst) + forest(cst, basis = z), d), "bart")
