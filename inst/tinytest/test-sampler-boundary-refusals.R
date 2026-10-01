# Arguments the sampler's methods refuse at the R and C boundary, and what a
# refusal leaves behind: the sampler a refused call reached must be the one it
# was before the call, mirrors included.

set.seed(1)
n <- 20L
x <- matrix(rnorm(n * 2L), n)
y <- rnorm(n)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  keepTrees = TRUE,
  updateState = FALSE,
  verbose = FALSE
)
sampler <- dbarts(x, y, control = control)

# --- run counts: a negative count wrapped to a huge size and reported draws
# the run never recorded; 0 + 0 returned NULL. Both refused as in 0.9-x, at
# the method and at the C entry beneath it
expect_error(
  sampler$run(-2L, 5L),
  "number of burn-in steps must be greater than or equal to 0"
)
expect_error(
  sampler$run(3L, -1L),
  "number of samples must be greater than or equal to 0"
)
expect_error(
  sampler$run(0L, 0L),
  "either number of burn-in or samples must be positive"
)
expect_error(sampler$run(2.5, 1L), "whole number")
expect_error(sampler$run(c(1L, 2L), 1L), "'numBurnIn' must be a single")
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    sampler$getPointer(),
    -1L,
    2L,
    NULL,
    NULL,
    TRUE
  ),
  "number of burn-in steps must be greater than or equal to 0"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    sampler$getPointer(),
    0L,
    0L,
    NULL,
    NULL,
    TRUE
  ),
  "either number of burn-in or samples must be positive"
)
# a refused run records nothing: the store still holds no draws
expect_error(sampler$predict(x), "holds no recorded draws")
r <- sampler$run(2L, 5L)
expect_true(all(is.finite(r$sigma)) && all(r$sigma > 1e-100))

# --- missing indices are refused by name, not by a bare missing condition
expect_error(
  sampler$getTrees(chainNums = NA_integer_),
  "'chainNums' contains missing values"
)
expect_error(
  sampler$getTrees(sampleNums = NA_integer_),
  "'sampleNums' contains missing values"
)
expect_error(
  sampler$getTrees(treeNums = NA_integer_),
  "'treeNums' contains missing values"
)
expect_error(
  sampler$setPredictor(x[, 1L], NA_integer_),
  "'column' contains missing values"
)

# --- non-finite swaps: one infinite value left sigma and every fit NaN for
# good, even after the value was put back; refused as creation refuses one,
# and the sampler is unchanged
expect_error(
  sampler$setResponse(replace(y, 1L, Inf)),
  "response contains non-finite values"
)
expect_error(
  sampler$setOffset(replace(numeric(n), 1L, -Inf)),
  "'offset' contains non-finite values"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    replace(y, 1L, Inf),
    FALSE,
    NULL
  ),
  "response contains non-finite values"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setOffset,
    sampler$getPointer(),
    replace(numeric(n), 2L, NaN),
    FALSE
  ),
  "offset contains non-finite values"
)
expect_identical(sampler$data@y, y)
expect_null(sampler$data@offset)
r <- sampler$run(0L, 2L)
expect_true(all(is.finite(r$sigma)) && all(is.finite(r$train)))

# --- updateScale is a single TRUE or FALSE: 1 rescaled the engine while the R
# side skipped its own rescaling steps
expect_error(
  sampler$setOffset(numeric(n), updateScale = 1),
  "'updateScale' must be TRUE or FALSE"
)
expect_error(
  sampler$setResponse(y, updateScale = NA),
  "'updateScale' must be TRUE or FALSE"
)

# --- a refused offset swap leaves the mirror a save and load re-creates from
set.seed(3)
nb <- 60L
xb <- matrix(runif(nb * 2L), nb)
zb <- rbinom(nb, 1L, 0.5)
yb <- xb[, 1L] + zb + rnorm(nb, sd = 0.2)
bcf <- dbarts(
  xb,
  yb,
  forests = list(forest(), forest(basis = ~ factor(zb), n.trees = 10L)),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    verbose = FALSE
  )
)
bcf$setOffset(rep(0.5, nb), updateScale = FALSE)
expect_error(
  bcf$setOffset(rep(5, nb), updateScale = NA),
  "'updateScale' must be TRUE or FALSE"
)
expect_identical(bcf$data@offset, rep(0.5, nb))
